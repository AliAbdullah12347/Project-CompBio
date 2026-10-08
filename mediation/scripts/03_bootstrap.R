#!/usr/bin/env Rscript
# ==============================================================================
# 03_bootstrap.R -- nonparametric bootstrap over SUBJECTS, B = 1000, for both
# mediator definitions.
#
# WHY THIS IS NOT OPTIONAL HERE
#
# 00_validate_med.R found the delta-method standard error to be about 5% below
# the truth on simulated data (95% intervals covering 93.0% and 94.4%), and
# 02b_diagnose.R found it about 10% below a direct bootstrap on the real data.
# The delta method is therefore mildly ANTICONSERVATIVE in this setting, so the
# bootstrap interval is the one to trust when the two disagree.
#
# The resampling unit is the SUBJECT, which is correct here because subjects are
# the independent units; each subject carries its exposure, its four or five
# balances, its covariates and its whole expression profile together.
#
# NDE is not stored: NDE = TE - NIE exactly, so it is reconstructed afterwards
# and 96 MB of memory is not spent twice.
# ==============================================================================

ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "mediation")); set.seed(481)
suppressPackageStartupMessages({library(edgeR); library(limma); library(parallel)})
source("scripts/med.R"); source("scripts/boot_common.R")
options(width = 150)
say <- function(...) cat(sprintf(...))

B      <- 1000L
NCORE  <- 4L
BLOCK  <- 50L                      # checkpoint every 50 replicates
VAR    <- list(lineage = "data/med_input.rds", proposal = "data/med_input_proposal.rds")
dir.create("results/boot", showWarnings = FALSE, recursive = TRUE)

prep <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))

for (vn in names(VAR)) {
  ckpt <- sprintf("results/boot/%s_boot.rds", vn)
  if (file.exists(ckpt)) { say("%s: checkpoint exists, skipping\n", vn); next }
  d   <- readRDS(VAR[[vn]])
  CNT <- prep$counts[d$gene_main, d$meta$title, drop = FALSE]
  storage.mode(CNT) <- "integer"
  M <- d$M; C <- d$C; x <- d$lithium; n <- nrow(M); G <- nrow(CNT)
  fit <- readRDS(sprintf("results/primary/%s_fit.rds", vn)); tb <- fit$tb
  # PRESPEC section 1.4 defines the headline estimand as the MEDIAN proportion
  # mediated among genes with a significant TOTAL effect. The gene set is fixed
  # from the primary analysis rather than re-selected inside each replicate, so
  # the bootstrap carries uncertainty in the EFFECTS, not in a moving selection
  # rule, which would be a different and much harder estimand.
  sigset <- which(tb$TE_q_BH < 0.05)
  say("\n== bootstrap: %s == %d genes x %d subjects, B = %d on %d cores\n",
      vn, G, n, B, NCORE)

  # the index sets are drawn ONCE, in the master, from the project seed, so the
  # run is reproducible regardless of how work lands on the workers
  set.seed(481)
  IDX <- lapply(seq_len(B), function(b) sample.int(n, n, replace = TRUE))

  cl <- make_med_cluster(NCORE, ROOT, c("CNT", "M", "C", "x"))
  on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)

  NIE <- matrix(NA_real_, B, G); TEm <- matrix(NA_real_, B, G)
  B1  <- matrix(NA_real_, B, ncol(M)); PMrep <- rep(NA_real_, B)
  nok <- 0L; t0 <- Sys.time()
  for (blk in seq_len(ceiling(B / BLOCK))) {
    rng <- ((blk - 1L) * BLOCK + 1L):min(blk * BLOCK, B)
    res <- parLapply(cl, IDX[rng], function(i)
      tryCatch(med_replicate(i, x[i], CNT, M, C), error = function(e) NULL))
    for (j in seq_along(rng)) {
      r <- res[[j]]; if (is.null(r)) next
      NIE[rng[j], ] <- r$NIE; TEm[rng[j], ] <- r$TE; B1[rng[j], ] <- r$b1
      PMrep[rng[j]] <- median(r$NIE[sigset] / r$TE[sigset])
      nok <- nok + 1L
    }
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    say("   %4d/%d done | %.1f min elapsed | ~%.1f min left\n",
        max(rng), B, el, el / max(rng) * (B - max(rng)))
  }
  stopCluster(cl)
  say("   %d of %d replicates usable\n", nok, B)

  ok <- !is.na(NIE[, 1])
  NIE <- NIE[ok, , drop = FALSE]; TEm <- TEm[ok, , drop = FALSE]; B1 <- B1[ok, , drop = FALSE]
  PMrep <- PMrep[ok]
  NDE <- TEm - NIE

  qs <- function(Mx) {
    list(se  = apply(Mx, 2, sd),
         lo  = apply(Mx, 2, quantile, 0.025, names = FALSE),
         hi  = apply(Mx, 2, quantile, 0.975, names = FALSE),
         med = apply(Mx, 2, median))
  }
  sN <- qs(NIE); sT <- qs(TEm); sD <- qs(NDE)

  # bias-corrected percentile: z0 from the share of replicates below the point
  # estimate (Efron 1987). Reported alongside, not instead of, the percentile.
  z0  <- qnorm(pmin(pmax(colMeans(sweep(NIE, 2, tb$NIE, "<")), 1 / nrow(NIE)),
                    1 - 1 / nrow(NIE)))
  al  <- c(0.025, 0.975)
  bcl <- sapply(seq_len(ncol(NIE)), function(j)
    quantile(NIE[, j], pnorm(2 * z0[j] + qnorm(al[1])), names = FALSE))
  bch <- sapply(seq_len(ncol(NIE)), function(j)
    quantile(NIE[, j], pnorm(2 * z0[j] + qnorm(al[2])), names = FALSE))

  out <- data.frame(
    gene = tb$gene,
    NIE = tb$NIE, NIE_se_delta = tb$NIE_se, NIE_se_boot = sN$se,
    NIE_lo = sN$lo, NIE_hi = sN$hi, NIE_lo_bc = bcl, NIE_hi_bc = bch,
    TE  = tb$TE,  TE_se_delta  = tb$TE_se,  TE_se_boot  = sT$se,
    TE_lo = sT$lo, TE_hi = sT$hi,
    NDE = tb$NDE, NDE_se_boot = sD$se, NDE_lo = sD$lo, NDE_hi = sD$hi,
    stringsAsFactors = FALSE)
  # bootstrap p-value: invert the percentile interval (Hall 1992) -- the
  # smallest alpha at which the interval excludes zero
  pb <- 2 * pmin(colMeans(NIE <= 0), colMeans(NIE >= 0))
  out$NIE_p_boot <- pmax(pb, 1 / nrow(NIE))
  out$NIE_q_boot_BH <- p.adjust(out$NIE_p_boot, "BH")
  out$NIE_sig_boot <- out$NIE_lo > 0 | out$NIE_hi < 0

  pm_obs <- median(tb$PM[sigset]); pm_ci <- quantile(PMrep, c(.025, .975), names = FALSE)
  say("\n   HEADLINE -- median proportion mediated over %d genes with significant TE\n",
      length(sigset))
  say("     point estimate              : %.4f\n", pm_obs)
  say("     bootstrap percentile 95%% CI : %.4f to %.4f\n", pm_ci[1], pm_ci[2])
  say("     bootstrap SE                : %.4f\n", sd(PMrep))
  say("     replicates above 0.50       : %.1f%%  (pre-declared 'substantial')\n",
      100 * mean(PMrep > 0.5))
  say("     replicates above 0          : %.1f%%\n\n", 100 * mean(PMrep > 0))
  say("   delta SE / bootstrap SE: median %.3f (IQR %.3f-%.3f)\n",
      median(out$NIE_se_delta / out$NIE_se_boot),
      quantile(out$NIE_se_delta / out$NIE_se_boot, .25),
      quantile(out$NIE_se_delta / out$NIE_se_boot, .75))
  say("   NIE excluding zero, uncorrected percentile CI: %d genes (%.1f%%)\n",
      sum(out$NIE_sig_boot), 100 * mean(out$NIE_sig_boot))
  say("   NIE at BH 5%% on bootstrap p: %d genes\n", sum(out$NIE_q_boot_BH < 0.05))
  say("   mediator b1 across replicates (mean, and share with t > 1.96):\n")
  for (k in seq_len(ncol(B1)))
    say("     %-6s mean %+.4f  sd %.4f\n", colnames(M)[k], mean(B1[, k]), sd(B1[, k]))

  write.csv(out, sprintf("results/boot/%s_boot_ci.csv", vn), row.names = FALSE)
  saveRDS(list(summary = out, b1 = B1, nrep = nrow(NIE),
               pm_rep = PMrep, pm_obs = pm_obs, pm_ci = pm_ci, pm_se = sd(PMrep),
               n_sig = length(sigset),
               delta_vs_boot = median(out$NIE_se_delta / out$NIE_se_boot)), ckpt)
  say("   wrote %s\n", ckpt)
}
cat("\n03_bootstrap complete.\n")
