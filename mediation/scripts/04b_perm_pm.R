#!/usr/bin/env Rscript
# ==============================================================================
# 04b_perm_pm.R -- a usable permutation null for the HEADLINE quantity.
#
# WHY 04's VERSION DOES NOT WORK
#
# 04_permutation.R calibrated the median proportion mediated by recomputing it
# inside each permuted replicate over that replicate's OWN significant genes.
# Under the null almost nothing is significant, so the statistic was defined in
# only 32 of 1000 replicates (28 of 1000 for the proposal tree). A null built on
# 3% of replicates, selected precisely because they happened to produce
# significant genes, is not a null worth quoting.
#
# WHAT THIS DOES INSTEAD
#
# Holds the gene set FIXED at the genes that were significant in the real
# analysis, and asks: when lithium is known to do nothing, what median
# proportion mediated do those same genes produce? That is the apples-to-apples
# comparison, it is defined in every replicate, and it answers the question that
# actually matters -- whether a median proportion mediated near 0.4 is more than
# a ratio of two noisy quantities will give you by itself.
#
# B = 300 rather than 1000: the quantity is a median over >1,000 genes and so is
# far more stable than a discovery count, and 300 resolves a one-sided p to
# 0.003, which is ample for this comparison.
# ==============================================================================

ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "mediation")); set.seed(4812)
suppressPackageStartupMessages({library(edgeR); library(limma); library(parallel)})
source("scripts/med.R"); source("scripts/boot_common.R")
options(width = 150); say <- function(...) cat(sprintf(...))

B <- 300L; NCORE <- 4L; BLOCK <- 50L
VAR <- list(lineage = "data/med_input.rds", proposal = "data/med_input_proposal.rds")
prep <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
OUT <- list()

for (vn in names(VAR)) {
  d   <- readRDS(VAR[[vn]])
  CNT <- prep$counts[d$gene_main, d$meta$title, drop = FALSE]
  storage.mode(CNT) <- "integer"
  M <- d$M; C <- d$C; x <- d$lithium
  fit <- readRDS(sprintf("results/primary/%s_fit.rds", vn)); tb <- fit$tb
  sigset <- which(tb$TE_q_BH < 0.05)
  obs <- median(tb$PM[sigset])
  say("\n== %s == fixed gene set: %d genes | observed median PM = %.4f\n",
      vn, length(sigset), obs)

  set.seed(4812)
  PRM <- lapply(seq_len(B), function(b) sample(x))
  cl <- make_med_cluster(NCORE, ROOT, c("CNT", "M", "C", "sigset"))
  on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
  pm <- rep(NA_real_, B); t0 <- Sys.time()
  for (blk in seq_len(ceiling(B / BLOCK))) {
    rng <- ((blk - 1L) * BLOCK + 1L):min(blk * BLOCK, B)
    res <- parLapply(cl, PRM[rng], function(xp)
      tryCatch({ r <- med_replicate(seq_len(ncol(CNT)), xp, CNT, M, C)
                 median(r$NIE[sigset] / r$TE[sigset]) }, error = function(e) NA_real_))
    pm[rng] <- unlist(res)
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    say("   %3d/%d | %.1f min | ~%.1f min left\n", max(rng), B, el,
        el / max(rng) * (B - max(rng)))
  }
  stopCluster(cl)
  pm <- pm[is.finite(pm)]
  p_two <- (1 + sum(abs(pm) >= abs(obs))) / (1 + length(pm))
  say("   null median PM: median %.4f | 2.5-97.5%% %.4f to %.4f | max |.| %.4f\n",
      median(pm), quantile(pm, .025), quantile(pm, .975), max(abs(pm)))
  say("   observed %.4f -> two-sided permutation p = %.4f  (%d usable replicates)\n",
      obs, p_two, length(pm))
  say("   -> %s\n", if (p_two < 0.05)
    "the mediated proportion is BEYOND what a ratio of noise produces" else
    "NOT distinguishable from a ratio of noise on this statistic")
  OUT[[vn]] <- list(obs = obs, null = pm, p = p_two, n_sig = length(sigset))
}
saveRDS(OUT, "results/perm/perm_pm_fixedset.rds")
write.csv(data.frame(tree = names(OUT),
                     observed = sapply(OUT, `[[`, "obs"),
                     null_median = sapply(OUT, function(o) median(o$null)),
                     null_lo = sapply(OUT, function(o) quantile(o$null, .025)),
                     null_hi = sapply(OUT, function(o) quantile(o$null, .975)),
                     perm_p = sapply(OUT, `[[`, "p")),
          "results/perm/perm_pm_fixedset.csv", row.names = FALSE)
cat("\n04b complete.\n")
