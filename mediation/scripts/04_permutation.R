#!/usr/bin/env Rscript
# ==============================================================================
# 04_permutation.R -- the permutation null, B = 1000, both mediator definitions.
#
# The proposal's flowchart asks for exactly this: "Redo everything on shuffled
# copies to get the equivalent of a permutation null and to compute the error
# bars."
#
# WHAT IS SHUFFLED, AND WHAT NULL THAT CREATES
#
# The lithium labels are permuted across subjects. Each subject keeps its own
# expression profile, its own balances and its own covariates; only who is
# called a lithium user changes. This breaks the exposure's link to BOTH the
# mediator and the outcome at once, so the null is "lithium does nothing",
# under which the true NDE, NIE and total effect are all zero.
#
# WHY THIS MATTERS MORE THAN USUAL HERE
#
# 02_effects.R produced a result that is unstable in a specific way: the two
# mediator definitions give per-gene NIE estimates correlating at 0.980, yet one
# yields 0 genes at BH 5% and the other 3,647. That is not a contradiction --
# Benjamini-Hochberg is a step-up procedure, so a diffuse shift in the whole
# p-value distribution can flip it between almost nothing and thousands. It does
# mean the DISCOVERY COUNT cannot be read at face value. The permutation null
# calibrates it: we compare the observed count against the distribution of
# counts obtained when lithium is known to do nothing.
#
# The permuted design changes, so voom is refitted against it, exactly as in the
# bootstrap.
# ==============================================================================

ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "mediation")); set.seed(481)
suppressPackageStartupMessages({library(edgeR); library(limma); library(parallel)})
source("scripts/med.R"); source("scripts/boot_common.R")
options(width = 150)
say <- function(...) cat(sprintf(...))

B <- 1000L; NCORE <- 4L; BLOCK <- 50L
VAR <- list(lineage = "data/med_input.rds", proposal = "data/med_input_proposal.rds")
dir.create("results/perm", showWarnings = FALSE, recursive = TRUE)
prep <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))

for (vn in names(VAR)) {
  ckpt <- sprintf("results/perm/%s_perm.rds", vn)
  if (file.exists(ckpt)) { say("%s: checkpoint exists, skipping\n", vn); next }
  d   <- readRDS(VAR[[vn]])
  CNT <- prep$counts[d$gene_main, d$meta$title, drop = FALSE]
  storage.mode(CNT) <- "integer"
  M <- d$M; C <- d$C; x <- d$lithium; n <- nrow(M); G <- nrow(CNT)
  fit <- readRDS(sprintf("results/primary/%s_fit.rds", vn)); tb <- fit$tb
  obs <- list(nNIE = sum(tb$NIE_q_BH < 0.05), nTE = sum(tb$TE_q_BH < 0.05),
              maxz = max(abs(tb$NIE_z)), medPM = median(tb$PM[tb$TE_q_BH < 0.05]))
  say("\n== permutation: %s == observed: NIE %d, TE %d, max|z| %.2f\n",
      vn, obs$nNIE, obs$nTE, obs$maxz)

  set.seed(481)
  PRM <- lapply(seq_len(B), function(b) sample(x))        # shuffled labels

  cl <- make_med_cluster(NCORE, ROOT, c("CNT", "M", "C"))
  on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)

  S <- data.frame(rep = seq_len(B), nNIE = NA_integer_, nTE = NA_integer_,
                  maxz = NA_real_, medPM = NA_real_, minp = NA_real_)
  B1T <- matrix(NA_real_, B, ncol(M))
  exceed <- integer(G)                      # per-gene |null NIE| >= |observed|
  nok <- 0L; t0 <- Sys.time()
  for (blk in seq_len(ceiling(B / BLOCK))) {
    rng <- ((blk - 1L) * BLOCK + 1L):min(blk * BLOCK, B)
    res <- parLapply(cl, PRM[rng], function(xp)
      tryCatch(med_replicate(seq_len(ncol(CNT)), xp, CNT, M, C, with_p = TRUE),
               error = function(e) NULL))
    for (j in seq_along(rng)) {
      r <- res[[j]]; if (is.null(r)) next
      i <- rng[j]
      qN <- p.adjust(r$NIE_p, "BH"); qT <- p.adjust(r$TE_p, "BH")
      S$nNIE[i] <- sum(qN < 0.05); S$nTE[i] <- sum(qT < 0.05)
      S$maxz[i] <- max(abs(qnorm(r$NIE_p / 2)))
      S$minp[i] <- min(r$NIE_p)
      sigT <- qT < 0.05
      S$medPM[i] <- if (any(sigT)) median((r$NIE / r$TE)[sigT]) else NA_real_
      B1T[i, ] <- r$b1_t
      exceed <- exceed + as.integer(abs(r$NIE) >= abs(tb$NIE))
      nok <- nok + 1L
    }
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    say("   %4d/%d | %.1f min elapsed | ~%.1f min left | null BH-sig so far: median %.0f, max %.0f\n",
        max(rng), B, el, el / max(rng) * (B - max(rng)),
        median(S$nNIE, na.rm = TRUE), max(S$nNIE, na.rm = TRUE))
  }
  stopCluster(cl)
  S <- S[!is.na(S$nNIE), ]; nrep <- nrow(S)
  say("   %d of %d permutations usable\n", nrep, B)

  # one-sided permutation p-values: signal can only ADD discoveries
  pval <- function(obsv, nullv) (1 + sum(nullv >= obsv)) / (1 + length(nullv))
  res <- data.frame(
    statistic = c("genes with NIE at BH 5%", "genes with TE at BH 5%",
                  "max |z| for NIE", "median proportion mediated"),
    observed  = c(obs$nNIE, obs$nTE, obs$maxz, obs$medPM),
    null_median = c(median(S$nNIE), median(S$nTE), median(S$maxz),
                    median(S$medPM, na.rm = TRUE)),
    null_p95  = c(quantile(S$nNIE, .95), quantile(S$nTE, .95),
                  quantile(S$maxz, .95), quantile(S$medPM, .95, na.rm = TRUE)),
    null_max  = c(max(S$nNIE), max(S$nTE), max(S$maxz), max(S$medPM, na.rm = TRUE)),
    perm_p    = c(pval(obs$nNIE, S$nNIE), pval(obs$nTE, S$nTE),
                  pval(obs$maxz, S$maxz), pval(obs$medPM, na.omit(S$medPM))),
    stringsAsFactors = FALSE)
  print(res, row.names = FALSE)

  # per-gene empirical p and a Storey-Tibshirani style empirical FDR
  gene_p <- (1 + exceed) / (1 + nrep)
  gq <- p.adjust(gene_p, "BH")
  say("   per-gene empirical p: min %.4f | genes at empirical BH 5%%: %d\n",
      min(gene_p), sum(gq < 0.05))
  say("   mediator b1 t-statistics under the null (|t| > 1.96 rate per balance):\n")
  for (k in seq_len(ncol(M)))
    say("     %-6s %.1f%%  (observed t = %.2f, perm p = %.4f)\n", colnames(M)[k],
        100 * mean(abs(B1T[seq_len(nrep), k]) > 1.96),
        fit$med_tab$t[k], pval(abs(fit$med_tab$t[k]), abs(B1T[seq_len(nrep), k])))

  write.csv(res, sprintf("results/perm/%s_perm_summary.csv", vn), row.names = FALSE)
  write.csv(data.frame(gene = tb$gene, NIE = tb$NIE, emp_p = gene_p, emp_q = gq),
            sprintf("results/perm/%s_perm_gene.csv", vn), row.names = FALSE)
  saveRDS(list(summary = res, per_rep = S, b1_t = B1T[seq_len(nrep), , drop = FALSE],
               gene_p = gene_p, gene_q = gq, nrep = nrep), ckpt)
  say("   wrote %s\n", ckpt)
}
cat("\n04_permutation complete.\n")
