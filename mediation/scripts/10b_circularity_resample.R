#!/usr/bin/env Rscript
# ==============================================================================
# 10b_circularity_resample.R -- intervals and nulls for the circularity sets.
#
# For each tree, from the SAME common fit as 10_circularity.R (TMM + voom on all
# 12,368 genes, refitted inside every replicate):
#
#   bootstrap   B = 1000 subjects-with-replacement. Each replicate gives the
#               median proportion mediated over each set's FIXED significant-TE
#               genes. Because C1, C2 and C3 come from the same replicate, the
#               differences C2 - C1 and C3 - C1 get paired intervals.
#   permutation B = 300 lithium shuffles, same statistic, same fixed gene sets
#               (the 04b design: defined in every replicate).
#
# Also keeps per-gene bootstrap intervals for all 12,368 genes, so marker genes
# get intervals too (the original 03 bootstrap covered non-markers only).
# ==============================================================================

ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "mediation")); set.seed(481)
suppressPackageStartupMessages({library(edgeR); library(limma); library(parallel)})
source("scripts/med.R"); source("scripts/boot_common.R")
options(width = 150); say <- function(...) cat(sprintf(...))

B_BOOT <- as.integer(Sys.getenv("B_BOOT", "1000")); B_PERM <- as.integer(Sys.getenv("B_PERM", "300"))
NCORE <- 4L; BLOCK <- as.integer(Sys.getenv("BLOCK", "50"))
TAG <- Sys.getenv("TAG", "")      # non-empty only for a dry run, so it never overwrites real output
VAR <- list(lineage = "data/med_input.rds", proposal = "data/med_input_proposal.rds")
prep <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))

run_blocks <- function(cl, tasks, fun, label) {
  out <- vector("list", length(tasks)); t0 <- Sys.time()
  for (blk in seq_len(ceiling(length(tasks) / BLOCK))) {
    rng <- ((blk - 1L) * BLOCK + 1L):min(blk * BLOCK, length(tasks))
    out[rng] <- parLapply(cl, tasks[rng], fun)
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    say("   %s %4d/%d | %.1f min | ~%.1f min left\n", label, max(rng), length(tasks), el,
        el / max(rng) * (length(tasks) - max(rng)))
  }
  out
}

for (vn in names(VAR)) {
  ckpt <- sprintf("results/circularity/%s_resample%s.rds", vn, TAG)
  if (file.exists(ckpt)) { say("%s: checkpoint exists, skipping\n", vn); next }
  d  <- readRDS(VAR[[vn]])
  cf <- readRDS(sprintf("results/circularity/%s_circularity.rds", vn))
  CNT <- prep$counts[cf$tb$gene, d$meta$title, drop = FALSE]
  storage.mode(CNT) <- "integer"
  M <- d$M; C <- d$C; x <- d$lithium; n <- length(x)
  SIG <- cf$sig
  obs <- sapply(SIG, function(s) median(cf$tb$PM[s]))
  say("\n== %s == fixed significant-TE sets: C1 %d, C2 %d, C3 %d | observed median PM %s\n",
      vn, length(SIG$C1), length(SIG$C2), length(SIG$C3),
      paste(sprintf("%.4f", obs), collapse = " / "))

  cl <- make_med_cluster(NCORE, ROOT, c("CNT", "M", "C", "x"))
  ## ---- bootstrap
  set.seed(481)
  IDX <- lapply(seq_len(B_BOOT), function(b) sample.int(n, n, replace = TRUE))
  rb <- run_blocks(cl, IDX, function(i)
    tryCatch(med_replicate(i, x[i], CNT, M, C), error = function(e) NULL), "boot")
  ok <- !vapply(rb, is.null, logical(1))
  NIE <- do.call(rbind, lapply(rb[ok], `[[`, "NIE"))
  TEm <- do.call(rbind, lapply(rb[ok], `[[`, "TE"))
  rm(rb); invisible(gc())
  PMb <- sapply(SIG, function(s) apply(NIE[, s, drop = FALSE] / TEm[, s, drop = FALSE], 1, median))
  ## ---- permutation
  set.seed(4812)
  PRM <- lapply(seq_len(B_PERM), function(b) sample(x))
  rp <- run_blocks(cl, PRM, function(xp)
    tryCatch({ r <- med_replicate(seq_len(ncol(CNT)), xp, CNT, M, C)
               list(NIE = r$NIE, TE = r$TE) }, error = function(e) NULL), "perm")
  stopCluster(cl)
  okp <- !vapply(rp, is.null, logical(1))
  PMp <- t(sapply(rp[okp], function(r) sapply(SIG, function(s) median(r$NIE[s] / r$TE[s]))))
  rm(rp); invisible(gc())

  ## ---- summaries
  ci <- function(v) quantile(v, c(.025, .975), names = FALSE)
  res <- data.frame(tree = vn, set = names(SIG), n_sig = lengths(SIG), observed = obs,
                    boot_lo = apply(PMb, 2, function(v) ci(v)[1]),
                    boot_hi = apply(PMb, 2, function(v) ci(v)[2]),
                    boot_se = apply(PMb, 2, sd),
                    null_median = apply(PMp, 2, median),
                    null_lo = apply(PMp, 2, function(v) ci(v)[1]),
                    null_hi = apply(PMp, 2, function(v) ci(v)[2]),
                    perm_p = sapply(seq_along(SIG), function(j)
                      (1 + sum(abs(PMp[, j]) >= abs(obs[j]))) / (1 + nrow(PMp))),
                    stringsAsFactors = FALSE)
  dif <- data.frame(tree = vn,
    contrast = c("C2 markers - C1 non-markers", "C3 all - C1 non-markers"),
    observed = c(obs["C2"] - obs["C1"], obs["C3"] - obs["C1"]),
    lo = c(ci(PMb[, "C2"] - PMb[, "C1"])[1], ci(PMb[, "C3"] - PMb[, "C1"])[1]),
    hi = c(ci(PMb[, "C2"] - PMb[, "C1"])[2], ci(PMb[, "C3"] - PMb[, "C1"])[2]),
    stringsAsFactors = FALSE)
  print(transform(res[, c("set", "n_sig", "observed", "boot_lo", "boot_hi", "null_median", "perm_p")],
                  observed = round(observed, 4), boot_lo = round(boot_lo, 4),
                  boot_hi = round(boot_hi, 4), null_median = round(null_median, 4),
                  perm_p = round(perm_p, 3)), row.names = FALSE)
  print(transform(dif, observed = round(observed, 4), lo = round(lo, 4), hi = round(hi, 4)),
        row.names = FALSE)

  # per-gene bootstrap intervals, all 12,368 genes
  pg <- data.frame(gene = cf$tb$gene, symbol = cf$tb$symbol, set = cf$tb$set,
                   NIE = cf$tb$NIE, NIE_se_boot = apply(NIE, 2, sd),
                   NIE_lo = apply(NIE, 2, quantile, .025, names = FALSE),
                   NIE_hi = apply(NIE, 2, quantile, .975, names = FALSE),
                   TE = cf$tb$TE, TE_lo = apply(TEm, 2, quantile, .025, names = FALSE),
                   TE_hi = apply(TEm, 2, quantile, .975, names = FALSE),
                   stringsAsFactors = FALSE)
  pg$NIE_p_boot <- pmax(2 * pmin(colMeans(NIE <= 0), colMeans(NIE >= 0)), 1 / nrow(NIE))
  write.csv(pg, sprintf("results/circularity/%s_pergene_bootstrap%s.csv", vn, TAG), row.names = FALSE)
  write.csv(res, sprintf("results/circularity/%s_pm_intervals%s.csv", vn, TAG), row.names = FALSE)
  write.csv(dif, sprintf("results/circularity/%s_pm_differences%s.csv", vn, TAG), row.names = FALSE)
  saveRDS(list(res = res, dif = dif, PMb = PMb, PMp = PMp, nboot = nrow(NIE), nperm = nrow(PMp)), ckpt)
  rm(NIE, TEm); invisible(gc())
}
cat("\n10b complete.\n")
