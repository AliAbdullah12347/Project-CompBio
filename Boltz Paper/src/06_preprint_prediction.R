#!/usr/bin/env Rscript
# The prediction analysis from the bioRxiv PREPRINT of Boltz et al. (10.1101/2023.05.24.542156 v1,
# 2023-05-25). It was removed before the AJHG publication, and is the "single 70/30 split" result our
# project's Arm 1 revisits. Preprint, Methods:
#   "We also used logistic regression to predict either case control status or lithium use (only in
#   BP cases) from cell type expression estimates after residualizing for 50 expression PCs.
#   Variable numbers of genes were included based on genes with most variance per cell type, using a
#   range of 100 to 1000 genes with an interval of 100. Covariates include age, sex, RNA RIN, RNA
#   concentration, and cell type proportion estimates. A random 70% of individuals were sampled to
#   use for training, and 30% for testing the prediction."
# Result: "gene expression does not provide additional predictive value over the cell type
# proportions for either case/control status or lithium use." Their Supplementary Table 2 (the
# numbers) is an image behind bioRxiv's bot wall and could not be retrieved, so there is no numeric
# target; the comparison is qualitative.
#
# Reproduced as written, including two things our Arm 1 must not do:
#   * residualisation and top-variance gene selection use all samples before the split (leakage);
#   * a single random split. Added here (not in the preprint): the same model over 100 random splits,
#     to show how far one split's AUC can move.
# Adaptations: no RNA concentration in GEO; the preprint reports no metric, so test-set AUC is used;
# top-variance genes are ranked within a pool of the 1,000 genes most variable in bulk (after the same
# PC residualisation), because bMIND on all 18,006 genes is ~13 h here. bMIND runs on that pool for
# real, cached per chunk.

suppressPackageStartupMessages(library(MIND))
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
NCORE <- 4
CHUNK <- 200
POOL <- 1000
K_GRID <- seq(100, 1000, by = 100)
N_REPEAT <- 100
SEED <- 20230525      # preprint posting date

x <- readRDS(file.path(root, "data/processed/log2tpm_pcs.rds"))
meta <- read.csv(file.path(root, "data/processed/sample_metadata.csv"), check.names = FALSE,
                 colClasses = c(title = "character"))
cells <- readLines(file.path(root, "data/processed/celltypes_selected.txt"))
ck <- file.path(root, "data/processed/bmind_chunks")
W <- as.matrix(meta[, cells]); W <- W / rowSums(W); rownames(W) <- meta$title

resid_pcs <- function(M) t(resid(lm(t(M) ~ x$pcs)))   # genes x samples, 50 PCs regressed out

# ---- bMIND on the high-variance pool ----------------------------------------------------------
bulk_var <- apply(resid_pcs(x$logx), 1, var)
pool <- names(sort(bulk_var, decreasing = TRUE))[1:POOL]
chunks <- split(pool, ceiling(seq_along(pool) / CHUNK))
cat(sprintf("bMIND on the %d most variable genes (%d chunks)\n", POOL, length(chunks)))
bind_genes <- function(arrs) {        # stack gene x cell x sample arrays along the gene axis
  out <- array(NA_real_, c(sum(sapply(arrs, function(a) dim(a)[1])), dim(arrs[[1]])[2:3]),
               dimnames = c(list(unlist(lapply(arrs, function(a) dimnames(a)[[1]]))), dimnames(arrs[[1]])[2:3]))
  start <- 0
  for (a in arrs) {
    out[start + seq_len(dim(a)[1]), , ] <- a
    start <- start + dim(a)[1]
  }
  out
}
# Checkpoints are keyed by a hash of their inputs and code (as in 04), so stale chunks are never reused.
fingerprint <- function(...) { f <- tempfile(); saveRDS(list(...), f); on.exit(unlink(f)); unname(tools::md5sum(f)) }
key <- fingerprint(x$logx[pool, ], W, CHUNK, as.character(packageVersion("MIND")))
kf <- file.path(ck, "ctspool.key")
if (!file.exists(kf) || readLines(kf) != key) {
  unlink(list.files(ck, "^ctspool_[0-9]+\\.rds$", full.names = TRUE))
  writeLines(key, kf)
}
A <- bind_genes(lapply(seq_along(chunks), function(i) {
  f <- file.path(ck, sprintf("ctspool_%03d.rds", i))
  if (!file.exists(f)) {             # resume: finished chunks are never recomputed
    t0 <- Sys.time()
    saveRDS(bMIND(x$logx[chunks[[i]], , drop = FALSE], frac = W, ncore = NCORE)$A, f)
    cat(sprintf("  pool chunk %d/%d: %.1f min\n", i, length(chunks), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    flush(stdout())
  }
  readRDS(f)
}))
stopifnot(identical(dimnames(A)[[3]], meta$title))

auc <- function(score, y) {           # Mann-Whitney AUC, ties counted half
  r <- rank(score); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
fit_auc <- function(d, train) {
  m <- suppressWarnings(glm(y ~ ., data = d[train, ], family = binomial))
  pr <- suppressWarnings(predict(m, newdata = d[!train, ], type = "link"))
  c(auc = auc(pr, d$y[!train]), rank_deficient = any(is.na(coef(m))), converged = m$converged)
}

outcomes <- list(lithium = list(rows = meta$case == 1, y = meta$lithium),
                 casecontrol = list(rows = rep(TRUE, nrow(meta)), y = meta$case))
single <- list(); repeated <- list()
for (oc in names(outcomes)) {
  r <- outcomes[[oc]]$rows
  base <- data.frame(y = outcomes[[oc]]$y[r], age = meta$age[r], sex = factor(meta$sex[r]),
                     rin = meta$rin[r], meta[r, cells])
  set.seed(SEED)
  train <- seq_len(sum(r)) %in% sample(sum(r), round(0.7 * sum(r)))
  cat(sprintf("\n=== %s: n = %d (train %d / test %d) ===\n", oc, sum(r), sum(train), sum(!train)))
  b <- fit_auc(base, train)
  cat(sprintf("proportions + covariates only: test AUC %.3f\n", b["auc"]))
  single[[length(single) + 1]] <- data.frame(outcome = oc, cell_type = "none (proportions only)", k = 0, t(b))
  for (ct in cells) {
    # "after residualizing for 50 expression PCs" -- on all samples, before the split, as described
    E <- resid_pcs(A[, ct, ])[, r]
    ranked <- names(sort(apply(E, 1, var), decreasing = TRUE))
    for (k in K_GRID) {
      d <- cbind(base, t(E[ranked[1:k], , drop = FALSE]))
      single[[length(single) + 1]] <- data.frame(outcome = oc, cell_type = ct, k = k, t(fit_auc(d, train)))
    }
    s <- do.call(rbind, single)
    s <- s[s$outcome == oc & s$cell_type == ct, ]
    cat(sprintf("  + %-28s AUC over k = 100..1000: %s\n", ct, paste(sprintf("%.2f", s$auc), collapse = " ")))
  }

  # Added check: the baseline and the k = 100 model over many random 70/30 splits.
  set.seed(SEED)
  splits <- replicate(N_REPEAT, seq_len(sum(r)) %in% sample(sum(r), round(0.7 * sum(r))))
  best_ct <- cells[1]
  E <- resid_pcs(A[, best_ct, ])[, r]
  d100 <- cbind(base, t(E[names(sort(apply(E, 1, var), decreasing = TRUE))[1:100], ]))
  rp <- t(apply(splits, 2, function(tr) c(base = fit_auc(base, tr)["auc"], k100 = fit_auc(d100, tr)["auc"])))
  colnames(rp) <- c("proportions_only", paste0("plus_100_", best_ct))
  repeated[[oc]] <- data.frame(outcome = oc, split = seq_len(N_REPEAT), rp, check.names = FALSE)
  q <- apply(rp, 2, quantile, c(0.025, 0.5, 0.975))
  cat(sprintf("Over %d random splits -- proportions only: median %.3f [95%% of splits %.3f-%.3f]; + 100 %s genes: median %.3f [%.3f-%.3f]\n",
              N_REPEAT, q[2, 1], q[1, 1], q[3, 1], best_ct, q[2, 2], q[1, 2], q[3, 2]))
}
write.csv(do.call(rbind, single), file.path(root, "results/preprint_prediction_single_split.csv"), row.names = FALSE)
write.csv(do.call(rbind, repeated), file.path(root, "results/preprint_prediction_repeated_splits.csv"), row.names = FALSE)
cat("\nWrote results/preprint_prediction_{single_split,repeated_splits}.csv\n")
