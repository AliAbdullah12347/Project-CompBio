#!/usr/bin/env Rscript
# bMIND cell-type expression and cell-type differential expression, Boltz et al. (2024) p. 324, 326, 333.
#
# Paper: "We log2-transformed the matrix of bulk TPM values before inputting these values into bMIND
# because the largest expression value was greater than 50 TPM. Using the cell-type proportions
# derived from CIBERSORTx in conjunction with these log-transformed bulk expression measures, we
# used bMIND to derive cell-type expression estimates, with flag np = TRUE." And: "For the
# cell-type-specific differential expression analysis, we used the bmind_de() function ... we also
# used the log2-transformed expression measures as inputs along with the first 50 expression PCs
# as covariates."
#
# How the paper maps onto MIND 0.3.3 (github.com/randel/MIND, the only public release):
#   * bMIND() has no `np` argument. With no single-cell prior supplied it builds a data-driven one
#     (profile = gene mean; covariance = var / sum(mean fraction^2)), the no-reference mode used
#     here. `np` exists only in bmind_de(), where np = TRUE means a non-informative prior.
#
# Compute budget (4-core laptop; measured): bMIND ~2.6 s/gene, bmind_de ~1.6 s/gene/contrast before
# its adaptive re-sampling (up to 1e6 draws for strong genes). All 18,006 genes would take ~30 h.
#   1. bMIND (random-effects MCMC, no closed form): run for real on a fixed random 2,000-gene subset.
#      Every downstream use (Fig 1B, Fig S2, Table S3) is a correlation or PCA across genes, which a
#      2,000-gene sample estimates closely.
#   2. bmind_de with np = TRUE and noRE = TRUE fits, per gene, the fixed-effects Bayesian regression
#          x ~ -1 + PCs + sum_k W_k:co + sum_k W_k:ca
#      under MCMCglmm's default priors: flat on coefficients (V = 1e10 I) and IW(V = 1, nu = 0), i.e.
#      p(sigma^2) proportional to 1/sigma^2, on the residual. Under those priors the posterior of each
#      contrast d_k = beta_k,ca - beta_k,co is exactly Student-t(n - p) centred on the OLS estimate,
#      so bmind_de's p = 2 min(P(d > 0), P(d < 0)) is the classical two-sided t-test p. MCMC only
#      approximates that limit. It is computed exactly for all genes here, and the real bmind_de()
#      is run on a 400-gene subset to confirm the two agree (results/bmind_de_mcmc_vs_exact.csv).

suppressPackageStartupMessages(library(MIND))
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
NCORE <- 4            # physical cores; hyperthreads gave no speed-up in benchmarking and cost RAM
CHUNK <- 200
N_CTS <- 2000         # genes for real bMIND
N_CHECK <- 400        # genes for the real-bmind_de-vs-exact check (first of the bMIND subset)
SEED <- 20240201      # Boltz publication date; only selects the gene subsets

ck <- file.path(root, "data/processed/bmind_chunks")
dir.create(ck, showWarnings = FALSE, recursive = TRUE)
x <- readRDS(file.path(root, "data/processed/log2tpm_pcs.rds"))
meta <- read.csv(file.path(root, "data/processed/sample_metadata.csv"), check.names = FALSE,
                 colClasses = c(title = "character"))
stopifnot(identical(colnames(x$logx), meta$title))
cells <- readLines(file.path(root, "data/processed/celltypes_selected.txt"))

# The eight selected types carry nearly all LM22 mass; renormalise so each row sums to 1, as bMIND's
# mixture model (bulk = sum_k w_k x_k) assumes.
W <- as.matrix(meta[, cells])
cat(sprintf("Selected-type share of LM22 mass: median %.3f, min %.3f\n", median(rowSums(W)), min(rowSums(W))))
W <- W / rowSums(W)
rownames(W) <- meta$title
set.seed(SEED)
sub_genes <- sample(rownames(x$logx), N_CTS)
check_genes <- sub_genes[1:N_CHECK]

# Checkpoints are keyed by a hash of everything that determines them -- the input rows, the function
# body and the MIND version -- so a changed input or edit can never silently reuse stale chunks.
fingerprint <- function(...) { f <- tempfile(); saveRDS(list(...), f); on.exit(unlink(f)); unname(tools::md5sum(f)) }
run_chunks <- function(tag, genes, fn, data) {
  chunks <- split(genes, ceiling(seq_along(genes) / CHUNK))
  key <- fingerprint(data, deparse(fn), CHUNK, as.character(packageVersion("MIND")))
  kf <- file.path(ck, paste0(tag, ".key"))
  if (!file.exists(kf) || readLines(kf) != key) {
    unlink(list.files(ck, paste0("^", tag, "_[0-9]+\\.rds$"), full.names = TRUE))
    writeLines(key, kf)
  }
  for (i in seq_along(chunks)) {
    f <- file.path(ck, sprintf("%s_%03d.rds", tag, i))
    if (file.exists(f)) next  # resume: finished chunks are never recomputed
    t0 <- Sys.time()
    saveRDS(fn(chunks[[i]]), f)
    cat(sprintf("  %s chunk %d/%d: %.1f min  [%s]\n", tag, i, length(chunks),
                as.numeric(difftime(Sys.time(), t0, units = "mins")), format(Sys.time(), "%H:%M")))
    flush(stdout())
  }
  lapply(seq_along(chunks), function(i) readRDS(file.path(ck, sprintf("%s_%03d.rds", tag, i))))
}

# Stack gene x cell x sample arrays along the gene axis (base R has no 3-D rbind).
bind_genes <- function(arrs) {
  out <- array(NA_real_, c(sum(sapply(arrs, function(a) dim(a)[1])), dim(arrs[[1]])[2:3]),
               dimnames = c(list(unlist(lapply(arrs, function(a) dimnames(a)[[1]]))), dimnames(arrs[[1]])[2:3]))
  start <- 0
  for (a in arrs) {
    out[start + seq_len(dim(a)[1]), , ] <- a
    start <- start + dim(a)[1]
  }
  out
}

# ---- 1. bmind_de, exact posterior (all genes) -----------------------------------------------------
bmind_de_exact <- function(X, Wr, y, pcs) {
  co <- as.numeric(y == 0); ca <- as.numeric(y == 1)
  D <- cbind(pcs, Wr * co, Wr * ca)            # column order mirrors bmind1_y's formula
  K <- ncol(Wr); P <- ncol(pcs)
  stopifnot(qr(D)$rank == ncol(D))
  XtXi <- solve(crossprod(D))
  B <- XtXi %*% crossprod(D, t(X))             # coefficients, (P + 2K) x genes, all genes at once
  res <- t(X) - D %*% B
  df <- nrow(D) - ncol(D)
  s2 <- colSums(res^2) / df
  pv <- sapply(seq_len(K), function(k) {
    cvec <- numeric(ncol(D)); cvec[P + K + k] <- 1; cvec[P + k] <- -1
    est <- drop(crossprod(cvec, B))
    se <- sqrt(s2 * drop(crossprod(cvec, XtXi %*% cvec)))
    2 * pt(-abs(est / se), df)
  })
  dimnames(pv) <- list(rownames(X), colnames(Wr))
  pv
}

contrasts <- list(casecontrol = list(rows = rep(TRUE, nrow(meta)), y = meta$case),
                  lithium = list(rows = meta$case == 1, y = meta$lithium))
for (tag in names(contrasts)) {
  r <- contrasts[[tag]]$rows
  y <- contrasts[[tag]]$y[r]
  cat(sprintf("\n=== bmind_de %s (n1 = %d, n0 = %d) ===\n", tag, sum(y == 1), sum(y == 0)))
  p <- bmind_de_exact(x$logx[, r], W[r, ], y, x$pcs[r, ])
  q <- apply(p, 2, p.adjust, method = "fdr")   # same BH-per-cell-type step as bmind_all()
  cat("DEGs at q < 0.05 per cell type:\n")
  print(colSums(q < 0.05))
  write.csv(data.frame(gene = rownames(p), p, check.names = FALSE),
            file.path(root, sprintf("results/bmind_de_%s_pvals.csv", tag)), row.names = FALSE)
  write.csv(data.frame(gene = rownames(q), q, check.names = FALSE),
            file.path(root, sprintf("results/bmind_de_%s_qvals.csv", tag)), row.names = FALSE)

  # The real bmind_de() on the check subset. Its p-values carry Monte Carlo error of roughly
  # sqrt(p / nsamp) and are floored at 1 / nsamp. Its adaptive re-sampling is capped at 1e4 draws
  # here (package default 1e6) because one gene at p ~ 1e-8 would otherwise cost about an hour;
  # the exact p is floored at the same 1e-4 for the comparison, so the cap cannot flatter it.
  mc <- do.call(rbind, run_chunks(paste0("demc_", tag), check_genes, function(g)
    bmind_de(x$logx[g, r, drop = FALSE], frac = W[r, ], y = y, covariate = x$pcs[r, ],
             covariate_bulk = colnames(x$pcs), np = TRUE, max_samp = 1e4, ncore = NCORE)$pval,
    data = list(x$logx[check_genes, r], W[r, ], y, x$pcs[r, ])))
  ex <- pmax(p[rownames(mc), ], 1e-4)
  colnames(mc) <- colnames(ex)  # bmind_all() strips the dots from cell-type names
  cmp <- data.frame(contrast = tag, cell_type = colnames(ex),
    n_genes = nrow(mc),
    cor_log10p = sapply(colnames(ex), function(k) cor(log10(mc[, k]), log10(ex[, k]))),
    median_abs_diff_p = sapply(colnames(ex), function(k) median(abs(mc[, k] - ex[, k]))),
    agree_p05 = sapply(colnames(ex), function(k) mean((mc[, k] < 0.05) == (ex[, k] < 0.05))))
  cat("Real bmind_de() vs exact posterior on the check subset:\n")
  print(format(cmp[, -1], digits = 3), row.names = FALSE)
  f <- file.path(root, "results/bmind_de_mcmc_vs_exact.csv")
  write.csv(rbind(if (tag != "casecontrol" && file.exists(f)) read.csv(f), cmp), f, row.names = FALSE)
}
# ---- 2. bMIND cell-type expression (real MCMC, 2,000 genes) --------------------------------------
cat(sprintf("\n=== bMIND cell-type expression: %d genes, %d cores ===\n", N_CTS, NCORE))
A <- bind_genes(run_chunks("cts", sub_genes, function(g) bMIND(x$logx[g, , drop = FALSE], frac = W, ncore = NCORE)$A,
                          data = list(x$logx[sub_genes, ], W)))
cat(sprintf("Estimated %d of %d genes (MCMCglmm errors drop the rest)\n", dim(A)[1], N_CTS))
saveRDS(list(A = A, frac = W), file.path(root, "data/processed/bmind_cts.rds"))

cat("\nWrote data/processed/bmind_cts.rds, results/bmind_de_*_{p,q}vals.csv, bmind_de_mcmc_vs_exact.csv\n")
