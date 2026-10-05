#!/usr/bin/env Rscript
# ==============================================================================
# 03_permutation_fdr.R -- the one correction here that assumes nothing about
# how genes depend on each other.
#
# WHY IT IS WORTH THREE HOURS
#
# Every method in 01 rests on an assumption about dependence:
#   BH      valid under independence or positive regression dependence (PRDS)
#   Storey  same, plus an estimated pi0
#   BY      valid under arbitrary dependence, but pays sum(1/i) ~ 9.9x at
#           m = 12,368, which is usually too much to accept
#   Holm    valid under anything, but controls FWER, a different and much
#           stricter goal
#
# Gene expression is co-expressed by construction -- that is the premise of the
# whole field, and Arm 3 of this project is about it. PRDS is PLAUSIBLE for
# expression data but it is an assumption, not a fact, and nothing in the data
# checks it.
#
# A permutation null needs no such assumption. Each permutation scrambles the
# group labels and leaves the expression matrix, and therefore every gene-gene
# correlation, exactly intact. The expected number of false positives at a
# threshold is then MEASURED rather than derived:
#
#     FDR(t) = pi0 * E_perm[ #{null p <= t} ] / max(1, #{observed p <= t})
#
# This is the Storey-Tibshirani plug-in, the same estimator SAM uses.
#
# WHAT IS PERMUTED, AND WHAT IS NOT
#
# Only the group column of the design matrix. All other covariates stay attached
# to their own subjects, so the permutation breaks the association between group
# and expression while preserving everything else -- including the covariate
# structure and the full correlation matrix of the outcome.
#
# For whole blood the entire voom + limma fit is repeated per permutation,
# because voom's precision weights depend on the design matrix and reusing them
# would leak the real grouping into the null.
#
# For cell type, bMIND is NOT repeated, and that is correct rather than a
# shortcut: 06 ran bMIND WITHOUT the phenotype, so the cell-type estimates do
# not depend on the labels in any way. A permutation would reproduce them
# bit-for-bit. The precision weights likewise depend only on the posterior SEs,
# so they are computed once. (Had the deconvolution used the phenotype -- as
# bmind_de does -- each permutation would have cost a 3-hour deconvolution and
# this analysis would be impossible.)
#
# B = 100. The quantity needed is a MEAN over permutations, not an extreme
# quantile, and a mean converges quickly; 100 gives a standard error on
# E[count] of about a tenth of its own standard deviation. v1 used 200 because
# it needed tail quantiles for a different question.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
HERE <- file.path(ROOT, "de_v2")
setwd(HERE); set.seed(481)
source("scripts/mtc.R")
RES <- "results"
PARTS <- file.path(RES, "perm_parts"); dir.create(PARTS, showWarnings = FALSE, recursive = TRUE)

p  <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
bp <- readRDS(file.path(ROOT, "de_analysis/results/bmind_profiles.rds"))
A <- bp$A; SE <- bp$SE; LINS <- bp$lineages
ALL <- readRDS(file.path(RES, "adjusted.rds"))

B <- 100
GRID <- c(1e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 0.01, 0.025, 0.05, 0.1)
cat(sprintf("== 03_permutation_fdr ==\nB=%d permutations | %d thresholds\n\n", B, length(GRID)))

design_for <- function(contrast, adjust = "raw") {
  s <- p$sel[[contrast]]
  mm <- p$meta[s$keep, , drop = FALSE]; mm$grp <- s$grp
  mm$plate <- droplevels(factor(mm$plate)); mm$sex <- droplevels(factor(mm$sex))
  if ("group" %in% s$cov) mm$group <- droplevels(factor(mm$group))
  cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]])) > 1, logical(1))]
  X <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse = " + "))), data = mm)
  if (adjust == "ilr") X <- cbind(X, p$ILR[mm$title, , drop = FALSE])
  ne <- nonEstimable(X); if (!is.null(ne)) X <- X[, setdiff(colnames(X), ne), drop = FALSE]
  list(X = X, samples = mm$title) }

make_weights <- function(S) {
  W <- 1 / (S^2); fin <- is.finite(W) & W > 0
  if (any(!fin)) W[!fin] <- min(W[fin]); W / mean(W) }

counts_at <- function(pv) vapply(GRID, function(t) sum(pv <= t, na.rm = TRUE), numeric(1))

run_perm <- function(key, contrast, level) {
  f <- file.path(PARTS, sprintf("%s.rds", key))
  if (file.exists(f)) { cat(sprintf("  %-12s from checkpoint\n", key)); return(readRDS(f)) }
  d <- design_for(contrast, "raw")
  t0 <- Sys.time()
  if (level == "WB") {
    cnt <- p$counts[, d$samples, drop = FALSE]
    dge <- calcNormFactors(DGEList(cnt), method = "TMM")
    NC <- t(vapply(seq_len(B), function(b) {
      Xp <- d$X; Xp[, "grpcase"] <- Xp[sample(nrow(Xp)), "grpcase"]
      v <- voom(dge, Xp)                       # re-run: weights depend on design
      tt <- topTable(eBayes(lmFit(v, Xp)), coef = "grpcase", number = Inf, sort.by = "none")
      if (b %% 25 == 0) { el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
        cat(sprintf("    %s %3d/%d  %.1f min, ETA %.1f\n", key, b, B, el, el/b*(B-b))); flush.console() }
      counts_at(tt$P.Value) }, numeric(length(GRID))))
  } else {
    Wl <- lapply(LINS, function(ct) make_weights(SE[, ct, d$samples])); names(Wl) <- LINS
    NC <- t(vapply(seq_len(B), function(b) {
      Xp <- d$X; Xp[, "grpcase"] <- Xp[sample(nrow(Xp)), "grpcase"]
      pv <- unlist(lapply(LINS, function(ct) {
        fit <- eBayes(lmFit(A[, ct, d$samples], Xp, weights = Wl[[ct]]), trend = TRUE)
        topTable(fit, coef = "grpcase", number = Inf, sort.by = "none")$P.Value }), use.names = FALSE)
      if (b %% 25 == 0) { el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
        cat(sprintf("    %s %3d/%d  %.1f min, ETA %.1f\n", key, b, B, el, el/b*(B-b))); flush.console() }
      counts_at(pv) }, numeric(length(GRID))))
  }
  colnames(NC) <- paste0("t_", GRID)
  saveRDS(NC, f)
  cat(sprintf("  %-12s done [%.1f min]\n", key, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  NC
}

JOBS <- data.frame(key = c("WB_LI_raw", "WB_BPD_raw", "CT_LI_raw", "CT_BPD_raw"),
                   contrast = c("LI", "BPD", "LI", "BPD"),
                   level = c("WB", "WB", "CT", "CT"), stringsAsFactors = FALSE)

OUT <- list(); QV <- list()
for (i in seq_len(nrow(JOBS))) {
  j <- JOBS[i, ]
  NC <- run_perm(j$key, j$contrast, j$level)
  obs <- ALL[[j$key]]$P.Value
  pi0 <- pi0_spline(obs)
  pf <- perm_fdr(obs, GRID, NC, pi0 = pi0)
  pf$analysis <- j$key; pf$pi0_used <- round(pi0, 4)
  # BH's implied FDR at the same thresholds, for a like-for-like comparison:
  # BH assumes E[#null <= t] = pi0 * m * t; the permutation measures it.
  # CARE WITH pi0 HERE. n_null_expected is the mean permutation count over ALL
  # m tests -- what would happen if every gene were null. The like-for-like
  # uniform comparison is therefore m * t, NOT pi0 * m * t. A first version
  # included pi0 on this line and so reported null inflation of 2-4x when the
  # true figure is 0.8-3.2x. pi0 belongs in the FDR formula above, where the
  # question is how many of the REAL nulls fall below t, and not in this
  # diagnostic.
  pf$n_null_if_uniform <- round(length(obs) * pf$threshold, 2)
  pf$null_inflation <- round(pf$n_null_expected / pmax(1e-9, pf$n_null_if_uniform), 3)
  OUT[[j$key]] <- pf
  QV[[j$key]] <- perm_qvalue(obs, pf)
  cat(sprintf("  %-12s pi0 %.3f | at p<=0.001: obs %d, null measured %.1f, null if uniform %.1f (x%.2f)\n",
              j$key, pi0, pf$n_observed[pf$threshold == 1e-3],
              pf$n_null_expected[pf$threshold == 1e-3],
              pf$n_null_if_uniform[pf$threshold == 1e-3],
              pf$null_inflation[pf$threshold == 1e-3]))
}
PF <- do.call(rbind, OUT)
write.csv(PF, file.path(RES, "permutation_fdr.csv"), row.names = FALSE)
saveRDS(QV, file.path(RES, "permutation_qvalues.rds"))

cat("\n=== how many genes does each correction keep, at 0.05? ===\n")
CMP <- do.call(rbind, lapply(JOBS$key, function(k) {
  f <- ALL[[k]]
  data.frame(analysis = k, n_tests = nrow(f),
             BH = sum(f$BH < 0.05), BY = sum(f$BY < 0.05),
             storey = sum(f$q_storey_spline < 0.05), holm = sum(f$holm < 0.05),
             permutation = sum(QV[[k]] < 0.05), stringsAsFactors = FALSE) }))
print(CMP, row.names = FALSE)
write.csv(CMP, file.path(RES, "correction_comparison.csv"), row.names = FALSE)
cat("\nwrote results/permutation_fdr.csv, correction_comparison.csv\n")
