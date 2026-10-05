#!/usr/bin/env Rscript
# ==============================================================================
# 07_celltype_limma.R -- cell-type differential expression, done in a design we
# control.
#
# Takes the per-sample, per-lineage expression estimates from 06 and tests them
# with limma. This replaces bmind_de(), which cannot answer the hypotheses
# because its p-values floor at 2e-4 and it discards the covariate argument.
#
# DECISIONS MADE HERE, AND WHY
#
# 1. PRECISION WEIGHTS FROM THE POSTERIOR SE. bMIND returns a posterior
#    standard error for every gene x lineage x sample. Treating its posterior
#    MEANS as if they were measured data would understate uncertainty and
#    inflate the DEG count -- which would bias CT_BPD toward "supported"
#    exactly where we want no thumb on the scale. Weighting each observation
#    by 1/SE^2 carries the deconvolution's own uncertainty into the test.
#
#    Non-finite weights are repaired rather than dropped, and the repair is
#    counted and reported: SE == 0 would give an infinite weight (capped at the
#    largest finite weight seen), and a missing or infinite SE would give a
#    missing weight (set to the smallest finite weight, i.e. treated as the
#    least informative observation). Silently passing NA weights to lmFit would
#    either error or quietly delete observations.
#
# 2. SAME DESIGN AS WHOLE BLOOD. Identical covariates, identical contrast
#    coding, identical FDR rule. The two levels then differ only in what is
#    being measured, which is the comparison CT_BPD actually makes.
#
# 3. BOTH ADJUSTMENTS. Unadjusted, and adjusted for the four ILR balances.
#    For cell-type estimates the adjusted version is the circularity control:
#    bMIND's posterior for a lineage depends on that sample's fraction for that
#    lineage, and both contrasts move those fractions.
#
# 4. The MEAN lineage fraction is reported beside every result, because power
#    to detect a within-lineage effect scales with how much of the mixture that
#    lineage occupies. A null in B cells at 1.5% of the mixture means something
#    very different from a null in granulocytes at 39%.
# ==============================================================================

suppressPackageStartupMessages({library(limma)})
HERE <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/de_analysis"
setwd(HERE); set.seed(481)
RES <- "results"
p <- readRDS("data/prep.rds")
bp <- readRDS(file.path(RES, "bmind_profiles.rds"))
A <- bp$A; SE <- bp$SE; LINS <- bp$lineages
NGENE <- dim(A)[1]
BAND_FEW <- 0.5; BAND_MANY <- 2.0
cat(sprintf("== 07_celltype_limma ==\nA: %s | lineages: %s\n",
            paste(dim(A), collapse = " x "), paste(LINS, collapse = ", ")))
cat(sprintf("A: %d non-finite of %d | SE: %d non-finite, %d exactly zero | SE range %.4g .. %.4g\n\n",
            sum(!is.finite(A)), length(A), sum(!is.finite(SE)), sum(SE == 0, na.rm = TRUE),
            min(SE[is.finite(SE)]), max(SE[is.finite(SE)])))

pi0_storey <- function(pv, l = 0.5) min(1, mean(pv > l, na.rm = TRUE) / (1 - l))
pi0_ci <- function(pv, B = 1000, l = 0.5) {
  n <- length(pv)
  quantile(vapply(seq_len(B), function(i) pi0_storey(pv[sample.int(n, n, TRUE)], l),
                  numeric(1)), c(.025, .975), names = FALSE) }

# 1/SE^2, with non-finite values repaired and counted (see note 1 above).
make_weights <- function(S) {
  W <- 1 / (S^2)
  fin <- is.finite(W) & W > 0
  n_hi <- sum(!is.finite(W) & !is.na(S) & S == 0)
  n_lo <- sum(!fin) - n_hi
  if (any(!fin)) {
    hi <- max(W[fin]); lo <- min(W[fin])
    W[!is.finite(W) & !is.na(S) & S == 0] <- hi
    W[!is.finite(W) | W <= 0] <- lo
  }
  list(W = W / mean(W), n_hi = n_hi, n_lo = n_lo)
}

design_for <- function(contrast, adjust) {
  s <- p$sel[[contrast]]
  mm <- p$meta[s$keep, , drop = FALSE]; mm$grp <- s$grp
  mm$plate <- droplevels(factor(mm$plate)); mm$sex <- droplevels(factor(mm$sex))
  if ("group" %in% s$cov) mm$group <- droplevels(factor(mm$group))
  cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]])) > 1, logical(1))]
  X <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse = " + "))), data = mm)
  if (adjust == "ilr") X <- cbind(X, p$ILR[mm$title, , drop = FALSE])
  ne <- nonEstimable(X); if (!is.null(ne)) X <- X[, setdiff(colnames(X), ne), drop = FALSE]
  list(X = X, samples = mm$title, meta = mm) }

PV <- list(); FC <- list(); ROWS <- list(); REPAIR <- list()
for (cn in c("LI", "BPD")) for (ad in c("raw", "ilr")) {
  d <- design_for(cn, ad)
  key <- sprintf("%s_%s", cn, ad)
  pm <- matrix(NA_real_, NGENE, length(LINS), dimnames = list(dimnames(A)[[1]], LINS))
  fm <- pm
  for (ct in LINS) {
    Y  <- A[, ct, d$samples]            # genes x samples
    wr <- make_weights(SE[, ct, d$samples])
    REPAIR[[length(REPAIR) + 1]] <- data.frame(run = key, lineage = ct,
      n_capped_high = wr$n_hi, n_set_low = wr$n_lo, n_cells = length(wr$W),
      stringsAsFactors = FALSE)
    fit <- eBayes(lmFit(Y, d$X, weights = wr$W), trend = TRUE)
    tt <- topTable(fit, coef = "grpcase", number = Inf, sort.by = "none")
    tt$gene <- rownames(tt)                      # row.names=FALSE below would lose these
    pm[rownames(tt), ct] <- tt$P.Value
    fm[rownames(tt), ct] <- tt$logFC
    q <- p.adjust(tt$P.Value, "BH"); ci <- pi0_ci(tt$P.Value)
    ROWS[[length(ROWS) + 1]] <- data.frame(
      contrast = cn, adjust = ad, lineage = ct,
      mean_frac = round(mean(p$lineage[, ct]), 4),
      n_genes = nrow(tt), deg = sum(q < 0.05, na.rm = TRUE),
      pct = round(100 * sum(q < 0.05, na.rm = TRUE) / nrow(tt), 3),
      min_p = signif(min(tt$P.Value, na.rm = TRUE), 3),
      pi0 = round(pi0_storey(tt$P.Value), 4),
      pi0_lo = round(ci[1], 4), pi0_hi = round(ci[2], 4), pi0_ci_incl_1 = ci[2] >= 1,
      max_abs_logFC = round(max(abs(tt$logFC)), 3), stringsAsFactors = FALSE)
    write.csv(tt[order(tt$P.Value), c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")],
              file.path(RES, sprintf("CT_%s_%s_toptable.csv", key, ct)), row.names = FALSE)
  }
  PV[[key]] <- pm; FC[[key]] <- fm
  cat(sprintf("%s: n=%d, min p across lineages = %.3g\n", key, length(d$samples), min(pm, na.rm = TRUE)))
}
PER <- do.call(rbind, ROWS); REP <- do.call(rbind, REPAIR)

cat("\n=== weight repair (how many posterior SEs needed fixing) ===\n")
print(aggregate(cbind(n_capped_high, n_set_low) ~ run, REP, sum), row.names = FALSE)

cat("\n=== per lineage (BH within lineage) ===\n")
print(PER[, c("contrast", "adjust", "lineage", "mean_frac", "deg", "pct", "min_p",
              "pi0", "pi0_lo", "pi0_hi")], row.names = FALSE)

## ---- pooled across lineages -------------------------------------------------
POOL <- do.call(rbind, lapply(names(PV), function(k) {
  P <- PV[[k]]; q <- matrix(p.adjust(as.vector(P), "BH"), nrow(P), dimnames = dimnames(P))
  hit <- rownames(P)[rowSums(q < 0.05, na.rm = TRUE) > 0]
  data.frame(run = k, n_tests = length(P), n_na = sum(is.na(P)),
             deg_tests = sum(q < 0.05, na.rm = TRUE),
             deg_unique_genes = length(hit),
             pct_of_genes = round(100 * length(hit) / nrow(P), 3), stringsAsFactors = FALSE) }))
cat("\n=== pooled (BH over all lineage x gene tests) ===\n"); print(POOL, row.names = FALSE)

## ---- decisions ---------------------------------------------------------------
wb <- read.csv(file.path(RES, "wholeblood_summary.csv"), row.names = 1)
cat("\n================== PRE-SPECIFIED DECISIONS ==================\n")
sLI <- PER[PER$contrast == "LI" & PER$adjust == "raw", ]
pLI <- POOL[POOL$run == "LI_raw", ]
vLI <- {
  if (pLI$pct_of_genes >= BAND_MANY) "REFUTED"
  else if (all(sLI$pct < BAND_FEW) && all(sLI$pi0_ci_incl_1)) "SUPPORTED"
  else "INCONCLUSIVE" }
cat(sprintf("\nCT_LI  (expect FEW)\n  per lineage: %s\n",
            paste(sprintf("%s=%.3f%%", sLI$lineage, sLI$pct), collapse = "  ")))
cat(sprintf("  all below %.1f%%: %s | all pi0 CIs include 1: %s\n",
            BAND_FEW, all(sLI$pct < BAND_FEW), all(sLI$pi0_ci_incl_1)))
cat(sprintf("  pooled %d genes (%.3f%%)\n  VERDICT: %s\n",
            pLI$deg_unique_genes, pLI$pct_of_genes, vLI))

pBP <- POOL[POOL$run == "BPD_raw", ]; wbp <- wb["WB_BPD_raw", "pct_of_genes"]
vBP <- {
  if (pBP$pct_of_genes >= BAND_FEW && wbp < BAND_FEW) "SUPPORTED"
  else if (pBP$pct_of_genes <= wbp) "REFUTED"
  else "INCONCLUSIVE" }
cat(sprintf("\nCT_BPD (expect MORE than whole blood)\n"))
cat(sprintf("  whole blood %.3f%% (%d genes) | cell type %.3f%% (%d genes)\n",
            wbp, wb["WB_BPD_raw", "deg_fdr05"], pBP$pct_of_genes, pBP$deg_unique_genes))
cat(sprintf("  VERDICT: %s\n", vBP))

DEC <- data.frame(hypothesis = c("CT_LI", "CT_BPD"), expect = c("few", "more than WB"),
                  pct = c(pLI$pct_of_genes, pBP$pct_of_genes),
                  verdict = c(vLI, vBP), stringsAsFactors = FALSE)
write.csv(PER,  file.path(RES, "celltype_limma_per_lineage.csv"), row.names = FALSE)
write.csv(POOL, file.path(RES, "celltype_limma_pooled.csv"), row.names = FALSE)
write.csv(DEC,  file.path(RES, "celltype_limma_decisions.csv"), row.names = FALSE)
write.csv(REP,  file.path(RES, "celltype_weight_repair.csv"), row.names = FALSE)
saveRDS(list(PV = PV, FC = FC), file.path(RES, "celltype_limma_pvalues.rds"))
cat("\nwrote results/celltype_limma_*.csv\n")
