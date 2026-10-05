#!/usr/bin/env Rscript
# ==============================================================================
# 03_wholeblood.R -- whole-blood differential expression, and the minimum
# detectable effect that makes a null interpretable.
#
# Covers hypotheses WB_LI (positive control) and WB_BPD, each fitted twice:
# unadjusted, and adjusted for the four ILR cell-composition balances.
#
# DECISIONS MADE HERE, AND WHY
#
# 1. voom + limma. It is the estimator the source paper used, so our numbers
#    are comparable to theirs, and it handles the mean-variance relationship of
#    counts properly. limma-trend and edgeR-QLF are run as pre-declared
#    sensitivity analyses rather than left as an open choice.
#
# 2. The gene set is FIXED at the 12,368 genes filtered once on all 474
#    samples. Percentages across contrasts are only comparable over a common
#    universe, and the hypotheses are stated as percentages.
#
# 3. pi0 gets TWO uncertainty estimates. The pre-specified one is a bootstrap
#    over genes, which is standard but assumes genes are independent -- they
#    are not, so its interval is too narrow. The permutation null (04) gives
#    the honest reference because each permutation preserves the real
#    gene-gene correlation structure. Both are reported; where they disagree,
#    the permutation one is the one to trust.
#
# 4. The minimum detectable effect is measured by spike-in, not computed from
#    a power formula. Formula-based power needs an assumed variance; the
#    spike-in uses the actual variance of this matrix at this sample size with
#    this design. It is computed on PERMUTED labels so no real signal
#    contaminates the calibration.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma); library(statmod)})
HERE <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/de_analysis"
setwd(HERE)
set.seed(481)
p <- readRDS("data/prep.rds")
RES <- "results"; dir.create(RES, showWarnings = FALSE, recursive = TRUE)

BAND_FEW  <- 0.5   # % of genes tested
BAND_MANY <- 2.0
NGENE <- nrow(p$counts)
cat(sprintf("== 03_wholeblood ==\ngenes tested (fixed): %d | bands: few <%.1f%% (%d genes), many >=%.1f%% (%d genes)\n\n",
            NGENE, BAND_FEW, round(BAND_FEW/100*NGENE), BAND_MANY, round(BAND_MANY/100*NGENE)))

## ---- helpers ---------------------------------------------------------------
pi0_storey <- function(pv, lambda = 0.5) min(1, mean(pv > lambda, na.rm = TRUE) / (1 - lambda))

pi0_boot_ci <- function(pv, B = 1000, lambda = 0.5) {
  n <- length(pv)
  v <- vapply(seq_len(B), function(i) pi0_storey(pv[sample.int(n, n, replace = TRUE)], lambda), numeric(1))
  quantile(v, c(0.025, 0.975), names = FALSE)
}

design_for <- function(contrast, adjust) {
  s <- p$sel[[contrast]]
  mm <- p$meta[s$keep, , drop = FALSE]
  mm$grp <- s$grp
  mm$plate <- droplevels(mm$plate); mm$sex <- droplevels(mm$sex)
  if ("group" %in% s$cov) mm$group <- droplevels(mm$group)
  cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]])) > 1, logical(1))]
  X <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse = " + "))), data = mm)
  if (adjust == "ilr") X <- cbind(X, p$ILR[mm$title, , drop = FALSE])
  ne <- nonEstimable(X); if (!is.null(ne)) X <- X[, setdiff(colnames(X), ne), drop = FALSE]
  list(X = X, samples = mm$title, meta = mm, dropped = setdiff(s$cov, cv))
}

fit_wb <- function(contrast, adjust, method = "voom", X_override = NULL) {
  d <- design_for(contrast, adjust)
  X <- if (is.null(X_override)) d$X else X_override
  cnt <- p$counts[, d$samples, drop = FALSE]
  dge <- calcNormFactors(DGEList(cnt), method = "TMM")
  if (method == "voom") {
    v <- voom(dge, X); f <- eBayes(lmFit(v, X))
    tt <- topTable(f, coef = "grpcase", number = Inf, sort.by = "none")
  } else if (method == "trend") {
    f <- eBayes(lmFit(cpm(dge, log = TRUE, prior.count = 3), X), trend = TRUE)
    tt <- topTable(f, coef = "grpcase", number = Inf, sort.by = "none")
  } else {
    dge <- estimateDisp(dge, X)
    r <- glmQLFTest(glmQLFit(dge, X), coef = "grpcase")
    tt <- topTags(r, n = Inf, sort.by = "none")$table
    names(tt)[names(tt) == "PValue"] <- "P.Value"; names(tt)[names(tt) == "FDR"] <- "adj.P.Val"
  }
  tt$gene <- rownames(tt)
  list(tt = tt, design = d, fit = if (method == "voom") f else NULL, dge = dge)
}

summarise <- function(tt, label) {
  n05 <- sum(tt$adj.P.Val < 0.05, na.rm = TRUE)
  pct <- 100 * n05 / nrow(tt)
  pi0 <- pi0_storey(tt$P.Value); ci <- pi0_boot_ci(tt$P.Value)
  data.frame(analysis = label, genes = nrow(tt),
             deg_fdr05 = n05, pct_of_genes = round(pct, 3),
             deg_fdr10 = sum(tt$adj.P.Val < 0.10, na.rm = TRUE),
             deg_bonf  = sum(p.adjust(tt$P.Value, "bonferroni") < 0.05, na.rm = TRUE),
             pi0 = round(pi0, 4), pi0_lo = round(ci[1], 4), pi0_hi = round(ci[2], 4),
             pi0_ci_includes_1 = ci[2] >= 1,
             max_abs_logFC = round(max(abs(tt$logFC)), 3),
             min_fdr = signif(min(tt$adj.P.Val, na.rm = TRUE), 3),
             stringsAsFactors = FALSE)
}

## ---- main fits -------------------------------------------------------------
cat("=== main fits ===\n")
FITS <- list(); S <- list()
for (cn in c("LI", "BPD")) for (ad in c("raw", "ilr")) {
  lab <- sprintf("WB_%s_%s", cn, ad)
  f <- fit_wb(cn, ad); FITS[[lab]] <- f
  S[[lab]] <- summarise(f$tt, lab)
  d <- f$design
  cat(sprintf("  %-12s n=%3d (%d vs %d) | dropped: %s\n", lab, nrow(d$meta),
              sum(d$meta$grp == "ctrl"), sum(d$meta$grp == "case"),
              if (length(d$dropped)) paste(d$dropped, collapse = ",") else "none"))
  write.csv(f$tt[order(f$tt$P.Value), ], file.path(RES, paste0(lab, "_toptable.csv")), row.names = FALSE)
}
SUM <- do.call(rbind, S)
cat("\n")
print(SUM[, c("analysis", "deg_fdr05", "pct_of_genes", "pi0", "pi0_lo", "pi0_hi",
              "pi0_ci_includes_1", "max_abs_logFC")], row.names = FALSE)

## ---- TREAT -----------------------------------------------------------------
cat("\n=== TREAT (formal test against a log-fold-change threshold) ===\n")
TR <- do.call(rbind, lapply(names(FITS), function(lab) {
  f <- FITS[[lab]]; if (is.null(f$fit)) return(NULL)
  do.call(rbind, lapply(c(0.1, 0.2), function(L) {
    tf <- treat(f$fit, lfc = L, trend = FALSE)
    tt <- topTreat(tf, coef = "grpcase", number = Inf, sort.by = "none")
    data.frame(analysis = lab, lfc = L, deg = sum(tt$adj.P.Val < 0.05, na.rm = TRUE),
               stringsAsFactors = FALSE) })) }))
print(TR, row.names = FALSE)

## ---- minimum detectable effect ---------------------------------------------
# Labels are PERMUTED first, so the matrix carries no real group difference and
# whatever is detected is attributable to the spike alone.
cat("\n=== minimum detectable effect (spike-in on permuted labels) ===\n")
EFF <- c(0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0)
NSPIKE <- 1200
mde_curve <- function(contrast) {
  d <- design_for(contrast, "raw")
  X <- d$X
  X[, "grpcase"] <- X[sample(nrow(X)), "grpcase"]          # break the real signal
  cnt <- p$counts[, d$samples, drop = FALSE]
  dge <- calcNormFactors(DGEList(cnt), method = "TMM")
  E <- cpm(dge, log = TRUE, prior.count = 3)
  spk <- sample(rownames(E), NSPIKE)
  grp <- X[, "grpcase"]
  do.call(rbind, lapply(EFF, function(L) {
    Y <- E; Y[spk, ] <- Y[spk, ] + matrix(rep(L * grp, each = NSPIKE), nrow = NSPIKE)
    tt <- topTable(eBayes(lmFit(Y, X), trend = TRUE), coef = "grpcase", number = Inf, sort.by = "none")
    det <- sum(tt[spk, "adj.P.Val"] < 0.05, na.rm = TRUE)
    data.frame(contrast = contrast, log2FC = L, detected = det,
               power = round(det / NSPIKE, 3), stringsAsFactors = FALSE) }))
}
MDE <- do.call(rbind, lapply(c("LI", "BPD"), mde_curve))
print(MDE, row.names = FALSE)
mde80 <- sapply(c("LI", "BPD"), function(cn) {
  m <- MDE[MDE$contrast == cn, ]
  if (all(m$power < 0.8)) return(NA_real_)
  approx(m$power, m$log2FC, xout = 0.8, ties = "ordered")$y })
cat(sprintf("\n  MDE at 80%% power:  LI = %.3f log2FC   BPD = %.3f log2FC\n", mde80["LI"], mde80["BPD"]))
cat(sprintf("  config meaningfulness floor = 0.5 log2FC -> null claims are %s\n",
            ifelse(all(mde80 <= 0.5, na.rm = TRUE), "INTERPRETABLE (MDE below the floor)",
                   "NOT interpretable at that floor")))

## ---- sensitivity: estimator -------------------------------------------------
cat("\n=== sensitivity: estimator ===\n")
SENS <- do.call(rbind, lapply(c("LI", "BPD"), function(cn)
  do.call(rbind, lapply(c("voom", "trend", "QLF"), function(mth) {
    tt <- fit_wb(cn, "raw", method = mth)$tt
    data.frame(contrast = cn, method = mth, deg = sum(tt$adj.P.Val < 0.05, na.rm = TRUE),
               pct = round(100 * sum(tt$adj.P.Val < 0.05, na.rm = TRUE) / nrow(tt), 3),
               pi0 = round(pi0_storey(tt$P.Value), 3), stringsAsFactors = FALSE) }))))
print(SENS, row.names = FALSE)

## ---- decision rules ---------------------------------------------------------
verdict <- function(pct, ci1, expect) {
  if (expect == "few") {
    if (pct < BAND_FEW && ci1) "SUPPORTED"
    else if (pct >= BAND_MANY) "REFUTED"
    else "INCONCLUSIVE"
  } else {
    if (pct >= BAND_MANY) "SUPPORTED"
    else if (pct < BAND_FEW) "REFUTED (pipeline failure)"
    else "INCONCLUSIVE"
  }
}
cat("\n================== PRE-SPECIFIED DECISIONS ==================\n")
DEC <- rbind(
  data.frame(hypothesis = "WB_LI (positive control)", expect = "many",
             pct = SUM["WB_LI_raw", "pct_of_genes"],
             pi0_ci_incl_1 = SUM["WB_LI_raw", "pi0_ci_includes_1"],
             verdict = verdict(SUM["WB_LI_raw", "pct_of_genes"], SUM["WB_LI_raw", "pi0_ci_includes_1"], "many")),
  data.frame(hypothesis = "WB_BPD", expect = "few",
             pct = SUM["WB_BPD_raw", "pct_of_genes"],
             pi0_ci_incl_1 = SUM["WB_BPD_raw", "pi0_ci_includes_1"],
             verdict = verdict(SUM["WB_BPD_raw", "pct_of_genes"], SUM["WB_BPD_raw", "pi0_ci_includes_1"], "few")))
print(DEC, row.names = FALSE)

write.csv(SUM, file.path(RES, "wholeblood_summary.csv"), row.names = FALSE)
write.csv(TR,  file.path(RES, "wholeblood_treat.csv"), row.names = FALSE)
write.csv(MDE, file.path(RES, "wholeblood_mde.csv"), row.names = FALSE)
write.csv(SENS, file.path(RES, "wholeblood_sensitivity_estimator.csv"), row.names = FALSE)
write.csv(DEC, file.path(RES, "wholeblood_decisions.csv"), row.names = FALSE)
saveRDS(list(SUM = SUM, MDE = MDE, mde80 = mde80, fits = lapply(FITS, `[[`, "tt")),
        file.path(RES, "wholeblood.rds"))
cat("\nwrote results/wholeblood_*.csv\n")
