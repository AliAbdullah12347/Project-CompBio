#!/usr/bin/env Rscript
# ==============================================================================
# 09_wholeblood_sensitivity.R -- the two pre-declared whole-blood sensitivity
# analyses that 03 did not cover.
#
# config.yaml declares four sensitivity analyses. 03 ran the estimator swap and
# 04 ran the permutation null. The remaining two are run here, so that running
# them is a commitment kept rather than a post-hoc search.
#
# ------------------------------------------------------------------------------
# A. TOBACCO MULTIPLE IMPUTATION, POOLED BY RUBIN'S RULES
#
# 20 of the 474 subjects had no recorded tobacco status and were multiply
# imputed, producing 20 completed columns. Analysing only column 1 would treat
# a guess as a measurement and understate uncertainty. Rubin's rules restore
# the missing-data uncertainty by adding the BETWEEN-imputation variance to the
# average WITHIN-imputation variance:
#
#     Qbar = mean_m(logFC_m)
#     Ubar = mean_m(SE_m^2)                 <- ordinary sampling variance
#     B    = var_m(logFC_m)                 <- extra variance from not knowing tobacco
#     T    = Ubar + (1 + 1/m) * B           <- total variance
#
# Degrees of freedom use the Barnard-Rubin small-sample correction, because the
# naive Rubin df can exceed the complete-data df, which is incoherent.
#
# WHAT WOULD FALSIFY THE MAIN RESULT: if the pooled DEG percentage moved across
# a band boundary, the headline verdict would depend on an imputation guess.
#
# ------------------------------------------------------------------------------
# B. DEPTH / RIN SCREENED SUBSET
#
# Krebs et al. state twice that no samples were removed for quality, so the
# main analysis keeps all 474 and records the flags only. But nine samples fall
# under 3M assigned reads and four have RIN < 5, and a reader is entitled to
# ask whether those drive the result. Re-fitting on the 461 that pass both
# screens answers that without changing the primary analysis.
#
# The GENE FILTER IS NOT RE-RUN. The hypotheses are percentages, and a
# percentage computed over a different gene universe is not comparable to the
# main result. This is a subject-level sensitivity check, not a new analysis.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})
HERE <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/de_analysis"
setwd(HERE); set.seed(481)
RES <- "results"
p <- readRDS("data/prep.rds")
NGENE <- nrow(p$counts)
BAND_FEW <- 0.5; BAND_MANY <- 2.0

pi0_storey <- function(pv, l = 0.5) min(1, mean(pv > l, na.rm = TRUE) / (1 - l))

# Build the design exactly as 03 does, but allow the subject set and the
# tobacco column to be swapped out.
design_for <- function(contrast, adjust, keep_extra = NULL, tob = NULL) {
  s <- p$sel[[contrast]]
  k <- s$keep; g <- s$grp
  if (!is.null(keep_extra)) { g <- g[keep_extra[k]]; k <- k & keep_extra }
  mm <- p$meta[k, , drop = FALSE]; mm$grp <- g
  if (!is.null(tob)) mm$tobacco <- tob[k]
  mm$plate <- droplevels(factor(mm$plate)); mm$sex <- droplevels(factor(mm$sex))
  if ("group" %in% s$cov) mm$group <- droplevels(factor(mm$group))
  cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]])) > 1, logical(1))]
  X <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse = " + "))), data = mm)
  if (adjust == "ilr") X <- cbind(X, p$ILR[mm$title, , drop = FALSE])
  ne <- nonEstimable(X); if (!is.null(ne)) X <- X[, setdiff(colnames(X), ne), drop = FALSE]
  list(X = X, samples = mm$title, meta = mm, dropped = setdiff(s$cov, cv))
}

fit_voom <- function(d) {
  dge <- calcNormFactors(DGEList(p$counts[, d$samples, drop = FALSE]), method = "TMM")
  eBayes(lmFit(voom(dge, d$X), d$X))
}

cat(sprintf("== 09_wholeblood_sensitivity ==\ngenes %d | bands few <%.1f%%  many >=%.1f%%\n",
            NGENE, BAND_FEW, BAND_MANY))

## ============================ A. Rubin pooling ===============================
M <- ncol(p$TOB)
cat(sprintf("\n=== A. tobacco multiple imputation: %d completed datasets, Rubin pooling ===\n", M))
cat(sprintf("    20 of 474 subjects imputed; smoker count ranges %d-%d across imputations\n",
            min(colSums(p$TOB == 1)), max(colSums(p$TOB == 1))))

rubin <- list(); PER_IMP <- list()
for (cn in c("LI", "BPD")) for (ad in c("raw", "ilr")) {
  lab <- sprintf("WB_%s_%s", cn, ad)
  Q <- NULL; U <- NULL; dfc <- NULL
  for (m in seq_len(M)) {
    d <- design_for(cn, ad, tob = p$TOB[, m])
    f <- fit_voom(d)
    Q <- cbind(Q, f$coefficients[, "grpcase"])
    U <- cbind(U, (f$stdev.unscaled[, "grpcase"] * sqrt(f$s2.post))^2)
    dfc <- f$df.total
    tt <- topTable(f, coef = "grpcase", number = Inf, sort.by = "none")
    PER_IMP[[length(PER_IMP) + 1]] <- data.frame(
      analysis = lab, imputation = m,
      deg = sum(tt$adj.P.Val < 0.05, na.rm = TRUE),
      pct = round(100 * sum(tt$adj.P.Val < 0.05, na.rm = TRUE) / nrow(tt), 3),
      pi0 = round(pi0_storey(tt$P.Value), 4), stringsAsFactors = FALSE)
  }
  Qbar <- rowMeans(Q); Ubar <- rowMeans(U)
  B    <- apply(Q, 1, var)
  Tv   <- Ubar + (1 + 1 / M) * B
  lam  <- pmin(1 - 1e-12, pmax(1e-12, (1 + 1 / M) * B / Tv))   # fraction of missing information
  df_old <- (M - 1) / lam^2
  df_obs <- (dfc + 1) / (dfc + 3) * dfc * (1 - lam)            # Barnard-Rubin
  df_adj <- df_old * df_obs / (df_old + df_obs)
  pv <- 2 * pt(-abs(Qbar / sqrt(Tv)), df = df_adj)
  qv <- p.adjust(pv, "BH")
  rubin[[lab]] <- data.frame(
    analysis = lab, m = M,
    deg_pooled = sum(qv < 0.05, na.rm = TRUE),
    pct_pooled = round(100 * sum(qv < 0.05, na.rm = TRUE) / length(qv), 3),
    pi0_pooled = round(pi0_storey(pv), 4),
    median_fmi = round(median(lam), 4),
    max_fmi = round(max(lam), 4),
    median_se_inflation = round(median(sqrt(Tv / Ubar)), 4),
    stringsAsFactors = FALSE)
  cat(sprintf("  %-12s pooled DEG %5d (%6.3f%%)  pi0 %.4f  median FMI %.4f  median SE inflation %.4fx\n",
              lab, rubin[[lab]]$deg_pooled, rubin[[lab]]$pct_pooled, rubin[[lab]]$pi0_pooled,
              rubin[[lab]]$median_fmi, rubin[[lab]]$median_se_inflation))
  flush.console()
}
RUB <- do.call(rbind, rubin); PI <- do.call(rbind, PER_IMP)

cat("\n  spread across the 20 single-imputation analyses (before pooling):\n")
for (lab in unique(PI$analysis)) {
  z <- PI[PI$analysis == lab, ]
  cat(sprintf("    %-12s DEG %d-%d | %%genes %.3f-%.3f | pi0 %.4f-%.4f\n",
              lab, min(z$deg), max(z$deg), min(z$pct), max(z$pct), min(z$pi0), max(z$pi0)))
}

## ============================ B. depth/RIN screen ============================
ok <- p$meta$depth_ok & p$meta$rin_ok
cat(sprintf("\n=== B. depth/RIN screened subset: %d of %d samples pass both ===\n", sum(ok), nrow(p$meta)))
SCR <- do.call(rbind, lapply(c("LI", "BPD"), function(cn) do.call(rbind, lapply(c("raw", "ilr"), function(ad) {
  lab <- sprintf("WB_%s_%s", cn, ad)
  d <- design_for(cn, ad, keep_extra = ok)
  tt <- topTable(fit_voom(d), coef = "grpcase", number = Inf, sort.by = "none")
  data.frame(analysis = lab, n_full = sum(p$sel[[cn]]$keep), n_screened = length(d$samples),
             deg = sum(tt$adj.P.Val < 0.05, na.rm = TRUE),
             pct = round(100 * sum(tt$adj.P.Val < 0.05, na.rm = TRUE) / nrow(tt), 3),
             pi0 = round(pi0_storey(tt$P.Value), 4), stringsAsFactors = FALSE) }))))
wb <- read.csv(file.path(RES, "wholeblood_summary.csv"), row.names = 1)
SCR$deg_main <- wb[SCR$analysis, "deg_fdr05"]
SCR$pct_main <- wb[SCR$analysis, "pct_of_genes"]
print(SCR, row.names = FALSE)

## ============================ verdict stability ==============================
band_of <- function(x) ifelse(x < BAND_FEW, "few", ifelse(x >= BAND_MANY, "many", "middle"))
STAB <- data.frame(
  analysis = RUB$analysis,
  pct_main = wb[RUB$analysis, "pct_of_genes"],
  pct_rubin = RUB$pct_pooled,
  pct_screened = SCR$pct[match(RUB$analysis, SCR$analysis)],
  stringsAsFactors = FALSE)
STAB$band_main <- band_of(STAB$pct_main)
STAB$band_rubin <- band_of(STAB$pct_rubin)
STAB$band_screened <- band_of(STAB$pct_screened)
STAB$all_agree <- STAB$band_main == STAB$band_rubin & STAB$band_main == STAB$band_screened
cat("\n=== does any verdict change? ===\n"); print(STAB, row.names = FALSE)
cat(sprintf("\n  every band identical across main / Rubin-pooled / screened: %s\n", all(STAB$all_agree)))

write.csv(RUB,  file.path(RES, "sens_tobacco_rubin.csv"), row.names = FALSE)
write.csv(PI,   file.path(RES, "sens_tobacco_per_imputation.csv"), row.names = FALSE)
write.csv(SCR,  file.path(RES, "sens_depthrin_screen.csv"), row.names = FALSE)
write.csv(STAB, file.path(RES, "sens_band_stability.csv"), row.names = FALSE)
cat("\nwrote results/sens_*.csv\n")
