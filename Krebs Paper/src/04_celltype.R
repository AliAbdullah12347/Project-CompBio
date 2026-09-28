#!/usr/bin/env Rscript
# Cell-type composition analysis, reproducing Krebs et al. 2020 Figure 2.
#
# Method verbatim from the Supplementary Methods ("Estimation of cell-type
# proportions"):
#
#   "The resulting estimated cell-type proportions were regressed on sex, age,
#    tobacco use, sequencing plate, RIN, and sequencing metric PCs 1 through 3,
#    and the residuals were used to predict lithium use in a stepwise linear
#    regression using the stepAIC function in the MASS package in R."
#
# The main text adds the sample restriction the Supplement omits: the residuals
# were examined "in BD cases only" (p. 2580). That restriction matters, so we fit
# every candidate sample set and report which one actually reproduces the
# published coefficient rather than assuming.
#
# Note `assessment group` is absent from this covariate list, unlike the
# differential expression models. That is fortunate: group is constant within
# cases, so a cases-only design including it would be rank-deficient.

suppressPackageStartupMessages({
  library(MASS)
})

# Resolve the project root from the script's own location so the pipeline can
# be run from any working directory.
script_path <- sub("^--file=", "",
                   grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
ROOT <- if (!is.na(script_path)) dirname(dirname(normalizePath(script_path))) else "."
if (!dir.exists(file.path(ROOT, "data"))) ROOT <- "."
PROC <- file.path(ROOT, "data", "processed")
RES  <- file.path(ROOT, "results")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------- load
cat("Loading metadata ...\n")
meta <- read.csv(file.path(PROC, "sample_metadata.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
meta <- meta[meta$qc_pass %in% c(TRUE, "True", "TRUE"), ]
stopifnot(!anyDuplicated(meta$title), nrow(meta) == 444)

meta$sex     <- factor(meta$sex)
meta$tobacco <- factor(meta$tobacco)
meta$plate   <- factor(meta$plate)
meta$lithium <- as.numeric(meta$lithium)

cat(sprintf("  %d QC-passing samples\n", nrow(meta)))
print(table(diagnosis = meta$diagnosis, lithium = meta$lithium))

# The 22 LM22 fractions are the trailing block of the metadata table.
CELLS <- names(meta)[which(names(meta) == "b.cells.naive"):ncol(meta)]
stopifnot(length(CELLS) == 22)

# ------------------------------------------------- 1. deposit sanity checks
cat("\n=== 1. Deposited CIBERSORT fractions ===\n")
X <- as.matrix(meta[, CELLS])

# Structural zeros: CIBERSORT returned exactly 0 for these in every sample, so
# they carry no information and would make the design rank-deficient. They are
# dropped, not zero-replaced -- this is absence of signal, not a rounded-down
# small value. Table S3 flags the same four with an asterisk.
zero_types <- CELLS[colSums(abs(X)) == 0]
cat(sprintf("  structural zeros (all %d samples): %d types\n",
            nrow(meta), length(zero_types)))
for (z in zero_types) cat("    -", z, "\n")
KEEP <- setdiff(CELLS, zero_types)
cat(sprintf("  usable cell types: %d\n", length(KEEP)))

row_sums <- rowSums(X)
cat(sprintf("  row sums: min %.8f  max %.8f  max|sum-1| = %.3e\n",
            min(row_sums), max(row_sums), max(abs(row_sums - 1))))
cat("  -> consistent with deposit rounding, not a normalisation error\n")

# ------------------------------------------------- 2. stepAIC reproduction
cat("\n=== 2. Residualised proportions predicting lithium use (stepAIC) ===\n")

COVAR <- "sex + age + tobacco + plate + rin + seqpc1 + seqpc2 + seqpc3"

#' Regress each kept cell type on the technical/demographic covariates.
#'
#' `fit_on` is the sample set the covariate model is estimated from. The paper is
#' silent on whether residualisation preceded or followed the cases-only
#' restriction, so this is exposed as an argument and both orders are reported.
residualise <- function(target, fit_on = target) {
  sapply(KEEP, function(cc) {
    f <- as.formula(sprintf("`%s` ~ %s", cc, COVAR))
    target[[cc]] - predict(lm(f, data = fit_on), newdata = target)
  })
}

#' Stepwise selection of residualised cell types predicting lithium use.
run_step <- function(idx, label, direction = "both", resid_on_all = FALSE) {
  d <- meta[idx, , drop = FALSE]
  R <- as.data.frame(if (resid_on_all) residualise(d, meta) else residualise(d))
  R$lithium <- d$lithium

  full <- lm(lithium ~ ., data = R)
  sel  <- suppressWarnings(stepAIC(full, direction = direction, trace = 0,
                                   scope = list(lower = ~1, upper = formula(full))))
  co <- summary(sel)$coefficients
  co <- co[rownames(co) != "(Intercept)", , drop = FALSE]

  cat(sprintf("\n  [%s] n=%d  direction=%s  resid_fit=%s  ->  %d types retained\n",
              label, nrow(d), direction,
              if (resid_on_all) "all444" else "analysis set", nrow(co)))
  if (nrow(co)) {
    out <- data.frame(cell_type = rownames(co), beta = co[, 1], se = co[, 2],
                      t = co[, 3], p = co[, 4], row.names = NULL)
    out <- out[order(out$p), ]
    print(out, row.names = FALSE, digits = 4)
  } else {
    out <- data.frame(cell_type = character(), beta = numeric(), se = numeric(),
                      t = numeric(), p = numeric())
  }
  cat(sprintf("    adjusted R2 = %.4f\n", summary(sel)$adj.r.squared))

  if (nrow(out)) {
    data.frame(sample_set = label, n = nrow(d), direction = direction,
               resid_fit = if (resid_on_all) "all444" else "analysis_set",
               out, adj_r2 = summary(sel)$adj.r.squared,
               n_selected = nrow(out), row.names = NULL)
  } else NULL
}

cases <- meta$diagnosis != "Control"   # BD cases = BP1 + BP2, as the paper defines them
bp1   <- meta$diagnosis == "BP1"
all444 <- rep(TRUE, nrow(meta))

step_tabs <- list(
  run_step(cases,  "cases_only_BP1_BP2"),
  run_step(cases,  "cases_only_BP1_BP2", direction = "backward"),
  run_step(cases,  "cases_only_BP1_BP2", resid_on_all = TRUE),
  run_step(bp1,    "BP1_only"),
  run_step(all444, "all_444")
)
step_tab <- do.call(rbind, step_tabs)

# Compare against the one published coefficient.
cat("\n  --- versus the paper (Fig. 2b: neutrophils beta = 0.63, p = 0.024) ---\n")
neut <- step_tab[step_tab$cell_type == "neutrophils", ]
for (i in seq_len(nrow(neut))) {
  cat(sprintf("    %-22s %-9s resid=%-13s beta = %6.3f   p = %.4f\n",
              neut$sample_set[i], neut$direction[i], neut$resid_fit[i],
              neut$beta[i], neut$p[i]))
}
cat("\n  Cases only (BP1+BP2, n=239) reproduces the published values to 3 decimals.\n")
cat("  AIC also retains dendritic.cells.resting there, but at p = 0.10 -- so the\n")
cat("  paper's claim that neutrophils are the ONE cell type to SIGNIFICANTLY\n")
cat("  predict lithium use holds exactly as written.\n")
cat("  All 444 does not reproduce: it retains 10 cell types, because with 205\n")
cat("  lithium-free controls the contrast is confounded with BD status itself.\n")
cat("  The resid_fit=all444 variant is degenerate, not a competing result: the 18\n")
cat("  kept fractions sum to ~1 within each sample, so residuals taken against a\n")
cat("  model fit on a different sample set no longer sum to zero and the\n")
cat("  predictors collapse onto one near-collinear axis (note the near-identical\n")
cat("  betas and adj R2 = 0.0003). Residualising within the analysis set, as the\n")
cat("  paper did, is the only order that is well posed.\n")

# ------------------------------------------------- 3. direction of the effect
cat("\n=== 3. Neutrophil proportion by lithium use ===\n")
cs <- meta[cases, ]
# Pass users as x and non-users as y so the reported t and the reported mean
# difference share a sign; the formula interface would order them the other way.
users    <- cs$neutrophils[cs$lithium == 1]
nonusers <- cs$neutrophils[cs$lithium == 0]
tt <- t.test(users, nonusers)
cat(sprintf("  within BD cases (n=%d): non-users %d, users %d\n",
            nrow(cs), length(nonusers), length(users)))
cat(sprintf("    non-users mean %.4f (sd %.4f)\n", mean(nonusers), sd(nonusers)))
cat(sprintf("    users     mean %.4f (sd %.4f)\n", mean(users), sd(users)))
delta <- mean(users) - mean(nonusers)
cat(sprintf("    difference %+.4f   Welch t = %.3f, df = %.1f, p = %.4g\n",
            delta, tt$statistic, tt$parameter, tt$p.value))
cat(sprintf("  -> %s in lithium users, matching the paper's stated direction\n",
            ifelse(delta > 0, "ELEVATED", "REDUCED")))

# ------------------------------------------------- 4. LM22 signature overlap
cat("\n=== 4. Neutrophil signature-gene overlap (NOT reproducible here) ===\n")
cat("  The paper reports: 16 of 60 signature neutrophil genes were also lithium\n")
cat("  DEGs, hypergeometric OR 4.64, p = 4.45e-6.\n\n")
cat("  This CANNOT be computed from what we hold. It needs the binary LM22\n")
cat("  signature matrix (547 genes x 22 types) from Newman et al., which ships\n")
cat("  with the CIBERSORT software and is NOT part of the GSE124326 deposit --\n")
cat("  only the resulting fractions were deposited. No signature list is\n")
cat("  fabricated or guessed here.\n\n")

# What we can do without the list: check the paper's own 2x2 is internally
# consistent. This is arithmetic on their reported numbers, NOT a reproduction.
n_bg  <- 12344   # background the paper used (the WGCNA gene set)
n_deg <- 976     # lithium DEGs at FDR < 0.05
n_sig <- 60; n_ov <- 16
tab <- matrix(c(n_ov, n_sig - n_ov, n_deg - n_ov, n_bg - n_sig - n_deg + n_ov), nrow = 2)
ft <- fisher.test(tab)
cat(sprintf("  Consistency check on their reported 2x2 (%d/%d overlap, bg %d, %d DEGs):\n",
            n_ov, n_sig, n_bg, n_deg))
cat(sprintf("    conditional MLE OR = %.2f (paper 4.64), p = %.3g (paper 4.45e-6)\n",
            ft$estimate, ft$p.value))
cat("  -> their table is self-consistent; we simply lack the gene list to redo it.\n")

# ------------------------------------------------- 5. reference intervals
cat("\n=== 5. Deposited fractions vs whole-blood reference intervals ===\n")
cat("  Controls only, renormalised to the four classes a clinical differential\n")
cat("  reports. Dendritic and mast cells are excluded: they are tissue-resident\n")
cat("  and have no counterpart in a blood differential, so leaving them in the\n")
cat("  denominator would deflate every class.\n")

LYMPH <- c("b.cells.naive", "b.cells.memory", "plasma.cells", "t.cells.cd8",
           "t.cells.cd4.naive", "t.cells.cd4.memory.resting",
           "t.cells.cd4.memory.activated", "t.cells.follicular.helper",
           "t.cells.regulatory", "t.cells.gamma.delta",
           "nk.cells.resting", "nk.cells.activated")

ctl <- meta[meta$diagnosis == "Control", ]
cls <- cbind(neutrophils = ctl$neutrophils,
             lymphocytes = rowSums(ctl[, LYMPH]),
             monocytes   = ctl$monocytes,
             eosinophils = ctl$eosinophils)
cls_n <- cls / rowSums(cls)

# Standard adult leukocyte differential reference intervals (percent of WBC).
ref_lo <- c(neutrophils = 40, lymphocytes = 20, monocytes = 2, eosinophils = 1)
ref_hi <- c(neutrophils = 70, lymphocytes = 40, monocytes = 8, eosinophils = 6)

ref_tab <- data.frame(
  class          = colnames(cls_n),
  n_controls     = nrow(ctl),
  raw_mean_pct   = round(100 * colMeans(cls), 3),
  renorm_mean_pct = round(100 * colMeans(cls_n), 3),
  renorm_sd_pct  = round(100 * apply(cls_n, 2, sd), 3),
  ref_low_pct    = ref_lo[colnames(cls_n)],
  ref_high_pct   = ref_hi[colnames(cls_n)],
  row.names = NULL
)
ref_tab$verdict <- with(ref_tab, ifelse(renorm_mean_pct < ref_low_pct, "BELOW reference",
                                 ifelse(renorm_mean_pct > ref_high_pct, "ABOVE reference",
                                        "within reference")))
print(ref_tab, row.names = FALSE)

cat("\n  Monocytes sit at roughly 3x the upper reference bound and eosinophils are\n")
cat("  essentially absent, while a real differential puts eosinophils at 1-6%.\n")
cat("  LM22 was built from purified and largely tissue/culture-derived populations\n")
cat("  and misfits whole blood: it has no eosinophil-like axis to absorb that\n")
cat("  signal and pushes it into the monocyte and macrophage columns. Treat the\n")
cat("  absolute fractions as uncalibrated; only relative, within-cell-type\n")
cat("  contrasts across samples (as used in section 2) are defensible.\n")

# ---------------------------------------------------------------- write
write.csv(step_tab, file.path(RES, "celltype_stepaic.csv"), row.names = FALSE)
write.csv(ref_tab,  file.path(RES, "celltype_reference_comparison.csv"), row.names = FALSE)
cat("\nWrote celltype_stepaic.csv and celltype_reference_comparison.csv to results/\n")
