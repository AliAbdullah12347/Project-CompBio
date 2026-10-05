#!/usr/bin/env Rscript
# ==============================================================================
# 95_eos_abundance_vs_induction.R -- abundance or induction?
#
# 94 showed a coherent eosinophil-axis signal (IL5RA, CLC, CYSLTR2, IL1RL1) that
# rises with lithium and that adjusting for the LM22 balances does not remove.
# There are two very different explanations:
#
#   ABUNDANCE  lithium raises circulating eosinophils/basophils. The genes go up
#              because there are more of those cells. This is COMPOSITION, and
#              the only reason it looks "direct" is that LM22 cannot measure
#              those cells. The project's mediation estimate is then an
#              UNDERestimate of the compositional share.
#
#   INDUCTION  lithium raises these transcripts inside cells. A genuine direct
#              effect, and the mediation estimate is right.
#
# They make different predictions about correlation structure:
#   * Under ABUNDANCE the marker genes must move together tightly across
#     samples, because one latent quantity (cell number) drives all of them,
#     and that latent factor should also track the LM22 eosinophil estimate --
#     noisy and tiny though it is, it is not nothing.
#   * Under INDUCTION they need not co-vary beyond the shared exposure.
#
# Benchmarking against the neutrophil panel gives the scale: neutrophils are a
# known-abundance-driven axis in this data, so whatever tightness they show is
# what an abundance signal looks like here.
# ==============================================================================

source("scripts/common.R")
p <- load_prep()
set.seed(481)

EOS  <- c("IL5RA", "CLC", "CYSLTR2", "IL1RL1")
NEUT <- c("FCGR3B", "CSF3R", "CEACAM3", "CXCR2", "S100A8", "S100A9", "S100A12")
TC   <- c("CD3E", "CD3D", "CD2", "IL7R", "LCK", "ZAP70")

ann <- read.csv(file.path(EXP_HOME, "runs", "just_composition", "per_gene.csv"),
                stringsAsFactors = FALSE)
sym2ens <- setNames(ann$gene, ann$symbol)

m <- p$meta[p$meta$dx == "BP1", ]
m$tob <- p$TOB[m$title, "tobacco_imp_01"]
covars <- c("age", "sex", "tob", "rin", "plate", "seqpc1", "seqpc2", "seqpc3")
m <- m[complete.cases(m[, c("lithium", covars)]), ]
m$plate <- droplevels(m$plate); m$sex <- droplevels(m$sex)

cnt <- p$counts[, m$title, drop = FALSE]
cnt <- cnt[filter_genes(cnt, 10, 0.90), , drop = FALSE]
lcpm <- cpm(calcNormFactors(DGEList(cnt), method = "TMM"), log = TRUE, prior.count = 3)

# Residualise on the covariates so the co-variation measured is not just shared
# technical or demographic structure.
X0 <- model.matrix(as.formula(paste("~", paste(covars, collapse = " + "))), data = m)
R <- t(residuals(lm.fit(X0, t(lcpm))))

panel_stats <- function(syms, label) {
  ids <- na.omit(sym2ens[syms]); ids <- ids[ids %in% rownames(R)]
  if (length(ids) < 2) return(NULL)
  C <- cor(t(R[ids, , drop = FALSE]))
  mean_r <- mean(C[upper.tri(C)])
  # First principal component = the latent axis the panel shares.
  pc <- prcomp(t(scale(t(R[ids, , drop = FALSE]))))
  frac1 <- (pc$sdev^2 / sum(pc$sdev^2))[1]
  score <- pc$rotation[, 1]
  if (mean(cor(score, t(R[ids, , drop = FALSE]))) < 0) score <- -score
  cat(sprintf("\n  %-12s n=%d  mean pairwise r = %+.3f  PC1 explains %.0f%% of panel variance\n",
              label, length(ids), mean_r, 100 * frac1))
  list(ids = ids, mean_r = mean_r, frac1 = frac1, score = score, label = label)
}

cat("=== co-variation within each marker panel (covariate-residualised) ===")
S <- Filter(Negate(is.null), list(panel_stats(EOS, "eosinophil"),
                                  panel_stats(NEUT, "neutrophil"),
                                  panel_stats(TC, "T cell")))

## ---- does each panel score track the deconvolution estimate? --------------
cat("\n=== panel score vs the LM22 estimate of that same cell type ===\n")
frac <- p$frac[m$title, ]
lin  <- p$lineage[m$title, ]
targets <- list(eosinophil = frac[, "Eosinophils"],
                neutrophil = frac[, "Neutrophils"],
                `T cell`   = lin[, "T"])
for (s in S) {
  tv <- targets[[s$label]]
  r <- suppressWarnings(cor(s$score, tv))
  rs <- suppressWarnings(cor(rank(s$score), rank(tv)))
  cat(sprintf("  %-12s vs LM22 %-12s  Pearson %+.3f  Spearman %+.3f   (LM22 mean %.4f, nonzero in %.0f%%)\n",
              s$label, s$label, r, rs, mean(tv), 100 * mean(tv > 0)))
}

## ---- does lithium move each panel score? ----------------------------------
cat("\n=== lithium effect on each panel score ===\n")
for (s in S) {
  f <- lm(s$score ~ m$lithium)
  cf <- summary(f)$coefficients["m$lithium", ]
  d <- (mean(s$score[m$lithium == 1]) - mean(s$score[m$lithium == 0])) / sd(s$score)
  cat(sprintf("  %-12s  beta %+.4f (p = %.3g)  Cohen d = %+.3f\n",
              s$label, cf[1], cf[4], d))
}

out <- data.frame(
  panel = vapply(S, function(s) s$label, character(1)),
  n_genes = vapply(S, function(s) length(s$ids), integer(1)),
  mean_pairwise_r = round(vapply(S, function(s) s$mean_r, numeric(1)), 3),
  pc1_var_explained = round(vapply(S, function(s) s$frac1, numeric(1)), 3),
  cor_with_lm22 = round(vapply(S, function(s)
    suppressWarnings(cor(s$score, targets[[s$label]])), numeric(1)), 3),
  stringsAsFactors = FALSE)
cat("\n=== verdict inputs ===\n"); print(out, row.names = FALSE)
cat("\n  A tightly co-varying eosinophil panel (mean r comparable to the\n")
cat("  neutrophil benchmark) points to ABUNDANCE -- one latent cell-number\n")
cat("  factor driving all four genes -- rather than independent induction.\n")

dir.create(file.path(EXP_HOME, "runs", "eos_abundance"), showWarnings = FALSE, recursive = TRUE)
write.csv(out, file.path(EXP_HOME, "runs", "eos_abundance", "summary.csv"), row.names = FALSE)
cat("\nwrote runs/eos_abundance/summary.csv\n")
