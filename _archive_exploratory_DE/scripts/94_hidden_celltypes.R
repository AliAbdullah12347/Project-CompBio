#!/usr/bin/env Rscript
# ==============================================================================
# 94_hidden_celltypes.R -- is the "direct effect" really composition LM22 cannot see?
#
# 92 found ~400 genes whose lithium association survives adjustment for the cell
# balances, and 93 showed the two strongest are IL5RA and CLC (galectin-10).
# Both are EOSINOPHIL/BASOPHIL genes. That matters because:
#
#   * LM22 has NO basophil cell type at all. Basophils are invisible to it.
#   * LM22's eosinophil estimate in this dataset is ~0.1% of leukocytes against
#     a clinical reference of 1-6%, i.e. it effectively cannot see them either.
#
# Adjusting expression for an estimated mixture can only remove composition in
# the cell types the estimate covers. If lithium shifts eosinophils or basophils,
# that shift is PURE composition but will appear as a "direct effect" simply
# because the mediator is blind to it.
#
# This distinguishes a real biological direct effect from a measurement gap, and
# the two have opposite implications for the project's central claim.
#
# Marker panels are defined from the literature and fixed BEFORE looking at the
# results, so this is a test rather than a description of what we already saw.
# ==============================================================================

source("scripts/common.R")

SETS <- list(
  # Eosinophil: granule proteins, lineage receptors, and the eosinophil/basophil
  # shared galectin CLC. IL5RA is the defining eosinophil/basophil growth-factor
  # receptor; SIGLEC8 and CCR3 are surface markers; EPX/PRG2/RNASE2/RNASE3 are
  # granule contents.
  eosinophil = c("IL5RA", "CLC", "SIGLEC8", "CCR3", "EPX", "PRG2", "RNASE2",
                 "RNASE3", "PRSS33", "ALOX15", "CYSLTR2", "IL1RL1"),
  # Basophil: FcepsilonRI beta chain, histidine decarboxylase, tryptase, and
  # the basophil-restricted transcription factor GATA2. NONE of these belong to
  # any LM22 cell type.
  basophil   = c("MS4A2", "HDC", "CPA3", "GATA2", "ENPP3", "IL3RA", "AKAP12"),
  # Positive control: LM22 sees neutrophils well, so these should be largely
  # ACCOUNTED FOR by the balance adjustment.
  neutrophil = c("FCGR3B", "CSF3R", "CEACAM3", "CXCR2", "S100A8", "S100A9",
                 "S100A12", "MMP8", "MMP9", "LCN2", "CAMP", "LTF", "PGLYRP1"),
  # Negative control: LM22 sees T cells well.
  tcell      = c("CD3E", "CD3D", "CD2", "IL7R", "CCR7", "LEF1", "TCF7", "CD27",
                 "CD28", "LCK", "ZAP70", "THEMIS")
)

d <- read.csv(file.path(EXP_HOME, "runs", "just_composition", "per_gene.csv"),
              stringsAsFactors = FALSE)
d <- d[!is.na(d$symbol), ]
direct <- d$fdr_direct < 0.05
cat(sprintf("%d genes annotated | %d with a direct effect at FDR 0.05 (%.1f%%)\n\n",
            nrow(d), sum(direct), 100 * mean(direct)))

res <- do.call(rbind, lapply(names(SETS), function(s) {
  g <- d[d$symbol %in% SETS[[s]], ]
  if (!nrow(g)) return(NULL)
  k <- sum(g$fdr_direct < 0.05); n <- nrow(g)
  # Hypergeometric: is this panel over-represented among direct-effect genes?
  ph <- stats::phyper(k - 1, sum(direct), nrow(d) - sum(direct), n, lower.tail = FALSE)
  data.frame(
    set = s, present = n, of = length(SETS[[s]]),
    mean_beta = mean(g$beta_lithium),
    mean_predicted_by_composition = mean(g$predicted_from_composition),
    mean_abs_residual = mean(abs(g$residual)),
    # The key ratio: what share of the marginal effect does composition explain?
    pct_explained = round(100 * mean(abs(g$predicted_from_composition)) /
                          mean(abs(g$beta_lithium)), 1),
    n_direct = k, pct_direct = round(100 * k / n, 1),
    hyper_p = signif(ph, 3), stringsAsFactors = FALSE)
}))
cat("=== does composition adjustment explain each marker panel? ===\n")
print(res, row.names = FALSE, digits = 3)
cat(sprintf("\nbackground rate of a direct effect: %.1f%%\n", 100 * mean(direct)))

cat("\n=== eosinophil / basophil genes, gene by gene ===\n")
eb <- d[d$symbol %in% c(SETS$eosinophil, SETS$basophil), ]
eb$panel <- ifelse(eb$symbol %in% SETS$eosinophil, "eos", "baso")
print(eb[order(-eb$beta_lithium),
         c("symbol", "panel", "beta_lithium", "predicted_from_composition",
           "residual", "fdr_marginal", "fdr_direct")], row.names = FALSE, digits = 3)

## ---- the decisive comparison ----------------------------------------------
# For a cell type LM22 SEES, the composition term should carry most of the
# marginal effect. For one it CANNOT see, it should carry almost none, and the
# effect should survive adjustment.
cat("\n=== decisive comparison ===\n")
vis <- res[res$set %in% c("neutrophil", "tcell"), ]
inv <- res[res$set %in% c("eosinophil", "basophil"), ]
cat(sprintf("  LM22 CAN see these  (neutrophil, T cell) : composition explains %.0f%% of |beta|, %.0f%% direct\n",
            mean(vis$pct_explained), mean(vis$pct_direct)))
cat(sprintf("  LM22 CANNOT see these (eos, basophil)    : composition explains %.0f%% of |beta|, %.0f%% direct\n",
            mean(inv$pct_explained), mean(inv$pct_direct)))
cat("\n  If the second row shows much less explained and many more direct hits,\n")
cat("  then part of the apparent DIRECT effect is composition in cell types the\n")
cat("  mediator is blind to -- not lithium acting inside cells.\n")

write.csv(res, file.path(EXP_HOME, "runs", "just_composition", "marker_panels.csv"),
          row.names = FALSE)
write.csv(eb, file.path(EXP_HOME, "runs", "just_composition", "eos_baso_genes.csv"),
          row.names = FALSE)
cat("\nwrote marker_panels.csv and eos_baso_genes.csv\n")
