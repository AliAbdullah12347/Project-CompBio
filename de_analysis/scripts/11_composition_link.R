#!/usr/bin/env Rscript
# ==============================================================================
# 11_composition_link.R -- is the lithium signature a CELL-COMPOSITION
# signature? Tested against our own data, not against prior literature.
#
# WHY THIS TEST
#
# The top lithium genes look like granulocyte genes. "Look like" is not
# evidence, and naming genes from memory is exactly the kind of unchecked claim
# this project is supposed to avoid. The checkable version of the claim is:
#
#     if the lithium signature is driven by a shift in the cell mixture, then
#     the genes it flags should be the genes whose expression tracks the
#     granulocyte fraction across subjects.
#
# That correlation is computed from our own matrix and our own fractions, over
# all 474 subjects, and it does not depend on any external annotation.
#
# A NOTE ON WHAT THIS CAN AND CANNOT SHOW. Both the fractions and the
# expression come from the same matrix, so a marker gene is correlated with its
# own cell type partly by construction. This is the circularity recorded in
# CLAUDE.md 4.2. It is why the test below is a DESCRIPTIVE check on the
# direction of the effect, not a formal mediation estimate -- that is Arm 2,
# and it has to exclude the LM22 marker genes to be meaningful.
# ==============================================================================

HERE <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/de_analysis"
setwd(HERE); set.seed(481)
RES <- "results"
p <- readRDS("data/prep.rds")
m <- read.delim("data/ens2sym.tsv", header = FALSE, col.names = c("ens", "sym"))
sym <- setNames(m$sym, m$ens)

E <- p$logtpm                       # genes x 474
gran <- p$lineage[colnames(E), "gran"]
cat(sprintf("== 11_composition_link ==\n%d genes x %d subjects | granulocyte fraction %.3f-%.3f\n",
            nrow(E), ncol(E), min(gran), max(gran)))

# Spearman, so a few extreme subjects cannot manufacture the association.
rg <- apply(E, 1, function(x) suppressWarnings(cor(x, gran, method = "spearman")))

li <- read.csv(file.path(RES, "WB_LI_raw_toptable.csv"), stringsAsFactors = FALSE)
bp <- read.csv(file.path(RES, "WB_BPD_raw_toptable.csv"), stringsAsFactors = FALSE)
li_deg <- li$gene[li$adj.P.Val < 0.05]
li_non <- li$gene[li$adj.P.Val >= 0.05]

cat("\n=== do the lithium genes track the granulocyte fraction? ===\n")
cat(sprintf("  lithium DEGs    (n=%5d): median |rho| = %.3f\n", length(li_deg), median(abs(rg[li_deg]))))
cat(sprintf("  non-DEGs        (n=%5d): median |rho| = %.3f\n", length(li_non), median(abs(rg[li_non]))))
w <- wilcox.test(abs(rg[li_deg]), abs(rg[li_non]))
cat(sprintf("  Wilcoxon p = %.3g\n", w$p.value))
cat(sprintf("  |rho| > 0.3 : %.1f%% of DEGs vs %.1f%% of non-DEGs\n",
            100 * mean(abs(rg[li_deg]) > 0.3), 100 * mean(abs(rg[li_non]) > 0.3)))

# Direction agreement: lithium raises granulocytes, so a gene that goes UP with
# lithium should also go UP with granulocyte fraction if composition drives it.
fc <- setNames(li$logFC, li$gene)
agree <- sign(fc[li_deg]) == sign(rg[li_deg])
cat(sprintf("\n  sign(logFC_lithium) == sign(rho_granulocyte) for %.1f%% of the %d DEGs\n",
            100 * mean(agree, na.rm = TRUE), length(li_deg)))
cat(sprintf("  (50%% would be chance; binomial p = %.3g)\n",
            binom.test(sum(agree, na.rm = TRUE), sum(!is.na(agree)))$p.value))

cat("\n=== top 15 lithium genes, with their granulocyte correlation ===\n")
t15 <- head(li[order(li$P.Value), ], 15)
t15$symbol <- ifelse(is.na(sym[t15$gene]), t15$gene, sym[t15$gene])
t15$rho_gran <- round(rg[t15$gene], 3)
print(t15[, c("symbol", "logFC", "adj.P.Val", "rho_gran")], row.names = FALSE, digits = 3)

cat("\n=== the 4 whole-blood BPD genes, same treatment ===\n")
b4 <- head(bp[order(bp$P.Value), ], 6)
b4$symbol <- ifelse(is.na(sym[b4$gene]), b4$gene, sym[b4$gene])
b4$rho_gran <- round(rg[b4$gene], 3)
print(b4[, c("symbol", "logFC", "adj.P.Val", "rho_gran")], row.names = FALSE, digits = 3)

## ---- how much does adjusting for composition shrink the effect? -------------
li_ilr <- read.csv(file.path(RES, "WB_LI_ilr_toptable.csv"), stringsAsFactors = FALSE)
fi <- setNames(li_ilr$logFC, li_ilr$gene)
shr <- 1 - abs(fi[li_deg]) / abs(fc[li_deg])
cat(sprintf("\n=== effect shrinkage after adjusting for the 4 ILR balances ===\n"))
cat(sprintf("  median shrinkage of |logFC| across the %d lithium DEGs: %.1f%%\n",
            length(li_deg), 100 * median(shr, na.rm = TRUE)))
cat(sprintf("  quartiles: %.1f%% / %.1f%% / %.1f%%\n",
            100 * quantile(shr, .25, na.rm = TRUE), 100 * median(shr, na.rm = TRUE),
            100 * quantile(shr, .75, na.rm = TRUE)))
cat(sprintf("  DEGs surviving adjustment: %d of %d (%.1f%%)\n",
            sum(li_ilr$adj.P.Val[match(li_deg, li_ilr$gene)] < 0.05, na.rm = TRUE),
            length(li_deg),
            100 * mean(li_ilr$adj.P.Val[match(li_deg, li_ilr$gene)] < 0.05, na.rm = TRUE)))
cat("\n  Shrinkage is a LOWER BOUND on the compositional share: the balances are\n")
cat("  estimated from the same matrix and are only 4 coordinates on 5 lineages,\n")
cat("  so they cannot absorb composition that LM22 never resolved.\n")

OUT <- data.frame(
  gene = li$gene, symbol = ifelse(is.na(sym[li$gene]), li$gene, sym[li$gene]),
  logFC_raw = li$logFC, adjP_raw = li$adj.P.Val,
  logFC_ilr = fi[li$gene], adjP_ilr = li_ilr$adj.P.Val[match(li$gene, li_ilr$gene)],
  rho_gran = round(rg[li$gene], 4), stringsAsFactors = FALSE)
write.csv(OUT[order(OUT$adjP_raw), ], file.path(RES, "composition_link_LI.csv"), row.names = FALSE)
write.csv(data.frame(set = c("LI_DEG", "LI_nonDEG"),
                     n = c(length(li_deg), length(li_non)),
                     median_abs_rho = c(median(abs(rg[li_deg])), median(abs(rg[li_non]))),
                     pct_abs_rho_gt_0.3 = c(100 * mean(abs(rg[li_deg]) > 0.3),
                                            100 * mean(abs(rg[li_non]) > 0.3))),
          file.path(RES, "composition_link_summary.csv"), row.names = FALSE)
cat("\nwrote results/composition_link_*.csv\n")
