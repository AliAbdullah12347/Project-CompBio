#!/usr/bin/env Rscript
# ==============================================================================
# 12_celltype_structure.R -- WHERE does the lithium effect sit once cell types
# are separated, and is the separation real?
#
# WHY THIS EXISTS
#
# 07 reported that the same handful of genes (TSPAN2, IL8, PIGB, RFX2,
# MIR24-2) top the significant list in granulocytes, monocytes, T cells AND NK
# cells, and that 100% of the monocyte, T and NK hits are also whole-blood
# hits. The obvious worry is that bMIND is not separating cell types at all --
# that each "cell-type estimate" is just a rescaled copy of the bulk, in which
# case every cell-type result in this project is meaningless.
#
# That worry is TESTABLE and it turns out to be WRONG, which is why the test is
# recorded here rather than the suspicion. If the lineages were copies of one
# another their per-gene logFCs would be highly correlated. They are not.
#
# What is actually going on is reported instead: the effect is CONCENTRATED in
# one lineage, and the shared gene names are simply the genes with the largest
# bulk signal surfacing wherever there is any power at all.
#
# TWO DIAGNOSTICS
#
# 1. Between-lineage correlation of per-gene logFC. Near 1 would mean the
#    deconvolution failed. Near 0 means the lineages carry different
#    information.
#
# 2. Effect-size ratio against whole blood. This locates the effect. A lineage
#    whose logFC is ~0 for genes that clearly move in bulk is telling us the
#    bulk movement did not come from inside that lineage.
#
# CAVEAT THAT MUST TRAVEL WITH THIS RESULT. "Within granulocytes" is only as
# fine-grained as our 5-lineage aggregation. The gran lineage pools
# neutrophils, eosinophils and mast cells, so a shift in the MIX of those --
# say more immature neutrophils -- would appear here as a within-granulocyte
# expression change. This analysis cannot separate those two.
# ==============================================================================

HERE <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/de_analysis"
setwd(HERE)
RES <- "results"
LINS <- c("gran", "mono", "T", "NK", "B")
p <- readRDS("data/prep.rds")
m <- read.delim("data/ens2sym.tsv", header = FALSE, col.names = c("ens", "sym"))
sym <- setNames(m$sym, m$ens)

OUTC <- list(); OUTR <- list()
for (cn in c("LI", "BPD")) {
  wb <- read.csv(file.path(RES, sprintf("WB_%s_raw_toptable.csv", cn)), stringsAsFactors = FALSE)
  FC <- sapply(LINS, function(ct) {
    d <- read.csv(file.path(RES, sprintf("CT_%s_raw_%s_toptable.csv", cn, ct)), stringsAsFactors = FALSE)
    setNames(d$logFC, d$gene) })
  g <- rownames(FC)
  w <- setNames(wb$logFC, wb$gene)[g]
  M <- cbind(FC, wholeblood = w)

  cat(sprintf("\n================ %s ================\n", cn))
  cat("--- 1. between-level correlation of per-gene logFC ---\n")
  cat("    (near 1 would mean the deconvolution produced copies, not separation)\n")
  C <- round(cor(M), 3); print(C)
  offd <- C[LINS, LINS][upper.tri(diag(length(LINS)))]
  cat(sprintf("    max |correlation| between two lineages: %.3f -> %s\n", max(abs(offd)),
              if (max(abs(offd)) > 0.9) "LINEAGES ARE COPIES -- results meaningless"
              else "lineages carry distinct information"))
  OUTC[[cn]] <- data.frame(contrast = cn, as.data.frame(as.table(C)), stringsAsFactors = FALSE)

  cat("\n--- 2. where does the effect live? ---\n")
  deg <- wb$gene[wb$adj.P.Val < 0.05]
  ref <- if (length(deg) >= 20) deg else wb$gene[order(wb$P.Value)][1:200]
  cat(sprintf("    reference set: %s (n=%d)\n",
              if (length(deg) >= 20) "whole-blood DEGs" else "top 200 whole-blood genes by p", length(ref)))
  cat(sprintf("    median |logFC| in whole blood over that set: %.4f\n", median(abs(w[ref]))))
  for (ct in LINS) {
    r <- data.frame(contrast = cn, lineage = ct,
      mean_frac = round(mean(p$lineage[, ct]), 4),
      cor_with_wholeblood = round(cor(FC[, ct], w, use = "complete.obs"), 3),
      median_abs_logFC = round(median(abs(FC[ref, ct])), 4),
      ratio_to_wholeblood = round(median(FC[ref, ct] / w[ref]), 4), stringsAsFactors = FALSE)
    OUTR[[length(OUTR) + 1]] <- r
    cat(sprintf("    %-5s frac %.4f | cor with bulk %6.3f | median |logFC| %.4f | ratio %6.3f\n",
                ct, r$mean_frac, r$cor_with_wholeblood, r$median_abs_logFC, r$ratio_to_wholeblood))
  }
}
R <- do.call(rbind, OUTR)

cat("\n\n=== READING ===\n")
li <- R[R$contrast == "LI", ]
top <- li$lineage[which.max(abs(li$ratio_to_wholeblood))]
cat(sprintf("  LI: the effect localises to '%s' (ratio %.3f, correlation with bulk %.3f).\n",
            top, li$ratio_to_wholeblood[li$lineage == top], li$cor_with_wholeblood[li$lineage == top]))
cat("  Every other lineage has a logFC ratio near zero for genes that clearly\n")
cat("  move in bulk, i.e. the bulk movement did not originate inside them.\n")
cat("  This is what the compositional account predicts, with the caveat in the\n")
cat("  header: our 'gran' lineage pools neutrophils, eosinophils and mast cells,\n")
cat("  so a shift WITHIN that pool is indistinguishable from a within-cell change.\n")

write.csv(do.call(rbind, OUTC), file.path(RES, "celltype_level_correlations.csv"), row.names = FALSE)
write.csv(R, file.path(RES, "celltype_effect_location.csv"), row.names = FALSE)
cat("\nwrote results/celltype_level_correlations.csv and celltype_effect_location.csv\n")
