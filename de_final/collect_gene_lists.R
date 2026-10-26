#!/usr/bin/env Rscript
# ==============================================================================
# collect_gene_lists.R -- every significant gene from every analysis, in one
# place, named and ordered, ready to drop into the master document.
#
# Four sources, kept separate because they answer different questions and are
# NOT subsets of one another:
#
#   1. Whole blood, lithium    1,426 genes at BH 0.05, of which 125 also pass
#                              Holm (family-wise error) -- the high-confidence set
#   2. Whole blood, bipolar    4 genes, all on the significance boundary
#   3. Cell type, lithium      161 gene x lineage pairs under per-lineage BH
#   4. Cell type, bipolar      none
#
# Plus the two rankings that do not depend on a significance threshold at all:
# the rank product across the five lineages, and the filter-stability ranking
# from the 31-filter sweep.
# ==============================================================================

ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "de_final")); options(width = 200)
RES <- "results"; OUT <- file.path(RES, "genelists"); dir.create(OUT, showWarnings = FALSE)

rd <- function(f) read.csv(gzfile(file.path(RES, f)), stringsAsFactors = FALSE)
wrap <- function(v, per = 10) {
  v <- sort(unique(v))
  paste(vapply(split(v, ceiling(seq_along(v)/per)), paste, character(1), collapse = ", "),
        collapse = "\n") }

cat("== collect_gene_lists ==\n")

## ---- 1. whole blood, lithium -------------------------------------------------
li <- rd("DE_WB_LI.csv.gz"); li <- li[li$significant, ]; li <- li[order(li$P), ]
V2 <- file.path(ROOT, "de_v2/results")
cand <- read.csv(file.path(V2, "candidate_genes_LI_significant.csv"), stringsAsFactors = FALSE)
holm <- cand$symbol[cand$holm < 0.05]
cat(sprintf("WB lithium  : %d at BH 0.05 | %d also pass Holm\n", nrow(li), length(holm)))
write.csv(li[, c("symbol","gene","rank","logFC","AveExpr","P","BH")],
          file.path(OUT, "WB_lithium_all_1426.csv"), row.names = FALSE)
writeLines(wrap(li$symbol), file.path(OUT, "WB_lithium_all_symbols.txt"))
writeLines(wrap(holm), file.path(OUT, "WB_lithium_holm_symbols.txt"))

cat("\n--- WB lithium, top 40 by p-value ---\n")
t40 <- head(li, 40)
print(data.frame(rank = t40$rank, symbol = t40$symbol, logFC = round(t40$logFC, 3),
                 fold = round(2^abs(t40$logFC), 2),
                 dir = ifelse(t40$logFC > 0, "up", "down"),
                 BH = signif(t40$BH, 3)), row.names = FALSE)

cat(sprintf("\n--- the %d Holm-significant genes (high confidence) ---\n", length(holm)))
cat(wrap(holm, 8), "\n")

## ---- 2. whole blood, bipolar -------------------------------------------------
bp <- rd("DE_WB_BPD.csv.gz"); bps <- bp[bp$significant, ]; bps <- bps[order(bps$P), ]
cat(sprintf("\n--- WB bipolar: all %d significant genes ---\n", nrow(bps)))
print(data.frame(rank = bps$rank, symbol = bps$symbol, logFC = round(bps$logFC, 3),
                 dir = ifelse(bps$logFC > 0, "up", "down"),
                 P = signif(bps$P, 3), BH = signif(bps$BH, 4)), row.names = FALSE)
write.csv(bps[, c("symbol","gene","rank","logFC","AveExpr","P","BH")],
          file.path(OUT, "WB_bipolar_all_4.csv"), row.names = FALSE)

## ---- 3. cell type, lithium ---------------------------------------------------
ct <- rd("DE_CT_LI.csv.gz"); cts <- ct[ct$significant, ]
cat(sprintf("\n--- cell type, lithium: %d gene x lineage pairs, %d distinct genes ---\n",
            nrow(cts), length(unique(cts$symbol))))
for (L in c("gran","mono","T","NK","B")) {
  z <- cts[cts$lineage == L, ]; z <- z[order(z$P), ]
  cat(sprintf("\n  %s  (%d genes)\n", L, nrow(z)))
  if (nrow(z)) cat("   ", wrap(z$symbol, 9), "\n") }
write.csv(cts[order(cts$lineage, cts$P),
              c("symbol","gene","lineage","rank","logFC","P","BH")],
          file.path(OUT, "CT_lithium_all_161.csv"), row.names = FALSE)

shared <- table(cts$symbol)
cat(sprintf("\n  genes significant in MORE THAN ONE lineage (%d of %d):\n",
            sum(shared > 1), length(shared)))
for (n in sort(unique(shared[shared > 1]), decreasing = TRUE))
  cat(sprintf("    in %d lineages: %s\n", n, paste(sort(names(shared)[shared == n]), collapse = ", ")))

## ---- 4. cell type, bipolar ---------------------------------------------------
cb <- rd("DE_CT_BPD.csv.gz")
cat(sprintf("\n--- cell type, bipolar: %d significant ---\n", sum(cb$significant)))
cat("  (the 10 closest to significance, for completeness)\n")
cbt <- head(cb[order(cb$P), ], 10)
print(data.frame(symbol = cbt$symbol, lineage = cbt$lineage, logFC = round(cbt$logFC, 3),
                 P = signif(cbt$P, 3), BH = signif(cbt$BH, 3)), row.names = FALSE)

## ---- 5. tier 1 ---------------------------------------------------------------
t1 <- read.csv(file.path(V2, "candidate_genes_tier1.csv"), stringsAsFactors = FALSE)
cat(sprintf("\n--- Tier 1: %d genes passing every test ---\n", nrow(t1)))
print(t1[order(t1$P), c("symbol","logFC","holm","filters_sig","filters_tested",
                        "rho_granulocyte","logFC_adj","pct_effect_lost","lineages")],
      row.names = FALSE)

## ---- 6. rank product ----------------------------------------------------------
for (CN in c("LI","BPD")) {
  rp <- rd(sprintf("followup_B_rankproduct_%s.csv.gz", CN))
  cat(sprintf("\n--- rank product across the 5 lineages, %s: top 30 ---\n", CN))
  print(head(rp[, c("symbol","rank_product","gran","mono","T","NK","B")], 30), row.names = FALSE)
  write.csv(head(rp, 200), file.path(OUT, sprintf("rankproduct_top200_%s.csv", CN)), row.names = FALSE) }

## ---- 7. filter-robust core ----------------------------------------------------
rb <- read.csv(file.path(ROOT, "de_v2/results/sweep/wb_gene_robustness.csv"), stringsAsFactors = FALSE)
rb <- rb[rb$analysis == "LI_raw", ]
core <- rb[rb$n_filters_sig == rb$n_filters_tested, ]
core <- core[order(core$median_rank), ]
cat(sprintf("\n--- filter-robust core: %d genes called under EVERY filter that retained them ---\n",
            nrow(core)))
cat("  top 40 by median rank:\n   ", wrap(head(core$symbol, 40), 8), "\n")
write.csv(core[, c("symbol","gene","n_filters_sig","n_filters_tested","best_rank",
                   "median_rank","median_logFC","min_BH")],
          file.path(OUT, "WB_lithium_filter_robust_core.csv"), row.names = FALSE)
writeLines(wrap(core$symbol), file.path(OUT, "WB_lithium_robust_core_symbols.txt"))

cat(sprintf("\nwrote %d files to results/genelists/\n", length(list.files(OUT))))
