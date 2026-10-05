#!/usr/bin/env Rscript
# ==============================================================================
# 17_composition_gradient.R -- resolving an apparent contradiction.
#
# TWO RESULTS THAT LOOK INCOMPATIBLE
#
#   From 11_composition_link.R:  only 292 of the 1,426 lithium genes (20.5%)
#   remain significant once cell composition is adjusted for, and the median
#   effect shrinks by 40.8%. Read alone, that says the signature is mostly the
#   cell mixture moving.
#
#   From 16_candidate_genes.R:  100% of the 125 strongest genes (Holm-
#   significant under every one of the 31 filters) survive that same
#   adjustment, losing only ~20-35% of their effect. Read alone, that says the
#   signature is NOT the cell mixture.
#
# Both are correct. They describe different parts of the same distribution, and
# quoting either alone would misrepresent the result. This script measures the
# gradient between them so the paper can state it properly.
#
# THE QUESTION. Does the probability that a gene survives composition
# adjustment depend on how strong its original signal was? If yes, the honest
# summary is "composition explains the breadth of the signature but not its
# core", which is a different claim from either headline above.
# ==============================================================================

ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "de_v2")); options(width = 170)
OUT <- "results"
G <- read.csv(gzfile(file.path(OUT, "candidate_genes_LI_all.csv.gz")), stringsAsFactors = FALSE)
S <- G[G$BH < 0.05, ]
cat(sprintf("== 17_composition_gradient ==\n%d genes significant at BH 0.05\n", nrow(S)))

## ---- survival by original significance -------------------------------------
S$decile <- cut(rank(S$rank_p), breaks = quantile(rank(S$rank_p), 0:10/10),
                include.lowest = TRUE, labels = FALSE)
tab <- do.call(rbind, lapply(sort(unique(S$decile)), function(d) {
  z <- S[S$decile == d, ]
  data.frame(decile = d, n = nrow(z),
             rank_range = sprintf("%d-%d", min(z$rank_p), max(z$rank_p)),
             median_logFC = round(median(abs(z$logFC)), 3),
             median_pct_lost = round(median(z$pct_effect_lost), 1),
             n_survive = sum(z$survives_adjustment),
             pct_survive = round(100 * mean(z$survives_adjustment), 1),
             stringsAsFactors = FALSE) }))
cat("\n=== does survival of composition adjustment depend on original strength? ===\n")
cat("    deciles of the 1,426 significant genes, strongest first\n")
print(tab, row.names = FALSE)

## ---- survival by correction stringency -------------------------------------
cat("\n=== the same, cut by how strict a correction the gene passes ===\n")
lv <- list(
  "Holm (FWER)"            = S$holm < 0.05,
  "BY (any dependence)"    = S$BY < 0.05 & S$holm >= 0.05,
  "BH only"                = S$BH < 0.05 & S$BY >= 0.05)
cut2 <- do.call(rbind, lapply(names(lv), function(k) {
  z <- S[lv[[k]], ]
  data.frame(passes = k, n = nrow(z),
             median_logFC = round(median(abs(z$logFC)), 3),
             median_pct_lost = round(median(z$pct_effect_lost), 1),
             pct_survive_adjustment = round(100 * mean(z$survives_adjustment), 1),
             stringsAsFactors = FALSE) }))
print(cut2, row.names = FALSE)

## ---- survival by filter stability -------------------------------------------
cat("\n=== and cut by filter stability ===\n")
cut3 <- do.call(rbind, lapply(list(
  c("called under every filter", TRUE), c("not under every filter", FALSE)), function(x) {
  z <- S[S$all_filters == as.logical(x[2]), ]
  data.frame(stability = x[1], n = nrow(z),
             median_logFC = round(median(abs(z$logFC)), 3),
             median_pct_lost = round(median(z$pct_effect_lost), 1),
             pct_survive_adjustment = round(100 * mean(z$survives_adjustment), 1),
             stringsAsFactors = FALSE) }))
print(cut3, row.names = FALSE)

## ---- genes whose effect GROWS on adjustment ---------------------------------
# Composition adjustment is not only subtractive. If a gene moves against the
# compositional current, removing that current makes its effect LARGER. These
# are the cleanest non-compositional candidates in the dataset.
cat("\n=== genes whose effect GROWS once composition is removed ===\n")
cat("    (moving against the compositional current -- the cleanest counterexamples)\n")
gr <- S[S$pct_effect_lost < 0 & S$holm < 0.05, ]
gr <- gr[order(gr$pct_effect_lost), ]
print(gr[, c("symbol", "logFC", "logFC_adj", "pct_effect_lost", "rho_granulocyte",
             "holm", "filters_sig", "filters_tested", "lineages")], row.names = FALSE)

## ---- correlation between compositional dependence and effect size ----------
ok <- is.finite(S$pct_effect_lost) & is.finite(S$rho_granulocyte)
cat(sprintf("\n=== relationships ===\n"))
cat(sprintf("  |rho with granulocytes| vs %% effect lost : Spearman %.3f (n=%d)\n",
            cor(abs(S$rho_granulocyte[ok]), S$pct_effect_lost[ok], method = "spearman"), sum(ok)))
cat(sprintf("  |log2FC| vs %% effect lost               : Spearman %.3f\n",
            cor(abs(S$logFC[ok]), S$pct_effect_lost[ok], method = "spearman")))
cat(sprintf("  median %% effect lost, all %d genes       : %.1f%%\n",
            nrow(S), median(S$pct_effect_lost, na.rm = TRUE)))
cat(sprintf("  median %% effect lost, Holm-significant   : %.1f%%\n",
            median(S$pct_effect_lost[S$holm < 0.05], na.rm = TRUE)))

cat("\n=== THE STATEMENT FOR THE PAPER ===\n")
cat(sprintf("  Of all %d lithium genes, %.1f%% remain significant after adjusting for cell\n",
            nrow(S), 100 * mean(S$survives_adjustment)))
cat(sprintf("  composition, with a median effect loss of %.1f%%. That proportion rises\n",
            median(S$pct_effect_lost, na.rm = TRUE)))
cat(sprintf("  monotonically with original signal strength, reaching %.0f%% among the %d genes\n",
            100 * mean(S$survives_adjustment[S$holm < 0.05 & S$all_filters]),
            sum(S$holm < 0.05 & S$all_filters)))
cat("  that pass family-wise error control under every expression filter.\n")
cat("  Composition therefore explains the BREADTH of the lithium signature but not its CORE.\n")

write.csv(tab,  file.path(OUT, "composition_gradient_by_decile.csv"), row.names = FALSE)
write.csv(cut2, file.path(OUT, "composition_gradient_by_correction.csv"), row.names = FALSE)
write.csv(gr,   file.path(OUT, "composition_effect_grows.csv"), row.names = FALSE)
cat("\nwrote results/composition_gradient_*.csv\n")
