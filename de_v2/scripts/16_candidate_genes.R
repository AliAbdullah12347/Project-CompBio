#!/usr/bin/env Rscript
# ==============================================================================
# 16_candidate_genes.R -- the genes worth writing about, with every piece of
# evidence attached to each one.
#
# WHY A COMPOSITE TABLE RATHER THAN A TOP-N LIST
#
# "The top 20 genes by p-value" is the wrong list for a discussion section. A
# p-value says a difference is unlikely under the null; it says nothing about
# whether the difference is large, whether it survives a stricter error rate,
# whether it depends on an arbitrary filter, or whether it is simply the cell
# mixture moving. Those are four separate questions and a gene can pass one and
# fail the others.
#
# So every gene gets all of it, joined into one row:
#
#   SIGNIFICANCE   BH, Storey, BY, Holm, Bonferroni at the baseline filter
#   EFFECT SIZE    log2FC, and whether it clears a formal TREAT threshold
#   STABILITY      how many of the 31 expression filters called it
#   COMPOSITION    Spearman rho with granulocyte fraction, and the log2FC after
#                  adjusting for cell composition -- the question this whole
#                  project exists to ask
#   LOCALISATION   which cell lineages it is called in, and how robustly
#
# TIERS. Genes are then sorted into evidence tiers, defined BEFORE looking at
# which genes land where. The tiers are about strength of evidence only; they
# are deliberately silent about biological interest, because a gene can be
# rock-solid and boring.
#
# THE TWO KINDS OF INTERESTING, AND WHY BOTH ARE REPORTED
#
# This project's thesis is that lithium's apparent effect on gene expression is
# largely the cell mixture shifting rather than cells changing behaviour. That
# makes two disjoint sets of genes worth discussing, for opposite reasons:
#
#   COMPOSITIONAL   strong effect, high correlation with granulocyte fraction,
#                   effect collapses once composition is adjusted for. These
#                   are the evidence FOR the thesis.
#   RESIDUAL        strong effect that SURVIVES composition adjustment. These
#                   are the counterexamples, and they are where any real
#                   within-cell lithium biology would have to live.
#
# A discussion that quotes only the first set is circular. Both are produced.
# ==============================================================================

suppressPackageStartupMessages({library(limma); library(edgeR)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
HERE <- file.path(ROOT, "de_v2")
setwd(HERE); set.seed(481)
OUT <- "results"; options(width = 190)

p <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
sym <- local({ m <- read.delim(file.path(ROOT, "de_analysis/data/ens2sym.tsv"),
                               header = FALSE, col.names = c("ens", "sym"))
               setNames(m$sym, m$ens) })

raw <- read.csv(gzfile(file.path(OUT, "tests_WB_LI_raw.csv.gz")), stringsAsFactors = FALSE)
ilr <- read.csv(gzfile(file.path(OUT, "tests_WB_LI_ilr.csv.gz")), stringsAsFactors = FALSE)
rob <- read.csv(file.path(OUT, "sweep/wb_gene_robustness.csv"), stringsAsFactors = FALSE)
ctr <- read.csv(file.path(OUT, "sweep/ct_gene_robustness.csv"), stringsAsFactors = FALSE)
cmp <- read.csv(file.path(ROOT, "de_analysis/results/composition_link_LI.csv"), stringsAsFactors = FALSE)
cat("== 16_candidate_genes ==\n")

## ---- TREAT: which genes clear a formal effect-size threshold ---------------
# wholeblood_treat.csv records only counts. The identities are needed here, so
# the test is re-run. TREAT is not "filter the DEG list by log2FC" -- that is a
# different and invalid procedure. It tests the null that the effect is BELOW
# the threshold, which is the question actually being asked.
design_for <- function(adjust) {
  s <- p$sel$LI
  mm <- p$meta[s$keep, , drop = FALSE]; mm$grp <- s$grp
  mm$plate <- droplevels(factor(mm$plate)); mm$sex <- droplevels(factor(mm$sex))
  cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]])) > 1, logical(1))]
  X <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse = " + "))), data = mm)
  if (adjust == "ilr") X <- cbind(X, p$ILR[mm$title, , drop = FALSE])
  list(X = X, samples = mm$title) }
d <- design_for("raw")
dge <- calcNormFactors(DGEList(p$counts[, d$samples, drop = FALSE]), method = "TMM")
fit <- eBayes(lmFit(voom(dge, d$X), d$X))
TR <- lapply(c(0.1, 0.2), function(L) {
  tt <- topTreat(treat(fit, lfc = L, trend = FALSE), coef = "grpcase", number = Inf, sort.by = "none")
  rownames(tt)[tt$adj.P.Val < 0.05] })
names(TR) <- c("treat_0.1", "treat_0.2")
cat(sprintf("TREAT: %d genes clear |log2FC| 0.1, %d clear 0.2\n", length(TR[[1]]), length(TR[[2]])))

## ---- assemble ---------------------------------------------------------------
G <- data.frame(gene = raw$gene, stringsAsFactors = FALSE)
G$symbol <- unname(ifelse(is.na(sym[G$gene]), G$gene, sym[G$gene]))
mi <- match(G$gene, raw$gene)
G$logFC      <- round(raw$logFC[mi], 4)
G$AveExpr    <- round(raw$AveExpr[mi], 3)
G$P          <- raw$P.Value[mi]
G$BH         <- raw$BH[mi]
G$storey     <- raw$q_storey_spline[mi]
G$BY         <- raw$BY[mi]
G$holm       <- raw$holm[mi]
G$bonferroni <- raw$bonferroni[mi]
G$rank_p     <- rank(raw$P.Value[mi], ties.method = "min")

mj <- match(G$gene, ilr$gene)
G$logFC_adj  <- round(ilr$logFC[mj], 4)
G$BH_adj     <- ilr$BH[mj]
G$survives_adjustment <- !is.na(G$BH_adj) & G$BH_adj < 0.05
# How much of the effect disappears when cell composition is accounted for.
G$pct_effect_lost <- round(100 * (1 - abs(G$logFC_adj) / abs(G$logFC)), 1)

mk <- match(G$gene, cmp$gene)
G$rho_granulocyte <- cmp$rho_gran[mk]
# Lithium raises granulocytes, so a purely compositional gene moves in the same
# direction as its granulocyte correlation.
G$direction_matches_composition <- sign(G$logFC) == sign(G$rho_granulocyte)

rb <- rob[rob$analysis == "LI_raw", ]
ml <- match(G$gene, rb$gene)
G$filters_sig    <- rb$n_filters_sig[ml]
G$filters_tested <- rb$n_filters_tested[ml]
G$best_rank      <- rb$best_rank[ml]
G$median_rank    <- rb$median_rank[ml]
G$all_filters    <- !is.na(ml) & rb$n_filters_sig[ml] == rb$n_filters_tested[ml]

G$treat_0.1 <- G$gene %in% TR[["treat_0.1"]]
G$treat_0.2 <- G$gene %in% TR[["treat_0.2"]]

ct <- ctr[ctr$analysis == "LI_raw", ]
agg <- do.call(rbind, lapply(split(ct, ct$gene), function(z) data.frame(
  gene = z$gene[1],
  lineages = paste(z$lineage[order(z$median_rank)], collapse = "/"),
  n_lineages = nrow(z),
  lineages_allsettings = sum(z$n_settings_sig == z$n_settings_tested),
  stringsAsFactors = FALSE)))
mm2 <- match(G$gene, agg$gene)
G$lineages <- ifelse(is.na(mm2), "", agg$lineages[mm2])
G$n_lineages <- ifelse(is.na(mm2), 0L, agg$n_lineages[mm2])
G$lineages_allsettings <- ifelse(is.na(mm2), 0L, agg$lineages_allsettings[mm2])

## ---- evidence tiers (defined before inspecting membership) ------------------
G$tier <- with(G, ifelse(
  holm < 0.05 & all_filters & treat_0.2 & n_lineages > 0, "1 - strongest",
  ifelse(holm < 0.05 & all_filters & treat_0.1, "2 - very strong",
  ifelse(BH < 0.05 & all_filters & treat_0.1, "3 - strong",
  ifelse(BH < 0.05 & all_filters, "4 - robust but small effect",
  ifelse(BH < 0.05, "5 - filter-dependent", "not significant"))))))

S <- G[G$BH < 0.05, ]
S <- S[order(S$rank_p), ]
cat("\n=== evidence tiers, whole-blood lithium ===\n")
print(as.data.frame(table(tier = S$tier)), row.names = FALSE)

cat("\n=== TIER 1: survive Holm, every filter, TREAT 0.2, and appear at cell-type level ===\n")
t1 <- S[S$tier == "1 - strongest", ]
print(t1[order(t1$rank_p), c("symbol", "logFC", "P", "holm", "filters_sig", "filters_tested",
                             "rho_granulocyte", "logFC_adj", "pct_effect_lost",
                             "survives_adjustment", "lineages")], row.names = FALSE)

## ---- the two discussion sets ------------------------------------------------
cat("\n=== COMPOSITIONAL: strongest genes whose effect collapses on adjustment ===\n")
cat("    (evidence FOR the thesis that lithium acts through the cell mixture)\n")
comp <- S[S$holm < 0.05 & S$all_filters & !is.na(S$rho_granulocyte) &
          abs(S$rho_granulocyte) > 0.3 & !S$survives_adjustment, ]
comp <- comp[order(-abs(comp$rho_granulocyte)), ]
print(head(comp[, c("symbol", "logFC", "rho_granulocyte", "logFC_adj", "pct_effect_lost",
                    "holm", "filters_sig", "filters_tested")], 15), row.names = FALSE)

cat("\n=== RESIDUAL: strongest genes that SURVIVE composition adjustment ===\n")
cat("    (counterexamples; where within-cell lithium biology would have to live)\n")
res <- S[S$holm < 0.05 & S$all_filters & S$survives_adjustment, ]
res <- res[order(res$BH_adj), ]
print(head(res[, c("symbol", "logFC", "logFC_adj", "pct_effect_lost", "BH_adj",
                   "rho_granulocyte", "lineages", "filters_sig", "filters_tested")], 20),
      row.names = FALSE)

cat(sprintf("\n  %d of %d Holm-significant, all-filter genes survive composition adjustment (%.1f%%)\n",
            nrow(res), sum(S$holm < 0.05 & S$all_filters),
            100 * nrow(res) / sum(S$holm < 0.05 & S$all_filters)))

## ---- bipolar, for completeness ----------------------------------------------
cat("\n=== BIPOLAR: every gene ever called, with why it fails ===\n")
braw <- read.csv(gzfile(file.path(OUT, "tests_WB_BPD_raw.csv.gz")), stringsAsFactors = FALSE)
brob <- rob[rob$analysis == "BPD_raw", ]
B <- braw[braw$BH < 0.05, ]
B$symbol <- unname(ifelse(is.na(sym[B$gene]), B$gene, sym[B$gene]))
mb <- match(B$gene, brob$gene)
B$filters_sig <- brob$n_filters_sig[mb]; B$filters_tested <- brob$n_filters_tested[mb]
B$holm_sig <- B$holm < 0.05; B$BY_sig <- B$BY < 0.05
print(B[order(B$P.Value), c("symbol", "logFC", "P.Value", "BH", "holm_sig", "BY_sig",
                            "filters_sig", "filters_tested")], row.names = FALSE)

write.csv(G[order(G$rank_p), ], gzfile(file.path(OUT, "candidate_genes_LI_all.csv.gz")), row.names = FALSE)
write.csv(S[, setdiff(names(S), "gene")], file.path(OUT, "candidate_genes_LI_significant.csv"), row.names = FALSE)
write.csv(t1, file.path(OUT, "candidate_genes_tier1.csv"), row.names = FALSE)
write.csv(comp, file.path(OUT, "candidate_genes_compositional.csv"), row.names = FALSE)
write.csv(res, file.path(OUT, "candidate_genes_residual.csv"), row.names = FALSE)
cat("\nwrote results/candidate_genes_*.csv\n")
