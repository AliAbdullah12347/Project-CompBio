#!/usr/bin/env Rscript
# ==============================================================================
# 02_hierarchical.R -- the multiplicity that a per-contrast correction misses.
#
# THE PROBLEM
#
# 01 corrects within each analysis. That is the usual practice and it is not
# enough here, because we did not run one analysis -- we ran four, and the
# cell-type ones are themselves five analyses each. Two separate inflations
# follow, and they need different machinery.
#
#   ACROSS LINEAGES. A cell-type contrast tests every gene five times, once per
#   lineage. Correcting inside each lineage and then reporting "any lineage"
#   gives the cell-type level five independent chances at significance, so it
#   would beat whole blood on multiplicity alone. Three treatments are reported:
#     (a) per-lineage BH     -- answers "is there signal in granulocytes?"
#     (b) pooled BH          -- answers "is there signal anywhere?", and is the
#                               only fair basis for comparing against whole blood
#     (c) Benjamini-Bogomolov -- answers (a) but pays for having looked at all
#                               five and reported only the interesting ones
#
#   ACROSS CONTRASTS. Four comparisons were run and the interesting ones get
#   written up. That selection is itself a multiple-testing problem, and nothing
#   in 01 accounts for it. Two treatments:
#     (d) global BH over all 148,416 tests -- the blunt instrument
#     (e) Benjamini-Bogomolov over the four families -- the right instrument
#
# WHY BENJAMINI-BOGOMOLOV RATHER THAN JUST POOLING EVERYTHING
#
# Pooling all 148,416 tests into one BH treats a granulocyte test and a
# whole-blood test as exchangeable members of one family. They are not: they
# have different power, different null behaviour, and answer different
# questions. Pooling lets a contrast with thousands of strong signals (WB_LI)
# effectively lend significance to a contrast with none (WB_BPD), because BH's
# threshold depends on the whole sorted vector. BB keeps the families separate
# and charges an explicit, interpretable price -- the fraction of families
# selected -- for having looked at all of them.
#
# Both are reported, because the difference between them is itself informative.
#
# NOTE ON DEPENDENCE. The four contrasts share subjects: the same 74 off-lithium
# bipolar I patients appear in WB_BPD and CT_BPD, and the same 226 bipolar I in
# WB_LI and CT_LI. The cell-type and whole-blood versions of one contrast are
# near-duplicates of each other. So these families are strongly positively
# dependent, which is the regime BH and BB are valid in, but it also means the
# four are nothing like four independent studies and must never be counted as
# such.
# ==============================================================================

ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
HERE <- file.path(ROOT, "de_v2")
setwd(HERE); set.seed(481)
source("scripts/mtc.R")
RES <- "results"
ALL <- readRDS(file.path(RES, "adjusted.rds"))
LINS <- c("gran", "mono", "T", "NK", "B")
Q <- 0.05
cat("== 02_hierarchical ==\n")

## ---- (a,b,c) within cell-type contrasts -------------------------------------
cat("\n--- across the five lineages, within each cell-type contrast ---\n")
CT <- list(); CTG <- list()
for (k in grep("^CT_", names(ALL), value = TRUE)) {
  f <- ALL[[k]]
  per <- do.call(rbind, lapply(LINS, function(ct) {
    z <- f[f$lineage == ct, ]
    data.frame(analysis = k, lineage = ct, n_tests = nrow(z),
               per_lineage_BH_05 = sum(p.adjust(z$P.Value, "BH") < Q),
               stringsAsFactors = FALSE) }))
  # pooled BH over all gene x lineage tests, collapsed to unique genes so the
  # denominator matches whole blood's
  qp <- p.adjust(f$P.Value, "BH")
  genes_hit <- unique(f$gene[qp < Q])
  ng <- length(unique(f$gene))
  # Benjamini-Bogomolov with lineages as families
  plist <- lapply(LINS, function(ct) f$P.Value[f$lineage == ct]); names(plist) <- LINS
  bb <- bb_hierarchical(plist, q = Q)
  per$BB_selected <- bb$selected
  per$BB_family_p <- signif(bb$family_p, 3)
  per$BB_05 <- vapply(LINS, function(ct)
    if (!bb$selected[ct]) 0L else sum(bb$adjusted[[ct]] < Q, na.rm = TRUE), integer(1))
  CT[[k]] <- per
  CTG[[k]] <- data.frame(analysis = k, n_tests = nrow(f), n_genes = ng,
                         pooled_BH_tests_05 = sum(qp < Q),
                         pooled_BH_unique_genes_05 = length(genes_hit),
                         BB_families_selected = bb$n_selected, BB_level2 = signif(bb$level2, 3),
                         BB_total_05 = sum(per$BB_05), stringsAsFactors = FALSE)
  cat(sprintf("  %-14s per-lineage BH: %s | pooled BH %d tests / %d genes | BB selected %d/5 -> %d\n",
              k, paste(sprintf("%s=%d", per$lineage, per$per_lineage_BH_05), collapse = " "),
              sum(qp < Q), length(genes_hit), bb$n_selected, sum(per$BB_05)))
}
PER <- do.call(rbind, CT); POOL <- do.call(rbind, CTG)
write.csv(PER,  file.path(RES, "celltype_per_lineage.csv"), row.names = FALSE)
write.csv(POOL, file.path(RES, "celltype_pooled.csv"), row.names = FALSE)

## ---- (d,e) across the four comparisons --------------------------------------
# Only the unadjusted ("raw") fits enter here. The ILR-adjusted fits are the
# same hypotheses asked a second way, not four more questions, so including
# them would double-count.
cat("\n--- across the four comparisons ---\n")
FAM <- c("WB_LI_raw", "WB_BPD_raw", "CT_LI_raw", "CT_BPD_raw")
plist <- lapply(FAM, function(k) ALL[[k]]$P.Value); names(plist) <- FAM
ntest <- vapply(plist, length, integer(1))
cat(sprintf("  tests per family: %s  (total %d)\n",
            paste(sprintf("%s=%d", FAM, ntest), collapse = " "), sum(ntest)))

# (d) global BH over everything
gp <- unlist(plist, use.names = FALSE)
gq <- p.adjust(gp, "BH")
idx <- rep(FAM, ntest)
glob <- vapply(FAM, function(k) sum(gq[idx == k] < Q), integer(1))

# (e) Benjamini-Bogomolov with the four contrasts as families
bb <- bb_hierarchical(plist, q = Q)
bbn <- vapply(FAM, function(k)
  if (!bb$selected[k]) 0L else sum(bb$adjusted[[k]] < Q, na.rm = TRUE), integer(1))

GL <- data.frame(
  family = FAM, n_tests = ntest,
  within_family_BH_05 = vapply(FAM, function(k) sum(p.adjust(plist[[k]], "BH") < Q), integer(1)),
  global_BH_05 = glob,
  BB_family_p = signif(bb$family_p, 3),
  BB_family_adj = signif(bb$family_adj, 3),
  BB_selected = bb$selected,
  BB_05 = bbn, stringsAsFactors = FALSE)
print(GL, row.names = FALSE)
cat(sprintf("\n  BB selected %d of %d families, so the within-family level becomes %.4f\n",
            bb$n_selected, bb$K, bb$level2))
cat("  (a family not selected reports nothing, which is the price of having looked at four)\n")
write.csv(GL, file.path(RES, "across_contrast_correction.csv"), row.names = FALSE)

cat("\nwrote results/celltype_*.csv and across_contrast_correction.csv\n")
