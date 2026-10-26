#!/usr/bin/env Rscript
# ==============================================================================
# 14_gene_lists.R -- every gene the mediation analyses flag, by name, in clearly
# defined categories. Feeds the master document; also written as CSV so the
# lists are machine-readable.
#
# THE CATEGORIES, AND WHAT EACH ONE MEANS
#
# The model is fitted gene by gene, so each gene gets its own total effect (TE),
# direct effect (NDE) and indirect effect through the cell mix (NIE). A gene is
# placed by asking, in order:
#
#   1  TOTAL       TE significant (BH 5% across the gene family). Lithium changes
#                  this gene's whole-blood expression, by whatever route. This is
#                  the solid list: it reproduces the DE arm (r = 0.9985) and beats
#                  its permutation null (p = 0.001).
#   2  DIRECT      TOTAL, and NDE also significant. The change survives
#                  accounting for the cell mix: the candidates for lithium acting
#                  INSIDE cells. Caveat: "direct" only means "not explained by
#                  this mediator" -- shifts between sub-types WITHIN a lineage
#                  are invisible to it and would show up here.
#   3  VIA_MIX     TOTAL, NDE not significant, and the bootstrap 95% interval of
#                  NIE excludes zero. The change is consistent with running
#                  through the cell mix and no direct part can be detected.
#                  NOT proof of a purely compositional effect: a non-significant
#                  NDE is absence of evidence, and the NIE interval is per gene,
#                  uncorrected for multiplicity.
#   4  UNRESOLVED  TOTAL, but neither of the above: lithium changes the gene and
#                  this study cannot say through which route.
#   5  OPPOSING    NDE significant but TE NOT: a direct effect that the indirect
#                  effect cancels ("inconsistent mediation", MacKinnon et al.
#                  2000). Lithium's net effect is hidden, not absent.
#   +  TOP_NIE     the 50 TOTAL genes with the largest |NIE| -- the genes whose
#                  lithium response is most driven by the cell mix
#   +  INTERACTION lithium x cell-mix interaction significant
#
# Per-gene significance of NIE on its own is NOT a category: the arm showed that
# count is unstable (0 on one tree, 3,647 on the other, and a permutation null
# that reaches 7,150). The proposal tree's NIE-significant genes are listed
# separately, as a consensus of three methods, with that warning attached.
# ==============================================================================

ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "mediation")); options(width = 160)
say <- function(...) cat(sprintf(...))
dir.create("results/genes", showWarnings = FALSE, recursive = TRUE)
e2s <- read.delim(file.path(ROOT, "de_analysis/data/ens2sym.tsv"), header = FALSE,
                  col.names = c("ens", "sym"), stringsAsFactors = FALSE)
sym <- setNames(e2s$sym, e2s$ens)

categorise <- function(vn) {
  tb <- readRDS(sprintf("results/primary/%s_fit.rds", vn))$tb
  bc <- read.csv(sprintf("results/boot/%s_boot_ci.csv", vn), stringsAsFactors = FALSE)
  pg <- read.csv(sprintf("results/perm/%s_perm_gene.csv", vn), stringsAsFactors = FALSE)
  stopifnot(identical(bc$gene, tb$gene), identical(pg$gene, tb$gene))
  g <- data.frame(gene = tb$gene, symbol = sym[tb$gene],
    TE = tb$TE, TE_q = tb$TE_q_BH, NDE = tb$NDE, NDE_q = tb$NDE_q_BH,
    NIE = tb$NIE, NIE_q_delta = tb$NIE_q_BH, NIE_lo_boot = bc$NIE_lo, NIE_hi_boot = bc$NIE_hi,
    NIE_q_boot = bc$NIE_q_boot_BH, NIE_q_perm = pg$emp_q,
    PM = ifelse(tb$TE_q_BH < 0.05, tb$PM, NA_real_), INT_q = tb$INT_q_BH,
    stringsAsFactors = FALSE)
  te  <- g$TE_q < 0.05; nde <- g$NDE_q < 0.05
  nie_ci <- g$NIE_lo_boot > 0 | g$NIE_hi_boot < 0
  g$category <- ifelse(te & nde, "DIRECT",
                ifelse(te & !nde & nie_ci, "VIA_MIX",
                ifelse(te, "UNRESOLVED",
                ifelse(nde, "OPPOSING", ""))))
  g$top_NIE <- FALSE
  o <- order(-abs(g$NIE) * te); g$top_NIE[o[seq_len(min(50, sum(te)))]] <- TRUE
  g$interaction <- g$INT_q < 0.05
  g$NIE_consensus <- g$NIE_q_delta < 0.05 & g$NIE_q_boot < 0.05 & g$NIE_q_perm < 0.05
  g$tree <- vn
  g
}
L <- categorise("lineage"); P <- categorise("proposal")
stopifnot(identical(L$gene, P$gene))

cnt <- function(g) c(TOTAL = sum(g$TE_q < 0.05), DIRECT = sum(g$category == "DIRECT"),
                     VIA_MIX = sum(g$category == "VIA_MIX"), UNRESOLVED = sum(g$category == "UNRESOLVED"),
                     OPPOSING = sum(g$category == "OPPOSING"), INTERACTION = sum(g$interaction),
                     NIE_delta = sum(g$NIE_q_delta < 0.05), NIE_boot = sum(g$NIE_q_boot < 0.05),
                     NIE_perm = sum(g$NIE_q_perm < 0.05), NIE_consensus = sum(g$NIE_consensus))
CT <- data.frame(category = names(cnt(L)), lineage = cnt(L), proposal = cnt(P),
                 both_trees = c(sum(L$TE_q < .05 & P$TE_q < .05),
                                sum(L$category == "DIRECT" & P$category == "DIRECT"),
                                sum(L$category == "VIA_MIX" & P$category == "VIA_MIX"),
                                sum(L$category == "UNRESOLVED" & P$category == "UNRESOLVED"),
                                sum(L$category == "OPPOSING" & P$category == "OPPOSING"),
                                sum(L$interaction & P$interaction),
                                sum(L$NIE_q_delta < .05 & P$NIE_q_delta < .05),
                                sum(L$NIE_q_boot < .05 & P$NIE_q_boot < .05),
                                sum(L$NIE_q_perm < .05 & P$NIE_q_perm < .05),
                                sum(L$NIE_consensus & P$NIE_consensus)),
                 stringsAsFactors = FALSE)
print(CT, row.names = FALSE)

## TOTAL genes split by direction, so the lists can be read biologically
for (g in list(L, P)) {
  te <- g$TE_q < 0.05
  say("  %s: TOTAL genes up %d / down %d | DIRECT up %d / down %d\n", g$tree[1],
      sum(te & g$TE > 0), sum(te & g$TE < 0),
      sum(g$category == "DIRECT" & g$TE > 0), sum(g$category == "DIRECT" & g$TE < 0))
}
## does the NIE share the sign of the total effect? (mediation in the same direction)
te <- L$TE_q < 0.05
say("  lineage TOTAL genes whose NIE has the same sign as TE: %.1f%%\n",
    100 * mean(sign(L$NIE[te]) == sign(L$TE[te])))

## LM22 marker genes, from the circularity common fit (lineage)
cf <- readRDS("results/circularity/lineage_circularity.rds")$tb
mk <- cf[cf$set == "marker", ]
mk$TE_q_own <- p.adjust(mk$TE_p, "BH"); mk$NDE_q_own <- p.adjust(mk$NDE_p, "BH")
acc <- read.csv("results/circularity/lm22_accounting.csv", stringsAsFactors = FALSE)
mk$top_celltype <- acc$top_celltype[match(toupper(mk$symbol), toupper(acc$symbol))]
MK <- data.frame(gene = mk$gene, symbol = mk$symbol, lm22_top_celltype = mk$top_celltype,
                 TE = mk$TE, TE_q = mk$TE_q_own, NDE = mk$NDE, NDE_q = mk$NDE_q_own,
                 NIE = mk$NIE, PM = ifelse(mk$TE_q_own < 0.05, mk$PM, NA_real_),
                 partialR2_cellmix = mk$partialR2_cellmix, stringsAsFactors = FALSE)
MK <- MK[order(MK$TE_q), ]
say("  LM22 markers with significant TE (own 350-gene family): %d\n", sum(MK$TE_q < 0.05))

ord <- function(g) g[order(g$TE_q, -abs(g$TE)), ]
write.csv(ord(L), "results/genes/genes_lineage_all.csv", row.names = FALSE)
write.csv(ord(P), "results/genes/genes_proposal_all.csv", row.names = FALSE)
write.csv(CT, "results/genes/category_counts.csv", row.names = FALSE)
write.csv(MK, "results/genes/lm22_markers_lineage.csv", row.names = FALSE)
saveRDS(list(lineage = L, proposal = P, counts = CT, markers = MK), "results/genes/gene_lists.rds")
cat("\n14_gene_lists complete.\n")
