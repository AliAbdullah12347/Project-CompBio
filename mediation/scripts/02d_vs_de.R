#!/usr/bin/env Rscript
# ==============================================================================
# 02d_vs_de.R -- EXTERNAL VALIDATION against the differential-expression arm.
#
# The mediation total effect is, by construction, the same quantity the DE arm
# estimated: lithium's effect on whole-blood expression. The two were computed
# by different code, different models (21 parameters here against ~12 there) and
# on different gene sets. If the mediation machinery is sound they must agree.
#
# This is the strongest check available that does not reuse any of this arm's
# own algebra, because the DE arm was written and validated months earlier and
# knows nothing about mediation.
#
# Exact agreement is NOT expected on counts: the mediation outcome model spends
# 9 extra degrees of freedom on the mediators and their interactions, so it has
# slightly less power than the DE model on the same data.
# ==============================================================================

ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "mediation")); options(width = 160)
say <- function(...) cat(sprintf(...))

de <- read.csv(gzfile(file.path(ROOT, "de_final/results/DE_WB_LI.csv.gz")),
               stringsAsFactors = FALSE)
f  <- readRDS("results/primary/lineage_fit.rds"); tb <- f$tb
m  <- match(tb$gene, de$gene)
stopifnot(!any(is.na(m)))

de_sig  <- intersect(de$gene[de$BH < 0.05], tb$gene)
med_sig <- tb$gene[tb$TE_q_BH < 0.05]
ov <- intersect(de_sig, med_sig)

say("== 02d: mediation total effect vs the DE arm ==\n\n")
say("  DE arm, all 12,368 genes, lithium DEGs at BH 5%%   : %d\n", sum(de$BH < 0.05))
say("  of those, inside our 12,018 non-LM22 primary set   : %d\n", length(de_sig))
say("  mediation TOTAL-EFFECT significant                 : %d\n", length(med_sig))
say("  overlap                                            : %d\n", length(ov))
say("  Jaccard                                            : %.3f\n",
    length(ov) / length(union(de_sig, med_sig)))
say("  share of DE genes recovered                        : %.1f%%\n",
    100 * length(ov) / length(de_sig))
say("\n  correlation of DE log2FC with mediation TE         : %.4f\n", cor(de$logFC[m], tb$TE))
say("  sign agreement                                      : %.1f%%\n",
    100 * mean(sign(de$logFC[m]) == sign(tb$TE)))
say("\n  for contrast, the DE log2FC against the two PARTS of that total:\n")
say("    vs natural DIRECT effect   : %.4f\n", cor(de$logFC[m], tb$NDE))
say("    vs natural INDIRECT effect : %.4f\n", cor(de$logFC[m], tb$NIE))
say("\n  verdict: %s\n", if (cor(de$logFC[m], tb$TE) > 0.99 &&
                            length(ov) / length(de_sig) > 0.8)
  "the mediation machinery reproduces the DE arm. The decomposition is splitting\n           a total effect that was independently established, not inventing one." else
  "**the two arms disagree more than they should -- investigate**")
stopifnot(cor(de$logFC[m], tb$TE) > 0.99)
writeLines(sprintf("corr %.4f | overlap %d/%d | jaccard %.3f | PASS",
                   cor(de$logFC[m], tb$TE), length(ov), length(de_sig),
                   length(ov) / length(union(de_sig, med_sig))),
           "results/de_crosscheck.txt")
