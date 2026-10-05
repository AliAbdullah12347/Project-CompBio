#!/usr/bin/env Rscript
# ==============================================================================
# baseline_relabel.R -- post-hoc reading of the one-subject sensitivity run.
#
# baseline_gap.R part B promotes one control at a time to a lithium-free case,
# emulating the single subject that Krebs count as a case and the deposit counts
# as a control, and refits the published specification each time. This script
# reads that output back and asks two things the loop itself did not record:
#
#   - how wide the spread is, and whether the published 976 sits inside it;
#   - whether the promoted subject's assessment group matters. It might: all
#     239 real cases are group A, so promoting a group-B control is the only
#     way a case can ever be group B, and that breaks the nesting of group
#     inside diagnosis that the rest of the design relies on. A group-A
#     promotion is the faithful emulation; a group-B one is a different
#     perturbation and should be read separately.
#
# No fitting. Writes only inside experimentation/runs/baseline/.
# ==============================================================================

setwd("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation")
source("scripts/common.R"); source("scripts/ref_krebs.R")

CACHE <- "C:/Users/hp/AppData/Local/Temp/claude/baseline_fits"
OUT   <- file.path(EXP_HOME, "runs", "baseline")
f <- file.path(OUT, "one_subject_sensitivity.csv")
if (!file.exists(f)) stop("run baseline_gap.R part B first")

rel  <- read.csv(f, colClasses = c(subject = "character"))
p    <- load_prep()
i    <- match(rel$subject, p$meta$title)
rel$group <- as.character(p$meta$group[i])
rel$sex   <- as.character(p$meta$sex[i])
rel$age   <- p$meta$age[i]
base <- n_deg(readRDS(file.path(CACHE, "spec_allcase2_444.rds"))$tt, 0.05)
rel$delta <- rel$DEG05 - base

cat(sprintf("== one-subject sensitivity of the published lithium specification ==\n"))
cat(sprintf("   unmodified fit: %d DEG at FDR .05   |   paper reports %d\n\n",
            base, KREBS_PAPER$li_ndeg))

summ <- function(d, lab) data.frame(subset = lab, n = nrow(d),
  min = min(d$DEG05), q25 = quantile(d$DEG05, .25), median = median(d$DEG05),
  q75 = quantile(d$DEG05, .75), max = max(d$DEG05), sd = round(sd(d$DEG05), 1),
  mean_delta = round(mean(d$delta), 1),
  n_below_base = sum(d$DEG05 < base), n_above_base = sum(d$DEG05 > base),
  covers_976 = min(d$DEG05) <= KREBS_PAPER$li_ndeg && max(d$DEG05) >= KREBS_PAPER$li_ndeg)
s <- summ(rel, "all promotions")
if (length(unique(rel$group)) > 1)
  s <- rbind(s, summ(rel[rel$group == "A", ], "group A only (faithful)"),
                summ(rel[rel$group == "B", ], "group B only (breaks nesting)"))
print(s, row.names = FALSE)

cat(sprintf("\n   range as a percentage of the unmodified count: %.1f%% to %.1f%%\n",
            100 * min(rel$DEG05) / base, 100 * max(rel$DEG05) / base))
cat(sprintf("   is the published 976 inside the observed range: %s\n",
            KREBS_PAPER$li_ndeg >= min(rel$DEG05) && KREBS_PAPER$li_ndeg <= max(rel$DEG05)))
cat(sprintf("   gap to explain (%d - %d) = %d, versus a one-subject sd of %.0f\n",
            base, KREBS_PAPER$li_ndeg, base - KREBS_PAPER$li_ndeg, sd(rel$DEG05)))

# ---------------------------------------------------------------------------
# why should one subject matter this much?
# ---------------------------------------------------------------------------
# A uniform rescaling of every p-value by 10% moves this fit's BH count by only
# ~96 genes, so a swing of several hundred is not a thresholding artifact: the
# test statistics themselves move. The candidate explanation is composition.
# The promoted subject joins an 87-person reference cell, and because cell
# mixture drives thousands of genes together, one atypical subject shifts the
# whole transcriptome coherently rather than gene by gene.
#
# Testable: the promoted subject's position on the myeloid/lymphoid balance
# should predict the direction of the change. A control with a lithium-like
# (myeloid-heavy) profile dilutes the contrast; a lymphoid-heavy one sharpens it.
ilr  <- p$ILR[rel$subject, , drop = FALSE]
lin  <- p$lineage[rel$subject, , drop = FALSE]
rel$b1_myeloid_lymphoid <- round(ilr[, "b1_myeloid_vs_lymphoid"], 4)
rel$neutrophil <- round(p$frac[rel$subject, "Neutrophils"], 4)
rel$gran <- round(lin[, "gran"], 4)

# reference point: the mean balance of the cell the subject is joining
refcell <- p$meta$title[p$meta$qc_pass & p$meta$dx != "Control" & p$meta$lithium == 0]
usecell <- p$meta$title[p$meta$qc_pass & p$meta$dx != "Control" & p$meta$lithium == 1]
mb_ref <- mean(p$ILR[refcell, "b1_myeloid_vs_lymphoid"])
mb_use <- mean(p$ILR[usecell, "b1_myeloid_vs_lymphoid"])
cat(sprintf("\n== composition test ==\n"))
cat(sprintf("   mean myeloid/lymphoid balance: non-using cases %.4f, lithium users %.4f\n",
            mb_ref, mb_use))
cat(sprintf("   (lithium users are the more myeloid group by %.4f)\n", mb_use - mb_ref))
ct <- cor.test(rel$b1_myeloid_lymphoid, rel$DEG05, method = "spearman", exact = FALSE)
cat(sprintf("   spearman(promoted subject's balance, DEG count) = %+.3f, p = %.3g, n = %d\n",
            ct$estimate, ct$p.value, nrow(rel)))
ct2 <- cor.test(rel$neutrophil, rel$DEG05, method = "spearman", exact = FALSE)
cat(sprintf("   spearman(promoted subject's neutrophil fraction, DEG count) = %+.3f, p = %.3g\n",
            ct2$estimate, ct2$p.value))
cat("   prediction: negative -- a myeloid-heavy (lithium-like) control dilutes the contrast.\n")

rel <- rel[order(rel$DEG05), ]
cat("\n   per-subject detail:\n")
print(rel[, c("subject", "group", "sex", "age", "n_genes", "DEG05", "delta",
              "b1_myeloid_lymphoid", "neutrophil")], row.names = FALSE)
comp <- data.frame(
  quantity = c("spearman balance vs DEG", "p", "spearman neutrophil vs DEG", "p",
               "mean balance non-using cases", "mean balance lithium users"),
  value = c(round(ct$estimate, 4), signif(ct$p.value, 3),
            round(ct2$estimate, 4), signif(ct2$p.value, 3),
            round(mb_ref, 4), round(mb_use, 4)))
write.csv(comp, file.path(OUT, "one_subject_composition.csv"), row.names = FALSE)
write.csv(s,   file.path(OUT, "one_subject_summary.csv"), row.names = FALSE)
write.csv(rel, file.path(OUT, "one_subject_sensitivity.csv"), row.names = FALSE)
cat("\n[baseline_relabel] done\n")
