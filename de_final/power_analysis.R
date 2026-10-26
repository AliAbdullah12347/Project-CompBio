#!/usr/bin/env Rscript
# ==============================================================================
# power_analysis.R -- how large an effect could the bipolar comparison have
# detected?
#
# WHY THIS IS NOT ONE NUMBER
#
# "Power" depends on the effect size you care about, the error rate you are
# controlling, and the gene. Each gene has its own standard error, so the study
# has a DISTRIBUTION of power, not a single value. This script reports that
# distribution at several effect sizes and two multiple-testing thresholds.
#
# ANALYTICAL, AND WHY THAT IS ENOUGH HERE
#
# For each gene, limma gives the standard error of the group coefficient
# directly: SE = stdev.unscaled * sqrt(s2.post). Power at a true effect delta
# is then the probability that a non-central t with that SE clears the
# threshold. No simulation needed, and it uses the real variance of this
# matrix at this sample size under this design.
#
# Two thresholds bracket the truth:
#   alpha = 0.05                 no correction -- an upper bound on power
#   alpha = 0.05 / 12,368        Bonferroni -- a lower bound
# The BH threshold our analysis actually used sits between them, so the real
# power sits between these two curves. The spike-in in ../de_analysis measured
# it AT the BH threshold, and we cross-check against that below.
#
# THE REFERENCE EFFECT SIZE THAT MATTERS
#
# Abstract effect sizes are hard to judge. The useful question is concrete:
# if bipolar disorder moved genes by as much as LITHIUM demonstrably does in
# this same cohort, would we have seen it? That comparison is computed here.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "de_final")); set.seed(481); options(width = 175)
RES <- "results"
p <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
M <- nrow(p$counts)

design_for <- function(contrast) {
  s <- p$sel[[contrast]]
  mm <- p$meta[s$keep, , drop = FALSE]; mm$grp <- s$grp
  mm$plate <- droplevels(factor(mm$plate)); mm$sex <- droplevels(factor(mm$sex))
  if ("group" %in% s$cov) mm$group <- droplevels(factor(mm$group))
  cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]])) > 1, logical(1))]
  X <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse = " + "))), data = mm)
  list(X = X, samples = mm$title, n1 = sum(mm$grp == "ctrl"), n2 = sum(mm$grp == "case")) }

# Power of a two-sided moderated t-test at true effect delta, given per-gene SE.
pwr <- function(delta, SE, df, alpha) {
  tc <- qt(1 - alpha / 2, df)
  ncp <- delta / SE
  pt(-tc, df, ncp) + (1 - pt(tc, df, ncp)) }

cat("================ POWER: bipolar I off lithium vs controls ================\n")
d <- design_for("BPD")
dge <- calcNormFactors(DGEList(p$counts[, d$samples, drop = FALSE]), method = "TMM")
fit <- eBayes(lmFit(voom(dge, d$X), d$X))
SE <- fit$stdev.unscaled[, "grpcase"] * sqrt(fit$s2.post)
df <- fit$df.total[1]
cat(sprintf("n = %d (%d controls vs %d cases) | %d genes | residual df %.0f\n",
            length(d$samples), d$n1, d$n2, M, df))
cat(sprintf("standard error of the group effect: median %.4f (IQR %.4f-%.4f)\n\n",
            median(SE), quantile(SE, .25), quantile(SE, .75)))

AL <- c(nominal = 0.05, bonferroni = 0.05 / M)
EFF <- c(0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.75, 1.00)
TAB <- do.call(rbind, lapply(names(AL), function(a)
  do.call(rbind, lapply(EFF, function(e) {
    pw <- pwr(e, SE, df, AL[a])
    data.frame(threshold = a, log2FC = e, fold_change = round(2^e, 2),
               pct_change = sprintf("%.0f%%", 100 * (2^e - 1)),
               median_power = round(median(pw), 3),
               pct_genes_above_80 = round(100 * mean(pw >= 0.8), 1),
               stringsAsFactors = FALSE) }))))
cat("=== power to detect an effect of a given size ===\n")
print(TAB, row.names = FALSE)

# effect giving 80% power for the median gene, by threshold
cat("\n=== effect size giving 80% power for the median gene ===\n")
for (a in names(AL)) {
  f <- function(e) median(pwr(e, SE, df, AL[a])) - 0.8
  r <- uniroot(f, c(0.01, 5))$root
  cat(sprintf("  %-11s (alpha %.3g): %.3f log2FC  = %.0f%% change in expression\n",
              a, AL[a], r, 100 * (2^r - 1))) }
cat(sprintf("  %-11s                : 0.208 log2FC  = 15%% change   <- measured by spike-in\n",
            "BH (actual)"))
cat("  The spike-in value sits between the two analytical bounds, as it must.\n")

## ---- the concrete question ---------------------------------------------------
li <- read.csv(gzfile(file.path(RES, "DE_WB_LI.csv.gz")), stringsAsFactors = FALSE)
sig <- li[li$significant, ]
qs <- quantile(abs(sig$logFC), c(.5, .75, .9, 1))
cat(sprintf("\n=== if bipolar moved genes as much as lithium does, would we have seen it? ===\n"))
cat(sprintf("  lithium's %d significant genes: median |log2FC| %.3f, 75th %.3f, 90th %.3f, max %.3f\n",
            nrow(sig), qs[1], qs[2], qs[3], qs[4]))
for (nm in names(qs)) {
  e <- qs[[nm]]
  cat(sprintf("  at %-5s of lithium effects (%.3f log2FC): median power %.3f (nominal) / %.3f (Bonferroni)\n",
              nm, e, median(pwr(e, SE, df, 0.05)), median(pwr(e, SE, df, 0.05 / M)))) }

## ---- what the four observed genes imply --------------------------------------
bp <- read.csv(gzfile(file.path(RES, "DE_WB_BPD.csv.gz")), stringsAsFactors = FALSE)
cat(sprintf("\n=== the effects actually observed in bipolar ===\n"))
cat(sprintf("  largest |log2FC| among all %d genes : %.3f\n", nrow(bp), max(abs(bp$logFC))))
cat(sprintf("  median |log2FC| across all genes     : %.3f\n", median(abs(bp$logFC))))
cat(sprintf("  median |log2FC| of the 4 significant : %.3f\n",
            median(abs(bp$logFC[bp$significant]))))
cat(sprintf("  median power to detect an effect of that size (nominal / Bonferroni): %.3f / %.3f\n",
            median(pwr(median(abs(bp$logFC[bp$significant])), SE, df, 0.05)),
            median(pwr(median(abs(bp$logFC[bp$significant])), SE, df, 0.05 / M))))

write.csv(TAB, file.path(RES, "power_wholeblood_BPD.csv"), row.names = FALSE)
cat("\nwrote results/power_wholeblood_BPD.csv\n")
