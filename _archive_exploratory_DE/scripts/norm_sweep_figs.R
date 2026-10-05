#!/usr/bin/env Rscript
# ==============================================================================
# norm_sweep_figs.R -- six panels carrying the argument of the sweep.
#
# Base graphics on purpose: no extra dependency, every element explicit.
# Figures go to results/, never alongside code.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})
source("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation/scripts/common.R")

OUTDIR <- file.path(EXP_HOME, "runs", "norm_sweep")
RES    <- file.path(EXP_HOME, "results"); dir.create(RES, showWarnings = FALSE)
p   <- load_prep()
OBJ <- readRDS(file.path(OUTDIR, "norm_sweep_objects.rds"))
NF  <- OBJ$NF
ps  <- read.csv(file.path(OUTDIR, "per_sample_normfactors.csv"), check.names = FALSE,
                colClasses = c(title = "character"))
lith <- ps$lithium
s    <- read.csv(file.path(OUTDIR, "summary.csv"))

png(file.path(RES, "norm_sweep.png"), width = 2000, height = 1300, res = 145)
par(mfrow = c(2, 3), mar = c(4.6, 4.4, 3.4, 1.2), cex.main = 1.0)

## A. the confound ------------------------------------------------------------
y <- log2(NF[, "TMM"])
boxplot(y ~ lith, names = c("non-user (74)", "lithium (152)"),
        col = c("grey88", "#b8d4e8"), outline = FALSE, xlab = "",
        ylab = "log2 TMM normalisation factor",
        main = "A. The scaling factor tracks the exposure")
set.seed(1); points(jitter(lith + 1, .45), y, pch = 16, cex = .45,
                    col = adjustcolor("grey25", .5))
segments(c(.7, 1.7), c(mean(y[lith == 0]), mean(y[lith == 1])), c(1.3, 2.3),
         col = "firebrick", lwd = 2.5)
mtext(sprintf("diff %.4f log2 units | partial r -0.244 | p = 3.1e-4 adjusted",
              mean(y[lith == 1]) - mean(y[lith == 0])), 3, .1, cex = .62)

## B. the mechanism -----------------------------------------------------------
plot(ps$gran, y, pch = 16, cex = .6,
     col = ifelse(lith == 1, adjustcolor("#1f6fb4", .7), adjustcolor("#c0392b", .7)),
     xlab = "granulocyte lineage fraction (CIBERSORTx)",
     ylab = "log2 TMM normalisation factor",
     main = "B. ...because granulocytes dominate the library")
abline(lm(y ~ ps$gran), lwd = 2)
legend("topright", c("lithium user", "non-user"), pch = 16, bty = "n", cex = .8,
       col = c("#1f6fb4", "#c0392b"))
mtext(sprintf("r = %.3f | top-50-gene share alone explains R2 = %.2f of the factor",
              cor(y, ps$gran), summary(lm(y ~ ps$top50))$r.squared), 3, .1, cex = .62)

## C. dropping normalisation adds one constant --------------------------------
a <- OBJ$fits_tt$voom_TMM;      a <- a[order(a$gene), ]
b <- OBJ$fits_tt$voom_none_CPM; b <- b[order(b$gene), ]
shift <- b$logFC - a$logFC
hist(shift, breaks = 80, col = "grey80", border = NA,
     xlab = "logFC(no normalisation) - logFC(TMM), per gene",
     main = "C. Dropping normalisation adds one constant")
abline(v = median(shift), col = "firebrick", lwd = 2)
mtext(sprintf("median %.4f, sd %.4f over %s genes -- a location shift, not new biology",
              median(shift), sd(shift), format(nrow(a), big.mark = ",")), 3, .1, cex = .62)

## D. the consequence for the gene list ---------------------------------------
sd_ <- s[order(s$deg_fdr05), ]
bp <- barplot(sd_$deg_fdr05, names.arg = NA, border = NA, ylim = c(0, 4700),
              col = ifelse(sd_$norm == "none", "#c0392b", "grey70"),
              ylab = "DEG at FDR 0.05", main = "D. ...and nearly triples the gene list")
text(bp, sd_$deg_fdr05 + 160, sd_$deg_fdr05, cex = .72)
text(bp, -190, sd_$config, srt = 35, adj = 1, xpd = NA, cex = .62)
abline(h = s$deg_fdr05[s$config == "voom_TMM"], lty = 3)
mtext("Jaccard of either red list against voom+TMM: 0.21", 3, .1, cex = .62)

## E. one scalar explains every method's gene list ----------------------------
plot(s$median_lfc_shift, log2(s$deg_up / s$deg_dn), pch = 21, bg = "#1f6fb4",
     cex = 1.5, xlab = "global logFC offset vs TMM (log2 units)",
     ylab = "log2( up-genes / down-genes )",
     main = "E. The offset determines the whole gene list")
abline(lm(log2(s$deg_up / s$deg_dn) ~ s$median_lfc_shift), lwd = 2, col = "grey40")
text(s$median_lfc_shift, log2(s$deg_up / s$deg_dn), sub("^(voom|trend)_", "", s$config),
     pos = c(4, 4, 4, 1, 4, 4, 2), cex = .6, offset = .6, xpd = NA)
abline(h = 0, lty = 3)
mtext(sprintf("Pearson r = %.4f across all seven configurations",
              cor(s$median_lfc_shift, log2(s$deg_up / s$deg_dn))), 3, .1, cex = .62)

## F. permutation null --------------------------------------------------------
pf <- file.path(OUTDIR, "permutation_null.csv")
if (file.exists(pf)) {
  pn <- read.csv(pf)
  bx <- split(pn$p_lt_05_frac, pn$config)
  boxplot(bx, col = c("grey70", "#c0392b", "grey85"), outline = TRUE, las = 1,
          names = sub("^voom_", "", names(bx)), ylab = "fraction of genes with p < 0.05",
          main = "F. Under the null every method is calibrated")
  abline(h = .05, col = "firebrick", lty = 2, lwd = 2)
  mtext(sprintf("%d permutations of the lithium column; nominal level 0.05 dashed",
                max(pn$perm)), 3, .1, cex = .62)
} else {
  plot.new(); title(main = "F. permutation null (pending)")
}
dev.off()
cat("wrote results/norm_sweep.png\n")
