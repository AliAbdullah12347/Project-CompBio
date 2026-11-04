#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Presentation figures for Arm 3. Reads only committed result CSVs — no model
# is refitted here, so a figure can never disagree with the table it came from.
#
#   Rscript network/make_figures.R        -> network/figures/*.png
#
# Sized for projection: large type, few words, one claim per figure.
#
# NOTE on figure 1. An earlier version put "technical 9.0%" and "age 5.7%" on
# the same bar chart as "composition 73.3%". Those are two different
# quantities — 9.0 and 5.7 are how much of the lithium *coefficient* those
# covariates explain away (README.md table 1), while 73.3 and 5.3 are shares of
# PC2's *variance* (table 2). Only the variance shares belong on that axis.
#
# Palette: validated categorical slots (CVD dE 9.1 protan, normal dE 22.9).
# ---------------------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(grDevices)))

R   <- "network/runs"
OUT <- "network/figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; AQUA <- "#1baf7a"
INK  <- "#0b0b0b"; INK2 <- "#52514e"; SURF <- "#ffffff"; GRID <- "#dcdbd6"

# pointsize 19 at res 150 makes the base font ~19pt on the slide, roughly
# double the previous figures. Every cex below is relative to that.
png_open <- function(f, w = 1800, h = 1050, ps = 16) {
  png(file.path(OUT, f), width = w, height = h, res = 150, bg = SURF, pointsize = ps)
  par(family = "sans", col.axis = INK2, col.lab = INK, fg = INK2,
      cex.axis = 1.0, cex.lab = 1.05, lend = "butt")
}
rd <- function(p) read.csv(file.path(R, p), check.names = FALSE, stringsAsFactors = FALSE)
# WGCNA appends two pseudo-modules: "gold" (a random gene set, the internal
# baseline) and "grey" (unassigned genes). Neither is a real module.
real <- function(m) !(m %in% c("gold", "grey"))

# title block: headline in ink, one supporting line in muted ink
head2 <- function(title, sub = NULL, at = NULL, l1 = 2.5, l2 = 1.0, cex = 1.25) {
  op <- par(xpd = NA); on.exit(par(op))
  mtext(title, side = 3, line = l1, adj = 0, outer = FALSE, cex = cex, font = 2, col = INK)
  if (!is.null(sub)) mtext(sub, side = 3, line = l2, adj = 0, cex = 0.78, col = INK2)
}

# ===========================================================================
# FIG 1 — the headline: the lithium network signal is cell composition
# ===========================================================================
t33 <- rd("trait-033/results_modtrait.csv"); t35 <- rd("trait-035/results_modtrait.csv")
li  <- unique(t33$module[t33$predictor == "lithium" & t33$padj < 0.05 & real(t33$module)])
ilr <- unique(t35$module[t35$padj < 0.05 & real(t35$module)])
n_mod <- length(unique(t33$module[real(t33$module)]))

png_open("fig1_lithium_is_composition.png", 1800, 1000)
layout(matrix(1:2, nrow = 1), widths = c(1.12, 1))

# -- A: share of PC2's variance each predictor explains on its own ----------
par(mar = c(7.4, 11.0, 6.0, 7.5))
v    <- c(5.3, 71.3, 73.3)
lab  <- c("Lithium status", "Cell fractions", "Cell composition")
cols <- c(ORANGE, BLUE, BLUE)
bp <- barplot(v, horiz = TRUE, col = cols, border = NA, xlim = c(0, 100),
              names.arg = rep("", 3), axes = FALSE, space = 0.55)
axis(1, at = seq(0, 100, 25), labels = paste0(seq(0, 100, 25), "%"), lwd = 0, lwd.ticks = 1)
text(-3, bp, lab, adj = 1, xpd = NA, cex = 1.0, col = INK)
text(v + 2.5, bp, sprintf("%.1f%%", v), adj = 0, xpd = NA, cex = 1.15, font = 2,
     col = ifelse(v > 50, BLUE, ORANGE))
head2("A   What is this axis made of?",
      "Variance explained, each predictor alone",
      l1 = 3.2, l2 = 1.6, cex = 1.2)
mtext("Adjusting for composition shrinks the lithium effect", side = 1, line = 4.0,
      adj = 0, cex = 0.76, col = INK2)
mtext(expression(paste("from ", beta, " = -0.91 (p = 0.0005) to ", beta, " = -0.24 (p = 0.081)")),
      side = 1, line = 5.2, adj = 0, cex = 0.76, col = INK2)

# -- B: the two module sets are nested, not merely overlapping --------------
par(mar = c(7.4, 2.0, 6.0, 2.0))
plot(0, 0, type = "n", xlim = c(0, 10), ylim = c(0, 10), axes = FALSE, xlab = "", ylab = "")
rect(0.3, 2.2, 9.7, 8.0, col = paste0(BLUE, "1F"), border = BLUE, lwd = 3)
rect(1.1, 3.0, 5.3, 7.2, col = paste0(ORANGE, "2E"), border = ORANGE, lwd = 3)
text(3.2, 5.6, length(li), cex = 2.6, font = 2, col = ORANGE)
text(3.2, 4.0, "track\nlithium", cex = 0.92, col = INK)
text(7.5, 5.6, length(ilr) - length(li), cex = 2.6, font = 2, col = BLUE)
text(7.5, 4.0, "composition\nonly", cex = 0.92, col = INK)
text(5.0, 8.7, sprintf("All %d track composition", length(ilr)), cex = 1.0, font = 2, col = BLUE, xpd = NA)
text(5.0, 1.15, "Lithium but NOT composition:  0", cex = 1.0, font = 2, col = ORANGE, xpd = NA)
head2("B   Lithium modules are a subset",
      "474 people, BH < 5%",
      l1 = 3.2, l2 = 1.6, cex = 1.2)
dev.off()

# ===========================================================================
# FIG 2 — the calibrated negative: controls vs unmedicated bipolar
# ===========================================================================
g <- rd("base-002-ceil/results_all_groups.csv")
ceil <- g[g$group == "ctrl_ceiling" & real(g$module), c("module", "medianRank")]
nol  <- g[g$group == "bp_nolith"    & real(g$module), c("module", "medianRank")]
m <- merge(ceil, nol, by = "module", suffixes = c("_ceil", "_nol"))
w <- suppressWarnings(wilcox.test(m$medianRank_ceil, m$medianRank_nol, paired = TRUE))
SD <- 18.2   # measured per-module shift sd at n = 74 (base-002-ceil)

png_open("fig2_calibrated_negative.png", 1180, 1150)
par(mar = c(6.6, 6.8, 7.4, 3.0))
lim <- c(0, max(m[, -1]) + 4)
plot(NA, xlim = lim, ylim = lim, xlab = "", ylab = "", axes = FALSE, asp = 1,
     xaxs = "i", yaxs = "i")
polygon(c(lim[1], lim[2], lim[2], lim[1]), c(lim[1] + SD, lim[2] + SD, lim[2] - SD, lim[1] - SD),
        col = GRID, border = NA)
abline(0, 1, col = INK2, lwd = 2.5)
points(m$medianRank_ceil, m$medianRank_nol, pch = 21, bg = paste0(BLUE, "CC"),
       col = SURF, cex = 1.9, lwd = 2)
axis(1, lwd = 0, lwd.ticks = 1); axis(2, lwd = 0, lwd.ticks = 1, las = 1)
mtext("Preservation rank in healthy controls", side = 1, line = 3.5, cex = 1.05, col = INK)
mtext("Rank in bipolar, off lithium", side = 2, line = 4.2, cex = 1.05, col = INK)
head2("No detectable difference",
      sprintf("%d modules, 74 people per group  ·  Wilcoxon p = %.2f", nrow(m), w$p.value),
      l1 = 4.1, l2 = 2.3, cex = 1.4)
mtext(sprintf("Grey band = what this measurement can resolve (+/- %.0f ranks)", SD),
      side = 3, line = 0.7, adj = 0, cex = 0.78, col = INK2)
text(lim[2] * 0.70, lim[2] * 0.28, "each dot is\none module", col = INK2, cex = 0.85)
dev.off()

# ===========================================================================
# FIG 3 — density survives, hub architecture does not
# ===========================================================================
p30 <- rd("pres-030/results_nolith.csv")
pr <- p30[real(p30$module), ]; pg <- p30[p30$module == "gold", ]

png_open("fig3_density_vs_connectivity.png", 1180, 1150)
par(mar = c(6.6, 6.8, 7.4, 3.0))
lim <- c(0, max(p30$Zdensity, p30$Zconnectivity, na.rm = TRUE) * 1.04)
plot(NA, xlim = lim, ylim = lim, xlab = "", ylab = "", axes = FALSE, asp = 1,
     xaxs = "i", yaxs = "i")
abline(0, 1, col = INK2, lwd = 2.5)
text(lim[2] * 0.58, lim[2] * 0.66, "equally preserved", col = INK2, cex = 0.85, srt = 45)
points(pr$Zdensity, pr$Zconnectivity, pch = 21, bg = paste0(BLUE, "CC"),
       col = SURF, cex = 1.9, lwd = 2)
points(pg$Zdensity, pg$Zconnectivity, pch = 23, bg = ORANGE, col = SURF, cex = 2.8, lwd = 2)
text(pg$Zdensity + 26, pg$Zconnectivity + 3, "random\ngene set", col = ORANGE, font = 2, cex = 0.88)
axis(1, lwd = 0, lwd.ticks = 1); axis(2, lwd = 0, lwd.ticks = 1, las = 1)
mtext("Do the module's genes still move together?", side = 1, line = 3.5, cex = 1.05, col = INK)
mtext("Is its hub structure intact?", side = 2, line = 4.2, cex = 1.05, col = INK)
head2("Modules stay coherent; their hubs do not",
      sprintf("%d of %d modules below the line: co-expression survives, hubs do not",
              sum(pr$Zdensity > pr$Zconnectivity), nrow(pr)),
      l1 = 4.1, l2 = 2.3, cex = 1.4)
mtext("Caution: the random gene set beats every real module here",
      side = 3, line = 0.7, adj = 0, cex = 0.78, col = ORANGE)
dev.off()

# ===========================================================================
# FIG 4 — modules carrying the lithium DEG signal
# ===========================================================================
a <- rd("ann-041/results_deg_overlap.csv")
w1 <- a[a$deg_list == "WB_LI" & a$padj < 0.05 & real(a$module), ]
w1 <- w1[order(w1$odds_ratio), ]
w1 <- tail(w1, 10)   # top 10 by enrichment — the rest crowd the axis

png_open("fig4_modules_carry_lithium_degs.png", 1700, 1050)
par(mar = c(6.4, 11.5, 6.8, 11.5))
xmax <- max(w1$odds_ratio) * 1.1
bp <- barplot(w1$odds_ratio, horiz = TRUE, col = BLUE, border = NA, axes = FALSE,
              names.arg = rep("", nrow(w1)), space = 0.45, xlim = c(0, xmax))
axis(1, lwd = 0, lwd.ticks = 1)
text(-xmax * 0.025, bp, w1$module, adj = 1, xpd = NA, cex = 0.95, col = INK)
text(w1$odds_ratio + xmax * 0.02, bp,
     sprintf("%.1fx   (%d genes, %.0f expected)", w1$odds_ratio, w1$n_overlap, w1$expected),
     adj = 0, xpd = NA, cex = 0.82, col = INK2)
head2("Which modules carry the lithium signal?",
      sprintf("Top 10 of %d enriched modules (Fisher exact, BH < 5%%)",
              nrow(a[a$deg_list == "WB_LI" & a$padj < 0.05 & real(a$module), ])),
      l1 = 3.5, l2 = 1.8, cex = 1.35)
mtext("Enrichment  (times more DE genes than expected by chance)", side = 1, line = 3.4,
      cex = 0.95, col = INK)
dev.off()

# ===========================================================================
# FIG 5 — the group-label permutation null
# ===========================================================================
pt <- rd("perm-031/results_permutation_test.csv")
pr5 <- pt[real(pt$module), ]; pgold <- pt[pt$module == "gold", ]
nd <- rd("perm-031/null_zsummary_alldraws.csv")
gold_null <- nd$Zsummary[nd$module == "gold"]; gold_obs <- pgold$obs_Zsummary_mean
below <- sum(pr5$obs_Zsummary_mean < pr5$null_Zsummary_mean)

png_open("fig5_permutation_null.png", 1800, 1000)
layout(matrix(1:2, nrow = 1), widths = c(1, 1))

# -- A: every module, observed vs null -------------------------------------
par(mar = c(6.6, 6.6, 7.0, 2.4))
lim <- c(0, max(pt$obs_Zsummary_mean, pt$null_Zsummary_mean, na.rm = TRUE) * 1.05)
plot(NA, xlim = lim, ylim = lim, xlab = "", ylab = "", axes = FALSE, asp = 1,
     xaxs = "i", yaxs = "i")
abline(0, 1, col = INK2, lwd = 2.5)
points(pr5$null_Zsummary_mean, pr5$obs_Zsummary_mean, pch = 21, bg = paste0(BLUE, "CC"),
       col = SURF, cex = 1.7, lwd = 1.8)
axis(1, lwd = 0, lwd.ticks = 1); axis(2, lwd = 0, lwd.ticks = 1, las = 1)
mtext("Expected by chance", side = 1, line = 3.4, cex = 1.0, col = INK)
mtext("Observed in lithium", side = 2, line = 3.9, cex = 1.0, col = INK)
head2("A   Across modules: a weak trend",
      sprintf("%d of %d (%.0f%%) below the line, vs 50%% by chance",
              below, nrow(pr5), 100 * below / nrow(pr5)),
      l1 = 3.5, l2 = 1.8, cex = 1.2)
mtext("Not a result on its own", side = 3, line = 0.6, adj = 0, cex = 0.78, col = INK2)

# -- B: the random-gene-set baseline, where the signal actually is ----------
par(mar = c(6.6, 3.0, 7.0, 2.4))
xr <- range(c(gold_null, gold_obs)) + c(-2.5, 2.5)
plot(NA, xlim = xr, ylim = c(0.2, 1.8), axes = FALSE, xlab = "", ylab = "")
set.seed(1); jit <- 1 + runif(length(gold_null), -0.26, 0.26)
points(gold_null, jit, pch = 21, bg = paste0(BLUE, "B0"), col = SURF, cex = 2.1, lwd = 1.8)
segments(gold_obs, 0.35, gold_obs, 1.52, col = ORANGE, lwd = 4)
text(gold_obs, 1.68, sprintf("observed  %.1f", gold_obs), col = ORANGE, font = 2, cex = 0.92)
axis(1, lwd = 0, lwd.ticks = 1)
mtext("Network preservation score", side = 1, line = 3.4, cex = 1.0, col = INK)
head2("B   The one clean signal",
      sprintf("%d of %d shuffled groups score above the real lithium group",
              sum(gold_null > gold_obs), length(gold_null)),
      l1 = 3.5, l2 = 1.8, cex = 1.2)
mtext("p = 0.05 is the floor at 20 shuffles", side = 3,
      line = 0.6, adj = 0, cex = 0.78, col = ORANGE)
dev.off()

cat(sprintf("wrote %d figures to %s\n", length(list.files(OUT, "\\.png$")), OUT))
for (f in sort(list.files(OUT, "\\.png$", full.names = TRUE)))
  cat(sprintf("  %-42s %5.0f KB\n", basename(f), file.size(f) / 1024))
cat(sprintf("\nfig1: %d lithium, %d composition, %d lithium-only, of %d modules\n",
            length(li), length(ilr), length(setdiff(li, ilr)), n_mod))
cat(sprintf("fig2: %d modules, Wilcoxon p = %.4f\n", nrow(m), w$p.value))
cat(sprintf("fig3: %d/%d density>connectivity\n", sum(pr$Zdensity > pr$Zconnectivity), nrow(pr)))
cat(sprintf("fig4: top %d of %d BH-sig modules shown\n", nrow(w1),
            nrow(a[a$deg_list == "WB_LI" & a$padj < 0.05 & real(a$module), ])))
cat(sprintf("fig5: %d/%d below diagonal; gold obs %.2f, %d/%d null draws higher\n",
            below, nrow(pr5), gold_obs, sum(gold_null > gold_obs), length(gold_null)))
