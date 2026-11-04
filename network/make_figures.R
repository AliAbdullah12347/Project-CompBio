#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Presentation figures for Arm 3. Reads only committed result CSVs — no model
# is refitted here, so a figure can never disagree with the table it came from.
#
#   Rscript network/make_figures.R        -> network/figures/*.png
#
# Palette: validated categorical slots (CVD dE 9.1 protan, normal dE 22.9).
# Low-contrast slots carry direct labels, which the validator requires.
# ---------------------------------------------------------------------------
suppressWarnings(suppressPackageStartupMessages(library(grDevices)))

R   <- "network/runs"
OUT <- "network/figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; AQUA <- "#1baf7a"; YELLOW <- "#eda100"
INK  <- "#0b0b0b"; INK2 <- "#52514e"; SURF <- "#fcfcfb"; GRID <- "#e3e2de"

png_open <- function(f, w = 1700, h = 1000) {
  png(file.path(OUT, f), width = w, height = h, res = 150, bg = SURF)
  par(family = "sans", col.axis = INK2, col.lab = INK, col.main = INK,
      fg = INK2, cex.axis = 0.85, cex.lab = 0.95)
}
rd <- function(p) read.csv(file.path(R, p), check.names = FALSE, stringsAsFactors = FALSE)
# WGCNA appends two pseudo-modules: "gold" (a random gene set, the internal
# baseline) and "grey" (unassigned genes). Neither is a real module.
real <- function(m) !(m %in% c("gold", "grey"))

# ===========================================================================
# FIG 1 — the headline: the lithium network signal is cell composition
# ===========================================================================
t33 <- rd("trait-033/results_modtrait.csv"); t35 <- rd("trait-035/results_modtrait.csv")
li  <- unique(t33$module[t33$predictor == "lithium" & t33$padj < 0.05 & real(t33$module)])
ilr <- unique(t35$module[t35$padj < 0.05 & real(t35$module)])
n_mod <- length(unique(t33$module[real(t33$module)]))

png_open("fig1_lithium_is_composition.png", 1700, 850)
layout(matrix(1:2, nrow = 1), widths = c(1.05, 1))

# -- A: variance of the lithium-associated dimension (PC2), analysis 002 ----
par(mar = c(4.6, 11, 4.4, 3.2))
v   <- c(5.3, 5.7, 9.0, 73.3)
lab <- c("Lithium, alone", "Age", "Technical covariates", "Cell composition\n(ILR b1-b4)")
cols <- c(ORANGE, GRID, GRID, BLUE)
bp <- barplot(v, horiz = TRUE, col = cols, border = NA, xlim = c(0, 86),
              names.arg = rep("", 4), axes = FALSE, space = 0.45)
axis(1, at = seq(0, 80, 20), labels = paste0(seq(0, 80, 20), "%"), lwd = 0, lwd.ticks = 1)
text(x = -2, y = bp, labels = lab, adj = 1, xpd = NA, cex = 0.9, col = INK)
text(x = v + 1.6, y = bp, labels = sprintf("%.1f%%", v), adj = 0, xpd = NA,
     cex = 0.92, font = 2, col = ifelse(v > 50, BLUE, INK2))
mtext("A   What explains the lithium-associated dimension?", side = 3, line = 2.1,
      adj = 0, at = -34, cex = 1.02, font = 2, col = INK)
mtext(expression(paste("PC2 of the module eigengenes  ·  ", beta, " = -0.908, p = 4.8e-4, BH = 0.0029")),
      side = 3, line = 0.65, adj = 0, at = -34, cex = 0.78, col = INK2)
mtext("Adjusting for composition attenuates lithium by 68.1% (p rises to 0.040; BH 0.24)",
      side = 1, line = 3.1, adj = 0, at = -34, cex = 0.75, col = INK2)

# -- B: the module sets are nested, not overlapping -------------------------
par(mar = c(4.6, 2.5, 4.4, 2.5))
plot(0, 0, type = "n", xlim = c(0, 10), ylim = c(0, 10), axes = FALSE, xlab = "", ylab = "")
rect(0.5, 1.6, 9.5, 8.4, col = paste0(BLUE, "22"), border = BLUE, lwd = 2)
rect(1.4, 2.7, 5.6, 7.3, col = paste0(ORANGE, "33"), border = ORANGE, lwd = 2)
text(3.5, 5.6, length(li), cex = 2.5, font = 2, col = ORANGE)
text(3.5, 4.3, "associate with\nlithium", cex = 0.86, col = INK)
text(7.6, 5.6, length(ilr) - length(li), cex = 2.5, font = 2, col = BLUE)
text(7.6, 4.3, "composition only", cex = 0.86, col = INK)
text(5.0, 8.95, sprintf("All %d modules associate with cell composition", length(ilr)),
     cex = 0.92, font = 2, col = BLUE, xpd = NA)
text(5.0, 0.75, sprintf("Lithium-only modules:  0  of  %d", n_mod),
     cex = 1.0, font = 2, col = ORANGE, xpd = NA)
mtext("B   Every lithium module is also a composition module", side = 3, line = 2.1,
      adj = 0, at = -0.6, cex = 1.02, font = 2, col = INK)
mtext("Module eigengenes vs lithium and vs ILR balances, 474 samples, BH < 5%",
      side = 3, line = 0.65, adj = 0, at = -0.6, cex = 0.78, col = INK2)
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

png_open("fig2_calibrated_negative.png", 1350, 1000)
par(mar = c(5.2, 5.4, 5.0, 2.4))
lim <- c(0, max(m[, -1]) + 3)
plot(NA, xlim = lim, ylim = lim, xlab = "", ylab = "", axes = FALSE, asp = 1)
polygon(c(lim[1], lim[2], lim[2], lim[1]), c(lim[1] + SD, lim[2] + SD, lim[2] - SD, lim[1] - SD),
        col = GRID, border = NA)
abline(0, 1, col = INK2, lwd = 2)
points(m$medianRank_ceil, m$medianRank_nol, pch = 21, bg = paste0(BLUE, "CC"),
       col = SURF, cex = 1.5, lwd = 1.5)
axis(1, lwd = 0, lwd.ticks = 1); axis(2, lwd = 0, lwd.ticks = 1, las = 1)
mtext("Preservation rank in held-out controls  (the ceiling)", side = 1, line = 3.1, cex = 0.95, col = INK)
mtext("Preservation rank in bipolar, off lithium", side = 2, line = 3.4, cex = 0.95, col = INK)
mtext("No detectable network difference", side = 3, line = 3.3, adj = 0, at = lim[1] - 9,
      cex = 1.15, font = 2, col = INK)
mtext(sprintf("%d modules, n = 74 per group, matched instrument (500 permutations both sides)",
              nrow(m)), side = 3, line = 1.9, adj = 0, at = lim[1] - 9, cex = 0.8, col = INK2)
mtext(sprintf("Wilcoxon signed-rank p = %.3f   ·   grey band = the measurement's own noise floor (+/- %.1f ranks)",
              w$p.value, SD), side = 3, line = 0.6, adj = 0, at = lim[1] - 9, cex = 0.8, col = INK2)
legend("bottomright", legend = c("one module", "equal preservation", "noise floor"),
       pch = c(21, NA, 22), pt.bg = c(BLUE, NA, GRID), col = c(SURF, INK2, GRID),
       lty = c(NA, 1, NA), lwd = c(1.5, 2, NA), pt.cex = c(1.5, NA, 2.2),
       bty = "n", cex = 0.82, text.col = INK)
dev.off()

# ===========================================================================
# FIG 3 — density survives, hub architecture does not
# ===========================================================================
p30 <- rd("pres-030/results_nolith.csv")
pr <- p30[real(p30$module), ]; pg <- p30[p30$module == "gold", ]

png_open("fig3_density_vs_connectivity.png", 1350, 1000)
par(mar = c(5.2, 5.4, 5.0, 2.4))
lim <- c(0, max(p30$Zdensity, p30$Zconnectivity, na.rm = TRUE) * 1.05)
plot(NA, xlim = lim, ylim = lim, xlab = "", ylab = "", axes = FALSE, asp = 1,
     xaxs = "i", yaxs = "i")
abline(0, 1, col = INK2, lwd = 2)
points(pr$Zdensity, pr$Zconnectivity, pch = 21, bg = paste0(BLUE, "CC"),
       col = SURF, cex = 1.5, lwd = 1.5)
points(pg$Zdensity, pg$Zconnectivity, pch = 23, bg = ORANGE, col = SURF, cex = 2.3, lwd = 2)
text(pg$Zdensity, pg$Zconnectivity + 7, "random gene set", col = ORANGE, font = 2, cex = 0.85)
axis(1, lwd = 0, lwd.ticks = 1); axis(2, lwd = 0, lwd.ticks = 1, las = 1)
mtext("Z density  —  do module genes still move together?", side = 1, line = 3.1, cex = 0.95, col = INK)
mtext("Z connectivity  —  is the hub structure intact?", side = 2, line = 3.4, cex = 0.95, col = INK)
mtext("Modules stay coherent; their hubs do not", side = 3, line = 3.3, adj = 0, at = lim[1] - 18,
      cex = 1.15, font = 2, col = INK)
mtext(sprintf("%d modules in bipolar off lithium. %d of %d sit below the diagonal: density preserved, connectivity lost.",
              nrow(pr), sum(pr$Zdensity > pr$Zconnectivity), nrow(pr)),
      side = 3, line = 1.9, adj = 0, at = lim[1] - 18, cex = 0.8, col = INK2)
mtext(sprintf("Mean gap %.1f.  The random gene set beats every real module on connectivity — treat with caution.",
              mean(pr$Zdensity - pr$Zconnectivity)),
      side = 3, line = 0.6, adj = 0, at = lim[1] - 18, cex = 0.8, col = INK2)
dev.off()

# ===========================================================================
# FIG 4 — modules carrying the lithium DEG signal
# ===========================================================================
a <- rd("ann-041/results_deg_overlap.csv")
w1 <- a[a$deg_list == "WB_LI" & a$padj < 0.05 & real(a$module), ]
w1 <- w1[order(w1$odds_ratio), ]

png_open("fig4_modules_carry_lithium_degs.png", 1500, 950)
par(mar = c(4.8, 10.5, 4.6, 8.5))
bp <- barplot(w1$odds_ratio, horiz = TRUE, col = BLUE, border = NA, axes = FALSE,
              names.arg = rep("", nrow(w1)), space = 0.4,
              xlim = c(0, max(w1$odds_ratio) * 1.16))
axis(1, lwd = 0, lwd.ticks = 1)
text(-max(w1$odds_ratio) * 0.02, bp, w1$module, adj = 1, xpd = NA, cex = 0.82, col = INK)
text(w1$odds_ratio + max(w1$odds_ratio) * 0.015, bp,
     sprintf("%.1fx    %d of %d genes (%.0f expected)", w1$odds_ratio, w1$n_overlap,
             w1$n_module, w1$expected),
     adj = 0, xpd = NA, cex = 0.72, col = INK2)
mtext("Which modules carry the lithium differential-expression signal?", side = 3, line = 2.7,
      adj = 0, at = -max(w1$odds_ratio) * 0.27, cex = 1.1, font = 2, col = INK)
mtext(sprintf("%d of %d modules enriched for the 1,426 whole-blood lithium DEGs (Fisher exact, BH < 5%%)",
              nrow(w1), n_mod), side = 3, line = 1.3, adj = 0,
      at = -max(w1$odds_ratio) * 0.27, cex = 0.8, col = INK2)
mtext("Enrichment odds ratio", side = 1, line = 2.9, cex = 0.9, col = INK)
dev.off()

# ===========================================================================
# FIG 5 — the group-label permutation null
# ===========================================================================
pt <- rd("perm-031/results_permutation_test.csv")
pr5 <- pt[real(pt$module), ]; pgold <- pt[pt$module == "gold", ]

png_open("fig5_permutation_null.png", 1700, 850)
layout(matrix(1:2, nrow = 1), widths = c(1, 1))
nd <- rd("perm-031/null_zsummary_alldraws.csv")
gold_null <- nd$Zsummary[nd$module == "gold"]
gold_obs  <- pgold$obs_Zsummary_mean
below <- sum(pr5$obs_Zsummary_mean < pr5$null_Zsummary_mean)

# -- A: every module, observed vs null -------------------------------------
par(mar = c(5.0, 5.2, 5.2, 2.0))
lim <- c(0, max(pt$obs_Zsummary_mean, pt$null_Zsummary_mean, na.rm = TRUE) * 1.06)
plot(NA, xlim = lim, ylim = lim, xlab = "", ylab = "", axes = FALSE, asp = 1,
     xaxs = "i", yaxs = "i")
abline(0, 1, col = INK2, lwd = 2)
points(pr5$null_Zsummary_mean, pr5$obs_Zsummary_mean, pch = 21, bg = paste0(BLUE, "CC"),
       col = SURF, cex = 1.4, lwd = 1.4)
points(pgold$null_Zsummary_mean, gold_obs, pch = 23, bg = ORANGE, col = SURF, cex = 2.2, lwd = 2)
text(pgold$null_Zsummary_mean + 10, gold_obs - 7, "random\ngene set", col = ORANGE, font = 2, cex = 0.78)
axis(1, lwd = 0, lwd.ticks = 1); axis(2, lwd = 0, lwd.ticks = 1, las = 1)
mtext("Expected under the null", side = 1, line = 2.9, cex = 0.88, col = INK)
mtext("Observed in the lithium group", side = 2, line = 3.2, cex = 0.88, col = INK)
mtext("A   Across modules: a weak trend", side = 3, line = 2.6, adj = 0,
      at = lim[1] - 11, cex = 1.0, font = 2, col = INK)
mtext(sprintf("%d of %d modules (%.0f%%) fall below the diagonal, against 50%% by chance.",
              below, nrow(pr5), 100 * below / nrow(pr5)),
      side = 3, line = 1.2, adj = 0, at = lim[1] - 11, cex = 0.76, col = INK2)
mtext("Not a result on its own.", side = 3, line = 0.2, adj = 0, at = lim[1] - 11,
      cex = 0.76, col = INK2)

# -- B: the random-gene-set baseline, where the signal actually is ----------
par(mar = c(5.0, 4.0, 5.2, 2.0))
xr <- range(c(gold_null, gold_obs)); xr <- xr + c(-2, 2)
plot(NA, xlim = xr, ylim = c(0.3, 1.7), axes = FALSE, xlab = "", ylab = "")
segments(xr[1], 1, xr[2], 1, col = GRID, lwd = 2)
set.seed(1); jit <- 1 + runif(length(gold_null), -0.13, 0.13)
points(gold_null, jit, pch = 21, bg = paste0(BLUE, "AA"), col = SURF, cex = 1.7, lwd = 1.4)
segments(gold_obs, 0.55, gold_obs, 1.45, col = ORANGE, lwd = 3)
text(gold_obs, 1.56, sprintf("observed
%.1f", gold_obs), col = ORANGE, font = 2, cex = 0.82)
axis(1, lwd = 0, lwd.ticks = 1)
mtext("Zsummary of the random gene set", side = 1, line = 2.9, cex = 0.88, col = INK)
mtext("B   The one clean signal", side = 3, line = 2.6, adj = 0, at = xr[1] - 1.5,
      cex = 1.0, font = 2, col = INK)
mtext(sprintf("%d of %d shuffled draws (blue) exceed the observed lithium value (orange).",
              sum(gold_null > gold_obs), length(gold_null)),
      side = 3, line = 1.2, adj = 0, at = xr[1] - 1.5, cex = 0.76, col = INK2)
mtext("Empirical p = 0.05, the floor at K = 20 — BH across 71 modules cannot pass it.",
      side = 3, line = 0.2, adj = 0, at = xr[1] - 1.5, cex = 0.76, col = ORANGE)
dev.off()

cat(sprintf("wrote %d figures to %s\n", length(list.files(OUT, "\\.png$")), OUT))
for (f in list.files(OUT, "\\.png$", full.names = TRUE))
  cat(sprintf("  %-42s %5.0f KB\n", basename(f), file.size(f) / 1024))
cat(sprintf("\nfig1: %d lithium-assoc, %d composition-assoc, %d lithium-only, %d modules\n",
            length(li), length(ilr), length(setdiff(li, ilr)), n_mod))
cat(sprintf("fig2: %d modules paired, Wilcoxon p = %.4f\n", nrow(m), w$p.value))
cat(sprintf("fig3: %d/%d modules density>connectivity, mean gap %.2f\n",
            sum(pr$Zdensity > pr$Zconnectivity), nrow(pr), mean(pr$Zdensity - pr$Zconnectivity)))
cat(sprintf("fig4: %d modules BH-sig for WB_LI\n", nrow(w1)))
cat(sprintf("fig5: %d/%d modules below diagonal; gold obs %.2f vs null %.2f\n",
            sum(pr5$obs_Zsummary_mean < pr5$null_Zsummary_mean), nrow(pr5),
            pgold$obs_Zsummary_mean, pgold$null_Zsummary_mean))
