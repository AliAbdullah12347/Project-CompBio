#!/usr/bin/env Rscript
# Figures for the prediction contests.
#
# Dot-and-interval, not bars. AUC's meaningful reference is 0.5, not 0, and a
# bar growing from zero would make 0.51 look like half a result instead of
# none. The chance line is drawn where it belongs and the bars are intervals.
#
# Each row also carries its own measured chance: the grey band runs from the
# mean of the permutation null to its 95th percentile. That band, not 0.5, is
# what a result has to clear -- small-n CV with thousands of features does not
# behave like a coin, and here the nulls sit slightly BELOW 0.5.
#
# Two series, blue and orange: 33.6 normal-vision separation and 24.7 under
# protanopia, so the pair survives colourblindness and greyscale print. Both
# are direct-labelled as well, so identity never rests on colour alone.
#
# A second copy of every figure goes to figures/for_docs at 300 dpi with the
# SAME geometry. ggplot lays text out in inches and does not reflow it, so
# re-saving at a smaller width runs the title off the edge; only the pixel
# density changes.
#
# Usage: Rscript 17_figures.R [summary_csv] [outdir]

suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
sumf <- if (length(args)) args[1] else "results/summary_auc.csv"
outd <- if (length(args) > 1) args[2] else "results/figures"
dir.create(file.path(outd, "for_docs"), showWarnings = FALSE, recursive = TRUE)

BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; NULLBAR <- "#d8d8d4"
INK <- "#1a1a19"; MUTED <- "#6b6b68"; GRID <- "#e6e6e3"; SURFACE <- "#fcfcfb"

base <- theme_minimal(base_size = 11) +
  theme(panel.background = element_rect(fill = SURFACE, colour = NA),
        plot.background  = element_rect(fill = SURFACE, colour = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_line(colour = GRID, linewidth = 0.3),
        axis.text  = element_text(colour = INK),
        axis.title = element_text(colour = MUTED, size = 9),
        plot.title    = element_text(colour = INK, face = "bold", size = 12),
        plot.subtitle = element_text(colour = MUTED, size = 9.5, lineheight = 1.25),
        strip.text = element_text(colour = INK, face = "bold", size = 9.5, lineheight = 1.15),
        plot.margin = margin(12, 16, 10, 12),
        legend.position = "none")

s0 <- read.csv(sumf, stringsAsFactors = FALSE)
LAB <- c(covariates  = "Covariates only",
         genes       = "Genes",
         genes_cov   = "Genes + covariates",
         genes_resid = "Genes, covariates\nregressed out")
MODNAME <- c(enet = "Elastic net", rf = "Random forest")

# Save at one geometry, twice: screen at 200 dpi, documents at 300.
save2 <- function(p, stem, w, h) {
  ggsave(file.path(outd, paste0(stem, ".png")), p, width = w, height = h, dpi = 200, bg = SURFACE)
  ggsave(file.path(outd, "for_docs", paste0(stem, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
}

for (md in c("enet", "rf")) {
  s <- s0[s0$model == md, ]
  if (!nrow(s)) next
  s$fs <- factor(LAB[s$feature_set], levels = rev(LAB))

  ## ---- Figure 1: the headline ------------------------------------------------
  f1 <- s[s$cohort %in% c("lithium", "groupA"), ]
  f1$contrast <- factor(ifelse(f1$cohort == "lithium",
                               "Lithium use\n152 users vs 74 non-users",
                               "Bipolar status\n74 unmedicated BP1 vs 112 controls"),
                        levels = c("Lithium use\n152 users vs 74 non-users",
                                   "Bipolar status\n74 unmedicated BP1 vs 112 controls"))
  p1 <- ggplot(f1, aes(auc, fs)) +
    geom_vline(xintercept = 0.5, colour = MUTED, linetype = "22", linewidth = 0.4) +
    geom_linerange(aes(xmin = null_mean, xmax = null_p95), linewidth = 6, colour = NULLBAR) +
    geom_linerange(aes(xmin = lo, xmax = hi, colour = contrast), linewidth = 2, alpha = 0.30) +
    geom_point(aes(colour = contrast), size = 2.8) +
    geom_text(aes(label = sprintf("%.2f", auc)), vjust = -1.4, size = 3.2, colour = INK) +
    facet_wrap(~contrast, nrow = 1) +
    scale_colour_manual(values = setNames(c(BLUE, ORANGE), levels(f1$contrast))) +
    scale_x_continuous(limits = c(0.42, 1.02), breaks = seq(0.5, 1.0, 0.1)) +
    scale_y_discrete(expand = expansion(add = c(0.45, 0.85))) +
    labs(title = "Whole blood predicts the drug, not the illness",
         subtitle = paste0(MODNAME[md],
           ". Pooled out-of-fold AUC, 25 repeats of 3-fold CV.\n",
           "Coloured bars are 95% CI over subjects. Grey band is measured chance,\n",
           "from label permutations: its mean to its 95th percentile."),
         x = "AUC", y = NULL) +
    base
  save2(p1, sprintf("fig1_%s", md), 8.4, 4.4)

  ## ---- Figure 2: why the control group decides the answer --------------------
  f2 <- s[s$cohort %in% c("groupA", "full", "groupA_plate"), ]
  CL <- c(full = "All controls\n74 cases vs 205",
          groupA = "Group A only\n74 vs 112",
          groupA_plate = "Group A, shared plates\n44 vs 109")
  f2$coh <- factor(CL[f2$cohort], levels = CL)
  p2 <- ggplot(f2, aes(auc, fs)) +
    geom_vline(xintercept = 0.5, colour = MUTED, linetype = "22", linewidth = 0.4) +
    geom_linerange(aes(xmin = null_mean, xmax = null_p95), linewidth = 6, colour = NULLBAR) +
    geom_linerange(aes(xmin = lo, xmax = hi, colour = cohort == "full"), linewidth = 2, alpha = 0.30) +
    geom_point(aes(colour = cohort == "full"), size = 2.8) +
    geom_text(aes(label = sprintf("%.2f", auc)), vjust = -1.4, size = 3.1, colour = INK) +
    facet_wrap(~coh, nrow = 1) +
    scale_colour_manual(values = c("TRUE" = ORANGE, "FALSE" = MUTED)) +
    scale_x_continuous(limits = c(0.42, 1.0), breaks = seq(0.5, 0.9, 0.2)) +
    scale_y_discrete(expand = expansion(add = c(0.45, 0.85))) +
    labs(title = "One cohort decision flips the bipolar result",
         subtitle = paste0(MODNAME[md],
           ". Assessment group B holds 109 controls and no cases, and is\n",
           "predictable from expression at AUC 0.94 -- so including it lets genes\n",
           "beat the covariates. Restricted to group A, that advantage disappears."),
         x = "AUC", y = NULL) +
    base
  save2(p2, sprintf("fig2_%s", md), 9.0, 4.4)
}

cat("wrote:\n")
for (f in list.files(outd, "png$", recursive = TRUE, full.names = TRUE))
  cat(sprintf("  %s\n", f))

## ---- Figure 3: Arm 1's primary question -------------------------------------
# Does a handful of compositional coordinates predict as well as the whole
# transcriptome? Only drawn when the cell-composition summary exists.
cf <- "results/summary_cells.csv"
if (file.exists(cf)) {
  sc <- read.csv(cf, stringsAsFactors = FALSE)
  # One covariate set per contrast, or the two lithium runs overplot each
  # other: `confounder` is lithium's primary (age + sex), `full` is the
  # case-control primary. confounder_tech is the sensitivity arm and is
  # reported in the table, not here.
  sc <- sc[(sc$cohort == "lithium" & sc$covset == "confounder") |
           (sc$cohort != "lithium" & sc$covset == "full"), ]
  sc <- sc[sc$model == "enet" & sc$feature_set %in%
             c("cells", "cells_cov", "genes", "genes_cells", "genes_decell"), ]
  L3 <- c(cells        = "Cell composition\n(4 ILR balances)",
          cells_cov    = "Cells + age, sex",
          genes        = "Genes\n(12,210)",
          genes_cells  = "Genes + cells",
          genes_decell = "Genes, cell mixture\nprojected out")
  sc$fs <- factor(L3[sc$feature_set], levels = rev(L3))
  sc$contrast <- factor(ifelse(sc$cohort == "lithium",
                               "Lithium use\n152 users vs 74 non-users",
                               "Bipolar status\n74 unmedicated BP1 vs 112 controls"),
                        levels = c("Lithium use\n152 users vs 74 non-users",
                                   "Bipolar status\n74 unmedicated BP1 vs 112 controls"))
  p3 <- ggplot(sc, aes(auc, fs)) +
    geom_vline(xintercept = 0.5, colour = MUTED, linetype = "22", linewidth = 0.4) +
    geom_linerange(aes(xmin = null_mean, xmax = null_p95), linewidth = 6, colour = NULLBAR) +
    geom_linerange(aes(xmin = lo, xmax = hi, colour = contrast), linewidth = 2, alpha = 0.30) +
    geom_point(aes(colour = contrast), size = 2.8) +
    geom_text(aes(label = sprintf("%.2f", auc)), vjust = -1.4, size = 3.2, colour = INK) +
    facet_wrap(~contrast, nrow = 1) +
    scale_colour_manual(values = setNames(c(BLUE, ORANGE), levels(sc$contrast))) +
    scale_x_continuous(limits = c(0.42, 1.02), breaks = seq(0.5, 1.0, 0.1)) +
    scale_y_discrete(expand = expansion(add = c(0.45, 0.85))) +
    labs(title = "A blood count is not a substitute for the transcriptome",
         subtitle = paste0("Elastic net. Cells fall 0.29 short of genes for lithium",
           " (CI -0.303 to -0.282),\nfar outside the pre-declared 0.05 equivalence margin.",
           " But projecting the\ncell mixture out of the genes costs only 0.05 -- the",
           " signal is not mostly mixture."),
         x = "AUC", y = NULL) +
    base
  save2(p3, "fig3_cells_vs_genes", 8.4, 4.8)
  cat("wrote fig3\n")
}
