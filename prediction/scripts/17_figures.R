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

## ---- Figure 4: the proposal's five inputs ------------------------------------
ff <- "results/summary_five_inputs.csv"
if (file.exists(ff)) {
  s5 <- read.csv(ff, stringsAsFactors = FALSE)
  s5$lab <- factor(s5$input, levels = rev(sort(unique(s5$input))))
  s5$mod <- factor(ifelse(s5$model == "enet", "Elastic net", "Random forest"),
                   levels = c("Elastic net", "Random forest"))
  p4 <- ggplot(s5, aes(auc, lab)) +
    geom_vline(xintercept = 0.5, colour = MUTED, linetype = "22", linewidth = 0.4) +
    geom_linerange(aes(xmin = 0.42, xmax = null_p95), linewidth = 6, colour = NULLBAR) +
    geom_linerange(aes(xmin = lo, xmax = hi, colour = mod), linewidth = 2, alpha = 0.30) +
    geom_point(aes(colour = mod), size = 2.8) +
    geom_text(aes(label = sprintf("%.2f", auc)), vjust = -1.4, size = 3.2, colour = INK) +
    facet_wrap(~mod, nrow = 1) +
    scale_colour_manual(values = c("Elastic net" = BLUE, "Random forest" = ORANGE)) +
    scale_x_continuous(limits = c(0.42, 1.02), breaks = seq(0.5, 1.0, 0.1)) +
    scale_y_discrete(expand = expansion(add = c(0.45, 0.85))) +
    labs(title = "Predicting lithium use: the five inputs we proposed",
         subtitle = paste0("152 users vs 74 non-users. Grey is measured chance, to the 95th\n",
           "percentile of the permutation null. Inputs 4 and 5 use 4,000 in-fold\n",
           "selected genes against 12,210 for inputs 1-3 -- a compute limit, not a choice."),
         x = "AUC", y = NULL) +
    base
  save2(p4, "fig4_five_inputs", 9.6, 4.8)
  cat("wrote fig4\n")
}

## ---- Figure 5: significance vs equivalence, side by side ---------------------
# One picture, two questions, and they are not the same question.
#   Does the interval cross the dotted zero line?  -> is it distinguishable
#     from the transcriptome at all?
#   Does the interval fit inside the shaded band?  -> is whatever difference
#     exists small enough to call the two interchangeable?
# An input can fail both, and at this sample size most do.
#
# Filled circle  = significantly worse (Holm p < 0.05).
# Hollow circle  = not distinguishable -- which means UNRESOLVED, not equal.
#
# Every number is a difference FROM the reference, so the reference has to be
# on the chart: it sits as the top row at zero with its own AUC, and the
# absolute AUC each difference lands on is printed beside every point. A
# reader should never have to go to another figure to learn what -0.306 means.
sf <- "results/significance_vs_genes.csv"
if (file.exists(sf)) {
  sg <- read.csv(sf, stringsAsFactors = FALSE)

  # Absolute AUCs: lithium from the five-input summary, bipolar from the
  # proposal summary. Both carry feature_set x model.
  abs_auc <- NULL
  if (file.exists("results/summary_five_inputs.csv")) {
    a <- read.csv("results/summary_five_inputs.csv", stringsAsFactors = FALSE)
    abs_auc <- rbind(abs_auc, data.frame(cohort = "lithium", feature_set = a$feature_set,
                                         model = a$model, abs = a$auc))
  }
  if (file.exists("results/summary_proposal.csv")) {
    b <- read.csv("results/summary_proposal.csv", stringsAsFactors = FALSE)
    b <- b[b$cohort == "groupA", ]
    abs_auc <- rbind(abs_auc, data.frame(cohort = "bipolar", feature_set = b$feature_set,
                                         model = b$model, abs = b$auc))
  }
  sg <- merge(sg, abs_auc, by = c("cohort", "feature_set", "model"), all.x = TRUE)

  REFLAB <- "1  genes (reference)"
  for (coh in unique(sg$cohort)) {
    x <- sg[sg$cohort == coh, ]
    ref <- abs_auc[abs_auc$cohort == coh & abs_auc$feature_set == "genes", ]
    # The reference itself, as a row at exactly zero.
    r <- data.frame(cohort = coh, feature_set = "genes", model = ref$model,
                    diff = 0, lo = 0, hi = 0, p = NA, p_holm = NA,
                    input = REFLAB, abs = ref$abs)
    x <- rbind(x[, names(r)], r)
    x$lab <- factor(x$input, levels = rev(sort(unique(x$input))))
    x$mod <- factor(ifelse(x$model == "enet", "Elastic net", "Random forest"),
                    levels = c("Elastic net", "Random forest"))
    x$sig <- factor(ifelse(x$input == REFLAB, "Reference: the full transcriptome",
                    ifelse(x$p_holm < 0.05, "Significantly worse than genes",
                                            "Not distinguishable (unresolved)")),
                    levels = c("Reference: the full transcriptome",
                               "Significantly worse than genes",
                               "Not distinguishable (unresolved)"))
    x$note <- ifelse(x$input == REFLAB, sprintf("AUC %.3f", x$abs),
                     sprintf("%+.3f   (AUC %.3f)", x$diff, x$abs))
    p5f <- ggplot(x, aes(diff, lab)) +
      annotate("rect", xmin = -0.05, xmax = 0.05, ymin = -Inf, ymax = Inf,
               fill = "#cfd8e3", alpha = 0.55) +
      geom_vline(xintercept = 0, colour = INK, linetype = "22", linewidth = 0.45) +
      geom_linerange(aes(xmin = lo, xmax = hi, colour = sig), linewidth = 1.9, alpha = 0.38) +
      geom_point(aes(colour = sig, fill = sig), size = 3.1, shape = 21, stroke = 1.3) +
      geom_text(aes(label = note), vjust = -1.5, size = 3.0, colour = INK) +
      facet_wrap(~mod, nrow = 1, scales = "free_x") +
      scale_colour_manual(values = c("Reference: the full transcriptome" = INK,
                                     "Significantly worse than genes" = ORANGE,
                                     "Not distinguishable (unresolved)" = BLUE),
                          drop = FALSE, name = NULL) +
      scale_fill_manual(values = c("Reference: the full transcriptome" = INK,
                                   "Significantly worse than genes" = ORANGE,
                                   "Not distinguishable (unresolved)" = "white"),
                        drop = FALSE, name = NULL) +
      scale_y_discrete(expand = expansion(add = c(0.55, 1.1))) +
      guides(colour = guide_legend(nrow = 1)) +
      labs(title = sprintf("%s: each input against the full transcriptome",
                           ifelse(coh == "lithium", "Lithium use", "Bipolar status")),
           subtitle = paste0("Everything is measured as a difference from input 1, the",
             " reference, whose own AUC is
the top row. Paired bootstrap over subjects",
             " (2,000), Holm-corrected.
",
             "FILLED = significantly worse.  HOLLOW = not distinguishable, which means the
",
             "data cannot tell them apart -- not that they are equal.
",
             "BLUE BAND = the +/-0.05 margin. An interval lying wholly inside it would mean
",
             "the two inputs are interchangeable. No interval here does."),
           x = "Difference in AUC from input 1 (genes)", y = NULL) +
      base +
      theme(legend.position = "top",
            legend.text = element_text(colour = INK, size = 8.5),
            legend.key.height = unit(10, "pt"),
            plot.subtitle = element_text(colour = MUTED, size = 9, lineheight = 1.3))
    save2(p5f, sprintf("fig5_vs_genes_%s", coh), 10.4, 5.8)
  }
  cat("wrote fig5
")
}
