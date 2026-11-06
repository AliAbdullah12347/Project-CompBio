#!/usr/bin/env Rscript
# Figures for the prediction contests.
#
# Dot-and-interval, not bars. AUC's meaningful reference is 0.5, not 0, and a
# bar growing from zero would make 0.51 look like half a result instead of
# none. The chance line is drawn where it belongs and the bars are intervals.
#
# Two series, blue and orange: 33.6 normal-vision separation and 24.7 under
# protanopia, so the pair survives colourblindness and greyscale print. Both
# are direct-labelled as well, so identity never rests on colour alone.
#
# Usage: Rscript 17_figures.R [summary_csv] [outdir]

suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
sumf <- if (length(args)) args[1] else "results/summary_auc.csv"
outd <- if (length(args) > 1) args[2] else "results/figures"
dir.create(outd, showWarnings = FALSE, recursive = TRUE)

BLUE <- "#2a78d6"; ORANGE <- "#eb6834"
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
        plot.subtitle = element_text(colour = MUTED, size = 9.5, lineheight = 1.2),
        plot.caption  = element_text(colour = MUTED, size = 8, hjust = 0),
        strip.text = element_text(colour = INK, face = "bold", size = 10),
        legend.position = "none")

s <- read.csv(sumf, stringsAsFactors = FALSE)
s <- s[s$model == "enet", ]
LAB <- c(covariates  = "Covariates only\n(age, sex, smoking, RIN, plate)",
         genes       = "Genes",
         genes_cov   = "Genes + covariates",
         genes_resid = "Genes, covariates\nregressed out")
s$fs <- factor(LAB[s$feature_set], levels = rev(LAB))

## ---- Figure 1: the headline -------------------------------------------------
f1 <- s[s$cohort %in% c("lithium", "groupA"), ]
f1$contrast <- ifelse(f1$cohort == "lithium", "Lithium use", "Bipolar status")
f1$contrast <- factor(f1$contrast, levels = c("Lithium use", "Bipolar status"))

p1 <- ggplot(f1, aes(auc, fs, colour = contrast)) +
  geom_vline(xintercept = 0.5, colour = MUTED, linetype = "22", linewidth = 0.4) +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 2, alpha = 0.28) +
  geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("%.2f", auc)), vjust = -1.25, size = 3.1,
            colour = INK, show.legend = FALSE) +
  facet_wrap(~contrast, nrow = 1) +
  scale_colour_manual(values = c("Lithium use" = BLUE, "Bipolar status" = ORANGE)) +
  scale_x_continuous(limits = c(0.42, 1.02), breaks = seq(0.5, 1.0, 0.1)) +
  labs(title = "Whole blood predicts the drug, not the illness",
       subtitle = paste("Pooled out-of-fold AUC, 25 repeats of 3-fold CV. Bars are 95% CI,",
                        "bootstrapped over subjects.\nDashed line is chance. Elastic net;",
                        "random forest gives the same ordering."),
       x = "AUC", y = NULL,
       caption = paste("Removing age, sex, smoking, RIN and plate from the genes leaves the lithium",
                       "signal intact and destroys the bipolar one.\nLithium: 152 users vs 74",
                       "non-users, all assessment group A. Bipolar: 74 unmedicated BP1 vs 112",
                       "group-A controls.")) +
  base
ggsave(file.path(outd, "fig1_lithium_vs_bipolar.png"), p1,
       width = 8.2, height = 4.3, dpi = 200, bg = SURFACE)

## ---- Figure 2: why the control group decides the answer ----------------------
f2 <- s[s$cohort %in% c("groupA", "full", "groupA_plate"), ]
CL <- c(full = "All controls\n(74 vs 205)",
        groupA = "Group A only\n(74 vs 112)",
        groupA_plate = "Group A, shared\nplates (44 vs 109)")
f2$coh <- factor(CL[f2$cohort], levels = CL)

p2 <- ggplot(f2, aes(auc, fs, colour = cohort == "full")) +
  geom_vline(xintercept = 0.5, colour = MUTED, linetype = "22", linewidth = 0.4) +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 2, alpha = 0.28) +
  geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("%.2f", auc)), vjust = -1.25, size = 3.0,
            colour = INK, show.legend = FALSE) +
  facet_wrap(~coh, nrow = 1) +
  scale_colour_manual(values = c("TRUE" = ORANGE, "FALSE" = MUTED)) +
  scale_x_continuous(limits = c(0.42, 1.0), breaks = seq(0.5, 0.9, 0.2)) +
  labs(title = "One cohort decision flips the bipolar result",
       subtitle = paste("Same pipeline, same models, same metric. Assessment group B holds 109",
                        "controls and no cases,\nand is predictable from expression at AUC 0.94 --",
                        "so including it lets genes beat the covariates."),
       x = "AUC", y = NULL,
       caption = paste("With all controls, genes beat covariates by +0.10 (CI excludes zero).",
                       "Restricted to group A, the\nsame comparison gives +0.003. The apparent",
                       "signal was recruitment structure.")) +
  base
ggsave(file.path(outd, "fig2_cohort_sensitivity.png"), p2,
       width = 8.6, height = 4.3, dpi = 200, bg = SURFACE)

cat(sprintf("wrote %s\n", file.path(outd, "fig1_lithium_vs_bipolar.png")))
cat(sprintf("wrote %s\n", file.path(outd, "fig2_cohort_sensitivity.png")))
