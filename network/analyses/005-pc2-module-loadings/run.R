#!/usr/bin/env Rscript
# ===========================================================================
# 005-pc2-module-loadings
#
# Which of the 33 modules drive PC2, the lithium-associated dimension?
#
# Analysis 002 showed PC2 of the 33-module eigengene matrix associates with
# lithium status (permutation p=0.0015) and is 73.3% cell composition.
# This analysis identifies which modules have the strongest PC2 loadings,
# giving the biological locus of the lithium-composition signal.
#
# Pre-specified threshold: |loading| > 0.3 (1.7× the expected null magnitude
# for an uncorrelated eigenvector of 33 variables: 1/sqrt(33) ≈ 0.174).
#
# Also reports per-module Pearson r with lithium status and ILR b1-b4.
# ===========================================================================
setwd("/Users/Aaylab/Project-CompBio")
suppressPackageStartupMessages({ library(WGCNA) })
options(stringsAsFactors = FALSE); enableWGCNAThreads()
OUT <- "network/analyses/005-pc2-module-loadings/results"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), paste0(...,collapse="")))
set.seed(20261014L)   # same seed as analysis 002

msg("Loading data ...")
source("network/00_load.R"); d <- load_project(cache = TRUE, verbose = FALSE)

# BP1 samples only (same as 001-004)
bp <- !is.na(d$groups) & d$groups %in% c("bp_nolith","bp_lith")
stopifnot(sum(bp) == 226L)
E   <- t(d$logtpm[, bp])       # 226 x 12368 (samples x genes)
md  <- d$meta[bp, ]
il  <- d$ilr[match(md$title, d$ilr$sample), c("b1","b2","b3","b4")]
lith <- as.numeric(md$lithium)

# Load module labels
msg("Computing module eigengenes (33 modules) ...")
lab <- readRDS("network/runs/mod-000-resolution/labels_p12_ds4.rds")
colors <- labels2colors(lab)
me  <- moduleEigengenes(E, colors = colors, excludeGrey = TRUE)$eigengenes
# me is 226 x 33 (rows=samples, cols=ME<module>)
msg(ncol(me), " module eigengenes")

# PCA — identical to analysis 002
msg("PCA ...")
pc <- prcomp(me, scale. = TRUE)

# PC2 loadings (which modules contribute to the lithium dimension?)
THRESHOLD <- 0.3  # pre-specified
loadings_pc2 <- pc$rotation[, 2]   # named vector: module -> loading
loadings_df  <- data.frame(
  module         = sub("^ME", "", names(loadings_pc2)),
  loading_pc2    = round(loadings_pc2, 4),
  abs_loading    = round(abs(loadings_pc2), 4),
  above_threshold = abs(loadings_pc2) > THRESHOLD,
  stringsAsFactors = FALSE
)
loadings_df <- loadings_df[order(-loadings_df$abs_loading), ]
msg(sum(loadings_df$above_threshold), " modules with |loading| > ", THRESHOLD)

cat("\n======= PC2 LOADINGS (top 15) =======\n")
print(head(loadings_df, 15), row.names = FALSE)
write.csv(loadings_df, file.path(OUT, "pc2_loadings.csv"), row.names = FALSE)

# Per-module Pearson r with lithium, b1-b4, and PC2 score
msg("Computing module-trait correlations ...")
pc2_scores <- pc$x[, 2]   # 226 PC2 scores (one per BP1 sample)

module_names <- colnames(me)
corr_df <- do.call(rbind, lapply(module_names, function(mn) {
  eg <- me[[mn]]
  data.frame(
    module       = sub("^ME", "", mn),
    r_lith       = round(cor(eg, lith, use = "pairwise"), 3),
    r_b1         = round(cor(eg, il$b1, use = "pairwise"), 3),
    r_b2         = round(cor(eg, il$b2, use = "pairwise"), 3),
    r_b3         = round(cor(eg, il$b3, use = "pairwise"), 3),
    r_b4         = round(cor(eg, il$b4, use = "pairwise"), 3),
    r_pc2        = round(cor(eg, pc2_scores, use = "pairwise"), 3),
    loading_pc2  = round(loadings_pc2[mn], 4),
    stringsAsFactors = FALSE
  )
}))
# Sort by |r_pc2| (highest correlation with the lithium dimension first)
corr_df <- corr_df[order(-abs(corr_df$r_pc2)), ]

cat("\n======= MODULE-TRAIT CORRELATIONS (top 15 by |r_pc2|) =======\n")
print(head(corr_df, 15), row.names = FALSE)
write.csv(corr_df, file.path(OUT, "module_trait_correlations.csv"), row.names = FALSE)

# Summarise: which modules above threshold associate most with lithium?
top_modules <- loadings_df$module[loadings_df$above_threshold]
if (length(top_modules) > 0) {
  cat(sprintf("\n== Top modules (|loading_pc2| > %.1f) ==\n", THRESHOLD))
  sub_corr <- corr_df[corr_df$module %in% top_modules, ]
  print(sub_corr[, c("module","r_lith","r_b1","r_b2","r_b3","r_b4","loading_pc2")],
        row.names = FALSE)
} else {
  cat("\nNo modules above threshold (", THRESHOLD, ").\n")
}

# PC2 variance explained
cat(sprintf("\nPC2 explains %.1f%% of ME variance. PC1: %.1f%%.\n",
    100*pc$sdev[2]^2/sum(pc$sdev^2),
    100*pc$sdev[1]^2/sum(pc$sdev^2)))

cat("\nANALYSIS_005_DONE\n")
