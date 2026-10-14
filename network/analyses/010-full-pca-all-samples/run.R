#!/usr/bin/env Rscript
# ===========================================================================
# 010-full-pca-all-samples
#
# Analyses 001-004 ran PCA on 226 BP samples only; PC2 captured the
# lithium-composition signal. Analysis 009 showed a network-wide downward
# shift in BP vs controls. This analysis runs PCA on ALL 460 samples to ask:
# does PC1 capture the BP disease effect?
#
# Pre-specified: PC1-PC5 variance explained; Pearson r with disease/lithium/
# ILR; permutation test K=1000 for PC1-disease association.
# ===========================================================================
setwd("/Users/Aaylab/Project-CompBio")
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(WGCNA) })
enableWGCNAThreads()
OUT <- "network/analyses/010-full-pca-all-samples/results"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), paste0(...,collapse="")))
set.seed(20261014L)
K_PERM <- 1000L

# Re-use sample_data from analysis 007 (460 samples, ILR)
msg("Loading sample data ...")
samp <- read.csv("network/analyses/007-darkturquoise-across-groups/results/sample_data.csv",
                 stringsAsFactors = FALSE)
samp$group <- factor(samp$group, levels = c("control","bp_nolith","bp_lith"))

# Load ME for all 460 samples from analysis 009
msg("Loading all-sample module eigengenes ...")
me_full <- read.csv("network/analyses/009-all-modules-across-groups/results/module_eigengenes_all.csv",
                    stringsAsFactors = FALSE)
stopifnot(nrow(me_full) == nrow(samp))
me_cols <- grep("^ME", names(me_full), value = TRUE)
me_mat  <- as.matrix(me_full[, me_cols])   # 460 x 33
rownames(me_mat) <- me_full$sample
msg(nrow(me_mat), " samples x ", ncol(me_mat), " modules")

# ------- PCA across all 460 samples ----------------------------------------
msg("PCA (460 samples × 33 modules) ...")
pc <- prcomp(me_mat, scale. = TRUE)
var_exp <- 100 * pc$sdev^2 / sum(pc$sdev^2)
cat("\n======= Variance explained (PC1-PC8) =======\n")
cat(sprintf("PC%d: %.2f%%\n", seq_along(var_exp[1:8]), var_exp[1:8]))
write.csv(data.frame(pc = 1:length(var_exp), var_pct = round(var_exp, 3)),
          file.path(OUT, "variance_explained.csv"), row.names = FALSE)

# ------- Group indicators ---------------------------------------------------
disease  <- as.integer(samp$group %in% c("bp_nolith","bp_lith"))  # 0=ctrl 1=BP
lithium  <- as.integer(samp$group == "bp_lith")                   # 0=ctrl+nolith 1=lith
grp_num  <- as.integer(samp$group) - 1L  # ctrl=0, nolith=1, lith=2

# ------- Pearson r of PC1-PC5 with group variables + ILR ------------------
msg("Pearson correlations with group variables ...")
pc_scores <- pc$x[, 1:8]
corr_results <- do.call(rbind, lapply(1:8, function(i) {
  sc <- pc_scores[, i]
  data.frame(
    pc       = i,
    var_pct  = round(var_exp[i], 2),
    r_disease = round(cor(sc, disease,  use = "pairwise"), 3),
    r_lithium = round(cor(sc, lithium,  use = "pairwise"), 3),
    r_grpnum  = round(cor(sc, grp_num,  use = "pairwise"), 3),
    r_b1      = round(cor(sc, samp$b1,  use = "pairwise"), 3),
    r_b2      = round(cor(sc, samp$b2,  use = "pairwise"), 3),
    r_b3      = round(cor(sc, samp$b3,  use = "pairwise"), 3),
    r_b4      = round(cor(sc, samp$b4,  use = "pairwise"), 3),
    stringsAsFactors = FALSE
  )
}))
cat("\n======= PC correlations with group and composition =======\n")
print(corr_results, row.names = FALSE)
write.csv(corr_results, file.path(OUT, "pc_correlations.csv"), row.names = FALSE)

# ------- Permutation test for PC1 × disease --------------------------------
msg("Permutation test (K=", K_PERM, ") for PC1 × disease ...")
obs_r <- as.numeric(cor(pc_scores[, 1], disease))
set.seed(20261014L + 10L)
perm_r <- replicate(K_PERM, {
  as.numeric(cor(pc_scores[, 1], sample(disease)))
})
perm_p <- mean(abs(perm_r) >= abs(obs_r))
cat(sprintf("\nPC1 r(disease) = %.4f; permutation p = %.5f (K=%d)\n",
            obs_r, perm_p, K_PERM))
write.csv(data.frame(observed_r = obs_r, perm_p = perm_p, K = K_PERM),
          file.path(OUT, "pc1_perm_test.csv"), row.names = FALSE)

# ------- PC1 scores by group (summary stats) --------------------------------
msg("PC1 by group ...")
pc1_grp <- do.call(rbind, lapply(levels(samp$group), function(g) {
  s <- pc_scores[samp$group == g, 1]
  data.frame(group=g, n=length(s), mean=round(mean(s),4), sd=round(sd(s),4),
             median=round(median(s),4))
}))
cat("\n======= PC1 scores by group =======\n")
print(pc1_grp, row.names = FALSE)
write.csv(pc1_grp, file.path(OUT, "pc1_by_group.csv"), row.names = FALSE)

# ------- Compare to BP-only PCA (analysis 002) -----------------------------
msg("For reference: BP-only PC1-PC4 (analysis 002) had Li perm p=0.0015 on PC2.")
cat("\nKey comparison:\n")
cat("  All-sample PCA: PC1 =", round(var_exp[1],1), "%, PC2 =", round(var_exp[2],1), "%\n")
cat("  BP-only PCA (002): PC1 explained majority; PC2 was lithium-associated\n")
cat("  Now including controls — PC1 should capture disease effect if BP globally lower.\n")

# PC loadings (top 10 modules by |loading| on PC1)
cat("\n======= PC1 top module loadings =======\n")
load1 <- sort(abs(pc$rotation[, 1]), decreasing = TRUE)
print(data.frame(module = sub("ME","",names(load1[1:10])),
                 loading = round(pc$rotation[names(load1[1:10]), 1], 4),
                 abs_loading = round(load1[1:10], 4)), row.names = FALSE)

cat("\nANALYSIS_010_DONE\n")
