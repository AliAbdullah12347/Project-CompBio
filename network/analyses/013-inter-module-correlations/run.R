#!/usr/bin/env Rscript
# ===========================================================================
# 013-inter-module-correlations
#
# Analysis 001 showed median inter-module |r| = 0.822. Finding 009 showed
# 31/33 modules have lower ME in BP vs controls. But does the correlation
# *structure* between modules change across groups, or only the levels?
#
# If the 33x33 correlation matrix is preserved across groups, BP disease
# reduces co-expression uniformly (a "volume knob") without reshuffling
# which genes are co-expressed. If the structure changes, the disease is
# genuinely reorganising gene-gene relationships.
#
# Pre-specified:
# 1. Compute 33x33 inter-module Pearson correlation matrix within each group.
# 2. Summarise: mean |r|, median |r|, fraction r>0.5 within each group.
# 3. Mantel-like: Pearson correlation between upper triangles of ctrl vs
#    nolith, ctrl vs lith, nolith vs lith correlation matrices.
# 4. Permutation test: shuffle sample group labels 1000 times, recompute
#    the ctrl vs nolith Mantel r. Compare to observed.
# ===========================================================================
setwd("/Users/Aaylab/Project-CompBio")
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(WGCNA) })
enableWGCNAThreads()
OUT <- "network/analyses/013-inter-module-correlations/results"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), paste0(...,collapse="")))
set.seed(20261014L)
K_PERM <- 1000L

msg("Loading data ...")
source("network/00_load.R"); d <- load_project(cache = TRUE, verbose = FALSE)
all_genes <- rownames(d$logtpm)
lab  <- readRDS("network/runs/mod-000-resolution/labels_p12_ds4.rds")
colors <- WGCNA::labels2colors(lab); names(colors) <- all_genes

samp <- read.csv("network/analyses/007-darkturquoise-across-groups/results/sample_data.csv",
                 stringsAsFactors = FALSE)
me_all <- read.csv("network/analyses/009-all-modules-across-groups/results/module_eigengenes_all.csv",
                   stringsAsFactors = FALSE)
stopifnot(nrow(me_all) == 460L)
me_cols <- grep("^ME", names(me_all), value = TRUE)
ME <- as.matrix(me_all[, me_cols])  # 460 x 33
grp <- samp$group

# ---- Group-specific ME matrices --------------------------------------------
ME_ctrl   <- ME[grp == "control",   ]  # 234 x 33
ME_nolith <- ME[grp == "bp_nolith", ]  #  74 x 33
ME_lith   <- ME[grp == "bp_lith",   ]  # 152 x 33
msg("n: ctrl=", nrow(ME_ctrl), " nolith=", nrow(ME_nolith), " lith=", nrow(ME_lith))

# ---- 33x33 inter-module correlations per group ----------------------------
compute_cor_matrix <- function(M) cor(M, use = "pairwise.complete.obs")
C_ctrl   <- compute_cor_matrix(ME_ctrl)
C_nolith <- compute_cor_matrix(ME_nolith)
C_lith   <- compute_cor_matrix(ME_lith)
C_all    <- compute_cor_matrix(ME)

# Save matrices
write.csv(as.data.frame(C_ctrl),   file.path(OUT, "cor_ctrl.csv"),   row.names = TRUE)
write.csv(as.data.frame(C_nolith), file.path(OUT, "cor_nolith.csv"), row.names = TRUE)
write.csv(as.data.frame(C_lith),   file.path(OUT, "cor_lith.csv"),   row.names = TRUE)

# ---- Summary statistics per group -----------------------------------------
upper_tri <- function(M) M[upper.tri(M, diag = FALSE)]

summ_stats <- function(M, label) {
  ut <- upper_tri(M)
  data.frame(group = label, n_pairs = length(ut),
             mean_abs_r = round(mean(abs(ut)), 4),
             median_abs_r = round(median(abs(ut)), 4),
             frac_r_gt05 = round(mean(abs(ut) > 0.5), 4),
             frac_r_gt08 = round(mean(abs(ut) > 0.8), 4),
             mean_r = round(mean(ut), 4))
}
summ <- rbind(summ_stats(C_ctrl,   "control"),
              summ_stats(C_nolith, "bp_nolith"),
              summ_stats(C_lith,   "bp_lith"),
              summ_stats(C_all,    "all_460"))
cat("\n======= Inter-module correlation summary =======\n")
print(summ, row.names = FALSE)
write.csv(summ, file.path(OUT, "inter_module_summary.csv"), row.names = FALSE)

# ---- Mantel-like: r between upper triangles of group correlation matrices --
ut_ctrl   <- upper_tri(C_ctrl)
ut_nolith <- upper_tri(C_nolith)
ut_lith   <- upper_tri(C_lith)

r_cn  <- as.numeric(cor(ut_ctrl, ut_nolith))
r_cl  <- as.numeric(cor(ut_ctrl, ut_lith))
r_nl  <- as.numeric(cor(ut_nolith, ut_lith))

cat("\n======= Mantel-like r (upper triangle correlation) =======\n")
cat(sprintf("ctrl vs nolith:   r = %.4f\n", r_cn))
cat(sprintf("ctrl vs lith:     r = %.4f\n", r_cl))
cat(sprintf("nolith vs lith:   r = %.4f\n", r_nl))

# ---- Permutation test for ctrl vs nolith Mantel r -------------------------
msg("Permutation test (K=", K_PERM, "): ctrl vs nolith structure ...")
all_idx <- seq_len(460L)
ctrl_n  <- sum(grp == "control")
nolith_n <- sum(grp == "bp_nolith")
set.seed(20261014L + 1L)
perm_r <- replicate(K_PERM, {
  perm_grp <- sample(all_idx)
  M_c <- ME[perm_grp[1:ctrl_n],   ]
  M_n <- ME[perm_grp[(ctrl_n+1):(ctrl_n+nolith_n)], ]
  as.numeric(cor(upper_tri(compute_cor_matrix(M_c)),
                 upper_tri(compute_cor_matrix(M_n))))
})
perm_p <- mean(perm_r <= r_cn)  # one-sided: does randomisation produce lower r?
# (lower observed r vs null = more structure change than expected by chance)
# Two-sided for reporting
perm_p2 <- mean(abs(perm_r - mean(perm_r)) >= abs(r_cn - mean(perm_r)))
cat(sprintf("\nPermutation p (one-sided, r_obs < random): %.5f\n", perm_p))
cat(sprintf("Permutation p (two-sided):                  %.5f\n", perm_p2))
cat(sprintf("Null mean: %.4f; SD: %.4f; obs: %.4f\n",
            mean(perm_r), sd(perm_r), r_cn))

perm_res <- data.frame(comparison = "ctrl_vs_nolith",
                       observed_r = r_cn,
                       null_mean = mean(perm_r),
                       null_sd = sd(perm_r),
                       perm_p_onesided = perm_p,
                       perm_p_twosided = perm_p2,
                       K = K_PERM)
write.csv(perm_res, file.path(OUT, "mantel_perm.csv"), row.names = FALSE)

# ---- Mantel summary table --------------------------------------------------
mantel_tbl <- data.frame(
  comparison = c("ctrl_vs_nolith","ctrl_vs_lith","nolith_vs_lith"),
  mantel_r   = round(c(r_cn, r_cl, r_nl), 4)
)
write.csv(mantel_tbl, file.path(OUT, "mantel_r.csv"), row.names = FALSE)

# ---- Modules with largest change in inter-module correlations --------------
# For each module, compare its average |r| with all other modules in ctrl vs nolith
msg("Per-module changes in inter-module correlations ...")
mod_names <- colnames(ME)
delta_r <- do.call(rbind, lapply(seq_along(mod_names), function(i) {
  mod <- mod_names[i]
  avg_ctrl   <- mean(abs(C_ctrl[i, -i]))
  avg_nolith <- mean(abs(C_nolith[i, -i]))
  avg_lith   <- mean(abs(C_lith[i, -i]))
  data.frame(module = sub("ME","",mod),
             mean_abs_r_ctrl = round(avg_ctrl, 4),
             mean_abs_r_nolith = round(avg_nolith, 4),
             mean_abs_r_lith = round(avg_lith, 4),
             delta_nolith_ctrl = round(avg_nolith - avg_ctrl, 4),
             delta_lith_ctrl = round(avg_lith - avg_ctrl, 4))
}))
delta_r <- delta_r[order(abs(delta_r$delta_nolith_ctrl), decreasing = TRUE), ]

cat("\n======= Per-module change in mean |r| with other modules (top 10) =======\n")
print(head(delta_r[, c("module","mean_abs_r_ctrl","mean_abs_r_nolith","delta_nolith_ctrl")], 10),
      row.names = FALSE)
write.csv(delta_r, file.path(OUT, "per_module_correlation_change.csv"), row.names = FALSE)

cat("\nANALYSIS_013_DONE\n")
