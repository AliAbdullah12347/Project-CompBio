#!/usr/bin/env Rscript
# ===========================================================================
# 014-effect-sizes
#
# Analysis 009 showed 31/33 modules significantly lower in BP vs controls
# (ctrl-nolith, Bonferroni<0.05, after composition adjustment). But the
# p-values reflect both effect size AND sample size (n=308 ctrl+nolith).
# This analysis computes standardised effect sizes for each module:
#
# 1. Cohen's d from the composition-adjusted model:
#    d = group_coefficient / residual_SD (from lm(ME ~ group + b1+b2+b3+b4))
#    (composition-adjusted Cohen's d, not raw)
#
# 2. Partial η² from the group term in ANOVA:
#    η² = SS_group / SS_total (analogous to R² for the group term alone)
#
# 3. Raw (unadjusted) mean difference and SD within each group.
#
# Pre-specified: same model as analysis 009.
# All 33 modules computed; ranked by ctrl-nolith |d|.
# ===========================================================================
setwd("/Users/Aaylab/Project-CompBio")
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(WGCNA) })
enableWGCNAThreads()
OUT <- "network/analyses/014-effect-sizes/results"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), paste0(...,collapse="")))
set.seed(20261014L)

msg("Loading data ...")
samp <- read.csv("network/analyses/007-darkturquoise-across-groups/results/sample_data.csv",
                 stringsAsFactors = FALSE)
me_all <- read.csv("network/analyses/009-all-modules-across-groups/results/module_eigengenes_all.csv",
                   stringsAsFactors = FALSE)
stopifnot(nrow(me_all) == 460L)
me_cols <- grep("^ME", names(me_all), value = TRUE)

samp$group <- factor(samp$group, levels = c("control","bp_nolith","bp_lith"))

# Load analysis 009 results for p-values
anova09 <- read.csv("network/analyses/009-all-modules-across-groups/results/all_modules_group_anova.csv",
                    stringsAsFactors = FALSE)

msg("Computing effect sizes for ", length(me_cols), " modules ...")
results <- do.call(rbind, lapply(me_cols, function(me_col) {
  mod <- sub("ME", "", me_col)
  y   <- me_all[[me_col]]

  # Fit composition-adjusted model (same as 009)
  df <- data.frame(ME = y, group = samp$group,
                   b1 = samp$b1, b2 = samp$b2, b3 = samp$b3, b4 = samp$b4)
  m  <- lm(ME ~ group + b1 + b2 + b3 + b4, data = df)

  # Residual SD
  res_sd <- sigma(m)

  # Group coefficients (vs control reference)
  cf <- coef(m)
  d_nolith <- cf["groupbp_nolith"] / res_sd   # Cohen's d, ctrl vs nolith
  d_lith   <- cf["groupbp_lith"]   / res_sd   # Cohen's d, ctrl vs lith

  # ANOVA: partial η² for group term
  av   <- anova(m)
  ss_g <- av["group", "Sum Sq"]
  ss_r <- av["Residuals", "Sum Sq"]
  eta2 <- ss_g / (ss_g + ss_r)   # partial η²

  # Raw means and SDs per group
  ctrl_me   <- y[samp$group == "control"]
  nolith_me <- y[samp$group == "bp_nolith"]
  lith_me   <- y[samp$group == "bp_lith"]
  pool_sd   <- sqrt(((length(ctrl_me)-1)*var(ctrl_me) +
                     (length(nolith_me)-1)*var(nolith_me)) /
                    (length(ctrl_me) + length(nolith_me) - 2))
  raw_d_cn  <- (mean(nolith_me) - mean(ctrl_me)) / pool_sd

  # Merge with 009 p-values
  a09 <- anova09[anova09$module == mod, ]

  data.frame(
    module            = mod,
    n_genes           = sum(WGCNA::labels2colors(
      readRDS("network/runs/mod-000-resolution/labels_p12_ds4.rds")) == mod),
    cohens_d_ctrl_nolith = round(d_nolith, 4),
    cohens_d_ctrl_lith   = round(d_lith,   4),
    raw_d_ctrl_nolith    = round(raw_d_cn,  4),
    partial_eta2         = round(eta2, 4),
    mean_ctrl   = round(mean(ctrl_me),   6),
    mean_nolith = round(mean(nolith_me), 6),
    mean_lith   = round(mean(lith_me),   6),
    sd_ctrl     = round(sd(ctrl_me),     6),
    p_bonf_ctrl_nolith = if (nrow(a09) > 0) round(a09$bonf_ctrl_nolith[1], 4) else NA,
    p_bonf_ctrl_lith   = if (nrow(a09) > 0) round(a09$bonf_ctrl_lith[1], 4) else NA,
    stringsAsFactors = FALSE
  )
}))

# Load module labels once for n_genes (fix: already loaded per row above — recompute cleanly)
lab <- readRDS("network/runs/mod-000-resolution/labels_p12_ds4.rds")
cols <- WGCNA::labels2colors(lab)
mod_sizes <- table(cols)
results$n_genes <- as.integer(mod_sizes[results$module])

results_cn <- results[order(results$cohens_d_ctrl_nolith), ]  # most negative first (BP lower)
write.csv(results_cn, file.path(OUT, "effect_sizes_all_modules.csv"), row.names = FALSE)

cat("\n======= All 33 modules: composition-adjusted Cohen's d (ctrl-nolith) =======\n")
cat("Positive d = nolith HIGHER than ctrl; negative = nolith LOWER\n\n")
print(results_cn[, c("module","n_genes","cohens_d_ctrl_nolith","raw_d_ctrl_nolith",
                     "partial_eta2","p_bonf_ctrl_nolith")], row.names = FALSE)

# Summary stats
sig <- results[!is.na(results$p_bonf_ctrl_nolith) & results$p_bonf_ctrl_nolith < 0.05, ]
cat(sprintf("\n%d / %d modules: ctrl-nolith sig (bonf<0.05)\n", nrow(sig), nrow(results)))
cat(sprintf("Median |d| (sig modules): %.3f\n", median(abs(sig$cohens_d_ctrl_nolith), na.rm=TRUE)))
cat(sprintf("Range |d| (sig modules): %.3f – %.3f\n",
            min(abs(sig$cohens_d_ctrl_nolith), na.rm=TRUE),
            max(abs(sig$cohens_d_ctrl_nolith), na.rm=TRUE)))
cat(sprintf("Median partial η²: %.4f\n", median(sig$partial_eta2, na.rm=TRUE)))

# Write summary stats
summ <- data.frame(
  n_sig_bonf05 = nrow(sig),
  pct_negative_d = round(100 * mean(sig$cohens_d_ctrl_nolith < 0), 1),
  median_abs_d = round(median(abs(sig$cohens_d_ctrl_nolith), na.rm=TRUE), 4),
  range_min_abs_d = round(min(abs(sig$cohens_d_ctrl_nolith), na.rm=TRUE), 4),
  range_max_abs_d = round(max(abs(sig$cohens_d_ctrl_nolith), na.rm=TRUE), 4),
  median_eta2 = round(median(sig$partial_eta2, na.rm=TRUE), 4)
)
write.csv(summ, file.path(OUT, "effect_size_summary.csv"), row.names = FALSE)

cat("\nANALYSIS_014_DONE\n")
