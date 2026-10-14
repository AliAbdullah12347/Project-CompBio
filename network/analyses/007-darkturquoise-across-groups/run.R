#!/usr/bin/env Rscript
# ===========================================================================
# 007-darkturquoise-across-groups
#
# Analyses 001-006 used only 226 BP samples. Darkturquoise is the ribosomal
# module that drives the lithium-composition PC2 (finding 006). We now extend
# to all 474 samples to ask:
#   (1) Does darkturquoise also track cell composition in controls?
#   (2) Is there a group difference after composition adjustment?
#
# Pre-specified tests: see config.json
# ===========================================================================
setwd("/Users/Aaylab/Project-CompBio")
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(WGCNA) })
enableWGCNAThreads()
OUT <- "network/analyses/007-darkturquoise-across-groups/results"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), paste0(...,collapse="")))
set.seed(20261014L)

msg("Loading project data ...")
source("network/00_load.R"); d <- load_project(cache = TRUE, verbose = FALSE)

all_genes <- rownames(d$logtpm)
msg(ncol(d$logtpm), " samples, ", nrow(d$logtpm), " genes")

# Module colours from mod-000-resolution (33-module ds4 decomposition)
msg("Loading module labels ...")
lab    <- readRDS("network/runs/mod-000-resolution/labels_p12_ds4.rds")
colors <- WGCNA::labels2colors(lab)
names(colors) <- all_genes
stopifnot(length(colors) == length(all_genes))
n_dt <- sum(colors == "darkturquoise")
msg("Darkturquoise: ", n_dt, " genes")

# ------- Compute darkturquoise ME across ALL 474 samples ------------------
msg("Computing darkturquoise ME for all 474 samples ...")
E_all <- t(d$logtpm)   # 474 x 12368
# moduleEigengenes needs the gene-coloured matrix
ME_all <- WGCNA::moduleEigengenes(E_all, colors = colors, excludeGrey = TRUE)$eigengenes
dt_ME <- ME_all[, "MEdarkturquoise"]   # named by rownames(E_all) = sample titles
stopifnot(length(dt_ME) == 474L)

# ------- Metadata ------------------------------------------------------------
md  <- d$meta                                     # 474 rows
grp <- d$groups                                   # ctrl / bp_nolith / bp_lith / NA
il  <- d$ilr[match(md$title, d$ilr$sample), ]   # ILR b1-b4

# Sanity: group counts
grp_tab <- table(grp, useNA = "ifany")
msg("Group counts: ", paste(names(grp_tab), grp_tab, sep="=", collapse="; "))

# Build analysis data frame (all samples with non-missing group and ILR)
ok <- !is.na(grp) & !is.na(il$b1)
msg(sum(ok), " samples with complete data (group + ILR)")
df <- data.frame(
  sample       = md$title[ok],
  group        = grp[ok],
  dt_ME        = dt_ME[ok],
  b1           = il$b1[ok],
  b2           = il$b2[ok],
  b3           = il$b3[ok],
  b4           = il$b4[ok],
  stringsAsFactors = FALSE
)
df$group <- factor(df$group, levels = c("control","bp_nolith","bp_lith"))

write.csv(df, file.path(OUT, "sample_data.csv"), row.names = FALSE)
msg("Sample data written (", nrow(df), " rows)")

# ------- (1) Composition correlation within each group -----------------------
msg("\n--- Pearson r(darkturquoise_ME, b2) by group ---")
r_by_grp <- vapply(levels(df$group), function(g) {
  s <- df[df$group == g, ]
  cor(s$dt_ME, s$b2, use = "pairwise.complete.obs")
}, numeric(1))
r_df <- data.frame(group = names(r_by_grp), r_dt_b2 = round(r_by_grp, 3))
print(r_df, row.names = FALSE)
write.csv(r_df, file.path(OUT, "composition_r_by_group.csv"), row.names = FALSE)

# Also r with b1-b4 across ALL samples
r_all_b <- vapply(c("b1","b2","b3","b4"), function(b) cor(df$dt_ME, df[[b]], use="pairwise"), numeric(1))
cat(sprintf("\nAll-sample r(dt_ME, b1..b4): %s\n",
    paste(names(r_all_b), round(r_all_b,3), sep="=", collapse="; ")))

# ------- (2) Unadjusted group ANOVA -----------------------------------------
msg("\n--- ANOVA: dt_ME ~ group (unadjusted) ---")
m0 <- lm(dt_ME ~ group, data = df)
a0 <- anova(m0)
print(a0)
write.csv(as.data.frame(a0), file.path(OUT, "anova_unadjusted.csv"), row.names = TRUE)

# Group means ± SD
grp_stats <- do.call(rbind, lapply(levels(df$group), function(g) {
  s <- df$dt_ME[df$group == g]
  data.frame(group=g, n=length(s), mean=round(mean(s),4), sd=round(sd(s),4))
}))
print(grp_stats, row.names=FALSE)
write.csv(grp_stats, file.path(OUT, "group_means.csv"), row.names = FALSE)

# ------- (3) Composition-adjusted group ANOVA --------------------------------
msg("\n--- ANOVA: dt_ME ~ group + b1 + b2 + b3 + b4 (composition-adjusted) ---")
m1 <- lm(dt_ME ~ group + b1 + b2 + b3 + b4, data = df)
a1 <- anova(m1)
print(a1)
write.csv(as.data.frame(a1), file.path(OUT, "anova_adjusted.csv"), row.names = TRUE)

# Effect size: partial eta-squared for group
ss_group <- a1["group","Sum Sq"]
ss_resid  <- a1["Residuals","Sum Sq"]
eta2_group_partial <- round(ss_group / (ss_group + ss_resid), 4)
cat(sprintf("\nPartial eta-sq for group (comp-adjusted): %.4f\n", eta2_group_partial))

# R² composition model alone (to see how much composition explains in dt_ME)
m_comp <- lm(dt_ME ~ b1 + b2 + b3 + b4, data = df)
r2_comp <- summary(m_comp)$r.squared
cat(sprintf("R2 of composition alone for dt_ME: %.4f (%.1f%%)\n", r2_comp, 100*r2_comp))

# ------- (4) Pairwise comparisons (Bonferroni 3) ----------------------------
msg("\n--- Pairwise t-tests (Bonferroni n=3) ---")
pairs <- list(
  c("control","bp_nolith"),
  c("control","bp_lith"),
  c("bp_nolith","bp_lith")
)
pwise <- do.call(rbind, lapply(pairs, function(p) {
  a <- df$dt_ME[df$group == p[1]]
  b_ <- df$dt_ME[df$group == p[2]]
  tt <- t.test(a, b_, var.equal = FALSE)
  data.frame(contrast=paste0(p,collapse=" vs "), t=round(tt$statistic,3),
             df=round(tt$parameter,1), p_raw=round(tt$p.value,5),
             p_bonf=round(min(1, tt$p.value * 3), 5),
             diff_means=round(diff(tt$estimate),4), stringsAsFactors=FALSE)
}))
print(pwise, row.names=FALSE)
write.csv(pwise, file.path(OUT, "pairwise_ttest.csv"), row.names = FALSE)

# Composition-adjusted pairwise: residuals from composition model
df$dt_resid <- resid(m_comp)
pwise_adj <- do.call(rbind, lapply(pairs, function(p) {
  a <- df$dt_resid[df$group == p[1]]
  b_ <- df$dt_resid[df$group == p[2]]
  tt <- t.test(a, b_, var.equal = FALSE)
  data.frame(contrast=paste0(p,collapse=" vs "), t=round(tt$statistic,3),
             df=round(tt$parameter,1), p_raw=round(tt$p.value,5),
             p_bonf=round(min(1, tt$p.value * 3), 5),
             diff_resid=round(diff(tt$estimate),4), stringsAsFactors=FALSE)
}))
msg("Composition-adjusted pairwise:")
print(pwise_adj, row.names=FALSE)
write.csv(pwise_adj, file.path(OUT, "pairwise_adjusted.csv"), row.names = FALSE)

cat("\nANALYSIS_007_DONE\n")
