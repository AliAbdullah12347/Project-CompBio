#!/usr/bin/env Rscript
# ===========================================================================
# 009-all-modules-across-groups
#
# Extends analyses 007-008 to all 33 modules. For each module eigengene:
# ANOVA (ME ~ group + b1+b2+b3+b4), BH across modules; pairwise by Bonferroni.
# ===========================================================================
setwd("/Users/Aaylab/Project-CompBio")
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(WGCNA) })
enableWGCNAThreads()
OUT <- "network/analyses/009-all-modules-across-groups/results"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), paste0(...,collapse="")))
set.seed(20261014L)

msg("Loading project data ...")
source("network/00_load.R"); d <- load_project(cache = TRUE, verbose = FALSE)
all_genes <- rownames(d$logtpm)
msg(ncol(d$logtpm), " samples, ", nrow(d$logtpm), " genes")

# Module colours
lab    <- readRDS("network/runs/mod-000-resolution/labels_p12_ds4.rds")
colors <- WGCNA::labels2colors(lab)
names(colors) <- all_genes
modules <- sort(unique(colors[colors != "grey"]))
msg(length(modules), " non-grey modules")

# Compute ME for ALL 474 samples
msg("Computing module eigengenes for all 474 samples ...")
E_all <- t(d$logtpm)
ME_all <- WGCNA::moduleEigengenes(E_all, colors = colors, excludeGrey = TRUE)$eigengenes
stopifnot(nrow(ME_all) == 474L)

# Metadata + ILR
md  <- d$meta
grp <- d$groups
il  <- d$ilr[match(md$title, d$ilr$sample), ]

# Filter: complete ILR + non-NA group
ok <- !is.na(grp) & !is.na(il$b1)
msg(sum(ok), " samples with complete data")
grp_f <- factor(grp[ok], levels = c("control","bp_nolith","bp_lith"))
ME_ok <- ME_all[ok, ]
il_ok <- il[ok, ]

write.csv(data.frame(sample=md$title[ok], group=grp_f,
                     ME_ok[,paste0("ME",modules)], il_ok[,c("b1","b2","b3","b4")]),
          file.path(OUT, "module_eigengenes_all.csv"), row.names=FALSE)

# ------- Per-module ANOVA: ME ~ group + b1+b2+b3+b4 -----------------------
msg("Running ANOVA for each module ...")
results <- do.call(rbind, lapply(modules, function(m) {
  y    <- ME_ok[[paste0("ME", m)]]
  df_m <- data.frame(y = y, group = grp_f, il_ok[, c("b1","b2","b3","b4")])
  m_add  <- lm(y ~ group + b1 + b2 + b3 + b4, data = df_m)
  m_null <- lm(y ~ b1 + b2 + b3 + b4, data = df_m)
  a      <- anova(m_null, m_add)
  r2_comp <- summary(m_null)$r.squared

  # Pairwise contrasts (Bonferroni n=3 within module)
  res <- coef(m_add)
  vc  <- vcov(m_add)
  b1m <- mean(df_m$b1); b2m <- mean(df_m$b2)
  b3m <- mean(df_m$b3); b4m <- mean(df_m$b4)
  # Contrast vectors (same structure as analysis 008 but without interaction term)
  # Intercept, groupbp_nolith, groupbp_lith, b1, b2, b3, b4
  c_ctrl   <- c(1, 0, 0, b1m, b2m, b3m, b4m)
  c_nolith <- c(1, 1, 0, b1m, b2m, b3m, b4m)
  c_lith   <- c(1, 0, 1, b1m, b2m, b3m, b4m)
  contrasts <- list(
    ctrl_nolith = c_nolith - c_ctrl,
    ctrl_lith   = c_lith   - c_ctrl,
    nolith_lith = c_lith   - c_nolith
  )
  praw <- vapply(contrasts, function(cc) {
    d_v <- as.numeric(res %*% cc)
    se  <- sqrt(as.numeric(t(cc) %*% vc %*% cc))
    t_v <- d_v / se
    2 * pt(-abs(t_v), df = m_add$df.residual)
  }, numeric(1))
  pbonf <- pmin(praw * 3, 1)

  data.frame(
    module         = m,
    F_group        = round(a$F[2], 3),
    p_group        = round(a$"Pr(>F)"[2], 6),
    R2_composition = round(r2_comp, 4),
    p_ctrl_nolith  = round(praw["ctrl_nolith"], 5),
    p_ctrl_lith    = round(praw["ctrl_lith"], 5),
    p_nolith_lith  = round(praw["nolith_lith"], 5),
    bonf_ctrl_nolith = round(pbonf["ctrl_nolith"], 5),
    bonf_ctrl_lith   = round(pbonf["ctrl_lith"], 5),
    bonf_nolith_lith = round(pbonf["nolith_lith"], 5),
    stringsAsFactors = FALSE
  )
}))
results$FDR_group <- p.adjust(results$p_group, "BH")
results <- results[order(results$p_group), ]

cat("\n======= Group ANOVA (sorted by p_group, top 15) =======\n")
show <- results[, c("module","F_group","p_group","FDR_group","bonf_ctrl_lith","bonf_nolith_lith")]
print(head(show, 15), row.names = FALSE)

sig_group <- results[!is.na(results$FDR_group) & results$FDR_group < 0.05, ]
cat(sprintf("\nModules with FDR < 0.05 for group: %d / %d\n", nrow(sig_group), nrow(results)))
if (nrow(sig_group) > 0) {
  cat("  Surviving: ", paste(sig_group$module, collapse=", "), "\n")
}

ctrl_lith_sig <- results[!is.na(results$bonf_ctrl_lith) & results$bonf_ctrl_lith < 0.05, ]
nolith_lith_sig <- results[!is.na(results$bonf_nolith_lith) & results$bonf_nolith_lith < 0.05, ]
cat(sprintf("Modules with bonf ctrl-lith < 0.05: %d — %s\n",
    nrow(ctrl_lith_sig), paste(ctrl_lith_sig$module, collapse=", ")))
cat(sprintf("Modules with bonf nolith-lith < 0.05: %d — %s\n",
    nrow(nolith_lith_sig), paste(nolith_lith_sig$module, collapse=", ")))

write.csv(results, file.path(OUT, "all_modules_group_anova.csv"), row.names = FALSE)
cat("\nANALYSIS_009_DONE\n")
