#!/usr/bin/env Rscript
# ===========================================================================
# samp-043  Eigengene-space structure across groups
#
# Spec: eigengene-space separation by group; outlier detection; hclust;
#       does any cluster track lithium?
#
# Design:
#   - Compute module eigengenes for all 474 samples from base-001 modules
#   - Hierarchical clustering of samples in eigengene space
#   - PCA in eigengene space
#   - Silhouette score for group separation
#   - Outlier detection: samples > 2.5 SD from group centroid
#   - Does sample clustering recover groups (ctrl/nolith/lith)?
#
# Dependency: requires base-001 reference cache.
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "samp-043"
stopifnot(id == "samp-043")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

run_dir <- file.path("network/runs", id)
cfg_dir <- "network/config"
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg_dir, showWarnings = FALSE, recursive = TRUE)

REF_CACHE_KEY <- "control_all_bicor_signed_p12"
REF_CACHE     <- file.path("network/cache", paste0("ref_", REF_CACHE_KEY, ".rds"))
SEED          <- 20261009L
OUTLIER_SD    <- 2.5

if (!file.exists(REF_CACHE)) {
  msg("DEPENDENCY MISSING: ", REF_CACHE, " — base-001 must complete first")
  quit(save = "no", status = 1L)
}

# ---- config ----------------------------------------------------------------
cfg_path <- file.path(cfg_dir, paste0(id, ".json"))
if (!file.exists(cfg_path)) {
  scl <- function(v) {
    if (is.logical(v)) return(if (isTRUE(v)) "true" else "false")
    if (is.integer(v) || is.numeric(v)) return(as.character(v[1]))
    paste0('"', gsub('"', '\\"', as.character(v[1]), fixed = TRUE), '"')
  }
  jsn <- function(nm, v) paste0('  "', nm, '": ', scl(v))
  sha <- tryCatch(trimws(system("git rev-parse --short HEAD 2>/dev/null", intern=TRUE)[1]),
                  error = function(e) "unknown")
  cfg <- list(id=id, family="samplestruct",
    spec="eigengene-space separation by group; outlier detection; hclust; does any cluster track lithium?",
    ref_cache=REF_CACHE_KEY,
    outlier_sd_threshold=OUTLIER_SD,
    cluster_method="complete linkage on euclidean distance in eigengene space",
    pca_components=10L,
    samples="all 474",
    seed=SEED,
    git_sha=if(length(sha)&&!is.na(sha))sha else "unknown",
    r_version=paste(R.version$major,R.version$minor,sep="."),
    wgcna_version=as.character(packageVersion("WGCNA")))
  lines <- mapply(jsn, names(cfg), cfg, SIMPLIFY=TRUE)
  writeLines(c("{", paste(lines, collapse=",\n"), "}"), cfg_path)
  msg("config written: ", cfg_path)
}

suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors = FALSE)

source("network/00_load.R")
d <- load_project()

ref_net      <- readRDS(REF_CACHE)
moduleColors <- ref_net$moduleColors
modules      <- sort(unique(moduleColors[moduleColors != "grey"]))
msg(length(modules), " modules (excl grey)")

# ---- module eigengenes for all 474 samples ---------------------------------
E   <- t(d$logtpm)
MEs <- moduleEigengenes(E, colors = moduleColors, excludeGrey = TRUE)$eigengenes
rownames(MEs) <- colnames(d$logtpm)
msg("eigengenes: ", nrow(MEs), " x ", ncol(MEs))
write.csv(MEs, file.path(run_dir, "module_eigengenes_474.csv"))

# group vector
grps <- as.character(d$groups)
grps[is.na(grps)] <- "other"

# ---- hierarchical clustering in eigengene space ----------------------------
set.seed(SEED)
dmat <- dist(MEs, method="euclidean")
hc   <- hclust(dmat, method="complete")
hc_k3 <- cutree(hc, k=3)
tab  <- table(hc_k3, grps)
msg("hclust k=3 vs groups:\n", paste(capture.output(print(tab)), collapse="\n"))
write.csv(data.frame(sample=rownames(MEs), hc_k3=hc_k3, group=grps),
          file.path(run_dir, "sample_hclust.csv"), row.names=FALSE)

# ---- PCA in eigengene space ------------------------------------------------
pc <- prcomp(MEs, center=TRUE, scale.=FALSE)
pc_var <- (pc$sdev^2 / sum(pc$sdev^2)) * 100
n_pc <- min(10L, ncol(MEs))
pc_df <- data.frame(sample=rownames(MEs), group=grps,
                    pc$x[, seq_len(n_pc), drop=FALSE])
write.csv(pc_df, file.path(run_dir, "eigengene_pca.csv"), row.names=FALSE)
write.csv(data.frame(PC=seq_len(n_pc), var_pct=pc_var[seq_len(n_pc)]),
          file.path(run_dir, "eigengene_pca_variance.csv"), row.names=FALSE)
msg("PC1 explains ", round(pc_var[1],1), "% | PC2 explains ", round(pc_var[2],1), "%")

# ---- group separation in PC1-PC2 -------------------------------------------
grp_levels <- c("control","bp_nolith","bp_lith")
grp_stats  <- do.call(rbind, lapply(grp_levels, function(g) {
  idx <- grps == g
  if (!any(idx)) return(NULL)
  data.frame(group=g, n=sum(idx),
             pc1_mean=mean(pc$x[idx,1]), pc1_sd=sd(pc$x[idx,1]),
             pc2_mean=mean(pc$x[idx,2]), pc2_sd=sd(pc$x[idx,2]))
}))
write.csv(grp_stats, file.path(run_dir, "group_centroids_pc12.csv"), row.names=FALSE)

# ---- outlier detection: distance from group centroid -----------------------
me_mat <- as.matrix(MEs)
outlier_rows <- list()
for (g in grp_levels) {
  idx <- grps == g
  if (!any(idx)) next
  cent  <- colMeans(me_mat[idx,, drop=FALSE])
  dists <- sqrt(rowSums((me_mat[idx,, drop=FALSE] - cent)^2))
  thr   <- mean(dists) + OUTLIER_SD * sd(dists)
  outs  <- rownames(me_mat)[idx][dists > thr]
  if (length(outs)) msg(length(outs), " outliers in ", g, ": ", paste(outs, collapse=","))
  outlier_rows[[g]] <- data.frame(sample=rownames(me_mat)[idx], group=g,
                                  dist_from_centroid=dists,
                                  is_outlier=dists > thr)
}
outlier_df <- do.call(rbind, outlier_rows)
write.csv(outlier_df, file.path(run_dir,"sample_outliers.csv"), row.names=FALSE)
msg(sum(outlier_df$is_outlier), " total outlier samples across all groups")

# ---- ANOVA per eigengene: does any module separate groups? -----------------
anova_res <- lapply(colnames(MEs), function(me_col) {
  fit <- aov(MEs[[me_col]] ~ factor(grps))
  pv  <- summary(fit)[[1]][["Pr(>F)"]][1]
  data.frame(module=sub("^ME","",me_col), anova_pval=pv)
})
anova_df <- do.call(rbind, anova_res)
anova_df$padj <- p.adjust(anova_df$anova_pval, method="BH")
write.csv(anova_df, file.path(run_dir,"eigengene_anova_by_group.csv"), row.names=FALSE)
msg(sum(anova_df$padj<0.05, na.rm=TRUE), " modules significantly differ by group (BH 5%)")

msg("=== ", id, " COMPLETE ===")
