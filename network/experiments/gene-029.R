#!/usr/bin/env Rscript
# ===========================================================================
# gene-029  Gene subset sensitivity: explicit non-DEGs
#
# Spec: ref=control;genes=explicit_non_DEGs
#
# Purpose:
#   Tests whether network disruption in DEG-excluded genes is also present,
#   indicating a broad effect rather than signal limited to DE genes.
#   If non-DEGs show similar preservation changes as the full set, the
#   differential expression is not the sole driver of module disruption.
#
# Design:
#   - Load significant_genes.csv from de_final/results/.
#   - Exclude ALL WB_LI and WB_BPD significant DEGs.
#   - Build reference modules on non-DEG genes, control group.
#   - Preservation in bp_nolith (500 perms) and bp_lith (20 draws x 50 perms).
#   - Otherwise identical to base-001 (power 12, bicor, signed).
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "gene-029"
stopifnot(id == "gene-029")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

run_dir   <- file.path("network/runs",   id)
cfg_dir   <- "network/config"
cache_dir <- "network/cache"
for (d_ in c(run_dir, cfg_dir, cache_dir))
  dir.create(d_, showWarnings = FALSE, recursive = TRUE)

POWER          <- 12L
POWER_MIN      <-  4L
NET_TYPE       <- "signed"
COR_FNC        <- "bicor"
MIN_MODULE     <- 30L
MERGE_THRESH   <- 0.25
DRAWS          <- 20L
PERMS_PER_DRAW <- 50L
NOLITH_PERMS   <- 500L
SUBSAMPLE_N    <- 74L
SEED_BASE      <- 20261009L + 29000L
REF_GROUP      <- "control"
CACHE_KEY      <- "control_nonDEGs_bicor_signed_p12"

# ---- load DEG file (before config so counts are known) ----------------------
DEG_FILE <- "de_final/results/significant_genes.csv"
if (!file.exists(DEG_FILE)) stop("DEG file not found: ", DEG_FILE)
msg("loading DEG file: ", DEG_FILE)
deg_df <- read.csv(DEG_FILE, stringsAsFactors = FALSE)
msg("DEG file loaded: ", nrow(deg_df), " rows, columns: ",
    paste(colnames(deg_df), collapse = ", "))

# Extract WB_LI and WB_BPD significant DEGs
WB_LI_DEGs  <- deg_df$gene[deg_df$comparison == "WB_LI"  & deg_df$significant == TRUE]
WB_BPD_DEGs <- deg_df$gene[deg_df$comparison == "WB_BPD" & deg_df$significant == TRUE]
all_DEGs     <- union(WB_LI_DEGs, WB_BPD_DEGs)
msg("WB_LI significant DEGs: ", length(WB_LI_DEGs))
msg("WB_BPD significant DEGs: ", length(WB_BPD_DEGs))
msg("Union of DEGs to exclude: ", length(all_DEGs))

suppressPackageStartupMessages({ library(WGCNA); library(fastcluster) })
enableWGCNAThreads()
options(stringsAsFactors = FALSE)

source("network/00_load.R")
d <- load_project()

# Non-DEGs: all genes in expression matrix NOT in any DEG set
keep_ens    <- setdiff(rownames(d$logtpm), all_DEGs)
n_genes     <- length(keep_ens)
n_non_deg   <- n_genes
n_deg_excluded <- length(intersect(all_DEGs, rownames(d$logtpm)))
msg("genes in expression matrix: ", nrow(d$logtpm))
msg("DEGs matched in expression set (excluded): ", n_deg_excluded)
msg("non-DEG genes kept: ", n_non_deg)
stopifnot(n_genes > 0)

# ---- config (written after data load so n_genes is known) -------------------
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
  cfg <- list(id=id, family="geneset",
    spec="ref=control;genes=explicit_non_DEGs",
    note="Non-DEG genes: tests whether network disruption in DEG-excluded genes is also present, indicating a broad effect rather than signal limited to DE genes",
    draws_spec=100L, draws_implemented=DRAWS,
    draws_revision_note="Same as base-001: revised to 20 (METHODS.md 2026-10-09)",
    deg_file=DEG_FILE,
    deg_comparisons_excluded="WB_LI,WB_BPD",
    deg_filter="significant==TRUE",
    n_WB_LI_DEGs=length(WB_LI_DEGs),
    n_WB_BPD_DEGs=length(WB_BPD_DEGs),
    n_deg_excluded=n_deg_excluded,
    n_non_deg=n_non_deg,
    n_genes=n_genes,
    corFnc=COR_FNC, networkType=NET_TYPE, power=POWER,
    power_source="WGCNA_FAQ_signed_n_gt_40",
    minModuleSize=MIN_MODULE, mergeCutHeight=MERGE_THRESH,
    maxBlockSize=15000L,
    draws_lith=DRAWS, perms_per_lith_draw=PERMS_PER_DRAW,
    nolith_perms=NOLITH_PERMS, subsample_n=SUBSAMPLE_N,
    seed_base=SEED_BASE, git_sha=if(length(sha)&&!is.na(sha))sha else "unknown",
    r_version=paste(R.version$major, R.version$minor, sep="."),
    wgcna_version=as.character(packageVersion("WGCNA")))
  lines <- mapply(jsn, names(cfg), cfg, SIMPLIFY = TRUE)
  writeLines(c("{", paste(lines, collapse = ",\n"), "}"), cfg_path)
  msg("config written: ", cfg_path)
}

# ---- helper to subset expression --------------------------------------------
expr_subset <- function(d, group, ens_keep, residualise = NULL) {
  E <- expr_for(d, group = group, residualise = residualise)
  E[, ens_keep[ens_keep %in% colnames(E)], drop = FALSE]
}

# ---- scale-free scan --------------------------------------------------------
r2_path <- file.path(run_dir, "sft_r2_curve.csv")
if (!file.exists(r2_path)) {
  msg("pickSoftThreshold scan on non-DEG gene set (", n_genes, " genes) ...")
  tmp_ref <- expr_subset(d, REF_GROUP, keep_ens)
  set.seed(SEED_BASE)
  sft <- pickSoftThreshold(tmp_ref, powerVector = 1:20,
                           networkType = NET_TYPE, corFnc = COR_FNC, verbose = 0)
  est <- sft$powerEstimate
  msg("powerEstimate=", if (is.na(est)) "NA" else est)
  r2df <- sft$fitIndices; r2df$power_used <- POWER; r2df$powerEstimate <- est
  r2df$n_genes <- n_genes
  write.csv(r2df, r2_path, row.names = FALSE)
  msg("R2 curve saved (max=", round(max(r2df$SFT.R.sq, na.rm = TRUE), 3), ")")
  rm(tmp_ref, sft, r2df); gc()
}

# ---- reference network ------------------------------------------------------
ref_cache <- file.path(cache_dir, paste0("ref_", CACHE_KEY, ".rds"))
if (file.exists(ref_cache)) {
  msg("loading cached reference: ", ref_cache)
  ref_net <- readRDS(ref_cache); moduleColors <- ref_net$moduleColors
  msg(ref_net$n_modules, " modules (non-DEG genes)")
} else {
  msg("building reference network (non-DEGs, power ", POWER, ") ...")
  tmp_ref <- expr_subset(d, REF_GROUP, keep_ens)
  set.seed(SEED_BASE)
  bwm <- blockwiseModules(tmp_ref, power = POWER,
    networkType = NET_TYPE, TOMType = NET_TYPE, corType = COR_FNC,
    minModuleSize = MIN_MODULE, mergeCutHeight = MERGE_THRESH,
    maxBlockSize = 15000L, numericLabels = FALSE, saveTOMs = FALSE, verbose = 3)
  moduleColors <- bwm$colors
  n_mods <- length(unique(moduleColors)) - 1L
  msg(n_mods, " modules detected")
  ref_net <- list(moduleColors = moduleColors, n_modules = n_mods,
                  module_sizes = table(moduleColors), power = POWER,
                  cache_key = CACHE_KEY, n_genes = n_genes,
                  n_non_deg = n_non_deg, n_deg_excluded = n_deg_excluded)
  saveRDS(ref_net, ref_cache)
  write.csv(data.frame(n_genes_total = nrow(d$logtpm),
                       n_deg_excluded = n_deg_excluded,
                       n_non_deg = n_non_deg),
            file.path(run_dir, "gene_selection_summary.csv"), row.names = FALSE)
  write.csv(as.data.frame(ref_net$module_sizes),
            file.path(run_dir, "reference_module_sizes.csv"), row.names = FALSE)
  rm(bwm, tmp_ref); gc()
}

refExpr <- expr_subset(d, REF_GROUP, keep_ens)

# Pre-compute all bp_lith residuals at once
lithAllFull <- expr_for(d, group = "bp_lith")
lithAll     <- lithAllFull[, keep_ens[keep_ens %in% colnames(lithAllFull)], drop = FALSE]
lith_sids   <- rownames(lithAll)
stopifnot(length(lith_sids) == 152L)
rm(lithAllFull); gc()

extract_mp <- function(mp) {
  obs  <- mp$preservation$observed[["ref.Ref"]][["inColumnsAlsoPresentIn.Test"]]
  zsco <- mp$preservation$Z[["ref.Ref"]][["inColumnsAlsoPresentIn.Test"]]
  if (is.null(obs) || is.null(zsco)) stop("unexpected modulePreservation structure")
  res <- data.frame(module = rownames(obs), Zsummary = zsco[, "Zsummary.pres"],
                    medianRank = obs[, "medianRank.pres"], stringsAsFactors = FALSE)
  res[res$module != "grey", ]
}

# ---- bp_nolith --------------------------------------------------------------
nolith_ckpt <- file.path(run_dir, "nolith_preservation.rds")
if (!file.exists(nolith_ckpt)) {
  msg("bp_nolith (non-DEGs, ", NOLITH_PERMS, " perms) ...")
  nolithExpr <- expr_subset(d, "bp_nolith", keep_ens)
  perm_file  <- file.path(run_dir, "permStats_nolith.RData")
  mp_nolith  <- modulePreservation(
    multiData  = list(Ref = list(data = refExpr), Test = list(data = nolithExpr)),
    multiColor = list(Ref = moduleColors), dataIsExpr = TRUE,
    networkType = NET_TYPE, corFnc = COR_FNC,
    nPermutations = NOLITH_PERMS, randomSeed = SEED_BASE + 1000L,
    savePermutedStatistics = TRUE, loadPermutedStatistics = FALSE,
    permutedStatisticsFile = perm_file, verbose = 1, indent = 2)
  nr <- extract_mp(mp_nolith); nr$group <- "bp_nolith"; nr$draw <- NA_integer_
  nr$nperms <- NOLITH_PERMS; nr$n_test <- nrow(nolithExpr)
  saveRDS(nr, nolith_ckpt)
  msg("bp_nolith done: Z[", round(min(nr$Zsummary), 2), ", ", round(max(nr$Zsummary), 2), "]")
  rm(nolithExpr, mp_nolith); gc()
} else { msg("bp_nolith exists, skipping") }

# ---- bp_lith draws ----------------------------------------------------------
for (draw in seq_len(DRAWS)) {
  ckpt <- file.path(run_dir, sprintf("lith_draw_%03d.rds", draw))
  if (file.exists(ckpt)) { msg("draw ", draw, " exists"); next }
  msg("draw ", draw, "/", DRAWS)
  set.seed(SEED_BASE + draw)
  sel_sids <- sample(lith_sids, SUBSAMPLE_N, replace = FALSE)
  testExpr <- lithAll[sel_sids, , drop = FALSE]
  perm_file <- file.path(run_dir, sprintf("permStats_lith_draw_%03d.RData", draw))
  mp_draw <- tryCatch(
    modulePreservation(
      multiData  = list(Ref = list(data = refExpr), Test = list(data = testExpr)),
      multiColor = list(Ref = moduleColors), dataIsExpr = TRUE,
      networkType = NET_TYPE, corFnc = COR_FNC,
      nPermutations = PERMS_PER_DRAW, randomSeed = SEED_BASE + draw,
      savePermutedStatistics = TRUE, loadPermutedStatistics = FALSE,
      permutedStatisticsFile = perm_file, verbose = 1, indent = 2),
    error = function(e) { msg("ERROR draw ", draw, ": ", conditionMessage(e)); NULL })
  if (is.null(mp_draw)) { rm(testExpr); gc(); next }
  dr <- extract_mp(mp_draw); dr$group <- "bp_lith"; dr$draw <- draw
  dr$nperms <- PERMS_PER_DRAW; dr$n_test <- SUBSAMPLE_N
  dr$samples <- paste(sel_sids, collapse = ",")
  saveRDS(dr, ckpt)
  msg("draw ", draw, " done: Z[", round(min(dr$Zsummary), 2), ", ", round(max(dr$Zsummary), 2), "]")
  rm(testExpr, mp_draw, dr); gc()
}

# ---- aggregate --------------------------------------------------------------
done_nolith <- file.exists(nolith_ckpt)
done_lith <- vapply(seq_len(DRAWS), function(i)
  file.exists(file.path(run_dir, sprintf("lith_draw_%03d.rds", i))), logical(1))
if (!done_nolith || sum(done_lith) < DRAWS) {
  msg("not complete: nolith=", done_nolith, " lith=", sum(done_lith), "/", DRAWS)
  quit(save = "no", status = 1L)
}
msg("aggregating ...")
nolith_res <- readRDS(nolith_ckpt)
write.csv(nolith_res[, c("group", "module", "Zsummary", "medianRank", "draw", "nperms", "n_test")],
          file.path(run_dir, "results_nolith.csv"), row.names = FALSE)
lith_all <- do.call(rbind, lapply(seq_len(DRAWS), function(i) {
  x <- readRDS(file.path(run_dir, sprintf("lith_draw_%03d.rds", i)))
  x[, c("group", "module", "Zsummary", "medianRank", "draw", "nperms", "n_test")]
}))
write.csv(lith_all, file.path(run_dir, "results_lith_alldraws.csv"), row.names = FALSE)
modules <- unique(lith_all$module)
lith_summary <- do.call(rbind, lapply(modules, function(m) {
  s <- lith_all[lith_all$module == m, ]
  data.frame(group = "bp_lith", module = m, n_draws = nrow(s),
             Zsummary_mean   = mean(s$Zsummary,   na.rm = TRUE),
             Zsummary_median = median(s$Zsummary,  na.rm = TRUE),
             Zsummary_sd     = sd(s$Zsummary,      na.rm = TRUE),
             medianRank_mean = mean(s$medianRank,  na.rm = TRUE),
             stringsAsFactors = FALSE)
}))
write.csv(lith_summary, file.path(run_dir, "results_lith_summary.csv"), row.names = FALSE)
msg("=== ", id, " COMPLETE ===")
