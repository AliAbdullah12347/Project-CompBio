#!/usr/bin/env Rscript
# ===========================================================================
# inp-024  voom-weighted expression as WGCNA input
#
# Spec: "ref=control;input=voom_weighted_from_counts;draws=100;perms_per_draw=50;
#        CHECKPOINT PER DRAW (spans sessions)"
#
# Design:
#   - Voom-normalise raw counts (DGEList -> calcNormFactors(TMM) -> voom).
#   - Use voom log2-CPM (v$E, transposed) as input to WGCNA; voom precision
#     weights are NOT carried into WGCNA (no approved WGCNA function supports
#     weighted correlation; see config note).
#   - Build reference modules on voom-normalised control.
#   - Test preservation in voom-normalised bp_nolith (500 perms, one run).
#   - Test preservation in voom-normalised bp_lith: pre-compute voom on all
#     152 samples; subsample residuals (rows) per draw; 20 draws x 50 perms.
#   - Draws=20 (revised from spec 100; METHODS.md 2026-10-09).
#   - CHECKPOINT PER DRAW.
#
# Limitation documented in config:
#   Voom precision weights account for mean-variance trend in count data but
#   are NOT used in WGCNA bicor. Benefit over inp-021 (log2TPM): TMM library-
#   size normalisation applied to raw counts rather than using pre-normalised
#   log2(TPM+1).
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "inp-024"
stopifnot(id == "inp-024")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

run_dir   <- file.path("network/runs",   id)
cfg_dir   <- "network/config"
cache_dir <- "network/cache"
for (d_ in c(run_dir, cfg_dir, cache_dir))
  dir.create(d_, showWarnings = FALSE, recursive = TRUE)

# ---- parameters (ALL pre-specified) ----------------------------------------
POWER          <- 12L
POWER_MIN      <-  4L
NET_TYPE       <- "signed"
COR_FNC        <- "bicor"
GENE_SET       <- "all"
REF_GROUP      <- "control"
MIN_MODULE     <- 30L
MERGE_THRESH   <- 0.25
DRAWS          <- 20L        # revised from spec 100 (METHODS.md 2026-10-09)
PERMS_PER_DRAW <- 50L
NOLITH_PERMS   <- 500L
SUBSAMPLE_N    <- 74L
SEED_BASE      <- 20261009L + 24000L
CACHE_KEY      <- "control_all_voom_bicor_signed_p12"

# ---- config ----------------------------------------------------------------
cfg_path <- file.path(cfg_dir, paste0(id, ".json"))
if (!file.exists(cfg_path)) {
  scl <- function(v) {
    if (is.logical(v)) return(if (isTRUE(v)) "true" else "false")
    if (is.integer(v) || is.numeric(v)) return(as.character(v[1]))
    paste0('"', gsub('"', '\\"', as.character(v[1]), fixed = TRUE), '"')
  }
  arr <- function(v) paste0("[", paste(vapply(v, scl, ""), collapse = ", "), "]")
  jsn <- function(nm, v) paste0('  "', nm, '": ', if (length(v) > 1) arr(v) else scl(v))
  sha <- tryCatch(trimws(system("git rev-parse --short HEAD 2>/dev/null", intern=TRUE)[1]),
                  error = function(e) "unknown")
  cfg <- list(
    id                  = id,
    family              = "input",
    spec                = "ref=control;input=voom_weighted_from_counts;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans sessions)",
    draws_spec          = 100L,
    draws_implemented   = DRAWS,
    draws_revision_note = "Revised to 20 for budget (METHODS.md 2026-10-09)",
    input_description   = "voom log2-CPM from raw counts: DGEList -> calcNormFactors(TMM) -> voom(). TMM corrects for library size differences across samples.",
    voom_weights_note   = "voom-normalized log2-CPM used as input to WGCNA. The voom precision weights (which account for mean-variance trend in count data) are NOT carried into the WGCNA correlation step -- no approved WGCNA function supports weighted correlation in this manner. The benefit is TMM normalization applied to raw counts rather than using pre-normalized log2(TPM+1).",
    voom_design_note    = "Voom normalization uses a minimal design (~1, intercept only) per group separately; no formula covariates applied here. This is purely for library-size normalisation prior to WGCNA, not for DE analysis.",
    ref_group           = REF_GROUP,
    test_groups         = c("bp_nolith", "bp_lith"),
    corFnc              = COR_FNC,
    networkType         = NET_TYPE,
    gene_set            = GENE_SET,
    power               = POWER,
    power_source        = "WGCNA_FAQ_signed_n_gt_40",
    minModuleSize       = MIN_MODULE,
    mergeCutHeight      = MERGE_THRESH,
    draws_lith          = DRAWS,
    perms_per_lith_draw = PERMS_PER_DRAW,
    nolith_perms        = NOLITH_PERMS,
    subsample_n         = SUBSAMPLE_N,
    bp_lith_subsampled  = TRUE,
    bp_nolith_subsampled = FALSE,
    stat                = c("Zsummary", "medianRank"),
    descriptive_only    = TRUE,
    seed_base           = SEED_BASE,
    cache_key           = CACHE_KEY,
    git_sha             = if (length(sha) && !is.na(sha)) sha else "unknown",
    r_version           = paste(R.version$major, R.version$minor, sep = "."),
    wgcna_version       = as.character(packageVersion("WGCNA")),
    limma_version       = as.character(packageVersion("limma")),
    edgeR_version       = as.character(packageVersion("edgeR"))
  )
  lines <- mapply(jsn, names(cfg), cfg, SIMPLIFY = TRUE)
  writeLines(c("{", paste(lines, collapse = ",\n"), "}"), cfg_path)
  msg("config written: ", cfg_path)
}

# ---- packages --------------------------------------------------------------
suppressPackageStartupMessages({
  library(WGCNA)
  library(fastcluster)
  library(limma)
  library(edgeR)
})
enableWGCNAThreads()
options(stringsAsFactors = FALSE)
msg("WGCNA ", as.character(packageVersion("WGCNA")),
    " | limma ", as.character(packageVersion("limma")),
    " | edgeR ", as.character(packageVersion("edgeR")), " loaded")

# ---- data ------------------------------------------------------------------
source("network/00_load.R")
d <- load_project()

# ---- voom normalisation function -------------------------------------------
# Returns a samples x genes matrix (same orientation as expr_for()).
# Uses a minimal design (~1) -- no formula covariates; this is for library-size
# normalisation only (voom's primary contribution is TMM + mean-variance trend).
voom_for <- function(d_obj, group) {
  keep       <- !is.na(d_obj$groups) & d_obj$groups == group
  counts_grp <- d_obj$counts[, keep, drop = FALSE]   # genes x samples
  y          <- DGEList(counts = counts_grp)
  y          <- calcNormFactors(y, method = "TMM")
  v          <- voom(y, design = NULL, plot = FALSE)
  # v$E is genes x samples; transpose to samples x genes for WGCNA
  t(v$E)
}

# ---- scale-free scan (voom-normalised reference) ---------------------------
r2_path <- file.path(run_dir, "sft_r2_curve.csv")
if (!file.exists(r2_path)) {
  msg("pickSoftThreshold scan on voom-normalised control expression ...")
  tmp_ref <- voom_for(d, group = REF_GROUP)
  set.seed(SEED_BASE)
  sft <- pickSoftThreshold(tmp_ref, powerVector = 1:20,
                            networkType = NET_TYPE, corFnc = COR_FNC, verbose = 0)
  est <- sft$powerEstimate
  msg("powerEstimate = ", if (is.na(est)) "NA" else est,
      " (informational; power fixed at ", POWER, ")")
  r2df <- sft$fitIndices
  r2df$power_used <- POWER; r2df$powerEstimate <- est
  write.csv(r2df, r2_path, row.names = FALSE)
  msg("R2 curve saved (max R2=", round(max(r2df$SFT.R.sq, na.rm = TRUE), 3), ")")
  rm(tmp_ref, sft, r2df); gc()
} else {
  msg("R2 curve checkpoint exists, skipping scan")
}

# ---- reference network (voom-normalised control, cached) -------------------
ref_cache <- file.path(cache_dir, paste0("ref_", CACHE_KEY, ".rds"))
if (file.exists(ref_cache)) {
  msg("loading cached reference modules: ", ref_cache)
  ref_net      <- readRDS(ref_cache)
  moduleColors <- ref_net$moduleColors
  msg(ref_net$n_modules, " modules loaded")
} else {
  msg("building reference network (voom-normalised, bicor, signed, power ", POWER, ") ...")
  tmp_ref <- voom_for(d, group = REF_GROUP)
  set.seed(SEED_BASE)
  bwm <- blockwiseModules(tmp_ref, power = POWER,
    networkType = NET_TYPE, TOMType = NET_TYPE, corType = COR_FNC,
    minModuleSize = MIN_MODULE, mergeCutHeight = MERGE_THRESH,
    maxBlockSize = 15000L, numericLabels = FALSE, saveTOMs = FALSE, verbose = 3)
  moduleColors <- bwm$colors
  n_mods       <- length(unique(moduleColors)) - 1L
  msg(n_mods, " modules detected (excluding grey); ",
      sum(moduleColors == "grey"), " unassigned genes")
  ref_net <- list(moduleColors = moduleColors, n_modules = n_mods,
                  module_sizes = table(moduleColors), power = POWER,
                  cache_key = CACHE_KEY, input_type = "voom_log2cpm")
  saveRDS(ref_net, ref_cache)
  write.csv(as.data.frame(ref_net$module_sizes),
            file.path(run_dir, "reference_module_sizes.csv"), row.names = FALSE)
  msg("reference modules cached: ", ref_cache)
  rm(bwm, tmp_ref); gc()
}

# refExpr for preservation calls
refExpr <- voom_for(d, group = REF_GROUP)   # 234 x 12368

# Pre-compute voom on ALL 152 bp_lith samples; subsample rows per draw.
# This is cheaper than re-running voom every draw and ensures within-group
# TMM factors are stable (computed once on the full bp_lith cohort).
msg("pre-computing voom on all 152 bp_lith samples ...")
lithAllVoom <- voom_for(d, group = "bp_lith")   # 152 x 12368
lith_sids   <- rownames(lithAllVoom)
stopifnot(length(lith_sids) == 152L)
msg("bp_lith voom done: matrix ", nrow(lithAllVoom), " x ", ncol(lithAllVoom))

# ---- helper ----------------------------------------------------------------
extract_mp <- function(mp) {
  obs  <- mp$preservation$observed[["ref.Ref"]][["inColumnsAlsoPresentIn.Test"]]
  zsco <- mp$preservation$Z[["ref.Ref"]][["inColumnsAlsoPresentIn.Test"]]
  if (is.null(obs) || is.null(zsco))
    stop("modulePreservation return structure unexpected")
  res <- data.frame(module = rownames(obs), Zsummary = zsco[, "Zsummary.pres"],
                    medianRank = obs[, "medianRank.pres"], stringsAsFactors = FALSE)
  res[res$module != "grey", ]
}

# ---- bp_nolith (one run, 500 perms) ----------------------------------------
nolith_ckpt <- file.path(run_dir, "nolith_preservation.rds")
if (!file.exists(nolith_ckpt)) {
  msg("bp_nolith preservation (voom-normalised, ", NOLITH_PERMS, " perms) ...")
  nolithExpr <- voom_for(d, group = "bp_nolith")   # 74 x 12368
  perm_file  <- file.path(run_dir, "permStats_nolith.RData")
  mp_nolith  <- modulePreservation(
    multiData = list(Ref = list(data = refExpr), Test = list(data = nolithExpr)),
    multiColor = list(Ref = moduleColors), dataIsExpr = TRUE,
    networkType = NET_TYPE, corFnc = COR_FNC,
    nPermutations = NOLITH_PERMS, randomSeed = SEED_BASE + 1000L,
    savePermutedStatistics = TRUE, loadPermutedStatistics = FALSE,
    permutedStatisticsFile = perm_file, verbose = 1, indent = 2)
  nr        <- extract_mp(mp_nolith)
  nr$group  <- "bp_nolith"; nr$draw  <- NA_integer_
  nr$nperms <- NOLITH_PERMS; nr$n_test <- nrow(nolithExpr)
  saveRDS(nr, nolith_ckpt)
  msg("bp_nolith done: Zsummary [", round(min(nr$Zsummary), 2),
      ", ", round(max(nr$Zsummary), 2), "]")
  rm(nolithExpr, mp_nolith); gc()
} else {
  msg("bp_nolith checkpoint exists, skipping")
}

# ---- bp_lith draws (20 x 50 perms, checkpoint per draw) --------------------
for (draw in seq_len(DRAWS)) {
  ckpt <- file.path(run_dir, sprintf("lith_draw_%03d.rds", draw))
  if (file.exists(ckpt)) { msg("draw ", draw, " exists, skipping"); next }
  msg("draw ", draw, "/", DRAWS, ": subsampling voom bp_lith ...")
  set.seed(SEED_BASE + draw)
  sel_sids <- sample(lith_sids, SUBSAMPLE_N, replace = FALSE)
  testExpr  <- lithAllVoom[sel_sids, , drop = FALSE]   # 74 x 12368
  perm_file <- file.path(run_dir, sprintf("permStats_lith_draw_%03d.RData", draw))
  mp_draw <- tryCatch(
    modulePreservation(
      multiData = list(Ref = list(data = refExpr), Test = list(data = testExpr)),
      multiColor = list(Ref = moduleColors), dataIsExpr = TRUE,
      networkType = NET_TYPE, corFnc = COR_FNC,
      nPermutations = PERMS_PER_DRAW, randomSeed = SEED_BASE + draw,
      savePermutedStatistics = TRUE, loadPermutedStatistics = FALSE,
      permutedStatisticsFile = perm_file, verbose = 1, indent = 2),
    error = function(e) { msg("ERROR draw ", draw, ": ", conditionMessage(e)); NULL })
  if (is.null(mp_draw)) { rm(testExpr); gc(); next }
  dr         <- extract_mp(mp_draw)
  dr$group   <- "bp_lith"; dr$draw <- draw
  dr$nperms  <- PERMS_PER_DRAW; dr$n_test <- SUBSAMPLE_N
  dr$samples <- paste(sel_sids, collapse = ",")
  saveRDS(dr, ckpt)
  msg("draw ", draw, " done: Zsummary [", round(min(dr$Zsummary), 2),
      ", ", round(max(dr$Zsummary), 2), "]")
  rm(testExpr, mp_draw, dr); gc()
}

# ---- aggregate (only after ALL complete) -----------------------------------
done_nolith <- file.exists(nolith_ckpt)
done_lith   <- vapply(seq_len(DRAWS), function(i)
  file.exists(file.path(run_dir, sprintf("lith_draw_%03d.rds", i))), logical(1))
if (!done_nolith || sum(done_lith) < DRAWS) {
  msg("not all complete: nolith=", done_nolith, " lith=", sum(done_lith), "/", DRAWS)
  quit(save = "no", status = 1L)
}
msg("aggregating results ...")

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
             Zsummary_median = median(s$Zsummary, na.rm = TRUE),
             Zsummary_sd     = sd(s$Zsummary,     na.rm = TRUE),
             medianRank_mean = mean(s$medianRank,  na.rm = TRUE),
             stringsAsFactors = FALSE)
}))
write.csv(lith_summary, file.path(run_dir, "results_lith_summary.csv"), row.names = FALSE)

msg("Results written to ", run_dir)
msg("=== ", id, " COMPLETE ===")
