#!/usr/bin/env Rscript
# ===========================================================================
# inp-022  KEY QUESTION: composition-residualised network
#
# Spec: ref=control;input=residualised_composition(ILR b1-b4) ** KEY QUESTION **
#       draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans sessions)
#
# Design:
#   - Residualise ALL expression data for cell composition (ILR b1-b4).
#   - Build reference modules on residualised control data.
#   - Test preservation in residualised bp_nolith and bp_lith.
#   - Draws=20 (same revision as base-001; METHODS.md 2026-10-09).
#
# Scientific purpose:
#   Compare Zsummary here vs. base-001 (raw expression). If low Zsummary in
#   base-001's bp_lith persists here, the disruption is NOT due to composition.
#   If Zsummary recovers here, composition was driving the base-001 finding.
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "inp-022"
stopifnot(id == "inp-022")

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
DRAWS          <- 20L       # revised from spec 100 (METHODS.md 2026-10-09)
PERMS_PER_DRAW <- 50L
NOLITH_PERMS   <- 500L
SUBSAMPLE_N    <- 74L
SEED_BASE      <- 20261009L
# Residualise for ILR b1-b4 (cell composition proxies)
RESID_FORMULA  <- ~ b1 + b2 + b3 + b4
RESID_NAME     <- "composition_ilr_b1b4"
CACHE_KEY      <- paste(REF_GROUP, GENE_SET, COR_FNC, NET_TYPE,
                        paste0("p", POWER), RESID_NAME, sep = "_")

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
    id = id, family = "input",
    spec = "ref=control;input=residualised_composition(ILR b1-b4) ** KEY QUESTION **;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans sessions)",
    draws_spec = 100L, draws_implemented = DRAWS,
    draws_revision_note = "Same as base-001: revised to 20 for budget (METHODS.md 2026-10-09)",
    residualise_formula = deparse(RESID_FORMULA),
    residualise_name    = RESID_NAME,
    residualise_note    = "ILR b1-b4 are log-ratio balances of 5 lineage cell fractions (gran,mono,T,NK,B); residualising removes the linear cell-composition effect from expression",
    ref_group = REF_GROUP, test_groups = c("bp_nolith","bp_lith"),
    corFnc = COR_FNC, networkType = NET_TYPE, gene_set = GENE_SET,
    power = POWER, power_source = "WGCNA_FAQ_signed_n_gt_40",
    minModuleSize = MIN_MODULE, mergeCutHeight = MERGE_THRESH,
    draws_lith = DRAWS, perms_per_lith_draw = PERMS_PER_DRAW,
    nolith_perms = NOLITH_PERMS, subsample_n = SUBSAMPLE_N,
    bp_lith_subsampled = TRUE, bp_nolith_subsampled = FALSE,
    stat = c("Zsummary","medianRank"), descriptive_only = TRUE,
    seed_base = SEED_BASE, cache_key = CACHE_KEY,
    git_sha   = if (length(sha) && !is.na(sha)) sha else "unknown",
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    wgcna_version = as.character(packageVersion("WGCNA"))
  )
  lines <- mapply(jsn, names(cfg), cfg, SIMPLIFY = TRUE)
  writeLines(c("{", paste(lines, collapse = ",\n"), "}"), cfg_path)
  msg("config written: ", cfg_path)
}

# ---- packages --------------------------------------------------------------
suppressPackageStartupMessages({
  library(WGCNA)
  library(fastcluster)
})
enableWGCNAThreads()
options(stringsAsFactors = FALSE)

# ---- data ------------------------------------------------------------------
source("network/00_load.R")
d <- load_project()

# ---- scale-free scan (residualised reference, record curve) ----------------
r2_path <- file.path(run_dir, "sft_r2_curve.csv")
if (!file.exists(r2_path)) {
  msg("pickSoftThreshold scan on residualised control expression ...")
  tmp_ref <- expr_for(d, group = REF_GROUP, residualise = RESID_FORMULA)
  set.seed(SEED_BASE)
  sft <- pickSoftThreshold(tmp_ref, powerVector = 1:20,
                           networkType = NET_TYPE, corFnc = COR_FNC, verbose = 0)
  est <- sft$powerEstimate
  msg("powerEstimate = ", if(is.na(est)) "NA" else est, " (expected 1 or NA)")
  r2df <- sft$fitIndices
  r2df$power_used <- POWER; r2df$powerEstimate <- est
  write.csv(r2df, r2_path, row.names = FALSE)
  msg("R2 curve saved (max R2=", round(max(r2df$SFT.R.sq, na.rm=TRUE), 3), ")")
  rm(tmp_ref, sft, r2df); gc()
}

# ---- reference network (residualised, cached) ------------------------------
ref_cache <- file.path(cache_dir, paste0("ref_", CACHE_KEY, ".rds"))
if (file.exists(ref_cache)) {
  msg("loading cached reference modules: ", ref_cache)
  ref_net      <- readRDS(ref_cache)
  moduleColors <- ref_net$moduleColors
  msg(ref_net$n_modules, " modules loaded")
} else {
  msg("building reference network (residualised bicor, signed, power ", POWER, ") ...")
  tmp_ref <- expr_for(d, group = REF_GROUP, residualise = RESID_FORMULA)
  set.seed(SEED_BASE)
  bwm <- blockwiseModules(tmp_ref, power = POWER,
    networkType = NET_TYPE, TOMType = NET_TYPE, corType = COR_FNC,
    minModuleSize = MIN_MODULE, mergeCutHeight = MERGE_THRESH,
    maxBlockSize = 15000L, numericLabels = FALSE, saveTOMs = FALSE, verbose = 3)
  moduleColors <- bwm$colors
  n_mods <- length(unique(moduleColors)) - 1L
  msg(n_mods, " modules detected (excluding grey)")
  ref_net <- list(moduleColors = moduleColors, n_modules = n_mods,
                  module_sizes = table(moduleColors), power = POWER,
                  cache_key = CACHE_KEY, residualise = deparse(RESID_FORMULA))
  saveRDS(ref_net, ref_cache)
  write.csv(as.data.frame(ref_net$module_sizes),
            file.path(run_dir, "reference_module_sizes.csv"), row.names = FALSE)
  rm(bwm, tmp_ref); gc()
}

# Residualised reference expression for preservation calls
refExpr <- expr_for(d, group = REF_GROUP, residualise = RESID_FORMULA)

# Pre-compute residualised bp_lith expression (all 152 samples) so we can
# subsample the residuals rather than re-running lm per draw.
lithAllRes <- expr_for(d, group = "bp_lith", residualise = RESID_FORMULA)  # 152 x 12368
lith_sids  <- rownames(lithAllRes)   # after residualise, rownames are sample IDs
stopifnot(length(lith_sids) == 152L)

# ---- helper ----------------------------------------------------------------
extract_mp <- function(mp) {
  obs  <- mp$preservation$observed[["ref.Ref"]][["inColumnsAlsoPresentIn.Test"]]
  zsco <- mp$preservation$Z[["ref.Ref"]][["inColumnsAlsoPresentIn.Test"]]
  if (is.null(obs) || is.null(zsco))
    stop("modulePreservation return structure unexpected")
  res <- data.frame(module = rownames(obs), Zsummary = zsco[,"Zsummary.pres"],
                    medianRank = obs[,"medianRank.pres"], stringsAsFactors = FALSE)
  res[res$module != "grey", ]
}

# ---- bp_nolith (one run) --------------------------------------------------
nolith_ckpt <- file.path(run_dir, "nolith_preservation.rds")
if (!file.exists(nolith_ckpt)) {
  msg("bp_nolith preservation (residualised, ", NOLITH_PERMS, " perms) ...")
  nolithExpr <- expr_for(d, group = "bp_nolith", residualise = RESID_FORMULA)
  perm_file  <- file.path(run_dir, "permStats_nolith.RData")
  mp_nolith  <- modulePreservation(
    multiData = list(Ref = list(data = refExpr), Test = list(data = nolithExpr)),
    multiColor = list(Ref = moduleColors), dataIsExpr = TRUE,
    networkType = NET_TYPE, corFnc = COR_FNC,
    nPermutations = NOLITH_PERMS, randomSeed = SEED_BASE + 1000L,
    savePermutedStatistics = TRUE, loadPermutedStatistics = FALSE,
    permutedStatisticsFile = perm_file, verbose = 1, indent = 2)
  nr           <- extract_mp(mp_nolith)
  nr$group     <- "bp_nolith"; nr$draw <- NA_integer_
  nr$nperms    <- NOLITH_PERMS; nr$n_test <- nrow(nolithExpr)
  saveRDS(nr, nolith_ckpt)
  msg("bp_nolith done: Zsummary [", round(min(nr$Zsummary),2),
      ", ", round(max(nr$Zsummary),2), "]")
  rm(nolithExpr, mp_nolith); gc()
} else {
  msg("bp_nolith checkpoint exists, skipping")
}

# ---- bp_lith draws ---------------------------------------------------------
for (draw in seq_len(DRAWS)) {
  ckpt <- file.path(run_dir, sprintf("lith_draw_%03d.rds", draw))
  if (file.exists(ckpt)) { msg("draw ", draw, " exists, skipping"); next }
  msg("draw ", draw, "/", DRAWS, ": subsampling residualised bp_lith ...")
  set.seed(SEED_BASE + draw)   # same seeds as base-001 for matched comparison
  sel_sids <- sample(lith_sids, SUBSAMPLE_N, replace = FALSE)
  testExpr <- lithAllRes[sel_sids, , drop = FALSE]   # 74 x 12368
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
  dr          <- extract_mp(mp_draw)
  dr$group    <- "bp_lith"; dr$draw <- draw
  dr$nperms   <- PERMS_PER_DRAW; dr$n_test <- SUBSAMPLE_N
  dr$samples  <- paste(sel_sids, collapse = ",")
  saveRDS(dr, ckpt)
  msg("draw ", draw, " done: Zsummary [", round(min(dr$Zsummary),2),
      ", ", round(max(dr$Zsummary),2), "]")
  rm(testExpr, mp_draw, dr); gc()
}

# ---- aggregate -------------------------------------------------------------
done_nolith <- file.exists(nolith_ckpt)
done_lith   <- vapply(seq_len(DRAWS), function(i)
  file.exists(file.path(run_dir, sprintf("lith_draw_%03d.rds", i))), logical(1))
if (!done_nolith || sum(done_lith) < DRAWS) {
  msg("not all complete: nolith=", done_nolith, " lith=", sum(done_lith), "/", DRAWS)
  quit(save = "no", status = 1L)
}
msg("aggregating results ...")
nolith_res <- readRDS(nolith_ckpt)
write.csv(nolith_res[, c("group","module","Zsummary","medianRank","draw","nperms","n_test")],
          file.path(run_dir, "results_nolith.csv"), row.names = FALSE)
lith_all <- do.call(rbind, lapply(seq_len(DRAWS), function(i) {
  x <- readRDS(file.path(run_dir, sprintf("lith_draw_%03d.rds", i)))
  x[, c("group","module","Zsummary","medianRank","draw","nperms","n_test")]
}))
write.csv(lith_all, file.path(run_dir, "results_lith_alldraws.csv"), row.names = FALSE)
modules <- unique(lith_all$module)
lith_summary <- do.call(rbind, lapply(modules, function(m) {
  s <- lith_all[lith_all$module == m, ]
  data.frame(group="bp_lith", module=m, n_draws=nrow(s),
             Zsummary_mean=mean(s$Zsummary,na.rm=T),
             Zsummary_median=median(s$Zsummary,na.rm=T),
             Zsummary_sd=sd(s$Zsummary,na.rm=T),
             medianRank_mean=mean(s$medianRank,na.rm=T),
             stringsAsFactors=FALSE)
}))
write.csv(lith_summary, file.path(run_dir, "results_lith_summary.csv"), row.names = FALSE)
msg("Results written to ", run_dir)
msg("=== ", id, " COMPLETE ===")
