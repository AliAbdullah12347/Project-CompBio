#!/usr/bin/env Rscript
# ===========================================================================
# base-002-ceil  Split-half ceiling for modulePreservation
#
# Spec: SPLIT-HALF CEILING. ref=160 randomly held-out controls;
#   test=the OTHER 74 controls (disjoint), plus bp_nolith(74) and
#   bp_lith(20 draws of 74) against the SAME reference.
#   corr=bicor;net=signed;genes=all;power=12;perms_per_draw=50;CHECKPOINT PER DRAW.
#   Gives the n=74 preservation CEILING: without it a Zsummary in bp_lith
#   cannot be distinguished from ordinary 74-sample noise.
#
# Why this matters:
#   base-001 uses a 234-sample reference and 74-sample test sets. A 234-sample
#   reference computes connectivity more precisely than any 74-sample test set
#   can reproduce, inflating apparent differences through precision alone.
#   The ceiling run uses a 160-sample reference vs. a DISJOINT 74 controls —
#   same population, matched test size — so any Zsummary for controls is the
#   upper bound expected under n=74 noise alone.
#
# Design (draws revised to 20 — same reasoning as base-001, METHODS.md 2026-10-09):
#   1. Partition 234 controls: 160 reference, 74 ceiling test (fixed seed).
#   2. Build reference modules on 160 controls.
#   3. Preservation in: 74-control ceiling, bp_nolith (74), 20 draws bp_lith.
#   4. Checkpoint per draw; aggregate after all done.
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "base-002-ceil"
stopifnot(id == "base-002-ceil")

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
MIN_MODULE     <- 30L
MERGE_THRESH   <- 0.25
DRAWS          <- 20L       # revised from spec 100 (METHODS.md 2026-10-09)
PERMS_PER_DRAW <- 50L
NOLITH_PERMS   <- 500L
SUBSAMPLE_N    <- 74L
N_REF_CTRL     <- 160L      # controls used to build reference (must be > SUBSAMPLE_N)
SEED_BASE      <- 20261009L
SEED_CTRL_SPLIT <- SEED_BASE + 2000L   # separate seed for control partition
CACHE_KEY      <- paste("ctrl160", GENE_SET, COR_FNC, NET_TYPE,
                        paste0("p", POWER), sep = "_")

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
  sha <- tryCatch(trimws(system("git rev-parse --short HEAD 2>/dev/null", intern = TRUE)[1]),
                  error = function(e) "unknown")
  cfg <- list(
    id = id, family = "baseline",
    spec = "SPLIT-HALF CEILING; ref=160 randomly held-out controls; test=74 disjoint controls+bp_nolith+bp_lith(20 draws)",
    draws_spec = 100L, draws_implemented = DRAWS,
    draws_revision_note = "Same as base-001: 100 draws x ~1962 s/draw = ~54 h; revised to 20",
    ref_controls_n    = N_REF_CTRL,
    ceiling_test_n    = SUBSAMPLE_N,
    seed_ctrl_split   = SEED_CTRL_SPLIT,
    corFnc = COR_FNC, networkType = NET_TYPE, gene_set = GENE_SET,
    power = POWER, power_source = "WGCNA_FAQ_signed_n_gt_40",
    minModuleSize = MIN_MODULE, mergeCutHeight = MERGE_THRESH,
    draws_lith = DRAWS, perms_per_lith_draw = PERMS_PER_DRAW,
    nolith_perms = NOLITH_PERMS,
    subsample_n = SUBSAMPLE_N,
    cache_key = CACHE_KEY,
    seed_base = SEED_BASE,
    git_sha   = if (length(sha <- tryCatch(trimws(system("git rev-parse --short HEAD 2>/dev/null", intern=TRUE)[1]), error=function(e)"unknown")) && !is.na(sha)) sha else "unknown",
    r_version    = paste(R.version$major, R.version$minor, sep = "."),
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

# ---- control partition (fixed seed, reproducible) --------------------------
ctrl_sids <- colnames(d$logtpm)[!is.na(d$groups) & d$groups == "control"]
stopifnot(length(ctrl_sids) == 234L)
set.seed(SEED_CTRL_SPLIT)
ref_ctrl  <- sample(ctrl_sids, N_REF_CTRL, replace = FALSE)  # 160 controls
ceil_ctrl <- setdiff(ctrl_sids, ref_ctrl)                     # 74 controls (disjoint)
stopifnot(length(ceil_ctrl) == (234L - N_REF_CTRL))
msg("Control split: ", length(ref_ctrl), " reference, ", length(ceil_ctrl), " ceiling test")

# Save the partition for reproducibility
part_path <- file.path(run_dir, "control_partition.csv")
if (!file.exists(part_path)) {
  write.csv(data.frame(sample_id = c(ref_ctrl, ceil_ctrl),
                       role      = c(rep("reference", length(ref_ctrl)),
                                     rep("ceiling_test", length(ceil_ctrl)))),
            part_path, row.names = FALSE)
}

# ---- scale-free scan (reference group, record curve) -----------------------
r2_path <- file.path(run_dir, "sft_r2_curve.csv")
if (!file.exists(r2_path)) {
  msg("pickSoftThreshold scan for 160-control reference ...")
  refExpr_tmp <- t(d$logtpm[, ref_ctrl])   # 160 x 12368
  set.seed(SEED_BASE)
  sft <- pickSoftThreshold(refExpr_tmp, powerVector = 1:20,
                           networkType = NET_TYPE, corFnc = COR_FNC, verbose = 0)
  est <- sft$powerEstimate
  msg("powerEstimate = ", if (is.na(est)) "NA" else est,
      " (expected 1 or NA: function failing)")
  r2df <- sft$fitIndices
  r2df$power_used <- POWER; r2df$powerEstimate <- est
  write.csv(r2df, r2_path, row.names = FALSE)
  msg("R2 curve saved (max R2=", round(max(r2df$SFT.R.sq, na.rm=TRUE), 3), ")")
  rm(refExpr_tmp, sft, r2df); gc()
}

# ---- reference network (160-control, cached separately from base-001) ------
ref_cache <- file.path(cache_dir, paste0("ref_", CACHE_KEY, ".rds"))
if (file.exists(ref_cache)) {
  msg("loading cached reference (160 controls): ", ref_cache)
  ref_net      <- readRDS(ref_cache)
  moduleColors <- ref_net$moduleColors
  msg(ref_net$n_modules, " modules (", sum(moduleColors != "grey"), " in non-grey)")
} else {
  msg("building reference network (160 controls, bicor, signed, power ", POWER, ") ...")
  refExpr_tmp <- t(d$logtpm[, ref_ctrl])   # 160 x 12368
  set.seed(SEED_BASE)
  bwm <- blockwiseModules(refExpr_tmp, power = POWER,
    networkType = NET_TYPE, TOMType = NET_TYPE, corType = COR_FNC,
    minModuleSize = MIN_MODULE, mergeCutHeight = MERGE_THRESH,
    maxBlockSize = 15000L, numericLabels = FALSE, saveTOMs = FALSE, verbose = 3)
  moduleColors <- bwm$colors
  n_mods <- length(unique(moduleColors)) - 1L
  msg("reference: ", n_mods, " modules")
  ref_net <- list(moduleColors = moduleColors, n_modules = n_mods,
                  module_sizes = table(moduleColors), power = POWER,
                  cache_key = CACHE_KEY, ref_ctrl_n = N_REF_CTRL)
  saveRDS(ref_net, ref_cache)
  write.csv(as.data.frame(ref_net$module_sizes),
            file.path(run_dir, "reference_module_sizes.csv"), row.names = FALSE)
  rm(bwm, refExpr_tmp); gc()
}

# refExpr for preservation calls
refExpr <- t(d$logtpm[, ref_ctrl])   # 160 x 12368

# ---- helper -----------------------------------------------------------------
extract_mp <- function(mp, ref_name = "Ref", test_name = "Test") {
  key_obs  <- paste0("ref.", ref_name)
  key_test <- paste0("inColumnsAlsoPresentIn.", test_name)
  obs  <- mp$preservation$observed[[key_obs]][[key_test]]
  zsco <- mp$preservation$Z[[key_obs]][[key_test]]
  if (is.null(obs) || is.null(zsco))
    stop("modulePreservation return structure unexpected for ", key_obs, " / ", key_test)
  res <- data.frame(module = rownames(obs), Zsummary = zsco[, "Zsummary.pres"],
                    medianRank = obs[, "medianRank.pres"], stringsAsFactors = FALSE)
  res[res$module != "grey", ]
}

# ---- ceiling: 74 disjoint controls (one run, many perms) -------------------
ceil_ckpt <- file.path(run_dir, "ceiling_ctrl_preservation.rds")
if (!file.exists(ceil_ckpt)) {
  msg("CEILING: 74 control preservation (", NOLITH_PERMS, " perms) ...")
  ceilExpr  <- t(d$logtpm[, ceil_ctrl])   # 74 x 12368
  perm_file <- file.path(run_dir, "permStats_ceiling_ctrl.RData")
  mp_ceil   <- modulePreservation(
    multiData = list(Ref = list(data = refExpr), Test = list(data = ceilExpr)),
    multiColor = list(Ref = moduleColors), dataIsExpr = TRUE,
    networkType = NET_TYPE, corFnc = COR_FNC,
    nPermutations = NOLITH_PERMS, randomSeed = SEED_BASE + 2001L,
    savePermutedStatistics = TRUE, loadPermutedStatistics = FALSE,
    permutedStatisticsFile = perm_file, verbose = 1, indent = 2)
  ceil_res            <- extract_mp(mp_ceil)
  ceil_res$group      <- "ctrl_ceiling"
  ceil_res$draw       <- NA_integer_
  ceil_res$nperms     <- NOLITH_PERMS
  ceil_res$n_test     <- nrow(ceilExpr)
  saveRDS(ceil_res, ceil_ckpt)
  msg("ceiling done: Zsummary range [", round(min(ceil_res$Zsummary), 2),
      ", ", round(max(ceil_res$Zsummary), 2), "]")
  rm(ceilExpr, mp_ceil); gc()
} else {
  msg("ceiling checkpoint exists, skipping")
}

# ---- bp_nolith preservation (one run, many perms) --------------------------
nolith_ckpt <- file.path(run_dir, "nolith_preservation.rds")
if (!file.exists(nolith_ckpt)) {
  msg("bp_nolith preservation (", NOLITH_PERMS, " perms) ...")
  nolithExpr <- expr_for(d, group = "bp_nolith")   # 74 x 12368
  perm_file  <- file.path(run_dir, "permStats_nolith.RData")
  mp_nolith  <- modulePreservation(
    multiData = list(Ref = list(data = refExpr), Test = list(data = nolithExpr)),
    multiColor = list(Ref = moduleColors), dataIsExpr = TRUE,
    networkType = NET_TYPE, corFnc = COR_FNC,
    nPermutations = NOLITH_PERMS, randomSeed = SEED_BASE + 3000L,
    savePermutedStatistics = TRUE, loadPermutedStatistics = FALSE,
    permutedStatisticsFile = perm_file, verbose = 1, indent = 2)
  nr             <- extract_mp(mp_nolith)
  nr$group       <- "bp_nolith"
  nr$draw        <- NA_integer_
  nr$nperms      <- NOLITH_PERMS
  nr$n_test      <- nrow(nolithExpr)
  saveRDS(nr, nolith_ckpt)
  msg("bp_nolith done: Zsummary range [", round(min(nr$Zsummary), 2),
      ", ", round(max(nr$Zsummary), 2), "]")
  rm(nolithExpr, mp_nolith); gc()
} else {
  msg("bp_nolith checkpoint exists, skipping")
}

# ---- bp_lith draws ---------------------------------------------------------
lith_sids <- colnames(d$logtpm)[!is.na(d$groups) & d$groups == "bp_lith"]
stopifnot(length(lith_sids) == 152L)

for (draw in seq_len(DRAWS)) {
  ckpt <- file.path(run_dir, sprintf("lith_draw_%03d.rds", draw))
  if (file.exists(ckpt)) { msg("draw ", draw, " exists, skipping"); next }
  msg("draw ", draw, "/", DRAWS, ": subsampling bp_lith ...")
  set.seed(SEED_BASE + draw)   # same seeds as base-001 for comparability
  sel_sids  <- sample(lith_sids, SUBSAMPLE_N, replace = FALSE)
  testExpr  <- t(d$logtpm[, sel_sids])
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
  dr             <- extract_mp(mp_draw)
  dr$group       <- "bp_lith"
  dr$draw        <- draw
  dr$nperms      <- PERMS_PER_DRAW
  dr$n_test      <- SUBSAMPLE_N
  dr$samples     <- paste(sel_sids, collapse = ",")
  saveRDS(dr, ckpt)
  msg("draw ", draw, " done: Zsummary range [", round(min(dr$Zsummary), 2),
      ", ", round(max(dr$Zsummary), 2), "]")
  rm(testExpr, mp_draw, dr); gc()
}

# ---- aggregate (only after all work done) ----------------------------------
done_ceil   <- file.exists(ceil_ckpt)
done_nolith <- file.exists(nolith_ckpt)
done_lith   <- vapply(seq_len(DRAWS), function(i)
  file.exists(file.path(run_dir, sprintf("lith_draw_%03d.rds", i))), logical(1))

if (!done_ceil || !done_nolith || sum(done_lith) < DRAWS) {
  msg("not all complete: ceiling=", done_ceil, " nolith=", done_nolith,
      " lith=", sum(done_lith), "/", DRAWS)
  quit(save = "no", status = 1L)
}

msg("aggregating ...")

# All groups combined
ceil_res   <- readRDS(ceil_ckpt)
nolith_res <- readRDS(nolith_ckpt)
lith_all   <- do.call(rbind, lapply(seq_len(DRAWS), function(i) {
  x <- readRDS(file.path(run_dir, sprintf("lith_draw_%03d.rds", i)))
  x[, c("group","module","Zsummary","medianRank","draw","nperms","n_test")]
}))

# Per-draw all in one file
all_oneshot <- rbind(
  ceil_res[, c("group","module","Zsummary","medianRank","draw","nperms","n_test")],
  nolith_res[, c("group","module","Zsummary","medianRank","draw","nperms","n_test")],
  lith_all
)
write.csv(all_oneshot, file.path(run_dir, "results_all_groups.csv"), row.names = FALSE)

# bp_lith summary
modules      <- unique(lith_all$module)
lith_summary <- do.call(rbind, lapply(modules, function(m) {
  s <- lith_all[lith_all$module == m, ]
  data.frame(group = "bp_lith", module = m, n_draws = nrow(s),
             Zsummary_mean = mean(s$Zsummary, na.rm=TRUE),
             Zsummary_median = median(s$Zsummary, na.rm=TRUE),
             Zsummary_sd = sd(s$Zsummary, na.rm=TRUE),
             medianRank_mean = mean(s$medianRank, na.rm=TRUE),
             stringsAsFactors = FALSE)
}))
write.csv(lith_summary, file.path(run_dir, "results_lith_summary.csv"), row.names = FALSE)

msg("Results written to ", run_dir)
msg("=== ", id, " COMPLETE ===")
