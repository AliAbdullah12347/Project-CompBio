#!/usr/bin/env Rscript
# ===========================================================================
# ref-004  Consensus reference network from all 3 groups
#
# Spec: "ref=consensus(all3);corr=bicor;net=signed;genes=all;draws=100;
#        perms_per_draw=50;CHECKPOINT PER DRAW"
#
# Design:
#   1. Build consensus network with WGCNA blockwiseConsensusModules() across
#      all 3 groups (control, bp_nolith, bp_lith subsample n=74 fixed draw).
#      Cache consensus moduleColors to network/cache/.
#   2. Test module preservation in each group using control expression as the
#      reference set (with consensus moduleColors applied):
#        - control (n=234, all samples): ONE run, 500 perms. CHECKPOINT.
#        - bp_nolith (n=74, all samples): ONE run, 500 perms. CHECKPOINT.
#        - bp_lith (n=152 -> subsample 74): 20 draws x 50 perms. CHECKPOINT PER DRAW.
#   3. Draws=20 (revised from spec 100; METHODS.md 2026-10-09).
#
# Scientific rationale:
#   Consensus modules are co-expressed across ALL groups simultaneously.
#   Testing preservation in each individual group then asks: are modules that
#   are globally co-expressed showing group-specific preservation differences?
#   (i.e., are there modules that dissolve in one group but survive in others?)
#
# Config note:
#   Consensus modules defined across all 3 groups using blockwiseConsensusModules().
#   modulePreservation uses control expression as the reference set but applies
#   consensus moduleColors. This tests: do modules that are co-expressed across
#   all groups show group-specific preservation differences?
#
# CACHE_KEY = "consensus_all3_bicor_signed_p12"
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "ref-004"
stopifnot(id == "ref-004")

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
DRAWS          <- 20L        # revised from spec 100 (METHODS.md 2026-10-09)
PERMS_PER_DRAW <- 50L
CTRL_PERMS     <- 500L       # control: one fixed run, many perms
NOLITH_PERMS   <- 500L       # bp_nolith: one fixed run, many perms
SUBSAMPLE_N    <- 74L
SEED_BASE      <- 20261009L + 4000L
# bp_lith fixed subsample seed for consensus modules input
LITH_CONSENSUS_SEED <- SEED_BASE + 999L
CACHE_KEY      <- "consensus_all3_bicor_signed_p12"

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
    id                   = id,
    family               = "reference",
    spec                 = "ref=consensus(all3);corr=bicor;net=signed;genes=all;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW",
    draws_spec           = 100L,
    draws_implemented    = DRAWS,
    draws_revision_note  = "Revised to 20 for budget (METHODS.md 2026-10-09)",
    consensus_groups     = c("control", "bp_nolith", "bp_lith"),
    consensus_note       = "blockwiseConsensusModules() applied across all 3 groups simultaneously. bp_lith subsampled to 74 samples for balance (LITH_CONSENSUS_SEED); one fixed subsample used for consensus definition.",
    lith_consensus_seed  = LITH_CONSENSUS_SEED,
    preservation_ref_set = "control",
    preservation_ref_note = "modulePreservation uses control expression (n=234) as the reference set with consensus moduleColors applied. Tests group-specific preservation of globally-defined modules.",
    corFnc               = COR_FNC,
    networkType          = NET_TYPE,
    gene_set             = GENE_SET,
    power                = POWER,
    power_source         = "WGCNA_FAQ_signed_n_gt_40",
    minModuleSize        = MIN_MODULE,
    mergeCutHeight       = MERGE_THRESH,
    ctrl_perms           = CTRL_PERMS,
    nolith_perms         = NOLITH_PERMS,
    draws_lith           = DRAWS,
    perms_per_lith_draw  = PERMS_PER_DRAW,
    subsample_n          = SUBSAMPLE_N,
    bp_lith_subsampled   = TRUE,
    bp_nolith_subsampled = FALSE,
    control_subsampled   = FALSE,
    stat                 = c("Zsummary", "medianRank"),
    descriptive_only     = TRUE,
    seed_base            = SEED_BASE,
    cache_key            = CACHE_KEY,
    dependency_check     = "none (no prior experiment required)",
    git_sha              = if (length(sha) && !is.na(sha)) sha else "unknown",
    r_version            = paste(R.version$major, R.version$minor, sep = "."),
    wgcna_version        = as.character(packageVersion("WGCNA"))
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
msg("WGCNA ", as.character(packageVersion("WGCNA")), " loaded; threads enabled")

# ---- data ------------------------------------------------------------------
source("network/00_load.R")
d <- load_project()

# ---- consensus network (blockwiseConsensusModules) -------------------------
consensus_cache <- file.path(cache_dir, paste0("ref_", CACHE_KEY, ".rds"))

if (file.exists(consensus_cache)) {
  msg("loading cached consensus modules: ", consensus_cache)
  con_net      <- readRDS(consensus_cache)
  moduleColors <- con_net$moduleColors
  msg(con_net$n_modules, " consensus modules loaded")
} else {
  msg("building consensus network across control + bp_nolith + bp_lith ...")

  # bp_lith subsample (fixed draw for consensus definition, not one of the 20)
  set.seed(LITH_CONSENSUS_SEED)
  lith_all_sids <- colnames(d$logtpm)[!is.na(d$groups) & d$groups == "bp_lith"]
  stopifnot(length(lith_all_sids) == 152L)
  lith_sel_sids <- sample(lith_all_sids, SUBSAMPLE_N, replace = FALSE)
  lith_subsample <- t(d$logtpm[, lith_sel_sids])   # 74 x 12368
  msg("bp_lith subsample for consensus: ", length(lith_sel_sids), " of 152 samples")

  multiExpr <- list(
    ctrl   = list(data = expr_for(d, group = "control")),    # 234 x 12368
    nolith = list(data = expr_for(d, group = "bp_nolith")),  # 74  x 12368
    lith   = list(data = lith_subsample)                     # 74  x 12368
  )
  msg("multiExpr: ctrl=", nrow(multiExpr$ctrl$data), " nolith=",
      nrow(multiExpr$nolith$data), " lith=", nrow(multiExpr$lith$data))

  set.seed(SEED_BASE)
  consensus_net <- blockwiseConsensusModules(
    multiExpr      = multiExpr,
    power          = POWER,
    networkType    = NET_TYPE,
    TOMType        = "signed",
    corType        = COR_FNC,
    minModuleSize  = MIN_MODULE,
    mergeCutHeight = MERGE_THRESH,
    maxBlockSize   = 15000L,
    numericLabels  = FALSE,
    verbose        = 3
  )

  moduleColors <- consensus_net$colors
  n_mods       <- length(unique(moduleColors)) - 1L
  msg("consensus network: ", n_mods, " modules (excluding grey); ",
      sum(moduleColors == "grey"), " unassigned genes")

  con_net <- list(
    moduleColors      = moduleColors,
    n_modules         = n_mods,
    module_sizes      = table(moduleColors),
    power             = POWER,
    cache_key         = CACHE_KEY,
    corFnc            = COR_FNC,
    networkType       = NET_TYPE,
    lith_sel_sids     = lith_sel_sids,
    lith_consensus_seed = LITH_CONSENSUS_SEED
  )
  saveRDS(con_net, consensus_cache)
  write.csv(as.data.frame(con_net$module_sizes),
            file.path(run_dir, "consensus_module_sizes.csv"), row.names = FALSE)
  msg("consensus modules cached: ", consensus_cache)
  rm(consensus_net, multiExpr, lith_subsample); gc()
}

# refExpr = control (n=234); used as reference set for ALL modulePreservation calls.
# Rationale: control is the largest, most stable group. Consensus moduleColors are
# applied to it. This lets us ask: do consensus modules dissolve differently in
# bp_nolith vs. bp_lith when measured against the healthy-control baseline?
refExpr <- expr_for(d, group = "control")   # 234 x 12368

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

# ---- control self-preservation (one run, 500 perms) ------------------------
# Tests how well consensus modules reproduce in the reference group itself.
# A Zsummary >> 10 here is expected and acts as a sanity check.
ctrl_ckpt <- file.path(run_dir, "ctrl_preservation.rds")
if (!file.exists(ctrl_ckpt)) {
  msg("control self-preservation (consensus modules, ", CTRL_PERMS, " perms) ...")
  # Use a held-out expression set: control itself (no separate test; same data).
  # This is explicitly a sanity check, not a discovery result.
  perm_file <- file.path(run_dir, "permStats_ctrl.RData")
  # For self-preservation, Ref and Test are both control -- Zsummary should be high.
  mp_ctrl <- modulePreservation(
    multiData = list(Ref  = list(data = refExpr),
                     Test = list(data = refExpr)),
    multiColor = list(Ref = moduleColors), dataIsExpr = TRUE,
    networkType = NET_TYPE, corFnc = COR_FNC,
    nPermutations = CTRL_PERMS, randomSeed = SEED_BASE + 500L,
    savePermutedStatistics = TRUE, loadPermutedStatistics = FALSE,
    permutedStatisticsFile = perm_file, verbose = 1, indent = 2)
  cr         <- extract_mp(mp_ctrl)
  cr$group   <- "control"; cr$draw <- NA_integer_
  cr$nperms  <- CTRL_PERMS; cr$n_test <- nrow(refExpr)
  cr$note    <- "self-preservation sanity check; ref=test=control with consensus modules"
  saveRDS(cr, ctrl_ckpt)
  msg("control self-preservation done: Zsummary [",
      round(min(cr$Zsummary), 2), ", ", round(max(cr$Zsummary), 2), "]")
  rm(mp_ctrl); gc()
} else {
  msg("control preservation checkpoint exists, skipping")
}

# ---- bp_nolith (one run, 500 perms) ----------------------------------------
nolith_ckpt <- file.path(run_dir, "nolith_preservation.rds")
if (!file.exists(nolith_ckpt)) {
  msg("bp_nolith preservation (consensus modules, ", NOLITH_PERMS, " perms) ...")
  nolithExpr <- expr_for(d, group = "bp_nolith")   # 74 x 12368
  perm_file  <- file.path(run_dir, "permStats_nolith.RData")
  mp_nolith  <- modulePreservation(
    multiData = list(Ref  = list(data = refExpr),
                     Test = list(data = nolithExpr)),
    multiColor = list(Ref = moduleColors), dataIsExpr = TRUE,
    networkType = NET_TYPE, corFnc = COR_FNC,
    nPermutations = NOLITH_PERMS, randomSeed = SEED_BASE + 1000L,
    savePermutedStatistics = TRUE, loadPermutedStatistics = FALSE,
    permutedStatisticsFile = perm_file, verbose = 1, indent = 2)
  nr         <- extract_mp(mp_nolith)
  nr$group   <- "bp_nolith"; nr$draw <- NA_integer_
  nr$nperms  <- NOLITH_PERMS; nr$n_test <- nrow(nolithExpr)
  saveRDS(nr, nolith_ckpt)
  msg("bp_nolith done: Zsummary [", round(min(nr$Zsummary), 2),
      ", ", round(max(nr$Zsummary), 2), "]")
  rm(nolithExpr, mp_nolith); gc()
} else {
  msg("bp_nolith checkpoint exists, skipping")
}

# ---- bp_lith draws (20 draws x 50 perms, checkpoint per draw) --------------
lith_all_sids <- colnames(d$logtpm)[!is.na(d$groups) & d$groups == "bp_lith"]
stopifnot(length(lith_all_sids) == 152L)

for (draw in seq_len(DRAWS)) {
  ckpt <- file.path(run_dir, sprintf("lith_draw_%03d.rds", draw))
  if (file.exists(ckpt)) { msg("lith draw ", draw, " exists, skipping"); next }
  msg("lith draw ", draw, "/", DRAWS, ": subsampling bp_lith ...")
  set.seed(SEED_BASE + draw)   # same seed offsets as base-001 for matched comparison
  sel_sids  <- sample(lith_all_sids, SUBSAMPLE_N, replace = FALSE)
  testExpr  <- t(d$logtpm[, sel_sids])   # 74 x 12368
  perm_file <- file.path(run_dir, sprintf("permStats_lith_draw_%03d.RData", draw))
  mp_draw <- tryCatch(
    modulePreservation(
      multiData = list(Ref  = list(data = refExpr),
                       Test = list(data = testExpr)),
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
  msg("lith draw ", draw, " done: Zsummary [", round(min(dr$Zsummary), 2),
      ", ", round(max(dr$Zsummary), 2), "]")
  rm(testExpr, mp_draw, dr); gc()
}

# ---- aggregate (only after ALL complete) -----------------------------------
done_ctrl   <- file.exists(ctrl_ckpt)
done_nolith <- file.exists(nolith_ckpt)
done_lith   <- vapply(seq_len(DRAWS), function(i)
  file.exists(file.path(run_dir, sprintf("lith_draw_%03d.rds", i))), logical(1))

if (!done_ctrl || !done_nolith || sum(done_lith) < DRAWS) {
  msg("not all complete: ctrl=", done_ctrl,
      " nolith=", done_nolith,
      " lith=", sum(done_lith), "/", DRAWS)
  quit(save = "no", status = 1L)
}

msg("aggregating results ...")

ctrl_res <- readRDS(ctrl_ckpt)
write.csv(ctrl_res[, c("group", "module", "Zsummary", "medianRank", "draw", "nperms", "n_test")],
          file.path(run_dir, "results_ctrl.csv"), row.names = FALSE)

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

msg("\n--- Consensus module preservation summary (descriptive only) ---")
msg("Control self-preservation (sanity check):")
cr_top <- ctrl_res[order(-ctrl_res$Zsummary), c("module", "Zsummary", "medianRank")]
print(head(cr_top, 5), row.names = FALSE)
msg("bp_nolith (top 5 modules by Zsummary):")
nr_top <- nolith_res[order(-nolith_res$Zsummary), c("module", "Zsummary", "medianRank")]
print(head(nr_top, 5), row.names = FALSE)
msg("bp_lith (top 5 modules by mean Zsummary across draws):")
ls_top <- lith_summary[order(-lith_summary$Zsummary_mean),
                        c("module", "Zsummary_mean", "Zsummary_sd", "medianRank_mean")]
print(head(ls_top, 5), row.names = FALSE)

msg("Results written to ", run_dir)
msg("=== ", id, " COMPLETE ===")
