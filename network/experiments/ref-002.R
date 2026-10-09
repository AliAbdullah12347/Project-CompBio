#!/usr/bin/env Rscript
# ===========================================================================
# ref-002  Alternative-reference network — bp_nolith as reference
#
# Spec (from queue):
#   ref=bp_nolith;corr=bicor;net=signed;genes=all;draws=100;perms_per_draw=50
#   CHECKPOINT PER DRAW (spans sessions)
#
# Design (draws revised to 20 — see network/METHODS.md entry 2026-10-09):
#   1. Build reference modules once on bp_nolith (n=74, all samples, no
#      subsampling for reference). Cache moduleColors to network/cache/.
#   2. bp_lith (n=152 -> subsample 74): 20 draws x 50 perms, checkpoint per draw.
#   3. control (n=234 -> subsample 74): 20 draws x 50 perms, checkpoint per draw.
#      Asks: do bipolar networks preserve in healthy controls?
#   4. Results aggregated only after all 20+20 draws complete; exits 0 then.
#
# Non-negotiable design points obeyed:
#   - Power 12 (WGCNA FAQ signed n>40; guard: reject < 4, fall back to 12)
#   - enableWGCNAThreads()
#   - permutedStatisticsFile set explicitly per draw (never default path)
#   - Zsummary and medianRank are DESCRIPTIVE ONLY: no p-values, no BH correction
#   - Ensembl ids unversioned (handled by 00_load.R)
#   - Config written before any computation
#   - R2 curve recorded (failure to reach R2>=0.80 is itself a result)
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "ref-002"
stopifnot(id == "ref-002")

# Use message() so progress reaches the log immediately (stderr is line-buffered).
now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

# ---------- 1. Directories ---------------------------------------------------
run_dir   <- file.path("network/runs",   id)
cfg_dir   <- "network/config"
cache_dir <- "network/cache"
for (d_ in c(run_dir, cfg_dir, cache_dir))
  dir.create(d_, showWarnings = FALSE, recursive = TRUE)

# ---------- 2. Parameters — ALL pre-specified BEFORE any computation ---------
POWER          <- 12L
POWER_MIN      <-  4L        # guard: reject auto-selected power < 4
NET_TYPE       <- "signed"
COR_FNC        <- "bicor"
GENE_SET       <- "all"
REF_GROUP      <- "bp_nolith"
MIN_MODULE     <- 30L
MERGE_THRESH   <- 0.25
DRAWS          <- 20L        # revised from spec 100 (see METHODS.md 2026-10-09)
PERMS_PER_DRAW <- 50L
NOLITH_PERMS   <- 500L       # kept for reference; not used here (ref IS nolith)
SUBSAMPLE_N    <- 74L
SEED_BASE      <- 20261009L + 2000L   # offset to avoid seed collision with base-001
CACHE_KEY      <- "bp_nolith_all_bicor_signed_p12"

# ---------- 3. Write config BEFORE any computation --------------------------
cfg_path <- file.path(cfg_dir, paste0(id, ".json"))
if (!file.exists(cfg_path)) {
  # Dependency-free JSON serializer for a flat named list.
  scl <- function(v) {
    if (is.logical(v)) return(if (isTRUE(v)) "true" else "false")
    if (is.integer(v) || is.numeric(v)) return(as.character(v[1]))
    paste0('"', gsub('"', '\\"', as.character(v[1]), fixed = TRUE), '"')
  }
  arr <- function(v) paste0("[", paste(vapply(v, scl, ""), collapse = ", "), "]")
  jsn <- function(nm, v) paste0('  "', nm, '": ', if (length(v) > 1) arr(v) else scl(v))
  sha <- tryCatch(
    trimws(system("git rev-parse --short HEAD 2>/dev/null", intern = TRUE)[1]),
    error = function(e) "unknown"
  )
  cfg_entries <- list(
    id                   = id,
    family               = "alt-reference",
    spec_from_queue      = "ref=bp_nolith;corr=bicor;net=signed;genes=all;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans sessions)",
    draws_spec           = 100L,
    draws_implemented    = DRAWS,
    draws_revision_note  = "Same as base-001: revised to 20 for budget (METHODS.md 2026-10-09)",
    ref_group            = REF_GROUP,
    test_groups          = c("bp_lith", "control"),
    corFnc               = COR_FNC,
    networkType          = NET_TYPE,
    power                = POWER,
    power_source         = "WGCNA_FAQ_signed_n_gt_40",
    minModuleSize        = MIN_MODULE,
    mergeCutHeight       = MERGE_THRESH,
    draws_lith           = DRAWS,
    draws_ctrl           = DRAWS,
    perms_per_draw       = PERMS_PER_DRAW,
    nolith_perms         = NOLITH_PERMS,
    seed_base            = SEED_BASE,
    cache_key            = CACHE_KEY,
    config_note          = "bp_nolith (n=74) as reference. Tests: bp_lith (20 draws of 74 from 152) and control (20 draws of 74 from 234). Asks: do bipolar networks preserve in healthy controls?",
    git_sha              = if (length(sha) && !is.na(sha)) sha else "unknown",
    r_version            = paste(R.version$major, R.version$minor, sep = "."),
    wgcna_version        = as.character(packageVersion("WGCNA"))
  )
  lines <- mapply(jsn, names(cfg_entries), cfg_entries, SIMPLIFY = TRUE)
  writeLines(c("{", paste(lines, collapse = ",\n"), "}"), cfg_path)
  msg("config written: ", cfg_path)
}

# ---------- 4. Packages ------------------------------------------------------
suppressPackageStartupMessages({
  library(WGCNA)
  library(fastcluster)
})
enableWGCNAThreads()
options(stringsAsFactors = FALSE)
msg("WGCNA ", as.character(packageVersion("WGCNA")), " loaded; threads enabled")

# ---------- 5. Data ----------------------------------------------------------
source("network/00_load.R")
d <- load_project()   # stops loudly if shape is wrong; verified facts checked inside

# ---------- 6. Scale-free topology scan (R2 curve is itself a result) --------
r2_path <- file.path(run_dir, "sft_r2_curve.csv")
if (!file.exists(r2_path)) {
  msg("pickSoftThreshold scan (powers 1:20) ...")
  set.seed(SEED_BASE)
  tmp_ref <- expr_for(d, group = REF_GROUP)   # 74 x 12368
  sft <- pickSoftThreshold(
    tmp_ref,
    powerVector = 1:20,
    networkType = NET_TYPE,
    corFnc      = COR_FNC,
    verbose     = 0
  )
  est <- sft$powerEstimate
  msg("pickSoftThreshold powerEstimate = ",
      if (is.na(est)) "NA" else as.character(est),
      " (informational only; power fixed at ", POWER, ")")
  if (is.na(est) || est < POWER_MIN) {
    msg("powerEstimate ", if (is.na(est)) "NA" else est,
        " < guard ", POWER_MIN, "; using pre-specified POWER=", POWER)
  }
  r2df               <- sft$fitIndices
  r2df$power_used    <- POWER
  r2df$powerEstimate <- est
  write.csv(r2df, r2_path, row.names = FALSE)
  msg("R2 curve saved (max signed-R2 = ",
      round(max(r2df$SFT.R.sq, na.rm = TRUE), 3), " at power ",
      r2df$Power[which.max(r2df$SFT.R.sq)], ")")
  rm(tmp_ref, sft, r2df); gc()
} else {
  msg("R2 curve checkpoint exists, skipping scan")
}

# ---------- 7. Reference network (module detection, cached) ------------------
ref_cache <- file.path(cache_dir, paste0("ref_", CACHE_KEY, ".rds"))
if (file.exists(ref_cache)) {
  msg("loading cached reference modules: ", ref_cache)
  ref_net      <- readRDS(ref_cache)
  moduleColors <- ref_net$moduleColors
  msg(ref_net$n_modules, " modules loaded (",
      sum(moduleColors != "grey"), " genes in non-grey modules)")
} else {
  msg("building reference network (bicor, signed, power ", POWER, ") on ",
      REF_GROUP, " (n=74, all samples) ...")
  tmp_ref <- expr_for(d, group = REF_GROUP)   # 74 x 12368
  set.seed(SEED_BASE)
  bwm <- blockwiseModules(
    tmp_ref,
    power          = POWER,
    networkType    = NET_TYPE,
    TOMType        = NET_TYPE,
    corType        = COR_FNC,
    minModuleSize  = MIN_MODULE,
    mergeCutHeight = MERGE_THRESH,
    maxBlockSize   = 15000L,   # keep all 12368 genes in one block
    numericLabels  = FALSE,
    saveTOMs       = FALSE,
    verbose        = 3
  )
  moduleColors <- bwm$colors
  n_mods <- length(unique(moduleColors)) - 1L   # subtract 1 for grey
  msg("reference network: ", n_mods, " modules detected (excluding grey); ",
      sum(moduleColors == "grey"), " unassigned genes")
  ref_net <- list(
    moduleColors = moduleColors,
    n_modules    = n_mods,
    module_sizes = table(moduleColors),
    power        = POWER,
    cache_key    = CACHE_KEY,
    corFnc       = COR_FNC,
    networkType  = NET_TYPE,
    ref_group    = REF_GROUP,
    gene_set     = GENE_SET
  )
  saveRDS(ref_net, ref_cache)
  msg("reference modules cached: ", ref_cache)
  # Module sizes as tracked CSV (small, useful for later rows)
  ms_path <- file.path(run_dir, "reference_module_sizes.csv")
  write.csv(as.data.frame(ref_net$module_sizes), ms_path, row.names = FALSE)
  rm(bwm, tmp_ref); gc()
}

# refExpr is needed by every modulePreservation call — build once here.
refExpr <- expr_for(d, group = REF_GROUP)   # 74 x 12368

# ---------- Helper: extract Zsummary + medianRank from modulePreservation ----
# Returns a data.frame with columns: module, Zsummary, medianRank.
# The grey module (unassigned genes) is excluded from the result.
extract_mp <- function(mp) {
  obs  <- mp$preservation$observed[["ref.Ref"]][["inColumnsAlsoPresentIn.Test"]]
  zsco <- mp$preservation$Z[["ref.Ref"]][["inColumnsAlsoPresentIn.Test"]]
  if (is.null(obs) || is.null(zsco)) {
    stop("modulePreservation return structure not as expected; ",
         "preservation$observed and/or preservation$Z absent for ref.Ref -> Test")
  }
  res <- data.frame(
    module     = rownames(obs),
    Zsummary   = zsco[, "Zsummary.pres"],
    medianRank = obs[,  "medianRank.pres"],
    stringsAsFactors = FALSE
  )
  res[res$module != "grey", ]
}

# ---------- 8. bp_lith draws (20 draws x 50 perms, checkpoint per draw) -----
# bp_lith has n=152; subsample to 74 each draw to match reference size.
lith_sids <- colnames(d$logtpm)[!is.na(d$groups) & d$groups == "bp_lith"]
stopifnot(length(lith_sids) == 152L)

for (draw in seq_len(DRAWS)) {
  ckpt <- file.path(run_dir, sprintf("lith_draw_%03d.rds", draw))
  if (file.exists(ckpt)) {
    msg("lith draw ", draw, "/", DRAWS, " checkpoint exists, skipping")
    next
  }
  msg("lith draw ", draw, "/", DRAWS,
      ": subsampling bp_lith (n=", SUBSAMPLE_N, " from 152) ...")
  set.seed(SEED_BASE + draw)
  sel_sids <- sample(lith_sids, SUBSAMPLE_N, replace = FALSE)
  testExpr <- t(d$logtpm[, sel_sids])   # 74 x 12368

  perm_file <- file.path(run_dir, sprintf("permStats_lith_draw_%03d.RData", draw))

  mp_draw <- tryCatch({
    modulePreservation(
      multiData = list(
        Ref  = list(data = refExpr),
        Test = list(data = testExpr)
      ),
      multiColor = list(Ref = moduleColors),
      dataIsExpr             = TRUE,
      networkType            = NET_TYPE,
      corFnc                 = COR_FNC,
      nPermutations          = PERMS_PER_DRAW,
      randomSeed             = SEED_BASE + draw,
      savePermutedStatistics = TRUE,
      loadPermutedStatistics = FALSE,
      permutedStatisticsFile = perm_file,
      verbose                = 1,
      indent                 = 2
    )
  }, error = function(e) {
    msg("ERROR in lith draw ", draw, ": ", conditionMessage(e))
    NULL
  })

  if (is.null(mp_draw)) {
    msg("lith draw ", draw, " failed; will retry in a later session if row is requeued")
    rm(testExpr); gc()
    next
  }

  draw_res         <- extract_mp(mp_draw)
  draw_res$group   <- "bp_lith"
  draw_res$draw    <- draw
  draw_res$nperms  <- PERMS_PER_DRAW
  draw_res$n_test  <- SUBSAMPLE_N
  draw_res$samples <- paste(sel_sids, collapse = ",")   # exact IDs for reproducibility
  saveRDS(draw_res, ckpt)
  msg("lith draw ", draw, "/", DRAWS, " done: Zsummary range [",
      round(min(draw_res$Zsummary), 2), ", ",
      round(max(draw_res$Zsummary), 2), "]")
  rm(testExpr, mp_draw, draw_res); gc()
}

# ---------- 9. control draws (20 draws x 50 perms, checkpoint per draw) -----
# control has n=234; subsample to 74 each draw to match reference size.
# Use SEED_BASE + draw + 1000 for ctrl draws to avoid seed collision with lith draws.
ctrl_sids <- colnames(d$logtpm)[!is.na(d$groups) & d$groups == "control"]
stopifnot(length(ctrl_sids) == 234L)

for (draw in seq_len(DRAWS)) {
  ckpt <- file.path(run_dir, sprintf("ctrl_draw_%03d.rds", draw))
  if (file.exists(ckpt)) {
    msg("ctrl draw ", draw, "/", DRAWS, " checkpoint exists, skipping")
    next
  }
  msg("ctrl draw ", draw, "/", DRAWS,
      ": subsampling control (n=", SUBSAMPLE_N, " from 234) ...")
  set.seed(SEED_BASE + draw + 1000L)
  sel_sids <- sample(ctrl_sids, SUBSAMPLE_N, replace = FALSE)
  testExpr <- t(d$logtpm[, sel_sids])   # 74 x 12368

  perm_file <- file.path(run_dir, sprintf("permStats_ctrl_draw_%03d.RData", draw))

  mp_draw <- tryCatch({
    modulePreservation(
      multiData = list(
        Ref  = list(data = refExpr),
        Test = list(data = testExpr)
      ),
      multiColor = list(Ref = moduleColors),
      dataIsExpr             = TRUE,
      networkType            = NET_TYPE,
      corFnc                 = COR_FNC,
      nPermutations          = PERMS_PER_DRAW,
      randomSeed             = SEED_BASE + draw + 1000L,
      savePermutedStatistics = TRUE,
      loadPermutedStatistics = FALSE,
      permutedStatisticsFile = perm_file,
      verbose                = 1,
      indent                 = 2
    )
  }, error = function(e) {
    msg("ERROR in ctrl draw ", draw, ": ", conditionMessage(e))
    NULL
  })

  if (is.null(mp_draw)) {
    msg("ctrl draw ", draw, " failed; will retry in a later session if row is requeued")
    rm(testExpr); gc()
    next
  }

  draw_res         <- extract_mp(mp_draw)
  draw_res$group   <- "control"
  draw_res$draw    <- draw
  draw_res$nperms  <- PERMS_PER_DRAW
  draw_res$n_test  <- SUBSAMPLE_N
  draw_res$samples <- paste(sel_sids, collapse = ",")   # exact IDs for reproducibility
  saveRDS(draw_res, ckpt)
  msg("ctrl draw ", draw, "/", DRAWS, " done: Zsummary range [",
      round(min(draw_res$Zsummary), 2), ", ",
      round(max(draw_res$Zsummary), 2), "]")
  rm(testExpr, mp_draw, draw_res); gc()
}

# ---------- 10. Aggregate (only after ALL draws complete) -------------------
done_lith <- vapply(seq_len(DRAWS), function(i)
  file.exists(file.path(run_dir, sprintf("lith_draw_%03d.rds", i))), logical(1))
done_ctrl <- vapply(seq_len(DRAWS), function(i)
  file.exists(file.path(run_dir, sprintf("ctrl_draw_%03d.rds", i))), logical(1))
n_lith_done <- sum(done_lith)
n_ctrl_done <- sum(done_ctrl)

if (n_lith_done < DRAWS || n_ctrl_done < DRAWS) {
  msg("not all work complete: lith draws done=", n_lith_done, "/", DRAWS,
      ", ctrl draws done=", n_ctrl_done, "/", DRAWS,
      "; exiting for resume in a later session")
  # Exit non-zero so the driver marks this PARTIAL, checks checkpoint progress,
  # refunds the attempt if new .rds files appeared, and requeues.
  quit(save = "no", status = 1L)
}

msg("all draws complete; aggregating results ...")

# --- bp_lith: per-draw data and per-module summary across draws ---
lith_all <- do.call(rbind, lapply(seq_len(DRAWS), function(i) {
  x <- readRDS(file.path(run_dir, sprintf("lith_draw_%03d.rds", i)))
  x[, c("group", "module", "Zsummary", "medianRank", "draw", "nperms", "n_test")]
}))
write.csv(lith_all, file.path(run_dir, "results_lith_alldraws.csv"), row.names = FALSE)

modules      <- unique(lith_all$module)
lith_summary <- do.call(rbind, lapply(modules, function(m) {
  s <- lith_all[lith_all$module == m, ]
  data.frame(
    group             = "bp_lith",
    module            = m,
    n_draws           = nrow(s),
    Zsummary_mean     = mean(s$Zsummary,   na.rm = TRUE),
    Zsummary_median   = median(s$Zsummary, na.rm = TRUE),
    Zsummary_sd       = sd(s$Zsummary,     na.rm = TRUE),
    Zsummary_min      = min(s$Zsummary,    na.rm = TRUE),
    Zsummary_max      = max(s$Zsummary,    na.rm = TRUE),
    medianRank_mean   = mean(s$medianRank,   na.rm = TRUE),
    medianRank_median = median(s$medianRank, na.rm = TRUE),
    medianRank_sd     = sd(s$medianRank,     na.rm = TRUE),
    stringsAsFactors  = FALSE
  )
}))
write.csv(lith_summary, file.path(run_dir, "results_lith_summary.csv"), row.names = FALSE)

# --- control: per-draw data and per-module summary across draws ---
ctrl_all <- do.call(rbind, lapply(seq_len(DRAWS), function(i) {
  x <- readRDS(file.path(run_dir, sprintf("ctrl_draw_%03d.rds", i)))
  x[, c("group", "module", "Zsummary", "medianRank", "draw", "nperms", "n_test")]
}))
write.csv(ctrl_all, file.path(run_dir, "results_ctrl_alldraws.csv"), row.names = FALSE)

modules      <- unique(ctrl_all$module)
ctrl_summary <- do.call(rbind, lapply(modules, function(m) {
  s <- ctrl_all[ctrl_all$module == m, ]
  data.frame(
    group             = "control",
    module            = m,
    n_draws           = nrow(s),
    Zsummary_mean     = mean(s$Zsummary,   na.rm = TRUE),
    Zsummary_median   = median(s$Zsummary, na.rm = TRUE),
    Zsummary_sd       = sd(s$Zsummary,     na.rm = TRUE),
    Zsummary_min      = min(s$Zsummary,    na.rm = TRUE),
    Zsummary_max      = max(s$Zsummary,    na.rm = TRUE),
    medianRank_mean   = mean(s$medianRank,   na.rm = TRUE),
    medianRank_median = median(s$medianRank, na.rm = TRUE),
    medianRank_sd     = sd(s$medianRank,     na.rm = TRUE),
    stringsAsFactors  = FALSE
  )
}))
write.csv(ctrl_summary, file.path(run_dir, "results_ctrl_summary.csv"), row.names = FALSE)

# Print tables for the log
msg("\n--- bp_lith Zsummary_mean across ", DRAWS, " draws (top 10, descriptive only) ---")
show_l <- lith_summary[order(-lith_summary$Zsummary_mean),
                       c("module", "Zsummary_mean", "Zsummary_sd", "medianRank_mean")]
print(head(show_l, 10), row.names = FALSE)

msg("\n--- control Zsummary_mean across ", DRAWS, " draws (top 10, descriptive only) ---")
show_c <- ctrl_summary[order(-ctrl_summary$Zsummary_mean),
                       c("module", "Zsummary_mean", "Zsummary_sd", "medianRank_mean")]
print(head(show_c, 10), row.names = FALSE)

msg("=== ", id, " COMPLETE ===")
# Exit 0: driver marks the row DONE.
