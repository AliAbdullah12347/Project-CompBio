#!/usr/bin/env Rscript
# ===========================================================================
# base-001  Baseline reference network — descriptive Zsummary + medianRank
#
# Spec (from queue):
#   ref=control; corr=bicor; net=signed; genes=all; draws=100; perms_per_draw=50
#   CHECKPOINT PER DRAW (spans sessions); stat=Zsummary+medianRank; DESCRIPTIVE ONLY
#
# Design (draws revised to 20 — see network/METHODS.md entry 2026-10-09):
#   1. Build reference modules once on control (n=234, bicor, signed, power 12).
#      Cache moduleColors to network/cache/ for reuse by later rows.
#   2. bp_nolith (n=74, never subsampled): ONE modulePreservation run, 500 perms.
#      This is equivalent to 100 draws × 50 perms on a fixed test set.
#   3. bp_lith (n=152 → subsample 74): 20 draws × 50 perms, checkpoint per draw.
#      Draws=20 matches the intended ~11 h budget given the actual ~1 962 s/draw
#      cost (ref TOM + test TOM + 50 perms, all inside modulePreservation).
#   4. Results aggregated only after all draws complete; exits 0 then.
#
# Non-negotiable design points obeyed:
#   - Power 12 (WGCNA FAQ signed n>40; guard: reject < 4, fall back to 12)
#   - enableWGCNAThreads()
#   - permutedStatisticsFile set explicitly per draw (never default path)
#   - Zsummary and medianRank are DESCRIPTIVE ONLY: no p-values, no BH correction
#   - Ensembl ids unversioned (handled by 00_load.R)
#   - Config written before any computation
#   - R2 curve recorded every time; failure to reach 0.80 is itself a result
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "base-001"
stopifnot(id == "base-001")

ts_str <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
log    <- function(...) { cat(sprintf("[%s] %s\n", ts_str(), paste0(..., collapse = ""))); flush.stdout() }

log("=== ", id, " starting ===")

# ---------- 1. Directories ---------------------------------------------------
run_dir   <- file.path("network/runs",   id)
cfg_dir   <- "network/config"
cache_dir <- "network/cache"
for (d_ in c(run_dir, cfg_dir, cache_dir))
  dir.create(d_, showWarnings = FALSE, recursive = TRUE)

# ---------- 2. Parameters — ALL pre-specified BEFORE any computation ---------
POWER          <- 12L
POWER_MIN      <-  4L       # guard: reject auto-selected power < 4
NET_TYPE       <- "signed"
COR_FNC        <- "bicor"
GENE_SET       <- "all"
REF_GROUP      <- "control"
MIN_MODULE     <- 30L
MERGE_THRESH   <- 0.25
DRAWS          <- 20L       # revised from spec 100 (see METHODS.md 2026-10-09)
PERMS_PER_DRAW <- 50L
NOLITH_PERMS   <- 500L      # bp_nolith: one run, many perms (no subsampling variance)
SUBSAMPLE_N    <- 74L
SEED_BASE      <- 20261009L
CACHE_KEY      <- paste(REF_GROUP, GENE_SET, COR_FNC, NET_TYPE,
                        paste0("p", POWER), sep = "_")

# ---------- 3. Write config BEFORE any computation --------------------------
cfg_path <- file.path(cfg_dir, paste0(id, ".json"))
if (!file.exists(cfg_path)) {
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
    family               = "baseline",
    spec_from_queue      = "ref=control;corr=bicor;net=signed;genes=all;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans sessions);stat=Zsummary+medianRank;DESCRIPTIVE ONLY",
    draws_spec           = 100L,
    draws_implemented    = DRAWS,
    draws_revision_note  = "100 draws x 1962 s/draw = 54 h; revised to 20 to fit budget (see METHODS.md 2026-10-09)",
    ref_group            = REF_GROUP,
    test_groups          = c("bp_nolith", "bp_lith"),
    corFnc               = COR_FNC,
    networkType          = NET_TYPE,
    gene_set             = GENE_SET,
    power                = POWER,
    power_source         = "WGCNA_FAQ_signed_n_gt_40",
    power_guard          = paste0("reject_lt_", POWER_MIN, "_fallback_", POWER),
    minModuleSize        = MIN_MODULE,
    mergeCutHeight       = MERGE_THRESH,
    draws_lith           = DRAWS,
    perms_per_lith_draw  = PERMS_PER_DRAW,
    nolith_perms         = NOLITH_PERMS,
    nolith_draws         = 1L,
    subsample_n          = SUBSAMPLE_N,
    bp_lith_subsampled   = TRUE,
    bp_nolith_subsampled = FALSE,
    stat                 = c("Zsummary", "medianRank"),
    descriptive_only     = TRUE,
    seed_base            = SEED_BASE,
    cache_key            = CACHE_KEY,
    git_sha              = if (length(sha) && !is.na(sha)) sha else "unknown",
    r_version            = paste(R.version$major, R.version$minor, sep = "."),
    wgcna_version        = as.character(packageVersion("WGCNA"))
  )
  lines <- mapply(jsn, names(cfg_entries), cfg_entries, SIMPLIFY = TRUE)
  writeLines(c("{", paste(lines, collapse = ",\n"), "}"), cfg_path)
  log("config written: ", cfg_path)
}

# ---------- 4. Packages ------------------------------------------------------
suppressPackageStartupMessages({
  library(WGCNA)
  library(fastcluster)
})
enableWGCNAThreads()
options(stringsAsFactors = FALSE)
log("WGCNA ", as.character(packageVersion("WGCNA")), " loaded; threads enabled")

# ---------- 5. Data ----------------------------------------------------------
source("network/00_load.R")
d <- load_project()   # stops if shape is wrong; verified facts are checked inside

# ---------- 6. Scale-free topology scan (R2 curve is itself a result) --------
r2_path <- file.path(run_dir, "sft_r2_curve.csv")
if (!file.exists(r2_path)) {
  log("pickSoftThreshold scan (powers 1:20) ...")
  set.seed(SEED_BASE)
  refExpr_tmp <- expr_for(d, group = REF_GROUP)   # 234 x 12368
  sft <- pickSoftThreshold(
    refExpr_tmp,
    powerVector = 1:20,
    networkType = NET_TYPE,
    corFnc      = COR_FNC,
    verbose     = 0
  )
  est <- sft$powerEstimate
  log("pickSoftThreshold powerEstimate = ",
      if (is.na(est)) "NA" else as.character(est),
      " (expected to be 1 or NA — function failing on this data, not a threshold)")
  if (is.na(est) || est < POWER_MIN) {
    log("powerEstimate ", if (is.na(est)) "NA" else est,
        " < guard ", POWER_MIN, "; using pre-specified POWER=", POWER)
  }
  r2df <- sft$fitIndices
  r2df$power_used    <- POWER
  r2df$power_spec    <- POWER
  r2df$powerEstimate <- est
  write.csv(r2df, r2_path, row.names = FALSE)
  log("R2 curve saved: ", r2_path, " (max signed-R2=",
      round(max(r2df$SFT.R.sq, na.rm = TRUE), 3), " at power ",
      r2df$Power[which.max(r2df$SFT.R.sq)], ")")
  rm(refExpr_tmp, sft, r2df); gc()
} else {
  log("R2 curve checkpoint exists, skipping scan")
}

# ---------- 7. Reference network (module detection, cached) ------------------
ref_cache <- file.path(cache_dir, paste0("ref_", CACHE_KEY, ".rds"))
if (file.exists(ref_cache)) {
  log("loading cached reference modules: ", ref_cache)
  ref_net      <- readRDS(ref_cache)
  moduleColors <- ref_net$moduleColors
  log(ref_net$n_modules, " modules loaded (",
      sum(moduleColors != "grey"), " genes in non-grey modules)")
} else {
  log("building reference network (bicor, signed, power ", POWER, ") ...")
  refExpr <- expr_for(d, group = REF_GROUP)   # 234 x 12368
  set.seed(SEED_BASE)
  bwm <- blockwiseModules(
    refExpr,
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
  n_mods <- length(unique(moduleColors)) - 1L   # -1 for grey
  log("reference network: ", n_mods, " modules (excl grey); ",
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
  log("reference modules cached: ", ref_cache)
  # Save module size table as a tracked CSV
  ms_path <- file.path(run_dir, "reference_module_sizes.csv")
  write.csv(as.data.frame(ref_net$module_sizes), ms_path, row.names = FALSE)
  rm(bwm, refExpr); gc()
}

refExpr <- expr_for(d, group = REF_GROUP)   # 234 x 12368 (needed by modulePreservation)

# Helper: extract Zsummary and medianRank from a modulePreservation result.
# Returns a data.frame with columns: module, Zsummary, medianRank.
# Grey module is excluded (it is the "background", not a real module).
extract_mp <- function(mp) {
  obs  <- mp$preservation$observed[["ref.Ref"]][["inColumnsAlsoPresentIn.Test"]]
  zsco <- mp$preservation$Z[["ref.Ref"]][["inColumnsAlsoPresentIn.Test"]]
  if (is.null(obs) || is.null(zsco)) stop("modulePreservation return structure unexpected")
  result <- data.frame(
    module     = rownames(obs),
    Zsummary   = zsco[, "Zsummary.pres"],
    medianRank = obs[, "medianRank.pres"],
    stringsAsFactors = FALSE
  )
  result[result$module != "grey", ]
}

# ---------- 8. bp_nolith preservation (one run, 500 perms) ------------------
nolith_ckpt <- file.path(run_dir, "nolith_preservation.rds")
if (!file.exists(nolith_ckpt)) {
  log("bp_nolith preservation (500 perms) ...")
  nolithExpr <- expr_for(d, group = "bp_nolith")   # 74 x 12368
  perm_file  <- file.path(run_dir, "permStats_nolith.RData")
  set.seed(SEED_BASE + 1000L)
  mp_nolith <- modulePreservation(
    multiData = list(
      Ref  = list(data = refExpr),
      Test = list(data = nolithExpr)
    ),
    multiColor = list(Ref = moduleColors),
    dataIsExpr           = TRUE,
    networkType          = NET_TYPE,
    corFnc               = COR_FNC,
    nPermutations        = NOLITH_PERMS,
    randomSeed           = SEED_BASE + 1000L,
    savePermutedStatistics = TRUE,
    loadPermutedStatistics = FALSE,
    permutedStatisticsFile = perm_file,
    verbose              = 1,
    indent               = 2
  )
  nolith_res <- extract_mp(mp_nolith)
  nolith_res$group      <- "bp_nolith"
  nolith_res$draw       <- NA_integer_
  nolith_res$nperms     <- NOLITH_PERMS
  nolith_res$n_test     <- nrow(nolithExpr)
  saveRDS(nolith_res, nolith_ckpt)
  log("bp_nolith done: ", nrow(nolith_res), " modules; Zsummary range [",
      round(min(nolith_res$Zsummary), 2), ", ",
      round(max(nolith_res$Zsummary), 2), "]")
  rm(nolithExpr, mp_nolith); gc()
} else {
  log("bp_nolith checkpoint exists, skipping")
}

# ---------- 9. bp_lith draws (20 draws × 50 perms, checkpoint per draw) ----
lith_sids <- colnames(d$logtpm)[!is.na(d$groups) & d$groups == "bp_lith"]
stopifnot(length(lith_sids) == 152L)

for (draw in seq_len(DRAWS)) {
  ckpt <- file.path(run_dir, sprintf("lith_draw_%03d.rds", draw))
  if (file.exists(ckpt)) {
    log("draw ", draw, "/", DRAWS, " checkpoint exists, skipping")
    next
  }
  log("draw ", draw, "/", DRAWS, ": subsampling bp_lith ...")
  set.seed(SEED_BASE + draw)
  sel_sids  <- sample(lith_sids, SUBSAMPLE_N, replace = FALSE)
  testExpr  <- t(d$logtpm[, sel_sids])    # 74 x 12368
  perm_file <- file.path(run_dir, sprintf("permStats_lith_draw_%03d.RData", draw))

  mp_draw <- tryCatch({
    set.seed(SEED_BASE + draw)
    modulePreservation(
      multiData = list(
        Ref  = list(data = refExpr),
        Test = list(data = testExpr)
      ),
      multiColor = list(Ref = moduleColors),
      dataIsExpr           = TRUE,
      networkType          = NET_TYPE,
      corFnc               = COR_FNC,
      nPermutations        = PERMS_PER_DRAW,
      randomSeed           = SEED_BASE + draw,
      savePermutedStatistics = TRUE,
      loadPermutedStatistics = FALSE,
      permutedStatisticsFile = perm_file,
      verbose              = 1,
      indent               = 2
    )
  }, error = function(e) {
    log("ERROR in draw ", draw, ": ", conditionMessage(e))
    NULL
  })

  if (is.null(mp_draw)) {
    log("draw ", draw, " failed; skipping (will retry in a later session if the row is requeued)")
    rm(testExpr); gc()
    next
  }

  draw_res <- extract_mp(mp_draw)
  draw_res$group    <- "bp_lith"
  draw_res$draw     <- draw
  draw_res$nperms   <- PERMS_PER_DRAW
  draw_res$n_test   <- SUBSAMPLE_N
  draw_res$samples  <- paste(sel_sids, collapse = ",")   # reproducibility: exact IDs used
  saveRDS(draw_res, ckpt)
  log("draw ", draw, "/", DRAWS, " done: Zsummary range [",
      round(min(draw_res$Zsummary), 2), ", ",
      round(max(draw_res$Zsummary), 2), "]")
  rm(testExpr, mp_draw, draw_res); gc()
}

# ---------- 10. Aggregate (only after ALL draws complete) -------------------
# Check whether all draws and nolith are done
nolith_ckpt <- file.path(run_dir, "nolith_preservation.rds")
done_draws  <- vapply(seq_len(DRAWS), function(i)
  file.exists(file.path(run_dir, sprintf("lith_draw_%03d.rds", i))), logical(1))
n_done      <- sum(done_draws)

if (!file.exists(nolith_ckpt) || n_done < DRAWS) {
  log("not all work complete: nolith=", file.exists(nolith_ckpt),
      ", lith draws done=", n_done, "/", DRAWS,
      "; exiting for resume in a later session")
  quit(status = 1)   # PARTIAL; driver will requeue if progress was made
}

log("all draws complete; aggregating results ...")

# bp_nolith
nolith_res  <- readRDS(nolith_ckpt)
nolith_res$samples <- NA_character_
nolith_out  <- nolith_res[, c("group", "module", "Zsummary", "medianRank",
                               "draw", "nperms", "n_test")]
write.csv(nolith_out,
          file.path(run_dir, "results_nolith.csv"), row.names = FALSE)

# bp_lith: per-draw data and summary across draws
lith_all <- do.call(rbind, lapply(seq_len(DRAWS), function(i) {
  x <- readRDS(file.path(run_dir, sprintf("lith_draw_%03d.rds", i)))
  x[, c("group", "module", "Zsummary", "medianRank", "draw", "nperms", "n_test")]
}))
write.csv(lith_all,
          file.path(run_dir, "results_lith_alldraws.csv"), row.names = FALSE)

# Per-module summary: mean, median, sd of Zsummary and medianRank across 20 draws
modules <- unique(lith_all$module)
lith_summary <- do.call(rbind, lapply(modules, function(m) {
  sub <- lith_all[lith_all$module == m, ]
  data.frame(
    group             = "bp_lith",
    module            = m,
    n_draws           = nrow(sub),
    Zsummary_mean     = mean(sub$Zsummary,   na.rm = TRUE),
    Zsummary_median   = median(sub$Zsummary, na.rm = TRUE),
    Zsummary_sd       = sd(sub$Zsummary,     na.rm = TRUE),
    Zsummary_min      = min(sub$Zsummary,    na.rm = TRUE),
    Zsummary_max      = max(sub$Zsummary,    na.rm = TRUE),
    medianRank_mean   = mean(sub$medianRank,   na.rm = TRUE),
    medianRank_median = median(sub$medianRank, na.rm = TRUE),
    medianRank_sd     = sd(sub$medianRank,     na.rm = TRUE),
    stringsAsFactors  = FALSE
  )
}))
write.csv(lith_summary,
          file.path(run_dir, "results_lith_summary.csv"), row.names = FALSE)

# Combined table for RESULTS.md: one row per module per group
nolith_sum <- nolith_out[, c("group", "module", "Zsummary", "medianRank")]
names(nolith_sum)[names(nolith_sum) == "Zsummary"]   <- "Zsummary_point"
names(nolith_sum)[names(nolith_sum) == "medianRank"] <- "medianRank_point"
nolith_sum$Zsummary_mean <- nolith_sum$Zsummary_point
nolith_sum$note <- paste0("single run, ", NOLITH_PERMS, " perms")

lith_comb <- data.frame(
  group              = lith_summary$group,
  module             = lith_summary$module,
  Zsummary_point     = lith_summary$Zsummary_median,
  medianRank_point   = lith_summary$medianRank_median,
  Zsummary_mean      = lith_summary$Zsummary_mean,
  note               = paste0(DRAWS, " draws x ", PERMS_PER_DRAW, " perms"),
  stringsAsFactors   = FALSE
)

combined <- rbind(
  nolith_sum[, intersect(names(nolith_sum), names(lith_comb))],
  lith_comb[, intersect(names(nolith_sum), names(lith_comb))]
)
write.csv(combined,
          file.path(run_dir, "results_combined.csv"), row.names = FALSE)

log("results written:")
log("  ", file.path(run_dir, "results_nolith.csv"))
log("  ", file.path(run_dir, "results_lith_alldraws.csv"))
log("  ", file.path(run_dir, "results_lith_summary.csv"))
log("  ", file.path(run_dir, "results_combined.csv"))

# Print a brief table of Zsummary statistics for the log
log("\n--- bp_nolith Zsummary (500 perms) ---")
nolith_show <- nolith_out[order(-nolith_out$Zsummary), c("module","Zsummary","medianRank")]
print(head(nolith_show, 10), row.names = FALSE)

log("\n--- bp_lith Zsummary_mean across ", DRAWS, " draws (top 10) ---")
lith_show <- lith_summary[order(-lith_summary$Zsummary_mean),
                           c("module","Zsummary_mean","Zsummary_sd","medianRank_mean")]
print(head(lith_show, 10), row.names = FALSE)

log("=== ", id, " COMPLETE ===")
