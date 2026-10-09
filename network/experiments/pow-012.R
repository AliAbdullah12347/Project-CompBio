#!/usr/bin/env Rscript
# ===========================================================================
# pow-012  Power=12 fixed sensitivity (independent reproducibility check)
#          — Zsummary + medianRank
#
# Spec (from queue):
#   ref=control;power=12 fixed;corr=bicor;net=signed;genes=all;draws=100;perms_per_draw=50
#   CHECKPOINT PER DRAW (spans sessions)
#
# Design (draws revised to 20 — see network/METHODS.md entry 2026-10-09):
#   1. Build reference modules once on control (n=234, bicor, signed, power 12).
#      Cache moduleColors to network/cache/ under a DISTINCT cache key from base-001
#      so that pow-012 builds its OWN reference independently (needed for a clean
#      reproducibility comparison — if both scripts share a cache, they trivially agree).
#   2. bp_nolith (n=74, never subsampled): ONE modulePreservation run, 500 perms.
#   3. bp_lith (n=152 -> subsample 74): 20 draws x 50 perms, checkpoint per draw.
#   4. Results aggregated only after all draws complete; exits 0 then.
#
# Sensitivity vs base-001:
#   Power=12 explicitly fixed (same as base-001's FAQ-recommended value). Builds a
#   separate independent network from the same parameters using an offset seed. If
#   results match base-001, confirms reproducibility. If they differ, seed differences
#   matter. No power guard is applied (power_source = "fixed_sensitivity"), matching
#   the convention established by pow-010 and pow-011.
#
# Non-negotiable design points obeyed:
#   - Power 12 (intentionally fixed; no guard applied — this IS the sensitivity test)
#   - enableWGCNAThreads()
#   - permutedStatisticsFile set explicitly per draw (never default path)
#   - Zsummary and medianRank are DESCRIPTIVE ONLY: no p-values, no BH correction
#   - Ensembl ids unversioned (handled by 00_load.R)
#   - Config written before any computation
#   - R2 curve recorded (failure to reach R2>=0.80 is itself a result)
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "pow-012"
stopifnot(id == "pow-012")

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
# NOTE: POWER is intentionally fixed at 12 for this sensitivity script.
# No power guard is applied — the point of this script IS the fixed-power test.
# SEED_BASE is offset by +12000 from the base-001 seed so draws are independent.
POWER          <- 12L
NET_TYPE       <- "signed"
COR_FNC        <- "bicor"
GENE_SET       <- "all"
REF_GROUP      <- "control"
MIN_MODULE     <- 30L
MERGE_THRESH   <- 0.25
DRAWS          <- 20L        # revised from spec 100 (see METHODS.md 2026-10-09)
PERMS_PER_DRAW <- 50L
NOLITH_PERMS   <- 500L       # bp_nolith: one run, many perms (no subsampling variance)
SUBSAMPLE_N    <- 74L
SEED_BASE      <- 20261009L + 12000L   # offset from base-001 for independence
CACHE_KEY      <- "control_all_bicor_signed_p12_fixed"
# Note: "_fixed" suffix distinguishes from base-001's "control_all_bicor_signed_p12"
# so pow-012 builds its OWN reference; comparison is then clean.

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
    family               = "sensitivity_power",
    spec_from_queue      = "ref=control;power=12 fixed;corr=bicor;net=signed;genes=all;draws=100;perms_per_draw=50;CHECKPOINT PER DRAW (spans sessions)",
    draws_spec           = 100L,
    draws_implemented    = DRAWS,
    draws_revision_note  = "Same as base-001: revised to 20 for budget (METHODS.md 2026-10-09)",
    ref_group            = REF_GROUP,
    test_groups          = c("bp_nolith", "bp_lith"),
    corFnc               = COR_FNC,
    networkType          = NET_TYPE,
    gene_set             = GENE_SET,
    power                = POWER,
    power_source         = "fixed_sensitivity",
    power_faq_note       = "FAQ recommended power for signed n>40 is 12. pow-012 explicitly fixes at 12 and builds an independent network; compare vs base-001 to confirm reproducibility",
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
    sensitivity_axis     = "power",
    sensitivity_note     = "Power=12 explicitly fixed (same as base-001's FAQ-recommended value); builds a separate independent network from the same parameters; if results match base-001, confirms reproducibility. If they differ, seed differences matter.",
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
# NOTE: R2 curve is saved for information but does NOT override POWER=12.
r2_path <- file.path(run_dir, "sft_r2_curve.csv")
if (!file.exists(r2_path)) {
  msg("pickSoftThreshold scan (powers 1:20) — informational only; POWER fixed at ", POWER, " ...")
  set.seed(SEED_BASE)
  tmp_ref <- expr_for(d, group = REF_GROUP)   # 234 x 12368
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
      " (informational only; this script uses fixed POWER=", POWER, ")")
  r2df               <- sft$fitIndices
  r2df$power_used    <- POWER
  r2df$powerEstimate <- est
  write.csv(r2df, r2_path, row.names = FALSE)
  msg("R2 curve saved (max signed-R2 = ",
      round(max(r2df$SFT.R.sq, na.rm = TRUE), 3), " at power ",
      r2df$Power[which.max(r2df$SFT.R.sq)], "; using fixed POWER=", POWER, ")")
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
  msg("building reference network (bicor, signed, power ", POWER, " fixed) ...")
  tmp_ref <- expr_for(d, group = REF_GROUP)   # 234 x 12368
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
  ms_path <- file.path(run_dir, "reference_module_sizes.csv")
  write.csv(as.data.frame(ref_net$module_sizes), ms_path, row.names = FALSE)
  rm(bwm, tmp_ref); gc()
}

# refExpr is needed by every modulePreservation call — build once here.
refExpr <- expr_for(d, group = REF_GROUP)   # 234 x 12368

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

# ---------- 8. bp_nolith preservation (one run, 500 perms) ------------------
nolith_ckpt <- file.path(run_dir, "nolith_preservation.rds")
if (!file.exists(nolith_ckpt)) {
  msg("bp_nolith preservation: 500 perms (fixed test set, no subsampling) ...")
  nolithExpr <- expr_for(d, group = "bp_nolith")   # 74 x 12368
  perm_file  <- file.path(run_dir, "permStats_nolith.RData")
  mp_nolith  <- modulePreservation(
    multiData = list(
      Ref  = list(data = refExpr),
      Test = list(data = nolithExpr)
    ),
    multiColor = list(Ref = moduleColors),
    dataIsExpr             = TRUE,
    networkType            = NET_TYPE,
    corFnc                 = COR_FNC,
    nPermutations          = NOLITH_PERMS,
    randomSeed             = SEED_BASE + 1000L,
    savePermutedStatistics = TRUE,
    loadPermutedStatistics = FALSE,
    permutedStatisticsFile = perm_file,
    verbose                = 1,
    indent                 = 2
  )
  nolith_res          <- extract_mp(mp_nolith)
  nolith_res$group    <- "bp_nolith"
  nolith_res$draw     <- NA_integer_
  nolith_res$nperms   <- NOLITH_PERMS
  nolith_res$n_test   <- nrow(nolithExpr)
  saveRDS(nolith_res, nolith_ckpt)
  msg("bp_nolith done: ", nrow(nolith_res), " modules; Zsummary range [",
      round(min(nolith_res$Zsummary), 2), ", ",
      round(max(nolith_res$Zsummary), 2), "]")
  rm(nolithExpr, mp_nolith); gc()
} else {
  msg("bp_nolith checkpoint exists, skipping")
}

# ---------- 9. bp_lith draws (20 draws x 50 perms, checkpoint per draw) -----
# bp_lith has n=152; subsample to 74 each draw to match bp_nolith.
lith_sids <- colnames(d$logtpm)[!is.na(d$groups) & d$groups == "bp_lith"]
stopifnot(length(lith_sids) == 152L)

for (draw in seq_len(DRAWS)) {
  ckpt <- file.path(run_dir, sprintf("lith_draw_%03d.rds", draw))
  if (file.exists(ckpt)) {
    msg("draw ", draw, "/", DRAWS, " checkpoint exists, skipping")
    next
  }
  msg("draw ", draw, "/", DRAWS, ": subsampling bp_lith (n=", SUBSAMPLE_N, " from 152) ...")
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
    msg("ERROR in draw ", draw, ": ", conditionMessage(e))
    NULL
  })

  if (is.null(mp_draw)) {
    msg("draw ", draw, " failed; will retry in a later session if row is requeued")
    rm(testExpr); gc()
    next
  }

  draw_res          <- extract_mp(mp_draw)
  draw_res$group    <- "bp_lith"
  draw_res$draw     <- draw
  draw_res$nperms   <- PERMS_PER_DRAW
  draw_res$n_test   <- SUBSAMPLE_N
  draw_res$samples  <- paste(sel_sids, collapse = ",")   # exact IDs for reproducibility
  saveRDS(draw_res, ckpt)
  msg("draw ", draw, "/", DRAWS, " done: Zsummary range [",
      round(min(draw_res$Zsummary), 2), ", ",
      round(max(draw_res$Zsummary), 2), "]")
  rm(testExpr, mp_draw, draw_res); gc()
}

# ---------- 10. Aggregate (only after ALL draws complete) -------------------
nolith_done <- file.exists(nolith_ckpt)
done_lith   <- vapply(seq_len(DRAWS), function(i)
  file.exists(file.path(run_dir, sprintf("lith_draw_%03d.rds", i))), logical(1))
n_lith_done <- sum(done_lith)

if (!nolith_done || n_lith_done < DRAWS) {
  msg("not all work complete: nolith=", nolith_done,
      ", lith draws done=", n_lith_done, "/", DRAWS,
      "; exiting for resume in a later session")
  quit(save = "no", status = 1L)
}

msg("all draws complete; aggregating results ...")

# bp_nolith — straightforward single-run result
nolith_res <- readRDS(nolith_ckpt)
nolith_out <- nolith_res[, c("group", "module", "Zsummary", "medianRank",
                              "draw", "nperms", "n_test")]
write.csv(nolith_out, file.path(run_dir, "results_nolith.csv"), row.names = FALSE)

# bp_lith — per-draw data and per-module summary across draws
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

# Combined result: one row per module per group
nolith_tidy <- data.frame(
  group            = nolith_out$group,
  module           = nolith_out$module,
  Zsummary_point   = nolith_out$Zsummary,    # single estimate (500 perms)
  medianRank_point = nolith_out$medianRank,
  Zsummary_mean    = nolith_out$Zsummary,
  note             = paste0("single run, ", NOLITH_PERMS, " perms"),
  stringsAsFactors = FALSE
)
lith_tidy <- data.frame(
  group            = lith_summary$group,
  module           = lith_summary$module,
  Zsummary_point   = lith_summary$Zsummary_median,   # median across 20 draws
  medianRank_point = lith_summary$medianRank_median,
  Zsummary_mean    = lith_summary$Zsummary_mean,
  note             = paste0(DRAWS, " draws x ", PERMS_PER_DRAW, " perms; point=median"),
  stringsAsFactors = FALSE
)
combined <- rbind(nolith_tidy, lith_tidy)
write.csv(combined, file.path(run_dir, "results_combined.csv"), row.names = FALSE)

# Print tables for the log
msg("\n--- bp_nolith Zsummary (", NOLITH_PERMS, " perms, descriptive only) ---")
show_n <- nolith_tidy[order(-nolith_tidy$Zsummary_point),
                      c("module", "Zsummary_point", "medianRank_point")]
print(head(show_n, 10), row.names = FALSE)

msg("\n--- bp_lith Zsummary_mean across ", DRAWS, " draws (top 10, descriptive only) ---")
show_l <- lith_summary[order(-lith_summary$Zsummary_mean),
                       c("module", "Zsummary_mean", "Zsummary_sd", "medianRank_mean")]
print(head(show_l, 10), row.names = FALSE)

msg("=== ", id, " COMPLETE ===")
# Exit 0: driver marks the row DONE.
