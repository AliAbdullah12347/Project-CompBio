#!/usr/bin/env Rscript
# ===========================================================================
# base-003-remerge  THE KEY RE-RUN. Preservation analysis at deepSplit=4, mch=0.
#
# Context: base-001 (v1) ran modulePreservation on a degenerate 4-module
# decomposition (mergeCutHeight=0.25 collapsed 33 modules to 3, largest 97.5%).
# base-001 (v2) was re-run with deepSplit=4, mergeCutHeight=0 (33 modules,
# largest 11.68%) during session 4 (2026-10-14). The v2 run IS this experiment.
#
# This script reads base-001 v2's completed results and writes them to the
# base-003-remerge run directory. It re-queues (exit 1) until base-001 v2
# is fully complete (all 20 lith draws + nolith).
#
# Why this design: base-001 v2 runs the actual modulePreservation computation.
# base-003-remerge provides a clean, named entry point for the "correct parameters"
# analysis while avoiding duplicating 2.2 h of computation.
#
# Queue spec: THE KEY RE-RUN. deepSplit=4, mergeCutHeight=0, draws=20, 500 nolith perms.
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "base-003-remerge"
stopifnot(id == "base-003-remerge")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))
msg("=== ", id, " starting ===")

# ---------- 1. Directories ---------------------------------------------------
run_dir <- file.path("network/runs", id)
cfg_dir <- "network/config"
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- 2. Parameters (pre-specified) ------------------------------------
DRAWS       <- 20L
REF_N       <- 234L
CACHE_KEY   <- "control_all_bicor_signed_p12_ds4_mch00"

# The source: base-001 v2 run at deepSplit=4, mergeCutHeight=0
SRC_DIR   <- "network/runs/base-001"
SRC_CACHE <- file.path("network/cache", paste0("ref_", CACHE_KEY, ".rds"))
SRC_V2    <- file.path(SRC_DIR, "v2_params.marker")

# ---------- 3. Write config BEFORE any computation --------------------------
cfg_path <- file.path(cfg_dir, paste0(id, ".json"))
if (!file.exists(cfg_path)) {
  sha <- tryCatch(
    trimws(system("git rev-parse --short HEAD 2>/dev/null", intern = TRUE)[1]),
    error = function(e) "unknown"
  )
  cfg <- paste0(
    '{\n',
    '  "id": "', id, '",\n',
    '  "family": "baseline",\n',
    '  "spec_note": "THE KEY RE-RUN at deepSplit=4, mergeCutHeight=0 (33 modules). ',
    'base-001 v2 (session 4) performs the actual modulePreservation; this script ',
    'reads its completed results.",\n',
    '  "source_experiment": "base-001_v2",\n',
    '  "source_dir": "', SRC_DIR, '",\n',
    '  "cache_key": "', CACHE_KEY, '",\n',
    '  "power": 12,\n',
    '  "deepSplit": 4,\n',
    '  "mergeCutHeight": 0.0,\n',
    '  "draws": ', DRAWS, ',\n',
    '  "perms_per_draw": 50,\n',
    '  "nolith_perms": 500,\n',
    '  "git_sha": "', if (!is.na(sha)) sha else "unknown", '",\n',
    '  "r_version": "', paste(R.version$major, R.version$minor, sep="."), '"\n',
    '}\n'
  )
  writeLines(cfg, cfg_path)
  msg("config written: ", cfg_path)
}

# ---------- 4. Check prerequisite: base-001 v2 must be complete --------------
msg("Checking base-001 v2 completion ...")

ref_ok   <- file.exists(SRC_CACHE)
v2_ok    <- file.exists(SRC_V2)
nolith_ok <- file.exists(file.path(SRC_DIR, "nolith_preservation.rds"))

done_lith <- vapply(seq_len(DRAWS), function(i)
  file.exists(file.path(SRC_DIR, sprintf("lith_draw_%03d.rds", i))), logical(1))
n_lith_done <- sum(done_lith)

msg("  ref cache: ", ref_ok)
msg("  v2 marker: ", v2_ok)
msg("  nolith:    ", nolith_ok)
msg("  lith draws done: ", n_lith_done, "/", DRAWS)

if (!ref_ok || !v2_ok || !nolith_ok || n_lith_done < DRAWS) {
  msg("base-001 v2 not yet complete; re-queuing (exit 1).")
  msg("Check back after base-001 finishes all ", DRAWS, " lith draws and nolith preservation.")
  quit(save = "no", status = 1L)
}

msg("base-001 v2 is complete. Reading results ...")

# ---------- 5. Read and re-report base-001 v2 results -----------------------
nolith_res <- readRDS(file.path(SRC_DIR, "nolith_preservation.rds"))
lith_all   <- do.call(rbind, lapply(seq_len(DRAWS), function(i) {
  x <- readRDS(file.path(SRC_DIR, sprintf("lith_draw_%03d.rds", i)))
  x[, intersect(c("group","module","Zsummary","medianRank","draw","nperms","n_test"),
                colnames(x))]
}))

# Nolith output table
nolith_out <- nolith_res[, c("group","module","Zsummary","medianRank","draw","nperms","n_test")]
write.csv(nolith_out, file.path(run_dir, "results_nolith.csv"), row.names = FALSE)

# Lith per-draw table
write.csv(lith_all, file.path(run_dir, "results_lith_alldraws.csv"), row.names = FALSE)

# Lith per-module summary
modules <- unique(lith_all$module)
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

# Gold-relative analysis: how many of 20 draws does each module beat gold?
gold_lith <- lith_all[lith_all$module == "gold", c("draw","Zsummary")]
colnames(gold_lith)[2] <- "gold_Z"
lith_with_gold <- merge(lith_all[lith_all$module != "gold", ],
                         gold_lith, by = "draw", all.x = TRUE)
beats_gold <- do.call(rbind, lapply(setdiff(modules, "gold"), function(m) {
  s <- lith_with_gold[lith_with_gold$module == m, ]
  data.frame(
    module       = m,
    n_draws_total = nrow(s),
    n_beats_gold  = sum(s$Zsummary > s$gold_Z, na.rm = TRUE),
    pct_beats     = round(100 * sum(s$Zsummary > s$gold_Z, na.rm = TRUE) / nrow(s), 1),
    stringsAsFactors = FALSE
  )
}))
# Sort by pct_beats descending
beats_gold <- beats_gold[order(-beats_gold$pct_beats), ]
write.csv(beats_gold, file.path(run_dir, "results_beats_gold.csv"), row.names = FALSE)

# Gold baseline stats for bp_nolith
nolith_gold_Z <- nolith_out$Zsummary[nolith_out$module == "gold"]
gold_nolith_note <- if (length(nolith_gold_Z) == 1)
  sprintf("gold Zsummary in bp_nolith = %.2f (random gene set)", nolith_gold_Z)
else "gold not found in nolith results"

# Combined (same format as base-001)
nolith_tidy <- data.frame(
  group            = nolith_out$group,
  module           = nolith_out$module,
  Zsummary_point   = nolith_out$Zsummary,
  medianRank_point = nolith_out$medianRank,
  Zsummary_mean    = nolith_out$Zsummary,
  note             = "single run, 500 perms",
  stringsAsFactors = FALSE
)
lith_tidy <- data.frame(
  group            = lith_summary$group,
  module           = lith_summary$module,
  Zsummary_point   = lith_summary$Zsummary_median,
  medianRank_point = lith_summary$medianRank_median,
  Zsummary_mean    = lith_summary$Zsummary_mean,
  note             = paste0("20 draws x 50 perms; point=median"),
  stringsAsFactors = FALSE
)
combined <- rbind(nolith_tidy, lith_tidy)
write.csv(combined, file.path(run_dir, "results_combined.csv"), row.names = FALSE)

# ---------- 6. Print summary for the log ------------------------------------
msg("\n--- bp_nolith Zsummary (500 perms, 33 modules, deepSplit=4 mch=0) ---")
msg(gold_nolith_note)
show_n <- nolith_out[order(-nolith_out$Zsummary), c("module","Zsummary","medianRank")]
print(head(show_n, 10), row.names = FALSE)

msg("\n--- bp_lith Zsummary across 20 draws (top 10) ---")
gold_lith_mean <- lith_summary$Zsummary_mean[lith_summary$module == "gold"]
msg(sprintf("gold mean across 20 draws = %.2f (sd=%.2f)",
    if (length(gold_lith_mean)) gold_lith_mean else NA,
    if (length(gold_lith_mean)) lith_summary$Zsummary_sd[lith_summary$module == "gold"] else NA))
show_l <- lith_summary[order(-lith_summary$Zsummary_mean),
                        c("module","Zsummary_mean","Zsummary_sd","medianRank_mean")]
print(head(show_l, 10), row.names = FALSE)

msg("\n--- Modules that beat gold: n_draws/20 ---")
print(head(beats_gold, 10), row.names = FALSE)

msg("=== ", id, " COMPLETE ===")
# Exit 0: driver marks the row DONE.
