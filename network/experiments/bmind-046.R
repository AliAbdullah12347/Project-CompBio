#!/usr/bin/env Rscript
# ===========================================================================
# bmind-046  REBUILD bMIND cell-type-specific expression profiles (~3.1h)
#
# Spec: REBUILD bMIND profiles: Rscript de_analysis/scripts/06_bmind_profiles.R
#       (~3.1h). CHECKPOINT. Own window.
#
# This script wraps the bMIND profile rebuild as a system() call. The wrapper:
#   1. Checks for a completion flag (bmind_profiles_done.flag) — skip if present.
#   2. Runs: Rscript de_analysis/scripts/06_bmind_profiles.R
#   3. Checks that output exists and writes the completion flag.
#
# IMPORTANT: run in its own terminal window (Own window) as it takes ~3.1h.
# bmind-047 (per-lineage networks) depends on this completing first.
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "bmind-046"
stopifnot(id == "bmind-046")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

run_dir <- file.path("network/runs", id)
cfg_dir <- "network/config"
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg_dir, showWarnings = FALSE, recursive = TRUE)

BMIND_SCRIPT <- "de_analysis/scripts/06_bmind_profiles.R"
DONE_FLAG    <- file.path(run_dir, "bmind_profiles_done.flag")

if (!file.exists(BMIND_SCRIPT)) stop("bMIND script not found: ", BMIND_SCRIPT)

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
  cfg <- list(id=id, family="bmind",
    spec="REBUILD bMIND profiles: Rscript de_analysis/scripts/06_bmind_profiles.R (~3.1h). CHECKPOINT. Own window.",
    bmind_script=BMIND_SCRIPT,
    done_flag=DONE_FLAG,
    expected_runtime_hours=3.1,
    note="wrapper script; runs de_analysis/scripts/06_bmind_profiles.R via system()",
    bmind047_dep="bmind-047 depends on this script completing (checks done_flag)",
    git_sha=if(length(sha)&&!is.na(sha))sha else "unknown",
    r_version=paste(R.version$major,R.version$minor,sep="."))
  lines <- mapply(jsn, names(cfg), cfg, SIMPLIFY=TRUE)
  writeLines(c("{", paste(lines, collapse=",\n"), "}"), cfg_path)
  msg("config written: ", cfg_path)
}

if (file.exists(DONE_FLAG)) {
  msg("CHECKPOINT: bMIND profiles already built (", DONE_FLAG, " exists)")
  msg("=== ", id, " COMPLETE (from checkpoint) ===")
  quit(save="no", status=0L)
}

msg("Starting bMIND profile rebuild (~3.1h) ...")
msg("Script: ", BMIND_SCRIPT)
msg("Working directory: ", getwd())

log_file <- file.path(run_dir, "bmind_profiles.log")
cmd      <- sprintf("Rscript %s > %s 2>&1", BMIND_SCRIPT, log_file)
msg("Running: ", cmd)
ret <- system(cmd, wait=TRUE)
msg("Return code: ", ret)

# Show last lines of log
if (file.exists(log_file)) {
  log_lines <- readLines(log_file)
  n <- min(20L, length(log_lines))
  msg("Last ", n, " lines of bMIND log:")
  for (ln in tail(log_lines, n)) message("  ", ln)
}

if (ret != 0) {
  msg("ERROR: bMIND script returned non-zero exit code ", ret)
  msg("Check log: ", log_file)
  quit(save="no", status=1L)
}

# ---- verify output exists ---------------------------------------------------
# 06_bmind_profiles.R writes output to de_analysis/data/ (typical bMIND path)
# Look for any new RDS/CSV that could be bMIND profiles
possible_outputs <- c(
  "de_analysis/data/bmind_profiles.rds",
  "de_analysis/data/bMIND_profiles.rds",
  "de_analysis/data/bmind_ct_expression.rds",
  "de_analysis/results/bmind_profiles.rds"
)
found <- Filter(file.exists, possible_outputs)
if (length(found) == 0) {
  msg("WARNING: no expected bMIND output found at standard paths")
  msg("Searched: ", paste(possible_outputs, collapse="; "))
  msg("Please verify output manually, then touch ", DONE_FLAG)
  quit(save="no", status=1L)
}
msg("bMIND output found: ", paste(found, collapse="; "))

writeLines(c(sprintf("completed: %s", now_()), paste("outputs:", paste(found, collapse="; "))),
           DONE_FLAG)
msg("flag written: ", DONE_FLAG)
msg("=== ", id, " COMPLETE ===")
