#!/usr/bin/env Rscript
# ===========================================================================
# pow-009  pickSoftThreshold R2 curves for all three groups
#
# Spec: pickSoftThreshold powerVector=1:20 per group; record scale-free R2
#       curves (curves ARE the result)
#
# This is fast (~5 min per group). The R2 curves are themselves a scientific
# finding: scale-free topology fails in this data (all groups, all powers 1:20),
# confirming the use of the FAQ-recommended power 12 for signed networks, n>40.
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "pow-009"
stopifnot(id == "pow-009")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

run_dir <- file.path("network/runs", id)
cfg_dir <- "network/config"
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg_dir, showWarnings = FALSE, recursive = TRUE)

# ---- parameters (pre-specified) -------------------------------------------
POWER_VEC   <- 1:20
NET_TYPE    <- "signed"
COR_FNC     <- "bicor"
POWER_USED  <- 12L          # what we actually use (FAQ recommendation for n>40 signed)
SEED_BASE   <- 20261009L
GROUPS      <- c("control", "bp_nolith", "bp_lith")
SUBSAMPLE_N <- 74L           # bp_lith subsampled to match other groups

# ---- config -----------------------------------------------------------------
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
    id           = id, family = "power",
    spec         = "pickSoftThreshold powerVector=1:20 per group; record scale-free R2 curves (curves ARE the result)",
    power_vector = paste(range(POWER_VEC), collapse = ":"),
    networkType  = NET_TYPE, corFnc = COR_FNC,
    power_used   = POWER_USED,
    power_rationale = "WGCNA FAQ signed n>40: use 12 when scale-free fit fails",
    groups       = GROUPS,
    bp_lith_subsample_n = SUBSAMPLE_N,
    seed_base    = SEED_BASE,
    git_sha      = if (length(sha) && !is.na(sha)) sha else "unknown",
    r_version    = paste(R.version$major, R.version$minor, sep = "."),
    wgcna_version = as.character(packageVersion("WGCNA"))
  )
  lines <- mapply(jsn, names(cfg), cfg, SIMPLIFY = TRUE)
  writeLines(c("{", paste(lines, collapse = ",\n"), "}"), cfg_path)
  msg("config written: ", cfg_path)
}

# ---- packages ---------------------------------------------------------------
suppressPackageStartupMessages(library(WGCNA))
enableWGCNAThreads()
options(stringsAsFactors = FALSE)

# ---- data -------------------------------------------------------------------
source("network/00_load.R")
d <- load_project()

# ---- results path -----------------------------------------------------------
out_path <- file.path(run_dir, "sft_r2_all_groups.csv")

if (file.exists(out_path)) {
  msg("result CSV exists, nothing to do: ", out_path)
  msg("=== ", id, " COMPLETE (already done) ===")
  quit(save = "no", status = 0L)
}

# ---- scan per group -------------------------------------------------------
all_rows <- list()
for (grp in GROUPS) {
  msg("pickSoftThreshold for group: ", grp)
  if (grp == "bp_lith") {
    # subsample to n=74 to match other groups
    lith_sids <- colnames(d$logtpm)[!is.na(d$groups) & d$groups == "bp_lith"]
    stopifnot(length(lith_sids) == 152L)
    set.seed(SEED_BASE + 9000L)   # distinct from draw seeds
    sel <- sample(lith_sids, SUBSAMPLE_N, replace = FALSE)
    E <- t(d$logtpm[, sel])   # 74 x 12368
  } else {
    E <- expr_for(d, group = grp)
  }
  msg("  n_samples=", nrow(E), " n_genes=", ncol(E))
  set.seed(SEED_BASE + match(grp, GROUPS))
  sft <- pickSoftThreshold(E, powerVector = POWER_VEC,
                           networkType = NET_TYPE, corFnc = COR_FNC, verbose = 0)
  est <- sft$powerEstimate
  msg("  powerEstimate=", if (is.na(est)) "NA" else est,
      " (expected 1 or NA: function failing on this data)")
  r2df          <- sft$fitIndices
  r2df$group    <- grp
  r2df$n_samp   <- nrow(E)
  r2df$powerEstimate <- est
  r2df$power_used    <- POWER_USED
  max_r2 <- max(r2df$SFT.R.sq, na.rm = TRUE)
  msg("  max signed-R2 = ", round(max_r2, 3),
      " at power ", r2df$Power[which.max(r2df$SFT.R.sq)],
      " (target: >=0.80; FAILS in this data)")
  all_rows[[grp]] <- r2df
  rm(E, sft, r2df); gc()
}

combined <- do.call(rbind, all_rows)
write.csv(combined, out_path, row.names = FALSE)
msg("R2 curves written: ", out_path)
msg("=== ", id, " COMPLETE ===")
