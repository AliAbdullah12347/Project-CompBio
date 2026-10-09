#!/usr/bin/env Rscript
# ===========================================================================
# bmind-047  Per-lineage co-expression networks from bMIND profiles
#
# Spec: per-lineage networks from bMIND profiles (requires bmind-046 DONE)
#
# bMIND deconvolution produces a cell-type-specific expression matrix:
#   T x G x N (T cell types, G genes, N samples), stored typically as a list
#   of T matrices each G x N.
#
# For each lineage (gran, mono, T, NK, B):
#   1. Extract lineage-specific expression matrix (samples x genes)
#   2. Build WGCNA reference modules on control samples
#   3. Test preservation in bp_nolith (full set) and bp_lith (20 draws x 50 perms)
#   4. Checkpoint per draw per lineage
#
# Design notes:
#   - Power=12 (same as base-001; WGCNA FAQ, signed, n>40)
#   - The lineage-specific matrices may be smaller (fewer expressed genes per
#     cell type); min module size stays 30 but fewer modules expected per lineage
#   - SEED_BASE offset from all other experiments
#
# Dependency: bmind-046 DONE FLAG must exist.
#             base-001 reference cache (for comparison).
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "bmind-047"
stopifnot(id == "bmind-047")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

run_dir   <- file.path("network/runs",  id)
cfg_dir   <- "network/config"
cache_dir <- "network/cache"
for (d_ in c(run_dir, cfg_dir, cache_dir))
  dir.create(d_, showWarnings = FALSE, recursive = TRUE)

BMIND_DONE_FLAG <- "network/runs/bmind-046/bmind_profiles_done.flag"
POWER           <- 12L
POWER_MIN       <-  4L
NET_TYPE        <- "signed"
COR_FNC         <- "bicor"
MIN_MODULE      <- 30L
MERGE_THRESH    <- 0.25
DRAWS           <- 20L
PERMS_PER_DRAW  <- 50L
NOLITH_PERMS    <- 200L   # fewer perms for lineage networks (4 modules → fast)
SUBSAMPLE_N     <- 74L
SEED_BASE       <- 20261009L + 47000L
LINEAGES        <- c("gran", "mono", "T", "NK", "B")

if (!file.exists(BMIND_DONE_FLAG)) {
  msg("DEPENDENCY MISSING: bmind-046 not done (", BMIND_DONE_FLAG, " absent)")
  msg("Exiting status=1 to requeue after bmind-046 completes.")
  quit(save = "no", status = 1L)
}

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
    spec="per-lineage networks from bMIND profiles (requires bmind-046 DONE)",
    lineages=paste(LINEAGES, collapse=";"),
    power=POWER, net_type=NET_TYPE, cor_fnc=COR_FNC,
    draws_lith=DRAWS, perms_per_draw=PERMS_PER_DRAW,
    nolith_perms=NOLITH_PERMS, subsample_n=SUBSAMPLE_N,
    seed_base=SEED_BASE,
    git_sha=if(length(sha)&&!is.na(sha))sha else "unknown",
    r_version=paste(R.version$major,R.version$minor,sep="."),
    wgcna_version=as.character(packageVersion("WGCNA")))
  lines <- mapply(jsn, names(cfg), cfg, SIMPLIFY=TRUE)
  writeLines(c("{", paste(lines, collapse=",\n"), "}"), cfg_path)
  msg("config written: ", cfg_path)
}

suppressPackageStartupMessages({ library(WGCNA); library(fastcluster) })
enableWGCNAThreads()
options(stringsAsFactors = FALSE)
source("network/00_load.R")
d <- load_project()

# ---- locate bMIND profiles ---------------------------------------------------
possible_paths <- c(
  "de_analysis/data/bmind_profiles.rds",
  "de_analysis/data/bMIND_profiles.rds",
  "de_analysis/data/bmind_ct_expression.rds",
  "de_analysis/results/bmind_profiles.rds"
)
bmind_path <- Filter(file.exists, possible_paths)[1]
if (is.na(bmind_path) || is.null(bmind_path)) {
  msg("Cannot find bMIND profiles at standard paths:")
  for (p in possible_paths) msg("  tried: ", p)
  quit(save="no", status=1L)
}
msg("loading bMIND profiles from: ", bmind_path)
bmind <- readRDS(bmind_path)

# bMIND output format: typically a list of cell-type matrices, each samples x genes
# Or a 3D array (T x G x N). Detect and normalise.
normalise_bmind <- function(bm, lineages) {
  if (is.list(bm)) {
    # Named list: each element is a matrix (samples x genes or genes x samples)
    nms <- names(bm)
    msg("bMIND is a list with ", length(nms), " elements: ", paste(head(nms,5), collapse=", "))
    return(bm)
  }
  if (is.array(bm) && length(dim(bm)) == 3) {
    # 3D array: dim = (cell_types, genes, samples) or (genes, samples, cell_types)
    msg("bMIND is 3D array: dim = ", paste(dim(bm), collapse=" x "))
    dn <- dimnames(bm)
    # Try to figure out which dimension is cell_types
    for (ax in 1:3) {
      if (length(dn[[ax]]) > 0 && any(lineages %in% dn[[ax]])) {
        msg("cell-type dimension = axis ", ax, ": ", paste(dn[[ax]], collapse=", "))
        out <- setNames(lapply(lineages[lineages %in% dn[[ax]]], function(lin) {
          slc <- switch(ax,
            `1` = bm[lin,,], `2` = bm[,lin,], `3` = bm[,,lin])
          t(slc)   # samples x genes
        }), lineages[lineages %in% dn[[ax]]])
        return(out)
      }
    }
  }
  stop("Unrecognised bMIND output format — inspect object manually")
}

bm_list <- tryCatch(normalise_bmind(bmind, LINEAGES),
                    error=function(e){msg("ERROR: ",conditionMessage(e));NULL})
if (is.null(bm_list)) quit(save="no", status=1L)
available_lineages <- intersect(LINEAGES, names(bm_list))
msg("available lineages: ", paste(available_lineages, collapse=", "))

extract_mp <- function(mp) {
  obs  <- mp$preservation$observed[["ref.Ref"]][["inColumnsAlsoPresentIn.Test"]]
  zsco <- mp$preservation$Z[["ref.Ref"]][["inColumnsAlsoPresentIn.Test"]]
  if (is.null(obs)||is.null(zsco)) stop("modulePreservation structure unexpected")
  res  <- data.frame(module=rownames(obs), Zsummary=zsco[,"Zsummary.pres"],
                     medianRank=obs[,"medianRank.pres"], stringsAsFactors=FALSE)
  res[res$module != "grey",]
}

all_lith_rows <- list()

for (lin in available_lineages) {
  msg("=== lineage: ", lin, " ===")
  lin_dir <- file.path(run_dir, lin)
  dir.create(lin_dir, showWarnings = FALSE)

  # Lineage-specific expression matrix: all samples x genes
  expr_lin <- bm_list[[lin]]
  msg("  expression dim: ", paste(dim(expr_lin), collapse=" x "))

  # Align sample order to d$groups
  sid_order <- colnames(d$logtpm)
  rn <- rownames(expr_lin)
  if (!is.null(rn) && any(rn %in% sid_order)) {
    expr_lin <- expr_lin[match(sid_order, rn), , drop=FALSE]
    rownames(expr_lin) <- sid_order
  }
  ctrl_idx   <- which(!is.na(d$groups) & d$groups == "control")
  nolith_idx <- which(!is.na(d$groups) & d$groups == "bp_nolith")
  lith_idx   <- which(!is.na(d$groups) & d$groups == "bp_lith")
  refExpr    <- expr_lin[ctrl_idx,  , drop=FALSE]
  nolithExpr <- expr_lin[nolith_idx, , drop=FALSE]
  lithAll    <- expr_lin[lith_idx,  , drop=FALSE]
  lith_sids  <- rownames(lithAll)

  # Reference network (per lineage)
  ref_cache_lin <- file.path(cache_dir, paste0("ref_", lin, "_bicor_signed_p12.rds"))
  if (!file.exists(ref_cache_lin)) {
    msg("  building ref network for ", lin, " (n=", nrow(refExpr), ") ...")
    bw <- tryCatch(
      blockwiseModules(refExpr, power=POWER, networkType=NET_TYPE, corType=COR_FNC,
                       TOMType="signed", minModuleSize=MIN_MODULE,
                       mergeCutHeight=MERGE_THRESH, numericLabels=FALSE,
                       saveTOMs=FALSE, verbose=0, nThreads=3L),
      error=function(e){msg("ERROR building modules for ",lin,": ",conditionMessage(e));NULL})
    if (is.null(bw)) next
    moduleColors_lin <- bw$colors
    n_mod <- length(setdiff(unique(moduleColors_lin),"grey"))
    msg("  ", n_mod, " modules in ", lin, " reference")
    saveRDS(list(moduleColors=moduleColors_lin, n_modules=n_mod, lineage=lin),
            ref_cache_lin)
  }
  ref_lin      <- readRDS(ref_cache_lin)
  moduleColors_lin <- ref_lin$moduleColors
  msg("  ref modules: ", ref_lin$n_modules, " (excl. grey)")

  # bp_nolith preservation
  nolith_ckpt <- file.path(lin_dir, "nolith_preservation.rds")
  if (!file.exists(nolith_ckpt)) {
    msg("  bp_nolith (", NOLITH_PERMS, " perms) ...")
    perm_file <- file.path(lin_dir, "permStats_nolith.RData")
    mp_nolith <- tryCatch(
      modulePreservation(
        multiData=list(Ref=list(data=refExpr), Test=list(data=nolithExpr)),
        multiColor=list(Ref=moduleColors_lin), dataIsExpr=TRUE,
        networkType=NET_TYPE, corFnc=COR_FNC,
        nPermutations=NOLITH_PERMS, randomSeed=SEED_BASE + which(LINEAGES==lin)*1000L,
        savePermutedStatistics=TRUE, loadPermutedStatistics=FALSE,
        permutedStatisticsFile=perm_file, verbose=1, indent=2),
      error=function(e){msg("ERROR: ",conditionMessage(e));NULL})
    if (!is.null(mp_nolith)) {
      nr <- extract_mp(mp_nolith); nr$lineage <- lin; nr$group <- "bp_nolith"
      saveRDS(nr, nolith_ckpt); msg("  nolith done")
      rm(mp_nolith); gc()
    }
  } else { msg("  nolith exists") }

  # bp_lith draws
  for (draw in seq_len(DRAWS)) {
    ckpt <- file.path(lin_dir, sprintf("lith_draw_%03d.rds", draw))
    if (file.exists(ckpt)) { msg("  draw ", draw, " exists"); next }
    set.seed(SEED_BASE + which(LINEAGES==lin)*1000L + draw)
    sel_sids <- sample(lith_sids, SUBSAMPLE_N, replace=FALSE)
    testExpr <- lithAll[match(sel_sids, lith_sids), , drop=FALSE]
    perm_file <- file.path(lin_dir, sprintf("permStats_lith_%03d.RData", draw))
    msg("  ", lin, " lith draw ", draw, "/", DRAWS)
    mp_draw <- tryCatch(
      modulePreservation(
        multiData=list(Ref=list(data=refExpr), Test=list(data=testExpr)),
        multiColor=list(Ref=moduleColors_lin), dataIsExpr=TRUE,
        networkType=NET_TYPE, corFnc=COR_FNC,
        nPermutations=PERMS_PER_DRAW, randomSeed=SEED_BASE + which(LINEAGES==lin)*1000L + draw,
        savePermutedStatistics=TRUE, loadPermutedStatistics=FALSE,
        permutedStatisticsFile=perm_file, verbose=1, indent=2),
      error=function(e){msg("ERROR draw ",draw,": ",conditionMessage(e));NULL})
    if (is.null(mp_draw)) { rm(testExpr); gc(); next }
    dr <- extract_mp(mp_draw); dr$lineage <- lin; dr$group <- "bp_lith"; dr$draw <- draw
    dr$n_test <- SUBSAMPLE_N
    saveRDS(dr, ckpt)
    rm(testExpr, mp_draw); gc()
  }
  rm(refExpr, nolithExpr, lithAll); gc()

  # Collect lith results for this lineage
  lith_files <- vapply(seq_len(DRAWS), function(i)
    file.path(lin_dir, sprintf("lith_draw_%03d.rds",i)), character(1))
  lith_done <- vapply(lith_files, file.exists, logical(1))
  if (any(lith_done)) {
    lin_rows <- do.call(rbind, lapply(lith_files[lith_done], readRDS))
    all_lith_rows <- c(all_lith_rows, list(lin_rows))
  }
}

if (length(all_lith_rows) > 0) {
  combined <- do.call(rbind, all_lith_rows)
  write.csv(combined, file.path(run_dir,"results_lith_alldraws.csv"), row.names=FALSE)
  msg("wrote lith results: ", nrow(combined), " rows")
}

# Completion check: all lineages, all draws, nolith
all_done <- vapply(available_lineages, function(lin) {
  ld <- file.path(run_dir, lin)
  nolith_ok <- file.exists(file.path(ld, "nolith_preservation.rds"))
  lith_ok   <- all(vapply(seq_len(DRAWS), function(i)
    file.exists(file.path(ld, sprintf("lith_draw_%03d.rds",i))), logical(1)))
  nolith_ok && lith_ok
}, logical(1))
if (!all(all_done)) {
  msg("incomplete lineages: ", paste(names(all_done)[!all_done], collapse=","))
  quit(save="no", status=1L)
}
msg("=== ", id, " COMPLETE ===")
