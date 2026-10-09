#!/usr/bin/env Rscript
# ===========================================================================
# pres-030  Module preservation: density vs connectivity components separately
#
# Spec: ref=control;report density and connectivity components SEPARATELY;
#       disagreement is informative
#
# Context: WGCNA's Zsummary combines density-based and connectivity-based
# statistics into a single composite. Disagreement between the two components
# is scientifically informative: a module can retain intra-module density
# (genes still co-expressed) but lose connectivity pattern (hub structure
# disrupted), or vice versa. This experiment reports them separately.
#
# Design:
#   - Same setup as base-001 (n=234 ctrl ref, bicor, signed, power 12)
#   - bp_nolith: 500 perms, full statistics
#   - bp_lith: 20 draws x 50 perms, checkpoint per draw
#   - Output: density, connectivity, Zsummary_density, Zsummary_connectivity
#     extracted per module per draw
#
# Dependency: requires base-001 reference cache.
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "pres-030"
stopifnot(id == "pres-030")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

run_dir <- file.path("network/runs", id)
cfg_dir <- "network/config"
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg_dir, showWarnings = FALSE, recursive = TRUE)

REF_CACHE_KEY  <- "control_all_bicor_signed_p12"
REF_CACHE      <- file.path("network/cache", paste0("ref_", REF_CACHE_KEY, ".rds"))
NET_TYPE       <- "signed"
COR_FNC        <- "bicor"
DRAWS          <- 20L
PERMS_PER_DRAW <- 50L
NOLITH_PERMS   <- 500L
SUBSAMPLE_N    <- 74L
SEED_BASE      <- 20261009L + 30000L   # offset from base-001

if (!file.exists(REF_CACHE)) {
  msg("DEPENDENCY MISSING: ", REF_CACHE, " — base-001 must complete first")
  quit(save = "no", status = 1L)
}

# ---- config ----------------------------------------------------------------
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
  cfg <- list(id=id, family="preserve",
    spec="ref=control;report density and connectivity components SEPARATELY; disagreement is informative",
    ref_cache=REF_CACHE_KEY,
    corFnc=COR_FNC, networkType=NET_TYPE, power=12L,
    draws_lith=DRAWS, perms_per_draw=PERMS_PER_DRAW,
    nolith_perms=NOLITH_PERMS, subsample_n=SUBSAMPLE_N,
    seed_base=SEED_BASE,
    output_note="Density-only and connectivity-only Z-scores extracted separately from Zsummary",
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

ref_net      <- readRDS(REF_CACHE)
moduleColors <- ref_net$moduleColors
msg(ref_net$n_modules, " modules from base-001")

refExpr <- expr_for(d, group = "control")

# Pre-compute bp_lith expression (all 152 samples)
lithAll   <- expr_for(d, group = "bp_lith")
lith_sids <- rownames(lithAll)
stopifnot(length(lith_sids) == 152L)

# ---- extraction helper — density + connectivity + Zsummary separately -------
extract_mp_full <- function(mp) {
  obs  <- mp$preservation$observed[["ref.Ref"]][["inColumnsAlsoPresentIn.Test"]]
  zsco <- mp$preservation$Z[["ref.Ref"]][["inColumnsAlsoPresentIn.Test"]]
  if (is.null(obs)||is.null(zsco)) stop("modulePreservation structure unexpected")
  res <- data.frame(
    module = rownames(obs),
    Zsummary      = zsco[,"Zsummary.pres"],
    Zdensity      = if ("Zdensity.pres" %in% colnames(zsco)) zsco[,"Zdensity.pres"] else NA_real_,
    Zconnectivity = if ("Zconnectivity.pres" %in% colnames(zsco)) zsco[,"Zconnectivity.pres"] else NA_real_,
    medianRank    = obs[,"medianRank.pres"],
    medianRankDens = if ("medianRankDens.pres" %in% colnames(obs)) obs[,"medianRankDens.pres"] else NA_real_,
    medianRankConn = if ("medianRankConn.pres" %in% colnames(obs)) obs[,"medianRankConn.pres"] else NA_real_,
    stringsAsFactors = FALSE)
  res[res$module != "grey",]
}

# ---- bp_nolith -------------------------------------------------------------
nolith_ckpt <- file.path(run_dir, "nolith_preservation.rds")
if (!file.exists(nolith_ckpt)) {
  msg("bp_nolith (", NOLITH_PERMS, " perms, density+connectivity separate) ...")
  nolithExpr <- expr_for(d, group = "bp_nolith")
  perm_file  <- file.path(run_dir, "permStats_nolith.RData")
  mp_nolith  <- modulePreservation(
    multiData=list(Ref=list(data=refExpr), Test=list(data=nolithExpr)),
    multiColor=list(Ref=moduleColors), dataIsExpr=TRUE,
    networkType=NET_TYPE, corFnc=COR_FNC,
    nPermutations=NOLITH_PERMS, randomSeed=SEED_BASE+1000L,
    savePermutedStatistics=TRUE, loadPermutedStatistics=FALSE,
    permutedStatisticsFile=perm_file, verbose=1, indent=2)
  nr <- extract_mp_full(mp_nolith); nr$group <- "bp_nolith"; nr$draw <- NA_integer_
  nr$nperms <- NOLITH_PERMS; nr$n_test <- nrow(nolithExpr)
  saveRDS(nr, nolith_ckpt)
  msg("bp_nolith done")
  rm(nolithExpr, mp_nolith); gc()
} else { msg("bp_nolith exists") }

# ---- bp_lith draws ---------------------------------------------------------
for (draw in seq_len(DRAWS)) {
  ckpt <- file.path(run_dir, sprintf("lith_draw_%03d.rds", draw))
  if (file.exists(ckpt)) { msg("draw ", draw, " exists"); next }
  msg("draw ", draw, "/", DRAWS)
  set.seed(SEED_BASE + draw)
  sel_sids <- sample(lith_sids, SUBSAMPLE_N, replace=FALSE)
  testExpr <- lithAll[sel_sids, , drop=FALSE]
  perm_file <- file.path(run_dir, sprintf("permStats_lith_draw_%03d.RData", draw))
  mp_draw <- tryCatch(
    modulePreservation(
      multiData=list(Ref=list(data=refExpr), Test=list(data=testExpr)),
      multiColor=list(Ref=moduleColors), dataIsExpr=TRUE,
      networkType=NET_TYPE, corFnc=COR_FNC,
      nPermutations=PERMS_PER_DRAW, randomSeed=SEED_BASE+draw,
      savePermutedStatistics=TRUE, loadPermutedStatistics=FALSE,
      permutedStatisticsFile=perm_file, verbose=1, indent=2),
    error=function(e){msg("ERROR draw ",draw,": ",conditionMessage(e));NULL})
  if (is.null(mp_draw)) { rm(testExpr); gc(); next }
  dr <- extract_mp_full(mp_draw); dr$group <- "bp_lith"; dr$draw <- draw
  dr$nperms <- PERMS_PER_DRAW; dr$n_test <- SUBSAMPLE_N
  dr$samples <- paste(sel_sids, collapse=",")
  saveRDS(dr, ckpt)
  msg("draw ", draw, " done")
  rm(testExpr, mp_draw, dr); gc()
}

done_nolith <- file.exists(nolith_ckpt)
done_lith   <- vapply(seq_len(DRAWS), function(i)
  file.exists(file.path(run_dir, sprintf("lith_draw_%03d.rds",i))), logical(1))
if (!done_nolith || sum(done_lith) < DRAWS) {
  msg("not complete: nolith=", done_nolith, " lith=", sum(done_lith), "/", DRAWS)
  quit(save="no", status=1L)
}
msg("aggregating ...")
cols <- c("group","module","Zsummary","Zdensity","Zconnectivity","medianRank",
          "medianRankDens","medianRankConn","draw","nperms","n_test")
nr <- readRDS(nolith_ckpt)
write.csv(nr[, cols[cols %in% colnames(nr)]], file.path(run_dir,"results_nolith.csv"), row.names=FALSE)
lith_all <- do.call(rbind, lapply(seq_len(DRAWS), function(i) {
  x <- readRDS(file.path(run_dir, sprintf("lith_draw_%03d.rds",i)))
  x[, cols[cols %in% colnames(x)]]
}))
write.csv(lith_all, file.path(run_dir,"results_lith_alldraws.csv"), row.names=FALSE)
modules <- unique(lith_all$module)
lith_sum <- do.call(rbind, lapply(modules, function(m) {
  s <- lith_all[lith_all$module==m,]
  data.frame(group="bp_lith",module=m,n_draws=nrow(s),
             Zsummary_mean=mean(s$Zsummary,na.rm=T),
             Zdensity_mean=mean(s$Zdensity,na.rm=T),
             Zconn_mean=mean(s$Zconnectivity,na.rm=T),
             stringsAsFactors=FALSE)
}))
write.csv(lith_sum, file.path(run_dir,"results_lith_summary.csv"), row.names=FALSE)
msg("=== ", id, " COMPLETE ===")
