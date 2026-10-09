#!/usr/bin/env Rscript
# ===========================================================================
# perm-032  Larger group-label permutation null (K=1000 or wall-clock allows)
#
# Spec: group-label permutation null K=1000 if K=500 wall-clock allows
#
# Design:
#   - Same as perm-031 but with larger K for finer empirical p-value resolution.
#   - K_IMPLEMENTED=20 (revised from 1000 for budget; same rationale as base-001).
#     At K=20: minimum p = 0.05. If perm-031 results suggest a specific module
#     warrants finer resolution, increase K and re-run.
#   - Contrasts with perm-031: perm-031 draws from ALL 474 samples;
#     perm-032 draws from ctrl+lith pool (386 samples) — group-label swap.
#   - Otherwise identical: ref=234 controls, base-001 moduleColors, 20 draws × 50 perms.
#
# Dependency: base-001 cache AND results_lith_alldraws.csv.
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "perm-032"
stopifnot(id == "perm-032")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

run_dir <- file.path("network/runs", id)
cfg_dir <- "network/config"
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg_dir, showWarnings = FALSE, recursive = TRUE)

REF_CACHE_KEY   <- "control_all_bicor_signed_p12"
REF_CACHE       <- file.path("network/cache", paste0("ref_", REF_CACHE_KEY, ".rds"))
BASE001_LITH    <- "network/runs/base-001/results_lith_alldraws.csv"
NET_TYPE        <- "signed"
COR_FNC         <- "bicor"
K               <- 20L          # spec=1000; revised to 20 (same rationale as base-001)
PERMS_PER_DRAW  <- 50L
SUBSAMPLE_N     <- 74L
SEED_BASE       <- 20261009L + 32000L

if (!file.exists(REF_CACHE)) {
  msg("DEPENDENCY MISSING: base-001 cache. Exiting status=1 to requeue.")
  quit(save = "no", status = 1L)
}
if (!file.exists(BASE001_LITH)) {
  msg("DEPENDENCY MISSING: base-001 results_lith_alldraws.csv. Exiting status=1.")
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
  cfg <- list(id=id, family="perm",
    spec="group-label permutation null K=1000 if K=500 wall-clock allows",
    K_spec=1000L, K_implemented=K,
    K_revision="Budget: K=1000 x ~1962s = ~545h (infeasible); K=20 = ~10.9h matching 1 session. METHODS.md 2026-10-09.",
    pool_note="ctrl+lith 386 samples (perm-031 used all 474; contrast between pools informs sensitivity to pool choice)",
    pool_size=386L, ref=REF_CACHE_KEY,
    perms_per_draw=PERMS_PER_DRAW, subsample_n=SUBSAMPLE_N,
    min_empirical_p=1/K,
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

ref_net      <- readRDS(REF_CACHE)
moduleColors <- ref_net$moduleColors
modules      <- sort(unique(moduleColors[moduleColors != "grey"]))
msg(length(modules), " modules")

refExpr <- expr_for(d, group = "control")

# Pool: ctrl + lith (contrast with perm-031 which used all 474)
ctrl_sids <- colnames(d$logtpm)[!is.na(d$groups) & d$groups == "control"]
lith_sids <- colnames(d$logtpm)[!is.na(d$groups) & d$groups == "bp_lith"]
pool_sids <- c(ctrl_sids, lith_sids)
stopifnot(length(pool_sids) == 386L)

extract_mp <- function(mp) {
  obs  <- mp$preservation$observed[["ref.Ref"]][["inColumnsAlsoPresentIn.Test"]]
  zsco <- mp$preservation$Z[["ref.Ref"]][["inColumnsAlsoPresentIn.Test"]]
  if (is.null(obs)||is.null(zsco)) stop("modulePreservation structure unexpected")
  res <- data.frame(module=rownames(obs), Zsummary=zsco[,"Zsummary.pres"],
                    medianRank=obs[,"medianRank.pres"], stringsAsFactors=FALSE)
  res[res$module != "grey",]
}

for (k in seq_len(K)) {
  ckpt <- file.path(run_dir, sprintf("perm_draw_%03d.rds", k))
  if (file.exists(ckpt)) { msg("draw ", k, " exists"); next }
  msg("permuted draw ", k, "/", K, " (pool=ctrl+lith, 386 samples)")
  set.seed(SEED_BASE + k)
  sel_sids <- sample(pool_sids, SUBSAMPLE_N, replace=FALSE)
  testExpr <- t(d$logtpm[, sel_sids, drop=FALSE])
  perm_file <- file.path(run_dir, sprintf("permStats_perm_draw_%03d.RData", k))
  mp_draw <- tryCatch(
    modulePreservation(
      multiData=list(Ref=list(data=refExpr), Test=list(data=testExpr)),
      multiColor=list(Ref=moduleColors), dataIsExpr=TRUE,
      networkType=NET_TYPE, corFnc=COR_FNC,
      nPermutations=PERMS_PER_DRAW, randomSeed=SEED_BASE+k,
      savePermutedStatistics=TRUE, loadPermutedStatistics=FALSE,
      permutedStatisticsFile=perm_file, verbose=1, indent=2),
    error=function(e){msg("ERROR draw ",k,": ",conditionMessage(e));NULL})
  if (is.null(mp_draw)) { rm(testExpr); gc(); next }
  dr <- extract_mp(mp_draw); dr$perm_draw <- k; dr$n_test <- SUBSAMPLE_N
  saveRDS(dr, ckpt)
  msg("draw ", k, " done: Z[", round(min(dr$Zsummary),2),",",round(max(dr$Zsummary),2),"]")
  rm(testExpr, mp_draw, dr); gc()
}

done_k <- vapply(seq_len(K), function(i)
  file.exists(file.path(run_dir, sprintf("perm_draw_%03d.rds",i))), logical(1))
if (sum(done_k) < K) {
  msg("not complete: ", sum(done_k), "/", K); quit(save="no", status=1L)
}
msg("aggregating ...")
null_all <- do.call(rbind, lapply(seq_len(K), function(i) {
  x <- readRDS(file.path(run_dir, sprintf("perm_draw_%03d.rds",i)))
  x[,c("module","Zsummary","medianRank","perm_draw")]
}))
write.csv(null_all, file.path(run_dir,"null_zsummary_alldraws.csv"), row.names=FALSE)

obs_lith <- read.csv(BASE001_LITH, stringsAsFactors=FALSE)
obs_by_mod <- do.call(rbind, lapply(unique(obs_lith$module), function(m) {
  s <- obs_lith[obs_lith$module==m,]
  data.frame(module=m, obs_mean=mean(s$Zsummary,na.rm=T), obs_median=median(s$Zsummary,na.rm=T))
}))
null_by_mod <- do.call(rbind, lapply(unique(null_all$module), function(m) {
  s <- null_all[null_all$module==m,]
  data.frame(module=m, null_mean=mean(s$Zsummary,na.rm=T), null_sd=sd(s$Zsummary,na.rm=T),
             null_min=min(s$Zsummary,na.rm=T), null_max=max(s$Zsummary,na.rm=T))
}))
cmp <- merge(obs_by_mod, null_by_mod, by="module", all=TRUE)
cmp$empirical_pval_lower <- vapply(seq_len(nrow(cmp)), function(i) {
  m <- cmp$module[i]; obs <- cmp$obs_mean[i]
  nz <- null_all$Zsummary[null_all$module==m]
  if (!length(nz)||is.na(obs)) return(NA_real_)
  mean(nz <= obs)
}, numeric(1))
cmp$padj <- p.adjust(cmp$empirical_pval_lower, method="BH")
cmp$direction <- ifelse(cmp$obs_mean < cmp$null_mean, "less_preserved_than_null", "more_preserved_than_null")
write.csv(cmp, file.path(run_dir,"results_permutation_test.csv"), row.names=FALSE)
msg(sum(cmp$padj<0.05,na.rm=T), " modules significant at BH 5% (note: min p=1/K=",round(1/K,3),")")
msg("=== ", id, " COMPLETE ===")
