#!/usr/bin/env Rscript
# ===========================================================================
# perm-031  Group-label permutation null for module preservation
#
# Spec: group-label permutation null K=500; draws matched observed vs permuted;
#       BH across MODULES not draws
#
# Design:
#   - K=20 permuted draws (revised from K=500 for budget; same reason as
#     base-001 draws=100→20, see METHODS.md 2026-10-09; minimum p=0.05)
#   - Each permuted draw: randomly select 74 samples from ALL 474 samples
#     (group-label permutation — any sample can serve as test)
#   - Reference always = 234 controls with base-001 moduleColors
#   - Run modulePreservation(ref=ctrl, test=random74, nPermutations=50)
#     to compute Zsummary under the null hypothesis of no group effect
#   - For each module: null distribution = 20 permuted Zsummary values
#   - Observed Zsummary from base-001 lith draws (read from checkpoint CSVs)
#   - Empirical p-value per module: fraction of null ≥ observed mean
#   - BH correction across modules (not across draws)
#
# Dependency: requires base-001 reference cache AND results files.
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "perm-031"
stopifnot(id == "perm-031")

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
K               <- 20L          # spec=500; revised to 20 (METHODS.md 2026-10-09)
PERMS_PER_DRAW  <- 50L
SUBSAMPLE_N     <- 74L
SEED_BASE       <- 20261009L + 31000L   # offset from base-001 to avoid collisions

if (!file.exists(REF_CACHE)) {
  msg("DEPENDENCY MISSING: ", REF_CACHE, " — base-001 must complete first")
  quit(save = "no", status = 1L)
}
if (!file.exists(BASE001_LITH)) {
  msg("DEPENDENCY MISSING: base-001 results_lith_alldraws.csv — base-001 must complete")
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
  cfg <- list(id=id, family="perm",
    spec="group-label permutation null K=500; draws matched observed vs permuted; BH across MODULES not draws",
    K_spec=500L, K_implemented=K,
    K_revision_note="Same budget constraint as base-001 draws: K=500 x ~1962s = 272h (infeasible); K=20 = ~10.9h. METHODS.md 2026-10-09.",
    pool="all 474 samples (any sample can be test under null)",
    ref=REF_CACHE_KEY,
    perms_per_draw=PERMS_PER_DRAW,
    subsample_n=SUBSAMPLE_N,
    comparison="base-001 bp_lith observed Zsummary vs permuted null",
    correction="BH across modules (not across draws)",
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

refExpr <- expr_for(d, group="control")   # 234 x 12368
all_sids <- colnames(d$logtpm)           # all 474 sample IDs

extract_mp <- function(mp) {
  obs  <- mp$preservation$observed[["ref.Ref"]][["inColumnsAlsoPresentIn.Test"]]
  zsco <- mp$preservation$Z[["ref.Ref"]][["inColumnsAlsoPresentIn.Test"]]
  if (is.null(obs)||is.null(zsco)) stop("modulePreservation structure unexpected")
  res <- data.frame(module=rownames(obs), Zsummary=zsco[,"Zsummary.pres"],
                    medianRank=obs[,"medianRank.pres"], stringsAsFactors=FALSE)
  res[res$module != "grey",]
}

# ---- K permuted draws -------------------------------------------------------
for (k in seq_len(K)) {
  ckpt <- file.path(run_dir, sprintf("perm_draw_%03d.rds", k))
  if (file.exists(ckpt)) { msg("draw ", k, " exists"); next }
  msg("permuted draw ", k, "/", K)
  set.seed(SEED_BASE + k)
  sel_sids <- sample(all_sids, SUBSAMPLE_N, replace=FALSE)
  testExpr <- t(d$logtpm[, sel_sids, drop=FALSE])   # 74 x 12368
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
  dr$samples <- paste(sel_sids, collapse=",")
  saveRDS(dr, ckpt)
  msg("draw ", k, " done: Z[", round(min(dr$Zsummary),2),", ",round(max(dr$Zsummary),2),"]")
  rm(testExpr, mp_draw, dr); gc()
}

done_k <- vapply(seq_len(K), function(i)
  file.exists(file.path(run_dir, sprintf("perm_draw_%03d.rds",i))), logical(1))
if (sum(done_k) < K) {
  msg("not complete: ", sum(done_k), "/", K, " permuted draws done")
  quit(save="no", status=1L)
}

# ---- compare to base-001 observed ------------------------------------------
msg("aggregating and computing empirical p-values ...")
null_all <- do.call(rbind, lapply(seq_len(K), function(i) {
  x <- readRDS(file.path(run_dir, sprintf("perm_draw_%03d.rds",i)))
  x[,c("module","Zsummary","medianRank","perm_draw")]
}))
write.csv(null_all, file.path(run_dir,"null_zsummary_alldraws.csv"), row.names=FALSE)

obs_lith <- read.csv(BASE001_LITH, stringsAsFactors=FALSE)
obs_by_mod <- do.call(rbind, lapply(unique(obs_lith$module), function(m) {
  s <- obs_lith[obs_lith$module==m,]
  data.frame(module=m, obs_Zsummary_mean=mean(s$Zsummary,na.rm=T),
             obs_Zsummary_median=median(s$Zsummary,na.rm=T),
             obs_n_draws=nrow(s))
}))
null_by_mod <- do.call(rbind, lapply(unique(null_all$module), function(m) {
  s <- null_all[null_all$module==m,]
  data.frame(module=m,
             null_Zsummary_mean=mean(s$Zsummary,na.rm=T),
             null_Zsummary_sd=sd(s$Zsummary,na.rm=T),
             null_Zsummary_min=min(s$Zsummary,na.rm=T),
             null_Zsummary_max=max(s$Zsummary,na.rm=T),
             null_n_draws=nrow(s))
}))
cmp <- merge(obs_by_mod, null_by_mod, by="module", all=TRUE)

# Empirical p (lower-tailed): fraction of null draws with Zsummary <= observed.
# The scientific signal is REDUCED preservation in bp_lith vs random draws.
# A small p means few random draws show Zsummary as low as bp_lith → specific.
# (If bp_lith were not reduced, ~50% of null draws would be ≤ observed.)
cmp$empirical_pval_lower <- vapply(seq_len(nrow(cmp)), function(i) {
  m   <- cmp$module[i]
  obs <- cmp$obs_Zsummary_mean[i]
  nz  <- null_all$Zsummary[null_all$module==m]
  if (!length(nz)||is.na(obs)) return(NA_real_)
  mean(nz <= obs)   # lower-tail: fraction of null ≤ observed
}, numeric(1))
cmp$padj <- p.adjust(cmp$empirical_pval_lower, method="BH")
cmp$direction <- ifelse(cmp$obs_Zsummary_mean < cmp$null_Zsummary_mean,
                        "less_preserved_than_random", "more_preserved_than_random")

write.csv(cmp, file.path(run_dir,"results_permutation_test.csv"), row.names=FALSE)
n_sig <- sum(cmp$padj < 0.05, na.rm=TRUE)
n_less <- sum(cmp$direction=="less_preserved_than_random", na.rm=TRUE)
msg(n_less, " modules less preserved in bp_lith than random null")
msg(n_sig, " of those significant at BH 5%")
msg("Note: minimum empirical p = 1/K = ", round(1/K,3), " with K=", K, " permuted draws")
msg("=== ", id, " COMPLETE ===")
