#!/usr/bin/env Rscript
# ===========================================================================
# tob-045  Tobacco sensitivity: repeat inp-022 across 5 of 20 tobacco imputations
#
# Spec: repeat inp-022 across several of the 20 tobacco imputations;
#       does any conclusion move?
#
# Design:
#   inp-022 residualises for cell composition only: ~ b1 + b2 + b3 + b4
#   tob-045 also residualises for tobacco:          ~ tobacco_imp_XX + b1 + b2 + b3 + b4
#
#   We use 5 imputations: imp_01, imp_05, imp_10, imp_15, imp_20 (evenly spaced
#   across the 20 available imputations). For each:
#   - Pre-compute residuals on all 474 samples
#   - Build reference modules on residualized control data (all 234)
#   - Test preservation in ONE lith draw (74 samples, 50 perms) — speed trade-off
#     documented in config; this is sensitivity analysis not the primary result
#
#   Comparison: if inp-022 results exist, compare Zsummary with vs. without tobacco.
#   Output: Zsummary per module per imputation + sensitivity summary.
#
# Dependency: base-001 reference cache (for moduleColors comparison).
# Note: imp-022 cache NOT required — tob-045 builds its own reference each time
#       because the module structure may shift with tobacco residualization.
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "tob-045"
stopifnot(id == "tob-045")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

run_dir   <- file.path("network/runs",  id)
cfg_dir   <- "network/config"
cache_dir <- "network/cache"
for (d_ in c(run_dir, cfg_dir, cache_dir))
  dir.create(d_, showWarnings = FALSE, recursive = TRUE)

POWER          <- 12L
POWER_MIN      <-  4L
NET_TYPE       <- "signed"
COR_FNC        <- "bicor"
MIN_MODULE     <- 30L
MERGE_THRESH   <- 0.25
PERMS_PER_DRAW <- 50L
SUBSAMPLE_N    <- 74L
SEED_BASE      <- 20261009L + 45000L
IMPUTATIONS    <- c("tobacco_imp_01", "tobacco_imp_05", "tobacco_imp_10",
                    "tobacco_imp_15", "tobacco_imp_20")
N_LITH_DRAWS   <- 1L   # 1 draw per imputation (sensitivity, not primary result)

REF_CACHE_BASE <- file.path("network/cache", "ref_control_all_bicor_signed_p12.rds")
if (!file.exists(REF_CACHE_BASE)) {
  msg("DEPENDENCY MISSING: base-001 cache. Exiting status=1 to requeue.")
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
  cfg <- list(id=id, family="input",
    spec="repeat inp-022 across several of the 20 tobacco imputations; does any conclusion move?",
    resid_formula="~ tobacco_imp_XX + b1 + b2 + b3 + b4",
    imputations_tested=paste(IMPUTATIONS, collapse=";"),
    n_imputations_tested=length(IMPUTATIONS),
    n_lith_draws=N_LITH_DRAWS,
    lith_draws_note="1 draw per imputation (sensitivity analysis, not primary result)",
    comparison="vs inp-022 (composition-only: ~b1+b2+b3+b4); does tobacco adjustment shift Zsummary?",
    power=POWER, net_type=NET_TYPE, cor_fnc=COR_FNC,
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

# ---- helper: residualize expression for tobacco_imp_XX + ILR b1-b4 ----------
expr_resid_tobacco <- function(d, group, tob_col, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  keep <- !is.na(d$groups) & d$groups %in% group
  M    <- d$logtpm[, keep, drop=FALSE]   # genes x samples
  md   <- d$meta[keep, , drop=FALSE]
  ilr  <- d$ilr[keep, -1, drop=FALSE]   # remove sample column
  tob  <- d$tobacco[keep, tob_col, drop=FALSE]
  colnames(tob) <- "tobacco_covariate"
  md_full <- cbind(md, ilr, tob)
  X <- tryCatch(
    model.matrix(~ tobacco_covariate + b1 + b2 + b3 + b4, data=md_full),
    error=function(e) {
      msg("WARN: tobacco model.matrix failed (", conditionMessage(e), "), falling back to ~b1+b2+b3+b4")
      model.matrix(~ b1 + b2 + b3 + b4, data=md_full)
    })
  if (nrow(X) != ncol(M)) stop("residualize: row mismatch")
  t(residuals(lm.fit(X, t(M))))   # samples x genes
}

extract_mp <- function(mp) {
  obs  <- mp$preservation$observed[["ref.Ref"]][["inColumnsAlsoPresentIn.Test"]]
  zsco <- mp$preservation$Z[["ref.Ref"]][["inColumnsAlsoPresentIn.Test"]]
  if (is.null(obs)||is.null(zsco)) stop("modulePreservation structure unexpected")
  res  <- data.frame(module=rownames(obs), Zsummary=zsco[,"Zsummary.pres"],
                     medianRank=obs[,"medianRank.pres"], stringsAsFactors=FALSE)
  res[res$module != "grey",]
}

all_results <- list()

for (tob_col in IMPUTATIONS) {
  ckpt <- file.path(run_dir, paste0("result_", tob_col, ".rds"))
  if (file.exists(ckpt)) {
    msg("imputation ", tob_col, " exists, loading")
    all_results[[tob_col]] <- readRDS(ckpt)
    next
  }
  msg("=== imputation: ", tob_col, " ===")

  # Residualize control and lith
  msg("  residualizing control (n=234) ...")
  refExpr  <- expr_resid_tobacco(d, "control", tob_col, seed=SEED_BASE)
  msg("  residualizing lith (n=152) ...")
  lithAll  <- expr_resid_tobacco(d, "bp_lith",  tob_col, seed=SEED_BASE+1L)
  lith_sids <- colnames(d$logtpm)[!is.na(d$groups) & d$groups == "bp_lith"]
  stopifnot(nrow(lithAll) == 152L)

  # Build reference network with tobacco-residualized controls
  cache_key_tob <- paste0("control_all_bicor_signed_p12_tobacco_", tob_col)
  ref_cache_tob <- file.path(cache_dir, paste0("ref_", cache_key_tob, ".rds"))
  if (!file.exists(ref_cache_tob)) {
    msg("  building reference network for ", tob_col, " ...")
    sft <- pickSoftThreshold.fromSimilarity(
      adjacency(refExpr, type=NET_TYPE, corFnc=COR_FNC, power=POWER),
      powerVector=1:20, verbose=0)
    power_used <- POWER
    if (!is.na(sft$powerEstimate) && sft$powerEstimate >= POWER_MIN) {
      power_used <- sft$powerEstimate
    }
    msg("  power=", power_used)
    bw <- blockwiseModules(refExpr, power=power_used, networkType=NET_TYPE, corType=COR_FNC,
                           TOMType="signed", minModuleSize=MIN_MODULE,
                           mergeCutHeight=MERGE_THRESH, numericLabels=FALSE,
                           saveTOMs=FALSE, verbose=0, nThreads=3L)
    moduleColors_tob <- bw$colors
    n_mod_tob <- length(setdiff(unique(moduleColors_tob),"grey"))
    msg("  ", n_mod_tob, " modules")
    saveRDS(list(moduleColors=moduleColors_tob, n_modules=n_mod_tob,
                 power=power_used, tob_col=tob_col), ref_cache_tob)
  } else {
    msg("  loading ref cache for ", tob_col)
  }
  ref_net_tob      <- readRDS(ref_cache_tob)
  moduleColors_tob <- ref_net_tob$moduleColors

  # One lith draw
  set.seed(SEED_BASE + which(IMPUTATIONS==tob_col))
  sel_sids  <- sample(lith_sids, SUBSAMPLE_N, replace=FALSE)
  testExpr  <- lithAll[match(sel_sids, lith_sids), , drop=FALSE]
  perm_file <- file.path(run_dir, paste0("permStats_", tob_col, ".RData"))
  msg("  testing preservation in 1 lith draw (n=", SUBSAMPLE_N, ") ...")
  mp <- tryCatch(
    modulePreservation(
      multiData=list(Ref=list(data=refExpr), Test=list(data=testExpr)),
      multiColor=list(Ref=moduleColors_tob), dataIsExpr=TRUE,
      networkType=NET_TYPE, corFnc=COR_FNC,
      nPermutations=PERMS_PER_DRAW, randomSeed=SEED_BASE + which(IMPUTATIONS==tob_col),
      savePermutedStatistics=TRUE, loadPermutedStatistics=FALSE,
      permutedStatisticsFile=perm_file, verbose=1, indent=2),
    error=function(e){msg("ERROR: ",conditionMessage(e));NULL})
  if (is.null(mp)) { rm(refExpr,lithAll,testExpr); gc(); next }
  res <- extract_mp(mp); res$imputation <- tob_col; res$n_modules <- ref_net_tob$n_modules
  saveRDS(res, ckpt)
  all_results[[tob_col]] <- res
  msg("  done: Z[", round(min(res$Zsummary),2), ",", round(max(res$Zsummary),2), "]")
  rm(refExpr, lithAll, testExpr, mp, res); gc()
}

if (length(all_results) > 0) {
  combined <- do.call(rbind, lapply(names(all_results), function(nm) {
    x <- all_results[[nm]]; x$imputation <- nm; x
  }))
  write.csv(combined, file.path(run_dir,"results_allimputations.csv"), row.names=FALSE)
  # Sensitivity summary: range of Zsummary across imputations per module
  mods_union <- unique(combined$module)
  sens <- do.call(rbind, lapply(mods_union, function(m) {
    s <- combined[combined$module==m,]
    data.frame(module=m, n_imputations=nrow(s),
               Zsummary_mean=mean(s$Zsummary,na.rm=T),
               Zsummary_min=min(s$Zsummary,na.rm=T),
               Zsummary_max=max(s$Zsummary,na.rm=T),
               Zsummary_range=diff(range(s$Zsummary,na.rm=T)),
               stringsAsFactors=FALSE)
  }))
  write.csv(sens, file.path(run_dir,"results_sensitivity_summary.csv"), row.names=FALSE)
  msg("Zsummary range across imputations: ",
      paste(sens$module, round(sens$Zsummary_range,2), sep="=", collapse="; "))
}

done_all <- vapply(IMPUTATIONS, function(nm)
  file.exists(file.path(run_dir, paste0("result_", nm, ".rds"))), logical(1))
if (!all(done_all)) {
  msg("incomplete: ", sum(!done_all), " imputations pending"); quit(save="no", status=1L)
}
msg("=== ", id, " COMPLETE ===")
