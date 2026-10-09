#!/usr/bin/env Rscript
# ===========================================================================
# gene-027  Circularity check: network without LM22 signature genes
#
# Spec: ref=control;genes=all_minus_547_LM22 ** CIRCULARITY CHECK **
#
# Purpose:
#   CIBERSORTx inferred cell fractions using the LM22 signature genes.
#   If module preservation results change substantially when those 547 genes
#   are excluded, the findings depended on the very genes that defined the
#   fractions — a circularity. Similar Zsummary without them means the result
#   is not an artifact of the fraction-estimation method.
#
# Design:
#   - Exclude LM22 genes (matched via d$sym Ensembl→HGNC, then intersected
#     with d$lm22 symbols).
#   - Build reference modules on non-LM22 genes, control group.
#   - Preservation in bp_nolith and bp_lith (draws=20).
#   - Otherwise identical to base-001 (power 12, bicor, signed).
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "gene-027"
stopifnot(id == "gene-027")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

run_dir   <- file.path("network/runs",   id)
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
DRAWS          <- 20L
PERMS_PER_DRAW <- 50L
NOLITH_PERMS   <- 500L
SUBSAMPLE_N    <- 74L
SEED_BASE      <- 20261009L
REF_GROUP      <- "control"

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
  cfg <- list(id=id, family="geneset",
    spec="ref=control;genes=all_minus_547_LM22 ** CIRCULARITY CHECK **",
    draws_spec=100L, draws_implemented=DRAWS,
    draws_revision_note="Same as base-001: revised to 20 (METHODS.md 2026-10-09)",
    lm22_exclusion="LM22 genes matched by HGNC symbol via d$sym; count determined at runtime",
    corFnc=COR_FNC, networkType=NET_TYPE, power=POWER,
    power_source="WGCNA_FAQ_signed_n_gt_40",
    minModuleSize=MIN_MODULE, mergeCutHeight=MERGE_THRESH,
    draws_lith=DRAWS, perms_per_lith_draw=PERMS_PER_DRAW,
    nolith_perms=NOLITH_PERMS, subsample_n=SUBSAMPLE_N,
    seed_base=SEED_BASE, git_sha=if(length(sha)&&!is.na(sha))sha else "unknown",
    r_version=paste(R.version$major,R.version$minor,sep="."),
    wgcna_version=as.character(packageVersion("WGCNA")))
  lines <- mapply(jsn, names(cfg), cfg, SIMPLIFY = TRUE)
  writeLines(c("{", paste(lines, collapse = ",\n"), "}"), cfg_path)
  msg("config written: ", cfg_path)
}

suppressPackageStartupMessages({ library(WGCNA); library(fastcluster) })
enableWGCNAThreads()
options(stringsAsFactors = FALSE)

source("network/00_load.R")
d <- load_project()

# ---- identify LM22 genes to exclude ----------------------------------------
lm22_sym <- d$lm22   # character vector of LM22 HGNC symbols
sym_map  <- d$sym    # named vector: Ensembl -> HGNC symbol
# Genes in our set whose HGNC symbol appears in LM22
lm22_ens <- names(sym_map)[sym_map %in% lm22_sym]
keep_ens <- setdiff(rownames(d$logtpm), lm22_ens)
n_lm22_excluded <- length(intersect(rownames(d$logtpm), lm22_ens))
n_kept           <- length(keep_ens)
msg("LM22 genes: ", length(lm22_sym), " symbols | matched to ",
    n_lm22_excluded, " Ensembl IDs in our set | kept: ", n_kept, " genes")
stopifnot(n_kept > 10000)   # sanity: should still have most genes

CACHE_KEY <- paste(REF_GROUP, paste0("noLM22_", n_kept), COR_FNC, NET_TYPE,
                   paste0("p", POWER), sep = "_")

# ---- helper to subset expression -------------------------------------------
expr_subset <- function(d, group, ens_keep, residualise = NULL) {
  E <- expr_for(d, group = group, residualise = residualise)   # samples x all_genes
  E[, ens_keep[ens_keep %in% colnames(E)], drop = FALSE]
}

# ---- scale-free scan -------------------------------------------------------
r2_path <- file.path(run_dir, "sft_r2_curve.csv")
if (!file.exists(r2_path)) {
  msg("pickSoftThreshold scan on non-LM22 gene set ...")
  tmp_ref <- expr_subset(d, REF_GROUP, keep_ens)
  set.seed(SEED_BASE)
  sft <- pickSoftThreshold(tmp_ref, powerVector=1:20,
                           networkType=NET_TYPE, corFnc=COR_FNC, verbose=0)
  est <- sft$powerEstimate
  msg("powerEstimate=", if(is.na(est))"NA" else est)
  r2df <- sft$fitIndices; r2df$power_used <- POWER; r2df$powerEstimate <- est
  r2df$n_genes <- n_kept
  write.csv(r2df, r2_path, row.names = FALSE)
  msg("R2 curve saved (max=", round(max(r2df$SFT.R.sq,na.rm=TRUE),3), ")")
  rm(tmp_ref, sft, r2df); gc()
}

# ---- reference network -----------------------------------------------------
ref_cache <- file.path(cache_dir, paste0("ref_", CACHE_KEY, ".rds"))
if (file.exists(ref_cache)) {
  msg("loading cached reference: ", ref_cache)
  ref_net <- readRDS(ref_cache); moduleColors <- ref_net$moduleColors
  msg(ref_net$n_modules, " modules (non-LM22 genes)")
} else {
  msg("building reference network (non-LM22, power ", POWER, ") ...")
  tmp_ref <- expr_subset(d, REF_GROUP, keep_ens)
  set.seed(SEED_BASE)
  bwm <- blockwiseModules(tmp_ref, power=POWER,
    networkType=NET_TYPE, TOMType=NET_TYPE, corType=COR_FNC,
    minModuleSize=MIN_MODULE, mergeCutHeight=MERGE_THRESH,
    maxBlockSize=15000L, numericLabels=FALSE, saveTOMs=FALSE, verbose=3)
  moduleColors <- bwm$colors
  n_mods <- length(unique(moduleColors)) - 1L
  msg(n_mods, " modules detected")
  ref_net <- list(moduleColors=moduleColors, n_modules=n_mods,
                  module_sizes=table(moduleColors), power=POWER, cache_key=CACHE_KEY,
                  n_genes_kept=n_kept, n_lm22_excluded=n_lm22_excluded)
  saveRDS(ref_net, ref_cache)
  write.csv(data.frame(n_genes_total=nrow(d$logtpm), n_lm22_excluded=n_lm22_excluded,
                       n_genes_kept=n_kept),
            file.path(run_dir, "gene_exclusion_summary.csv"), row.names=FALSE)
  write.csv(as.data.frame(ref_net$module_sizes),
            file.path(run_dir, "reference_module_sizes.csv"), row.names=FALSE)
  rm(bwm, tmp_ref); gc()
}

refExpr <- expr_subset(d, REF_GROUP, keep_ens)

# Pre-compute all bp_lith residuals at once
lithAllFull <- expr_for(d, group = "bp_lith")   # 152 x all_genes
lithAll     <- lithAllFull[, keep_ens[keep_ens %in% colnames(lithAllFull)], drop=FALSE]
lith_sids   <- rownames(lithAll)
stopifnot(length(lith_sids) == 152L)
rm(lithAllFull); gc()

extract_mp <- function(mp) {
  obs  <- mp$preservation$observed[["ref.Ref"]][["inColumnsAlsoPresentIn.Test"]]
  zsco <- mp$preservation$Z[["ref.Ref"]][["inColumnsAlsoPresentIn.Test"]]
  if (is.null(obs)||is.null(zsco)) stop("unexpected modulePreservation structure")
  res <- data.frame(module=rownames(obs), Zsummary=zsco[,"Zsummary.pres"],
                    medianRank=obs[,"medianRank.pres"], stringsAsFactors=FALSE)
  res[res$module != "grey",]
}

# ---- bp_nolith -------------------------------------------------------------
nolith_ckpt <- file.path(run_dir, "nolith_preservation.rds")
if (!file.exists(nolith_ckpt)) {
  msg("bp_nolith (non-LM22, ", NOLITH_PERMS, " perms) ...")
  nolithExpr <- expr_subset(d, "bp_nolith", keep_ens)
  perm_file  <- file.path(run_dir, "permStats_nolith.RData")
  mp_nolith  <- modulePreservation(
    multiData=list(Ref=list(data=refExpr), Test=list(data=nolithExpr)),
    multiColor=list(Ref=moduleColors), dataIsExpr=TRUE,
    networkType=NET_TYPE, corFnc=COR_FNC,
    nPermutations=NOLITH_PERMS, randomSeed=SEED_BASE+1000L,
    savePermutedStatistics=TRUE, loadPermutedStatistics=FALSE,
    permutedStatisticsFile=perm_file, verbose=1, indent=2)
  nr <- extract_mp(mp_nolith); nr$group <- "bp_nolith"; nr$draw <- NA_integer_
  nr$nperms <- NOLITH_PERMS; nr$n_test <- nrow(nolithExpr)
  saveRDS(nr, nolith_ckpt)
  msg("bp_nolith done: Z[", round(min(nr$Zsummary),2), ", ", round(max(nr$Zsummary),2), "]")
  rm(nolithExpr, mp_nolith); gc()
} else { msg("bp_nolith exists, skipping") }

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
  dr <- extract_mp(mp_draw); dr$group <- "bp_lith"; dr$draw <- draw
  dr$nperms <- PERMS_PER_DRAW; dr$n_test <- SUBSAMPLE_N
  dr$samples <- paste(sel_sids, collapse=",")
  saveRDS(dr, ckpt)
  msg("draw ", draw, " done: Z[", round(min(dr$Zsummary),2), ", ", round(max(dr$Zsummary),2), "]")
  rm(testExpr, mp_draw, dr); gc()
}

# ---- aggregate -------------------------------------------------------------
done_nolith <- file.exists(nolith_ckpt)
done_lith <- vapply(seq_len(DRAWS), function(i)
  file.exists(file.path(run_dir, sprintf("lith_draw_%03d.rds",i))), logical(1))
if (!done_nolith || sum(done_lith) < DRAWS) {
  msg("not complete: nolith=", done_nolith, " lith=", sum(done_lith), "/", DRAWS)
  quit(save="no", status=1L)
}
msg("aggregating ...")
nolith_res <- readRDS(nolith_ckpt)
write.csv(nolith_res[,c("group","module","Zsummary","medianRank","draw","nperms","n_test")],
          file.path(run_dir,"results_nolith.csv"), row.names=FALSE)
lith_all <- do.call(rbind, lapply(seq_len(DRAWS), function(i) {
  x <- readRDS(file.path(run_dir, sprintf("lith_draw_%03d.rds",i)))
  x[,c("group","module","Zsummary","medianRank","draw","nperms","n_test")]
}))
write.csv(lith_all, file.path(run_dir,"results_lith_alldraws.csv"), row.names=FALSE)
modules <- unique(lith_all$module)
lith_summary <- do.call(rbind, lapply(modules, function(m) {
  s <- lith_all[lith_all$module==m,]
  data.frame(group="bp_lith",module=m,n_draws=nrow(s),
             Zsummary_mean=mean(s$Zsummary,na.rm=T),
             Zsummary_median=median(s$Zsummary,na.rm=T),
             Zsummary_sd=sd(s$Zsummary,na.rm=T),
             medianRank_mean=mean(s$medianRank,na.rm=T),
             stringsAsFactors=FALSE)
}))
write.csv(lith_summary, file.path(run_dir,"results_lith_summary.csv"), row.names=FALSE)
msg("=== ", id, " COMPLETE ===")
