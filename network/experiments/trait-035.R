#!/usr/bin/env Rscript
# ===========================================================================
# trait-035  Module eigengene ~ ILR balances b1-b4
#
# Spec: eigengene ~ ILR balances b1-b4
#
# ILR balances b1-b4 are uncorrelated log-ratio coordinates of the 5 lineage
# fractions. Unlike raw fractions, ILR coords have no sum-to-1 constraint and
# can all be included without singularity.
#
# Dependency: requires base-001 reference cache.
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "trait-035"
stopifnot(id == "trait-035")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

run_dir <- file.path("network/runs", id)
cfg_dir <- "network/config"
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg_dir, showWarnings = FALSE, recursive = TRUE)

REF_CACHE_KEY <- "control_all_bicor_signed_p12"
REF_CACHE     <- file.path("network/cache", paste0("ref_", REF_CACHE_KEY, ".rds"))
SEED          <- 20261009L
ILR_COLS      <- c("b1","b2","b3","b4")

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
  cfg <- list(id=id, family="modtrait",
    spec="eigengene ~ ILR balances b1-b4",
    ref_cache=REF_CACHE_KEY,
    formula="eigengene ~ b1 + b2 + b3 + b4",
    ilr_note="b1-b4 are ILR coordinates of 5-lineage fractions; uncorrelated, no sum-to-1 constraint",
    samples="all 474",
    correction="BH across modules, per ILR coordinate",
    seed=SEED,
    git_sha=if(length(sha)&&!is.na(sha))sha else "unknown",
    r_version=paste(R.version$major,R.version$minor,sep="."),
    wgcna_version=as.character(packageVersion("WGCNA")))
  lines <- mapply(jsn, names(cfg), cfg, SIMPLIFY=TRUE)
  writeLines(c("{", paste(lines, collapse=",\n"), "}"), cfg_path)
  msg("config written: ", cfg_path)
}

suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors = FALSE)

source("network/00_load.R")
d <- load_project()

ref_net      <- readRDS(REF_CACHE)
moduleColors <- ref_net$moduleColors
modules      <- sort(unique(moduleColors))
modules      <- modules[modules != "grey"]
msg(length(modules), " modules")

E  <- t(d$logtpm)
MEs <- moduleEigengenes(E, colors = moduleColors, excludeGrey = TRUE)$eigengenes
rownames(MEs) <- colnames(d$logtpm)
msg("eigengenes: ", nrow(MEs), " x ", ncol(MEs))

# ILR: 474 x 4+1 (sample col + b1..b4)
ilr <- d$ilr
stopifnot(all(ilr$sample == colnames(d$logtpm)))
ilr_mat <- ilr[, ILR_COLS[ILR_COLS %in% colnames(ilr)], drop=FALSE]
df <- cbind(MEs, ilr_mat)

results <- list()
for (mod in modules) {
  me_col <- paste0("ME", mod)
  if (!me_col %in% colnames(df)) me_col <- paste0("me", mod)
  if (!me_col %in% colnames(df)) next
  form <- as.formula(paste0("`", me_col, "` ~ ", paste(ILR_COLS, collapse=" + ")))
  fit <- lm(form, data=df)
  cs  <- summary(fit)$coefficients
  for (pred in ILR_COLS) {
    if (!pred %in% rownames(cs)) next
    results[[length(results)+1]] <- data.frame(
      module=mod, predictor=pred, coef=cs[pred,"Estimate"],
      tstat=cs[pred,"t value"], pval=cs[pred,"Pr(>|t|)"],
      stringsAsFactors=FALSE)
  }
}
res <- do.call(rbind, results)
res$padj <- NA_real_
for (pred in ILR_COLS) {
  idx <- res$predictor == pred
  if (sum(idx) > 0) res$padj[idx] <- p.adjust(res$pval[idx], method="BH")
}
write.csv(res, file.path(run_dir,"results_modtrait.csv"), row.names=FALSE)
n_sig <- sum(res$padj < 0.05, na.rm=TRUE)
msg(n_sig, " module-ILR pairs significant at BH 5%")
msg("=== ", id, " COMPLETE ===")
