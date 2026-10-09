#!/usr/bin/env Rscript
# ===========================================================================
# ann-041  Module overlap with lithium + bipolar DEG lists
#
# Spec: module overlap with lithium + bipolar DEG lists (de_final/results/);
#       Fisher exact, BH across modules
#
# For each base-001 module: Fisher exact test for over/under-representation of
# significant DEGs from de_final/results/significant_genes.csv.
# DEG lists used:
#   WB_LI  (1426 genes) — lithium effect in whole blood
#   WB_BPD (4 genes)    — bipolar disorder effect
#   CT_LI  (161 genes)  — cell-type specific lithium effect
#
# Background: 12,368 filtered genes in base-001.
# BH correction across modules within each DEG list.
#
# Dependency: requires base-001 reference cache.
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "ann-041"
stopifnot(id == "ann-041")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

run_dir <- file.path("network/runs", id)
cfg_dir <- "network/config"
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg_dir, showWarnings = FALSE, recursive = TRUE)

REF_CACHE_KEY <- "control_all_bicor_signed_p12"
REF_CACHE     <- file.path("network/cache", paste0("ref_", REF_CACHE_KEY, ".rds"))
DEG_FILE      <- "de_final/results/significant_genes.csv"
SEED          <- 20261009L

if (!file.exists(REF_CACHE)) {
  msg("DEPENDENCY MISSING: ", REF_CACHE, " — base-001 must complete first")
  quit(save = "no", status = 1L)
}
if (!file.exists(DEG_FILE)) stop("DEG file not found: ", DEG_FILE)

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
  cfg <- list(id=id, family="annot",
    spec="module overlap with lithium + bipolar DEG lists (de_final/results/); Fisher exact, BH across modules",
    deg_file=DEG_FILE, ref_cache=REF_CACHE_KEY,
    deg_lists="WB_LI;WB_BPD;CT_LI",
    background="12368 filtered genes",
    test="Fisher exact (one-sided, over-representation)",
    correction="BH across modules within each DEG list",
    seed=SEED,
    git_sha=if(length(sha)&&!is.na(sha))sha else "unknown",
    r_version=paste(R.version$major,R.version$minor,sep="."))
  lines <- mapply(jsn, names(cfg), cfg, SIMPLIFY=TRUE)
  writeLines(c("{", paste(lines, collapse=",\n"), "}"), cfg_path)
  msg("config written: ", cfg_path)
}

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(WGCNA))

source("network/00_load.R")
d <- load_project()

ref_net      <- readRDS(REF_CACHE)
moduleColors <- ref_net$moduleColors
all_genes    <- rownames(d$logtpm)
N            <- length(all_genes)
modules      <- sort(unique(moduleColors[moduleColors != "grey"]))
msg(length(modules), " modules | background N=", N)

deg_df   <- read.csv(DEG_FILE, stringsAsFactors=FALSE)
deg_lists <- split(deg_df$gene[deg_df$significant], deg_df$comparison[deg_df$significant])
deg_lists <- lapply(deg_lists, function(g) intersect(g, all_genes))   # filter to background
for (nm in names(deg_lists)) msg("DEG list ", nm, ": ", length(deg_lists[[nm]]), " genes in background")

# ---- Fisher exact per module per DEG list ----------------------------------
results <- list()
for (list_nm in names(deg_lists)) {
  deg_set <- deg_lists[[list_nm]]
  K <- length(deg_set)   # total DEGs in background
  for (mod in modules) {
    in_module  <- all_genes[moduleColors == mod]
    n_mod      <- length(in_module)
    # 2x2: module x DEG
    #         | DEG   | non-DEG
    # module  | a     | b=n_mod-a
    # non-mod | c-a   | N-n_mod-K+a
    a <- sum(in_module %in% deg_set)
    mat <- matrix(c(a, n_mod-a, K-a, N-n_mod-K+a), nrow=2)
    ft  <- fisher.test(mat, alternative="greater")   # over-representation
    results[[length(results)+1]] <- data.frame(
      deg_list=list_nm, module=mod,
      n_module=n_mod, n_deg_total=K,
      n_overlap=a, expected=n_mod*K/N,
      odds_ratio=ft$estimate, pval=ft$p.value,
      stringsAsFactors=FALSE)
  }
}
res <- do.call(rbind, results)
res$padj <- NA_real_
for (list_nm in unique(res$deg_list)) {
  idx <- res$deg_list == list_nm
  res$padj[idx] <- p.adjust(res$pval[idx], method="BH")
}
write.csv(res, file.path(run_dir,"results_deg_overlap.csv"), row.names=FALSE)
n_sig <- sum(res$padj < 0.05, na.rm=TRUE)
msg(n_sig, " significant module-DEG overlaps at BH 5%")

# ---- also write top overlapping modules per list ---------------------------
for (list_nm in names(deg_lists)) {
  sub <- res[res$deg_list==list_nm,]
  sub <- sub[order(sub$pval),]
  msg(list_nm, " top modules: ",
      paste(head(sub$module[sub$padj<0.05],5), collapse=","), " (padj<0.05)")
}
msg("=== ", id, " COMPLETE ===")
