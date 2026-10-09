#!/usr/bin/env Rscript
# ===========================================================================
# ann-042  Module overlap with LM22 markers and cell-type DE results
#
# Spec: module overlap with LM22 markers and cell-type DE results; Fisher exact, BH
#
# Two overlap analyses:
#   A. LM22 signature genes (547 genes used by CIBERSORTx) vs base-001 modules
#   B. Cell-type specific DE genes (de_final/results/DE_CT_LI.csv.gz) per lineage
#      vs base-001 modules
#
# Background: 12,368 filtered genes in base-001.
# BH correction across modules within each gene set.
#
# Dependency: requires base-001 reference cache.
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "ann-042"
stopifnot(id == "ann-042")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

run_dir <- file.path("network/runs", id)
cfg_dir <- "network/config"
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg_dir, showWarnings = FALSE, recursive = TRUE)

REF_CACHE_KEY <- "control_all_bicor_signed_p12"
REF_CACHE     <- file.path("network/cache", paste0("ref_", REF_CACHE_KEY, ".rds"))
CT_LI_FILE    <- "de_final/results/DE_CT_LI.csv.gz"
CT_BPD_FILE   <- "de_final/results/DE_CT_BPD.csv.gz"
SEED          <- 20261009L

if (!file.exists(REF_CACHE)) {
  msg("DEPENDENCY MISSING: ", REF_CACHE, " — base-001 must complete first")
  quit(save = "no", status = 1L)
}
if (!file.exists(CT_LI_FILE)) stop("CT DE file not found: ", CT_LI_FILE)

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
    spec="module overlap with LM22 markers and cell-type DE results; Fisher exact, BH",
    lm22_source="data/LM22.txt (via d$lm22 from 00_load.R, matched by HGNC symbol)",
    ct_de_source=CT_LI_FILE,
    ref_cache=REF_CACHE_KEY,
    background="12368 filtered genes",
    test="Fisher exact (one-sided, over-representation)",
    correction="BH across modules within each gene set",
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
sym_map      <- d$sym   # Ensembl -> HGNC
msg(length(modules), " modules | background N=", N)

# ---- gene sets to test -----------------------------------------------------
gene_sets <- list()

# A. LM22 signature genes (matched via symbol)
lm22_sym <- d$lm22
lm22_ens <- names(sym_map)[sym_map %in% lm22_sym]
lm22_ens_bg <- intersect(lm22_ens, all_genes)
gene_sets[["LM22_all"]] <- lm22_ens_bg
msg("LM22: ", length(lm22_sym), " signature symbols | ",
    length(lm22_ens_bg), " matched Ensembl IDs in background")

# B. Cell-type specific DE genes (LI), per lineage
if (file.exists(CT_LI_FILE)) {
  ct_li <- read.csv(CT_LI_FILE, stringsAsFactors=FALSE)
  lineages <- unique(ct_li$lineage)
  for (lin in lineages) {
    sub <- ct_li[ct_li$lineage==lin & ct_li$significant, ]
    ens <- intersect(sub$gene, all_genes)
    if (length(ens) > 0) {
      gene_sets[[paste0("CT_LI_", lin)]] <- ens
      msg("CT_LI_", lin, ": ", length(ens), " genes")
    }
  }
}

# C. Cell-type BPD DE genes if available
if (file.exists(CT_BPD_FILE)) {
  ct_bpd <- read.csv(CT_BPD_FILE, stringsAsFactors=FALSE)
  lineages <- unique(ct_bpd$lineage)
  for (lin in lineages) {
    sub <- ct_bpd[ct_bpd$lineage==lin & ct_bpd$significant, ]
    ens <- intersect(sub$gene, all_genes)
    if (length(ens) > 0) {
      gene_sets[[paste0("CT_BPD_", lin)]] <- ens
      msg("CT_BPD_", lin, ": ", length(ens), " genes")
    }
  }
}

msg(length(gene_sets), " gene sets to test")

# ---- Fisher exact per module per gene set ----------------------------------
fisher_rows <- function(gs_name, gene_set) {
  K <- length(gene_set)
  lapply(modules, function(mod) {
    in_module <- all_genes[moduleColors == mod]
    n_mod <- length(in_module)
    a <- sum(in_module %in% gene_set)
    mat <- matrix(c(a, n_mod-a, K-a, N-n_mod-K+a), nrow=2)
    ft  <- fisher.test(mat, alternative="greater")
    data.frame(gene_set=gs_name, module=mod,
               n_module=n_mod, n_geneset=K,
               n_overlap=a, expected=round(n_mod*K/N,2),
               odds_ratio=ft$estimate, pval=ft$p.value,
               stringsAsFactors=FALSE)
  })
}
all_rows <- unlist(lapply(names(gene_sets), function(gs) fisher_rows(gs, gene_sets[[gs]])),
                   recursive=FALSE)
res <- do.call(rbind, all_rows)
res$padj <- NA_real_
for (gs in unique(res$gene_set)) {
  idx <- res$gene_set == gs
  res$padj[idx] <- p.adjust(res$pval[idx], method="BH")
}
write.csv(res, file.path(run_dir,"results_geneset_overlap.csv"), row.names=FALSE)
n_sig <- sum(res$padj < 0.05, na.rm=TRUE)
msg(n_sig, " significant module-geneset overlaps at BH 5%")
for (gs in names(gene_sets)) {
  sub <- res[res$gene_set==gs,]
  sub <- sub[order(sub$pval),]
  top <- head(sub$module[sub$padj<0.05], 5)
  if (length(top)) msg(gs, " top modules: ", paste(top, collapse=","))
}
msg("=== ", id, " COMPLETE ===")
