#!/usr/bin/env Rscript
# ===========================================================================
# filt-044  Filter robustness: which modules survive all 31 expression filters?
#
# Spec: reuse 31 filters from de_v2/scripts/13_filter_sweep_wb.R;
#       which modules survive all?
#
# The 31 filters are:
#   30 count-based: >C counts in >=P% of samples,
#                   C in {0,1,5,10,20,50} x P in {50,75,90,99,100}%
#    1 edgeR:       filterByExpr()
#   baseline filter = >10 in >=90% (c10_p90) = 12,368 genes in base-001.
#
# For each filter we compute:
#   - which genes from each base-001 module survive
#   - n_surviving and fraction_surviving per module
#   - "survives all" = gene is in the intersection of all 31 filter gene sets
#
# No DE analysis, no module re-computation — purely descriptive.
# Dependency: base-001 reference cache.
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "filt-044"
stopifnot(id == "filt-044")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

run_dir <- file.path("network/runs", id)
cfg_dir <- "network/config"
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg_dir, showWarnings = FALSE, recursive = TRUE)

REF_CACHE_KEY <- "control_all_bicor_signed_p12"
REF_CACHE     <- file.path("network/cache", paste0("ref_", REF_CACHE_KEY, ".rds"))
COUNTS_FILE   <- "data/cohort_474/counts_474.tsv.gz"

if (!file.exists(REF_CACHE)) {
  msg("DEPENDENCY MISSING: base-001 cache. Exiting status=1 to requeue.")
  quit(save = "no", status = 1L)
}
if (!file.exists(COUNTS_FILE)) stop("counts file not found: ", COUNTS_FILE)

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
    spec="reuse 31 filters from de_v2/scripts/13_filter_sweep_wb.R; which modules survive all?",
    filter_grid="30 count-based (C in {0,1,5,10,20,50} x P in {50,75,90,99,100}%) + 1 edgeR filterByExpr",
    baseline_filter="c10_p90: >10 in >=90% of samples = 12,368 genes",
    analysis="gene survival rate per module per filter (no DE, no module re-computation)",
    ref_cache=REF_CACHE_KEY,
    counts_file=COUNTS_FILE,
    git_sha=if(length(sha)&&!is.na(sha))sha else "unknown",
    r_version=paste(R.version$major,R.version$minor,sep="."))
  lines <- mapply(jsn, names(cfg), cfg, SIMPLIFY=TRUE)
  writeLines(c("{", paste(lines, collapse=",\n"), "}"), cfg_path)
  msg("config written: ", cfg_path)
}

suppressPackageStartupMessages({ library(WGCNA); library(edgeR) })
options(stringsAsFactors = FALSE)
source("network/00_load.R")
d <- load_project()

ref_net      <- readRDS(REF_CACHE)
moduleColors <- ref_net$moduleColors
all_genes    <- rownames(d$logtpm)
modules      <- sort(unique(moduleColors[moduleColors != "grey"]))
module_genes <- setNames(lapply(modules, function(m) all_genes[moduleColors == m]), modules)
msg(length(modules), " modules | background N=", length(all_genes))

msg("reading raw counts ...")
cnt <- as.matrix(read.delim(COUNTS_FILE, row.names=1, check.names=FALSE))
# strip version suffixes (ENSG00000000003.10 -> ENSG00000000003)
rownames(cnt) <- vapply(strsplit(rownames(cnt), ".", fixed=TRUE), `[[`, character(1), 1)
if (anyDuplicated(rownames(cnt))) {
  cnt <- cnt[!duplicated(rownames(cnt)), , drop=FALSE]
}
N_samples <- ncol(cnt)
msg("raw counts: ", nrow(cnt), " x ", N_samples, " (after version strip)")

# ---- build 31 filter gene sets -----------------------------------------------
MC <- c(0L, 1L, 5L, 10L, 20L, 50L)
MP <- c(0.50, 0.75, 0.90, 0.99, 1.00)
GRID <- expand.grid(min_count=MC, min_prop=MP, KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
GRID$fid <- sprintf("c%d_p%02d", GRID$min_count, round(GRID$min_prop * 100))
GRID$rule <- sprintf(">%d in >=%g%%", GRID$min_count, GRID$min_prop * 100)
GRID$is_baseline <- GRID$fid == "c10_p90"

# edgeR filterByExpr
tryCatch({
  dge_tmp <- DGEList(cnt)
  fbe     <- rownames(cnt)[filterByExpr(dge_tmp)]
  GRID    <- rbind(GRID, data.frame(min_count=NA_integer_, min_prop=NA_real_,
                                    fid="edgeR_fbe", rule="edgeR filterByExpr",
                                    is_baseline=FALSE, stringsAsFactors=FALSE))
  GRID$fbe_genes <- NA
  # store edgeR gene list separately
  filter_genes <- lapply(seq_len(nrow(GRID)-1), function(i)
    rownames(cnt)[rowSums(cnt > GRID$min_count[i]) >= ceiling(GRID$min_prop[i] * N_samples)])
  filter_genes[[nrow(GRID)]] <- fbe
  names(filter_genes) <- GRID$fid
}, error=function(e) {
  msg("edgeR filterByExpr failed: ", conditionMessage(e))
  filter_genes <<- lapply(seq_len(nrow(GRID)), function(i)
    rownames(cnt)[rowSums(cnt > GRID$min_count[i]) >= ceiling(GRID$min_prop[i] * N_samples)])
  names(filter_genes) <<- GRID$fid
})
if (!exists("filter_genes")) {
  filter_genes <- lapply(seq_len(nrow(GRID)), function(i)
    rownames(cnt)[rowSums(cnt > GRID$min_count[i]) >= ceiling(GRID$min_prop[i] * N_samples)])
  names(filter_genes) <- GRID$fid
}

n_genes_per_filter <- vapply(filter_genes, length, integer(1))
msg("filter gene counts: min=", min(n_genes_per_filter), " max=", max(n_genes_per_filter),
    " | baseline c10_p90=", n_genes_per_filter["c10_p90"])
msg("31 filter gene sets built")

# Genes that survive ALL 31 filters
genes_all_filters <- Reduce(intersect, filter_genes)
msg("genes surviving ALL filters: ", length(genes_all_filters))

# ---- per-module survival table -----------------------------------------------
rows <- list()
for (fid in names(filter_genes)) {
  fg <- filter_genes[[fid]]
  for (mod in modules) {
    mg <- module_genes[[mod]]
    n_surv <- sum(mg %in% fg)
    rows[[length(rows)+1]] <- data.frame(
      filter=fid, module=mod,
      n_module=length(mg), n_in_filter=length(fg),
      n_surviving=n_surv,
      frac_surviving=n_surv/length(mg),
      stringsAsFactors=FALSE)
  }
}
surv_table <- do.call(rbind, rows)
write.csv(surv_table, file.path(run_dir,"results_filter_survival.csv"), row.names=FALSE)

# ---- module-level summary: fraction of genes that survive EACH filter --------
mod_summary <- do.call(rbind, lapply(modules, function(mod) {
  mg <- module_genes[[mod]]
  n_survive_all <- sum(mg %in% genes_all_filters)
  frac_survive_all <- n_survive_all / length(mg)
  min_frac <- min(surv_table$frac_surviving[surv_table$module==mod])
  max_frac <- max(surv_table$frac_surviving[surv_table$module==mod])
  data.frame(module=mod, n_module=length(mg),
             n_survive_all_filters=n_survive_all,
             frac_survive_all_filters=frac_survive_all,
             min_frac_any_filter=min_frac,
             max_frac_any_filter=max_frac,
             stringsAsFactors=FALSE)
}))
write.csv(mod_summary, file.path(run_dir,"results_module_survival_summary.csv"), row.names=FALSE)

# ---- per-gene table: which genes survive all 31 filters, per module ----------
gene_rows <- list()
for (mod in modules) {
  mg <- module_genes[[mod]]
  for (g in mg) {
    n_filters_survived <- sum(vapply(filter_genes, function(fg) g %in% fg, logical(1)))
    gene_rows[[length(gene_rows)+1]] <- data.frame(
      module=mod, gene=g, sym=d$sym[g],
      n_filters_survived=n_filters_survived,
      survives_all=g %in% genes_all_filters,
      stringsAsFactors=FALSE)
  }
}
gene_table <- do.call(rbind, gene_rows)
write.csv(gene_table, file.path(run_dir,"results_gene_survival.csv"), row.names=FALSE)

# ---- print summary -----------------------------------------------------------
msg("=== Filter survival summary ===")
for (i in seq_len(nrow(mod_summary))) {
  m <- mod_summary[i,]
  msg(m$module, ": ", m$n_survive_all_filters, "/", m$n_module,
      " (", round(m$frac_survive_all_filters*100,1), "%) survive all filters;",
      " worst filter=", round(m$min_frac_any_filter*100,1), "%")
}
msg("=== ", id, " COMPLETE ===")
