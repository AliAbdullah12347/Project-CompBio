#!/usr/bin/env Rscript
# ===========================================================================
# cent-040  Degree, betweenness, eigenvector centrality; between-group differences
#
# Spec (from queue):
#   degree, betweenness, eigenvector centrality; between-group differences
#
# Design:
#   For each group (control, bp_nolith, bp_lith ONE fixed subsample):
#     1. Compute adjacency: A = bicor-based, signed, power 12
#        WGCNA: adjacency(E, power=12, type="signed", corFnc="bicor")
#     2. Convert to igraph: graph_from_adjacency_matrix(A, weighted=TRUE)
#     3. Centralities:
#          - degree:       strength(g)                 (weighted degree)
#          - eigenvector:  eigen_centrality(g)$vector  (weighted)
#          - betweenness:  threshold A > 0.05 first (sparse);
#                          betweenness(g_thresh, normalized=TRUE)
#     4. bp_lith: ONE fixed subsample (seed=SEED_BASE) for efficiency
#        (betweenness is O(n^3) on 12368 nodes — avoid 20 draws)
#     5. Compare centralities between groups per module: Wilcoxon test,
#        BH correction across modules.
#
# CHECKPOINT per group: save centrality_<group>.csv before comparing.
# Skip computation if checkpoint exists.
#
# Non-negotiable design points obeyed:
#   - betweenness thresholded at A > 0.05 (sparse; O(n^3) prohibitive otherwise)
#   - SEED_BASE = 20261009L + 40000L
#   - Config written before any computation
#   - igraph already approved; no new packages beyond approved list
#   - Dependency check: REF_CACHE must exist; quit(status=1) if missing
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "cent-040"
stopifnot(id == "cent-040")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

# ---------- 1. Directories ---------------------------------------------------
run_dir <- file.path("network/runs", id)
cfg_dir <- "network/config"
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- 2. Parameters — ALL pre-specified BEFORE any computation ---------
SEED_BASE        <- 20261009L + 40000L
SUBSAMPLE_N      <- 74L
POWER            <- 12L
NET_TYPE         <- "signed"
COR_FNC          <- "bicor"
ADJ_THRESH       <- 0.05     # betweenness: zero out edges with A <= 0.05
REF_CACHE        <- "network/cache/ref_control_all_bicor_signed_p12.rds"

# ---------- 3. Dependency check ----------------------------------------------
if (!file.exists(REF_CACHE)) {
  msg("DEPENDENCY: base-001 cache missing: ", REF_CACHE)
  quit(save = "no", status = 1L)
}

# ---------- 4. Write config BEFORE any computation ---------------------------
cfg_path <- file.path(cfg_dir, paste0(id, ".json"))
if (!file.exists(cfg_path)) {
  scl <- function(v) {
    if (is.logical(v)) return(if (isTRUE(v)) "true" else "false")
    if (is.integer(v) || is.numeric(v)) return(as.character(v[1]))
    paste0('"', gsub('"', '\\"', as.character(v[1]), fixed = TRUE), '"')
  }
  arr <- function(v) paste0("[", paste(vapply(v, scl, ""), collapse = ", "), "]")
  jsn <- function(nm, v) paste0('  "', nm, '": ', if (length(v) > 1) arr(v) else scl(v))
  sha <- tryCatch(
    trimws(system("git rev-parse --short HEAD 2>/dev/null", intern = TRUE)[1]),
    error = function(e) "unknown"
  )
  cfg_entries <- list(
    id                    = id,
    family                = "central",
    spec                  = "degree, betweenness, eigenvector centrality; between-group differences",
    ref_cache             = REF_CACHE,
    groups                = c("control", "bp_nolith", "bp_lith"),
    adjacency_type        = NET_TYPE,
    adjacency_corFnc      = COR_FNC,
    power                 = POWER,
    adjacency_formula     = "((1 + bicor(E)) / 2)^12 via WGCNA adjacency()",
    bp_lith_draws         = 1L,
    bp_lith_note          = "ONE fixed subsample (SEED_BASE) — betweenness O(n^3) prohibits 20 draws",
    subsample_n_lith      = SUBSAMPLE_N,
    betweenness_threshold = ADJ_THRESH,
    betweenness_note      = paste0("edges with A <= ", ADJ_THRESH, " zeroed before betweenness (12368 nodes, O(n^3) cost)"),
    degree_metric         = "strength(g) — weighted degree = sum of adjacency weights",
    eigenvector_metric    = "eigen_centrality(g, weights=E(g)$weight)$vector",
    betweenness_metric    = "betweenness(g_thresh, normalized=TRUE) on thresholded graph",
    comparison_test       = "Wilcoxon rank-sum (two-sided), BH across modules",
    checkpoint_per_group  = TRUE,
    seed_base             = SEED_BASE,
    git_sha               = if (length(sha) && !is.na(sha)) sha else "unknown",
    r_version             = paste(R.version$major, R.version$minor, sep = "."),
    wgcna_version         = as.character(packageVersion("WGCNA")),
    igraph_version        = as.character(packageVersion("igraph"))
  )
  lines <- mapply(jsn, names(cfg_entries), cfg_entries, SIMPLIFY = TRUE)
  writeLines(c("{", paste(lines, collapse = ",\n"), "}"), cfg_path)
  msg("config written: ", cfg_path)
}

# ---------- 5. Packages ------------------------------------------------------
suppressPackageStartupMessages({
  library(WGCNA)
  library(fastcluster)
  library(igraph)
})
enableWGCNAThreads()
options(stringsAsFactors = FALSE)
msg("WGCNA ", as.character(packageVersion("WGCNA")),
    " | igraph ", as.character(packageVersion("igraph")),
    " loaded; threads enabled")

# ---------- 6. Data ----------------------------------------------------------
source("network/00_load.R")
d <- load_project()

ref_net      <- readRDS(REF_CACHE)
moduleColors <- ref_net$moduleColors
modules      <- sort(unique(moduleColors[moduleColors != "grey"]))
genes        <- rownames(d$logtpm)
msg(length(modules), " non-grey modules; ", length(genes), " genes")

# ---------- 7. Helper: compute all centrality measures for one expr matrix ----
# E: samples x genes (WGCNA orientation)
# Returns data.frame: gene, degree, eigenvector, betweenness, module
compute_centrality <- function(E, group_label) {
  msg(group_label, ": computing adjacency (",
      nrow(E), " samples x ", ncol(E), " genes) ...")
  A <- adjacency(E, power = POWER, type = NET_TYPE, corFnc = COR_FNC)
  # A is genes x genes; diagonal = 1, set to 0 before igraph
  diag(A) <- 0

  msg(group_label, ": building igraph (full weighted graph) ...")
  g <- graph_from_adjacency_matrix(A, mode = "undirected",
                                   weighted = TRUE, diag = FALSE)

  msg(group_label, ": degree centrality (weighted) ...")
  deg <- strength(g)

  msg(group_label, ": eigenvector centrality (weighted) ...")
  eig_res <- tryCatch(
    eigen_centrality(g, weights = E(g)$weight)$vector,
    error = function(e) {
      msg("WARNING: eigen_centrality failed: ", conditionMessage(e), "; using NA")
      setNames(rep(NA_real_, vcount(g)), V(g)$name)
    }
  )

  # Betweenness: threshold at ADJ_THRESH to get sparse graph
  msg(group_label, ": thresholding A at ", ADJ_THRESH,
      " for betweenness (sparse graph) ...")
  A_thresh <- A
  A_thresh[A_thresh <= ADJ_THRESH] <- 0
  n_edges_sparse <- sum(A_thresh > 0) / 2
  msg(group_label, ": sparse graph has ", n_edges_sparse, " edges (threshold A > ",
      ADJ_THRESH, ")")
  g_thresh <- graph_from_adjacency_matrix(A_thresh, mode = "undirected",
                                          weighted = TRUE, diag = FALSE)
  msg(group_label, ": betweenness (normalized) on sparse graph ...")
  btw <- tryCatch(
    betweenness(g_thresh, normalized = TRUE),
    error = function(e) {
      msg("WARNING: betweenness failed: ", conditionMessage(e), "; using NA")
      setNames(rep(NA_real_, vcount(g_thresh)), V(g_thresh)$name)
    }
  )

  # Assemble output: align to gene order
  data.frame(
    gene        = genes,
    degree      = deg[genes],
    eigenvector = eig_res[genes],
    betweenness = btw[genes],
    module      = moduleColors,
    group       = group_label,
    stringsAsFactors = FALSE,
    row.names   = NULL
  )
}

# ---------- 8. control centrality --------------------------------------------
ctrl_ckpt <- file.path(run_dir, "centrality_control.csv")
if (file.exists(ctrl_ckpt)) {
  msg("centrality_control.csv checkpoint exists; loading ...")
  cent_ctrl <- read.csv(ctrl_ckpt, stringsAsFactors = FALSE)
} else {
  E_ctrl    <- expr_for(d, group = "control")   # 234 x 12368
  cent_ctrl <- compute_centrality(E_ctrl, "control")
  write.csv(cent_ctrl, ctrl_ckpt, row.names = FALSE)
  msg("centrality_control.csv written: ", nrow(cent_ctrl), " genes")
  rm(E_ctrl); gc()
}

# ---------- 9. bp_nolith centrality ------------------------------------------
nolith_ckpt <- file.path(run_dir, "centrality_nolith.csv")
if (file.exists(nolith_ckpt)) {
  msg("centrality_nolith.csv checkpoint exists; loading ...")
  cent_nolith <- read.csv(nolith_ckpt, stringsAsFactors = FALSE)
} else {
  E_nolith    <- expr_for(d, group = "bp_nolith")   # 74 x 12368
  cent_nolith <- compute_centrality(E_nolith, "bp_nolith")
  write.csv(cent_nolith, nolith_ckpt, row.names = FALSE)
  msg("centrality_nolith.csv written: ", nrow(cent_nolith), " genes")
  rm(E_nolith); gc()
}

# ---------- 10. bp_lith centrality (ONE fixed subsample) ---------------------
lith_ckpt <- file.path(run_dir, "centrality_lith.csv")
if (file.exists(lith_ckpt)) {
  msg("centrality_lith.csv checkpoint exists; loading ...")
  cent_lith <- read.csv(lith_ckpt, stringsAsFactors = FALSE)
} else {
  lith_sids <- colnames(d$logtpm)[!is.na(d$groups) & d$groups == "bp_lith"]
  stopifnot(length(lith_sids) == 152L)
  set.seed(SEED_BASE)
  sel_sids  <- sample(lith_sids, SUBSAMPLE_N, replace = FALSE)
  msg("bp_lith: using fixed subsample of ", length(sel_sids), " samples (seed=", SEED_BASE, ")")
  E_lith    <- t(d$logtpm[, sel_sids, drop = FALSE])   # 74 x 12368
  cent_lith <- compute_centrality(E_lith, "bp_lith")
  write.csv(cent_lith, lith_ckpt, row.names = FALSE)
  msg("centrality_lith.csv written: ", nrow(cent_lith), " genes")
  rm(E_lith); gc()
}

# ---------- 11. Per-module comparison: Wilcoxon test, BH correction ----------
msg("computing per-module centrality comparisons ...")
metrics <- c("degree", "eigenvector", "betweenness")

comparison_rows <- do.call(rbind, lapply(modules, function(m) {
  # Genes in module (same set for all groups; module labels from ref)
  in_mod_ctrl   <- cent_ctrl$module == m
  in_mod_nolith <- cent_nolith$module == m
  in_mod_lith   <- cent_lith$module == m
  if (!any(in_mod_ctrl)) return(NULL)

  do.call(rbind, lapply(metrics, function(met) {
    v_ctrl   <- cent_ctrl[[met]][in_mod_ctrl]
    v_nolith <- cent_nolith[[met]][in_mod_nolith]
    v_lith   <- cent_lith[[met]][in_mod_lith]

    # ctrl vs nolith
    wt_nolith <- tryCatch(
      wilcox.test(v_ctrl, v_nolith, exact = FALSE),
      error = function(e) list(statistic = NA, p.value = NA)
    )
    # ctrl vs lith
    wt_lith <- tryCatch(
      wilcox.test(v_ctrl, v_lith, exact = FALSE),
      error = function(e) list(statistic = NA, p.value = NA)
    )

    data.frame(
      module              = m,
      metric              = met,
      n_genes             = sum(in_mod_ctrl),
      mean_ctrl           = mean(v_ctrl,   na.rm = TRUE),
      mean_nolith         = mean(v_nolith, na.rm = TRUE),
      mean_lith           = mean(v_lith,   na.rm = TRUE),
      W_ctrl_nolith       = as.numeric(wt_nolith$statistic),
      pval_ctrl_nolith    = as.numeric(wt_nolith$p.value),
      W_ctrl_lith         = as.numeric(wt_lith$statistic),
      pval_ctrl_lith      = as.numeric(wt_lith$p.value),
      stringsAsFactors    = FALSE
    )
  }))
}))

# BH correction: separately per metric (across modules)
for (met in metrics) {
  idx_nolith <- comparison_rows$metric == met
  idx_lith   <- comparison_rows$metric == met
  if (sum(idx_nolith) > 0) {
    comparison_rows$padj_ctrl_nolith[idx_nolith] <-
      p.adjust(comparison_rows$pval_ctrl_nolith[idx_nolith], method = "BH")
    comparison_rows$padj_ctrl_lith[idx_lith] <-
      p.adjust(comparison_rows$pval_ctrl_lith[idx_lith], method = "BH")
  }
}
comparison_rows$padj_ctrl_nolith[is.na(comparison_rows$padj_ctrl_nolith)] <- NA_real_
comparison_rows$padj_ctrl_lith[is.na(comparison_rows$padj_ctrl_lith)]     <- NA_real_

write.csv(comparison_rows, file.path(run_dir, "centrality_comparison.csv"),
          row.names = FALSE)
msg("centrality_comparison.csv written: ", nrow(comparison_rows), " module-metric rows")

# Summary
for (met in metrics) {
  sub <- comparison_rows[comparison_rows$metric == met, ]
  n_sig_lith   <- sum(sub$padj_ctrl_lith   < 0.05, na.rm = TRUE)
  n_sig_nolith <- sum(sub$padj_ctrl_nolith < 0.05, na.rm = TRUE)
  msg(met, ": ", n_sig_lith, " modules differ ctrl vs lith (BH 5%); ",
      n_sig_nolith, " modules differ ctrl vs nolith (BH 5%)")
}

msg("=== ", id, " COMPLETE ===")
