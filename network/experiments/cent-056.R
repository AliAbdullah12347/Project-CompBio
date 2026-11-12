#!/usr/bin/env Rscript
# ===========================================================================
# cent-056  The control-split null cent-040 needs
#
# WHAT cent-040 REPORTS AND WHY IT CANNOT BE READ
#
# cent-040's matched comparison -- bp_nolith against bp_lith, both n=74 --
# says 52 of 71 modules differ on weighted degree and 58 of 71 on eigenvector
# centrality, while betweenness on a thresholded graph says 0 of 71.
#
# That is almost certainly one global difference counted once per module. The
# lithium graph is 11-22% sparser (three independent subsamples agree), 46 of
# the 52 significant modules move in the same direction, and a Wilcoxon over
# ~50 genes drawn from a single correlation matrix treats them as independent
# when they all shift together as the matrix's overall level shifts. The test
# cannot separate "this module is rewired" from "the whole graph is thinner".
#
# THE NULL THAT SETTLES IT
#
# Split the 234 controls into two disjoint halves of 74 and run the identical
# comparison. Both halves come from one group, so by construction NO module
# differs between them for any biological reason. Whatever count the
# procedure returns is what it returns on noise, at the same n, the same gene
# set and the same module definitions.
#
#   if the null count is near 52  -> cent-040's number is the procedure
#                                    talking, and the module counts must be
#                                    withdrawn.
#   if the null count is near 0   -> the lithium comparison is detecting
#                                    something a same-group split does not,
#                                    and the question becomes what.
#
# ONE SPLIT IS ONE DRAW. Each pod runs its own split under its own seed, and
# the counts are pooled into a distribution. The observed value is read
# against it with the add-one estimator, like every other null in this arm.
#
# COST. Two adjacency builds per split at ~15 min each, plus betweenness on
# the thresholded graphs at ~20 min each: about 70 minutes per split. Running
# one split per pod makes that the wall clock for the whole thing.
#
# SPLIT_ID selects which split this process runs (1..N_SPLITS);
# SPLIT_AGGREGATE=1 skips the computation and pools whatever exists.
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "cent-056"
stopifnot(id == "cent-056")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

run_dir <- file.path("network/runs", id)
cfg_dir <- "network/config"
for (d_ in c(run_dir, cfg_dir)) dir.create(d_, showWarnings = FALSE, recursive = TRUE)

# Same network parameters as cent-040, or the comparison is not the same one.
POWER       <- 12L
NET_TYPE    <- "signed"
COR_FNC     <- "bicor"
ADJ_THRESH  <- 0.05
HALF_N      <- 74L            # matches bp_nolith and cent-040's lith subsample
# 20, not the 10 this started with. At K=10 the add-one floor is 1/11 = 0.091,
# so a count that beats every null draw still cannot reach 0.05 -- the same
# K-floor that made perm-031, perm-032 and hub-037 unreadable. Degree came
# back at 52 against a null median of 36.5 and a max of 50, which is exactly
# the case the extra draws are for. A split costs about an hour on its own pod.
N_SPLITS    <- 20L
SEED_BASE   <- 20261009L + 56000L
METRICS     <- c("degree", "eigenvector", "betweenness")

# SLICE_ID is what mkjobs.py exports to each pod; one split per pod.
SPLIT_ID  <- as.integer(Sys.getenv("SLICE_ID", Sys.getenv("SPLIT_ID", "0")))
AGG_ONLY  <- Sys.getenv("SPLIT_AGGREGATE", "0") == "1"

suppressPackageStartupMessages({ library(WGCNA); library(igraph) })
enableWGCNAThreads()
options(stringsAsFactors = FALSE)
source("network/00_load.R")
d <- load_project(cache = TRUE, verbose = FALSE)

REF <- "network/cache/ref_control_all_bicor_signed_p12_ds4_mch00_offsetrm.rds"
if (!ref_cache_ok(REF)) stop("verified reference missing or unstamped: ", REF)
ref <- readRDS(REF)
moduleColors <- ref$moduleColors
genes <- colnames(expr_for(d, group = "control"))
if (is.null(names(moduleColors))) names(moduleColors) <- genes
modules <- setdiff(unique(moduleColors), "grey")
msg("reference: ", length(modules), " non-grey modules, ", length(genes), " genes")

ctrl_ids <- colnames(d$logtpm)[!is.na(d$groups) & d$groups == "control"]
msg("controls available: ", length(ctrl_ids), "; each split uses ",
    2 * HALF_N, " of them")
stopifnot(length(ctrl_ids) >= 2 * HALF_N)

## ---- config -----------------------------------------------------------------
cfg_path <- file.path(cfg_dir, paste0(id, ".json"))
if (!file.exists(cfg_path)) {
  sha <- tryCatch(trimws(system("git rev-parse --short HEAD", intern = TRUE)[1]),
                  error = function(e) "unknown", warning = function(w) "unknown")
  writeLines(c("{",
    paste0('  "id": "', id, '",'),
    '  "family": "centrality",',
    '  "spec": "control-split null for cent-040: two disjoint halves of the control group, identical comparison",',
    paste0('  "power": ', POWER, ', "net_type": "', NET_TYPE, '", "cor_fnc": "', COR_FNC, '",'),
    paste0('  "adj_threshold": ', ADJ_THRESH, ','),
    paste0('  "half_n": ', HALF_N, ', "n_splits": ', N_SPLITS, ','),
    '  "observed_from": "cent-040 padj_nolith_lith: degree 52/71, eigenvector 58/71, betweenness 0/71",',
    '  "estimator": "add-one (1 + sum(null >= obs)) / (1 + K)",',
    paste0('  "seed_base": ', SEED_BASE, ','),
    paste0('  "git_sha": "', if (length(sha) && !is.na(sha)) sha else "unknown", '",'),
    paste0('  "r_version": "', paste(R.version$major, R.version$minor, sep = "."), '",'),
    paste0('  "wgcna_version": "', as.character(packageVersion("WGCNA")), '"'),
    "}"), cfg_path)
  msg("config written: ", cfg_path)
}

## ---- centrality, copied from cent-040 so the comparison is identical --------
compute_centrality <- function(E, label) {
  msg("  ", label, ": adjacency (", nrow(E), " x ", ncol(E), ") ...")
  A <- adjacency(E, power = POWER, type = NET_TYPE, corFnc = COR_FNC)
  diag(A) <- 0
  g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)

  deg <- strength(g)
  eig <- tryCatch(eigen_centrality(g, weights = igraph::E(g)$weight)$vector,
                  error = function(e) { msg("  WARNING eigen: ", conditionMessage(e))
                                        setNames(rep(NA_real_, vcount(g)), V(g)$name) })

  A[A <= ADJ_THRESH] <- 0
  n_edges <- sum(A > 0) / 2
  msg("  ", label, ": sparse graph ", n_edges, " edges (A > ", ADJ_THRESH, ")")
  gt  <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
  rm(A); gc(FALSE)
  btw <- tryCatch(betweenness(gt, normalized = TRUE),
                  error = function(e) { msg("  WARNING betweenness: ", conditionMessage(e))
                                        setNames(rep(NA_real_, vcount(gt)), V(gt)$name) })
  rm(g, gt); gc(FALSE)
  list(tab = data.frame(gene = genes, degree = deg[genes], eigenvector = eig[genes],
                        betweenness = btw[genes], module = unname(moduleColors[genes])),
       n_edges = n_edges)
}

# cent-040's per-module test, verbatim: Wilcoxon on the genes of each module,
# BH across modules within each metric.
compare_halves <- function(ta, tb) {
  rows <- do.call(rbind, lapply(modules, function(m) {
    ia <- ta$module == m; ib <- tb$module == m
    if (!any(ia)) return(NULL)
    do.call(rbind, lapply(METRICS, function(met) {
      va <- ta[[met]][ia]; vb <- tb[[met]][ib]
      w <- tryCatch(wilcox.test(va, vb, exact = FALSE),
                    error = function(e) list(p.value = NA_real_))
      data.frame(module = m, metric = met, n_genes = sum(ia),
                 mean_a = mean(va, na.rm = TRUE), mean_b = mean(vb, na.rm = TRUE),
                 pval = as.numeric(w$p.value))
    }))
  }))
  rows$padj <- ave(rows$pval, rows$metric, FUN = function(x) p.adjust(x, "BH"))
  rows
}

## ---- run this pod's split ---------------------------------------------------
if (!AGG_ONLY && SPLIT_ID >= 1L && SPLIT_ID <= N_SPLITS) {
  ck <- file.path(run_dir, sprintf("split_%02d.rds", SPLIT_ID))
  if (file.exists(ck)) {
    msg("split ", SPLIT_ID, " checkpoint exists, skipping")
  } else {
    set.seed(SEED_BASE + SPLIT_ID)
    picked <- sample(ctrl_ids, 2 * HALF_N, replace = FALSE)
    ids_a  <- picked[1:HALF_N]; ids_b <- picked[(HALF_N + 1):(2 * HALF_N)]
    stopifnot(length(intersect(ids_a, ids_b)) == 0)
    msg("split ", SPLIT_ID, ": halves of ", HALF_N, ", disjoint, drawn from ",
        length(ctrl_ids), " controls")

    Ea <- expr_samples(d, ids_a); Eb <- expr_samples(d, ids_b)
    ra <- compute_centrality(Ea, sprintf("split%02d-A", SPLIT_ID)); rm(Ea); gc(FALSE)
    rb <- compute_centrality(Eb, sprintf("split%02d-B", SPLIT_ID)); rm(Eb); gc(FALSE)

    cmp <- compare_halves(ra$tab, rb$tab)
    counts <- vapply(METRICS, function(met)
      sum(cmp$padj[cmp$metric == met] < 0.05, na.rm = TRUE), integer(1))
    saveRDS(list(split = SPLIT_ID, counts = counts, n_modules = length(modules),
                 n_edges_a = ra$n_edges, n_edges_b = rb$n_edges, cmp = cmp), ck)
    msg("split ", SPLIT_ID, " done: ",
        paste(sprintf("%s=%d", METRICS, counts), collapse = "  "),
        "  | edges ", ra$n_edges, " vs ", rb$n_edges)
  }
}

## ---- aggregate --------------------------------------------------------------
fs <- sort(list.files(run_dir, "^split_[0-9]+[.]rds$", full.names = TRUE))
if (!length(fs)) { msg("no splits yet; nothing to aggregate"); quit(save = "no", status = 0L) }
parts <- lapply(fs, readRDS)
null_tab <- do.call(rbind, lapply(parts, function(x)
  data.frame(split = x$split, as.list(x$counts),
             n_edges_a = x$n_edges_a, n_edges_b = x$n_edges_b)))
null_tab <- null_tab[order(null_tab$split), ]
write.csv(null_tab, file.path(run_dir, "null_split_counts.csv"), row.names = FALSE)

# cent-040's matched comparison, the thing being calibrated.
OBS <- c(degree = 52L, eigenvector = 58L, betweenness = 0L)
obs_file <- "network/runs/cent-040/centrality_comparison.csv"
if (file.exists(obs_file)) {
  cc <- read.csv(obs_file)
  if ("padj_nolith_lith" %in% names(cc)) {
    OBS <- vapply(METRICS, function(met)
      sum(cc$padj_nolith_lith[cc$metric == met] < 0.05, na.rm = TRUE), integer(1))
    msg("observed counts read from cent-040: ",
        paste(sprintf("%s=%d", METRICS, OBS), collapse = "  "))
  }
}

K <- nrow(null_tab)
res <- do.call(rbind, lapply(METRICS, function(met) {
  nv  <- null_tab[[met]]
  obs <- OBS[[met]]
  data.frame(metric = met, observed = obs, n_modules = parts[[1]]$n_modules,
             K = K, null_mean = mean(nv), null_sd = sd(nv),
             null_median = median(nv), null_min = min(nv), null_max = max(nv),
             n_null_ge_obs = sum(nv >= obs),
             empirical_p_addone = (1 + sum(nv >= obs)) / (1 + K))
}))
write.csv(res, file.path(run_dir, "results_control_split_null.csv"), row.names = FALSE)

msg("")
msg("=== control-split null, K=", K, " splits ===")
for (i in seq_len(nrow(res))) {
  r <- res[i, ]
  msg(sprintf("  %-12s observed %2d/%d | null median %4.1f (min %d, max %d) | p = %.3f",
              r$metric, r$observed, r$n_modules, r$null_median, r$null_min,
              r$null_max, r$empirical_p_addone))
}
msg("")
msg("A null median near the observed count means the procedure returns that")
msg("many 'differing' modules on two halves of ONE group, and cent-040's")
msg("counts say nothing about lithium.")

if (K < N_SPLITS) {
  msg("NOTE: ", K, " of ", N_SPLITS, " splits present; partial but valid.")
  quit(save = "no", status = 1L)
}
msg("=== ", id, " COMPLETE ===")
