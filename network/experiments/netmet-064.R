#!/usr/bin/env Rscript
# ===========================================================================
# netmet-064  Diameter, path length and global connectivity, on the same
#             draws as netmet-060 and netmet-062
#
# WHY THIS EXISTS SEPARATELY
#
# Diameter and characteristic path length were part of the original request
# and netmet-060 did not compute them. They are split into their own row
# rather than folded back in because netmet-060 and netmet-062 have already
# run: re-running them to add two numbers would cost six hours of cluster
# time and would redraw nothing else.
#
# THE DRAWS ARE THE SAME DRAWS. The seed, the group pools and the sampling
# sequence are copied from netmet-060 exactly, so draw k here holds the same
# 296 subjects as draw k there and the two result sets join per draw. An
# assertion below checks that against netmet-060's own saved selections
# whenever those files are present.
#
# WHY SAMPLED SOURCES RATHER THAN EXACT DIAMETER
#
# Exact diameter is all-pairs shortest paths: 12,368 nodes and 1.0-3.9M
# edges, so roughly 5e10 operations per graph, times four groups, ten draws
# and two variants. Instead each graph gets BFS from N_SRC random sources:
#
#   diameter_lb   max eccentricity over the sampled sources. A LOWER BOUND
#                 on the true diameter, never an estimate of it -- any real
#                 source can only be further from something, not nearer.
#   apl           mean finite hop distance over sampled pairs. This one is
#                 unbiased: the sources are uniform, so the pairs are too.
#   efficiency    mean of 1/d over sampled pairs, which unlike apl stays
#                 defined when a pair is unreachable.
#
# Distances are HOP COUNTS on the thresholded graph, not weighted paths.
# A co-expression weight is a similarity, so a weighted shortest path needs
# an arbitrary similarity-to-distance map (1/w, -log w) and the answer moves
# with the choice. Hop count on a graph thresholded at a pre-specified 0.05
# has no such free parameter.
#
# WHAT ELSE IS RECORDED, because it is nearly free once the graph is built
# and because at this density it is what diameter mostly reflects anyway:
# edge count, density, component structure and the global clustering
# coefficient.
# ===========================================================================

VARIANTS <- list(
  "netmet-064" = NULL,                                            # raw, as netmet-060
  "netmet-065" = ~ age + sex + tobacco + rin + plate + seqpc1 + seqpc2 + seqpc3)

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "netmet-064"
stopifnot(id %in% names(VARIANTS))
RESID <- VARIANTS[[id]]

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))
msg("=== ", id, " starting ===")

run_dir <- file.path("network/runs", id)
cfg_dir <- "network/config"
for (d_ in c(run_dir, cfg_dir)) dir.create(d_, showWarnings = FALSE, recursive = TRUE)

## ---- parameters, copied from netmet-060 where they are shared ---------------
SEED_BASE  <- 20261009L + 60000L   # identical to netmet-060: same draws
N_DRAWS    <- 10L
N_MATCH    <- 74L
POWER      <- 12L
NET_TYPE   <- "signed"
COR_FNC    <- "bicor"
BTW_THRESH <- 0.05                 # the same threshold netmet-060 used
N_SRC      <- 300L                 # BFS sources per graph
REF_CACHE  <- "network/cache/ref_control_all_bicor_signed_p12_ds4_mch00_offsetrm.rds"

SLICE_ID <- suppressWarnings(as.integer(Sys.getenv("SLICE_ID", "0")))
if (is.na(SLICE_ID)) SLICE_ID <- 0L
SUBSET <- suppressWarnings(as.integer(Sys.getenv("NETMET_SUBSET", "0")))
if (is.na(SUBSET)) SUBSET <- 0L
msg("SLICE_ID=", SLICE_ID, " (0 = every draw)")

if (!file.exists(REF_CACHE)) { msg("DEPENDENCY: missing ", REF_CACHE); quit(save = "no", status = 1L) }

suppressPackageStartupMessages({ library(WGCNA); library(igraph) })
enableWGCNAThreads(); options(stringsAsFactors = FALSE)
source("network/00_load.R")
d <- load_project(cache = TRUE, verbose = FALSE)
genes <- rownames(d$logtpm)

if (SUBSET > 0L) {
  msg("*** NETMET_SUBSET=", SUBSET, " -- TEST RUN, results are not valid ***")
  set.seed(1L)
  genes <- sort(sample(genes, min(SUBSET, length(genes))))
  run_dir <- paste0(run_dir, "_test")
  dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
}

sid     <- colnames(d$logtpm)
ctrl_id <- sid[!is.na(d$groups) & d$groups == "control"]
nol_id  <- sid[!is.na(d$groups) & d$groups == "bp_nolith"]
lit_id  <- sid[!is.na(d$groups) & d$groups == "bp_lith"]
stopifnot(length(ctrl_id) >= 2L * N_MATCH, length(nol_id) == N_MATCH,
          length(lit_id) >= N_MATCH)

if (is.null(RESID)) {
  Emat <- function(ids) expr_samples(d, ids)[, genes, drop = FALSE]
} else {
  msg("residualising on ", deparse(RESID), " over all ", ncol(d$logtpm), " samples ...")
  ERES <- expr_for(d, group = NULL, residualise = RESID)
  Emat <- function(ids) ERES[ids, genes, drop = FALSE]
}

if (SUBSET == 0L) {
  cfg_path <- file.path(cfg_dir, paste0(id, ".json"))
  if (!file.exists(cfg_path)) {
    sha <- tryCatch(trimws(system("git rev-parse --short HEAD", intern = TRUE)[1]),
                    error = function(e) "unknown", warning = function(w) "unknown")
    writeLines(c("{",
      paste0('  "id": "', id, '", "family": "netmetrics",'),
      '  "spec": "diameter lower bound, characteristic path length, efficiency and connectivity on netmet-060 draws",',
      paste0('  "n_draws": ', N_DRAWS, ', "n_matched": ', N_MATCH, ', "n_sources": ', N_SRC, ','),
      '  "distance": "unweighted hop count on the graph thresholded at adjacency > 0.05",',
      '  "diameter": "LOWER BOUND: max eccentricity over sampled sources, never an estimate",',
      paste0('  "residualise": ', if (is.null(RESID)) '"none"' else paste0('"', deparse(RESID), '"'), ','),
      paste0('  "seed_base": ', SEED_BASE, ', "shares_draws_with": "netmet-060",'),
      paste0('  "git_sha": "', sha, '", "r_version": "',
             paste(R.version$major, R.version$minor, sep = "."), '"'),
      "}"), cfg_path)
    msg("config written: ", cfg_path)
  }
}

## ---- one graph --------------------------------------------------------------
graph_stats <- function(Ex, label) {
  msg("  [", label, "] adjacency ", nrow(Ex), " x ", ncol(Ex), " ...")
  A <- adjacency(Ex, power = POWER, type = NET_TYPE, corFnc = COR_FNC)
  diag(A) <- 0
  A[A <= BTW_THRESH] <- 0
  g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
  rm(A); gc(verbose = FALSE)

  nv <- vcount(g); ne <- ecount(g)
  cmp <- components(g)
  giant <- max(cmp$csize)
  msg("  [", label, "] ", nv, " nodes, ", ne, " edges, ", cmp$no,
      " component(s), giant ", giant)

  # BFS from a fixed sample of sources. Seeded off the label so each group of a
  # draw uses a different source set -- sharing one set across groups would make
  # the four numbers agree more than the graphs do.
  set.seed(SEED_BASE + nchar(label) + ne %% 1000L)
  src <- sample.int(nv, min(N_SRC, nv))
  D <- igraph::distances(g, v = src, mode = "all", weights = NA)
  fin <- is.finite(D) & D > 0
  ecc <- apply(D, 1, function(r) { r <- r[is.finite(r) & r > 0]; if (length(r)) max(r) else NA_real_ })

  out <- data.frame(
    group        = label,
    n_nodes      = nv,
    n_edges      = ne,
    density      = 2 * ne / (nv * (nv - 1)),
    mean_degree  = 2 * ne / nv,
    n_components = cmp$no,
    giant_frac   = giant / nv,
    n_isolated   = sum(degree(g) == 0),
    diameter_lb  = suppressWarnings(max(ecc, na.rm = TRUE)),
    mean_ecc     = mean(ecc, na.rm = TRUE),
    apl          = mean(D[fin]),
    efficiency   = mean(1 / D[fin]),
    p_unreachable = mean(!is.finite(D[row(D) != col(D)])),
    transitivity = transitivity(g, type = "global"),
    n_sources    = length(src))
  rm(g, D); gc(verbose = FALSE)
  out
}

## ---- one draw ---------------------------------------------------------------
do_draw <- function(k) {
  out <- file.path(run_dir, sprintf("draw_%02d.csv", k))
  if (file.exists(out)) { msg("draw ", k, " exists; skipping"); return(invisible(NULL)) }
  msg("--- draw ", k, " of ", N_DRAWS, " ---")

  # Sampling copied verbatim from netmet-060 so the draws coincide.
  set.seed(SEED_BASE + k)
  cs  <- sample(ctrl_id, 2L * N_MATCH, replace = FALSE)
  sel <- list(ctrlA  = cs[seq_len(N_MATCH)],
              ctrlB  = cs[N_MATCH + seq_len(N_MATCH)],
              nolith = nol_id,
              lith   = sample(lit_id, N_MATCH, replace = FALSE))

  # If netmet-060's draw is on disk, prove the selections match rather than
  # assuming the copied code still agrees with it.
  prior <- file.path("network/runs/netmet-060", sprintf("draw_%02d.rds", k))
  if (SUBSET == 0L && file.exists(prior)) {
    P <- readRDS(prior)$sel
    same <- all(vapply(names(sel), function(g) identical(sort(sel[[g]]), sort(P[[g]])), TRUE))
    if (!same) stop("draw ", k, " does not match netmet-060's selection; the sampling has drifted")
    msg("  draw ", k, " selection verified identical to netmet-060")
  }

  R <- do.call(rbind, lapply(names(sel), function(g) graph_stats(Emat(sel[[g]]), g)))
  R$draw <- k
  R$variant <- if (is.null(RESID)) "raw" else "resid"
  write.csv(R, out, row.names = FALSE)
  msg("draw ", k, " written: ", out)
}

if (SLICE_ID > 0L) {
  if (SLICE_ID > N_DRAWS) msg("SLICE_ID ", SLICE_ID, " > N_DRAWS; nothing to do") else do_draw(SLICE_ID)
  msg("=== ", id, " slice ", SLICE_ID, " COMPLETE ===")
  quit(save = "no", status = 0L)
}

for (k in seq_len(N_DRAWS)) do_draw(k)

## ---- summarise whatever exists ----------------------------------------------
f <- list.files(run_dir, pattern = "^draw_\\d+\\.csv$", full.names = TRUE)
if (length(f)) {
  A <- do.call(rbind, lapply(f, read.csv))
  write.csv(A, file.path(run_dir, "all_draws.csv"), row.names = FALSE)
  msg("")
  msg("=== mean over ", length(unique(A$draw)), " draw(s) ===")
  cols <- c("n_edges", "mean_degree", "diameter_lb", "apl", "efficiency",
            "transitivity", "giant_frac")
  agg <- aggregate(A[, cols], by = list(group = A$group), FUN = mean)
  for (i in seq_len(nrow(agg)))
    msg(sprintf("  %-7s edges %9.0f  meandeg %6.1f  diam>=%3.0f  apl %5.3f  eff %5.3f  C %5.3f  giant %5.3f",
        agg$group[i], agg$n_edges[i], agg$mean_degree[i], agg$diameter_lb[i],
        agg$apl[i], agg$efficiency[i], agg$transitivity[i], agg$giant_frac[i]))
  msg("")
  msg("  Paired tests use the same within-kind vs cross-diagnosis split as")
  msg("  netmet-061; run network/experiments/netmet-066.R to score them.")
}
msg("=== ", id, " COMPLETE ===")
