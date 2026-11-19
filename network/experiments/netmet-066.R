#!/usr/bin/env Rscript
# ===========================================================================
# netmet-066  Redo the lithium graph-shape comparison with a design that can
#             actually support it
#
# WHY netmet-064 CANNOT CARRY ITS OWN RESULT
#
# It reported that lithium thins the graph: 24% fewer edges, longer paths,
# wider diameter, p <= 0.008. Three things are wrong with how that was tested.
#
# 1. nolith HAS NO SAMPLING VARIANCE. It is all 74 bipolar-I subjects off
#    lithium, so the same people every draw and literally the same number in
#    every row: 2,732,350 edges, sd 0. A paired test against a constant asks
#    only whether 10 draws of lith fall below one fixed number. It cannot ask
#    whether that number was itself high by chance, and with one sample there
#    is no way to find out.
#
# 2. EDGE COUNT AT A FIXED THRESHOLD IS A VERY NOISY STATISTIC. Two disjoint
#    halves of the control group -- same people, matched covariates, nothing
#    to find -- differ by a mean of 22.8% and range from -22% to +53%. The
#    standard deviation of that difference is 563,573 edges. The lithium
#    difference is 659,007, which is 1.17 times that. A threshold turns a
#    small shift in the correlation distribution into a large change in count.
#
# 3. THE DIAMETER CLAIM DOES NOT SURVIVE AT ALL. Against the 20 control
#    samples rather than against fixed nolith, lith's diameter is 7.9 against
#    the controls' 8.3 and p = 0.22. The "lithium widens the diameter" result
#    was really "nolith happens to be narrow", and nolith is one sample.
#
# WHAT THIS ROW CHANGES
#
# EVERY GROUP IS SUBSAMPLED, so every group has a sampling distribution, and
# every group gets its own internal yardstick. n = 37 is the largest size at
# which bp_nolith (74) can supply two DISJOINT subsamples:
#
#     ctrlA  ctrlB     two disjoint halves of the control pool
#     nolA   nolB      two disjoint halves of bp_nolith   <- the missing null
#     litA   litB      two disjoint halves of bp_lith
#
# Within-group pairs say what the measurement does when nothing differs, for
# EACH group separately rather than only for controls. Cross-group pairs are
# then read against them.
#
# THRESHOLD-FREE MEASURES ARE PRIMARY. Mean adjacency and the quantiles of the
# correlation distribution use no cutoff, so they cannot be moved by one. Edge
# count is still recorded, at four thresholds rather than one, so the
# threshold sensitivity is visible rather than hidden.
#
# n = 37 is small and the networks are correspondingly noisy. That is the
# price of a valid null, and it is paid by every group equally.
# ===========================================================================

VARIANTS <- list("netmet-066" = NULL,
                 "netmet-067" = ~ age + sex + tobacco + rin + plate + seqpc1 + seqpc2 + seqpc3)

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "netmet-066"
stopifnot(id %in% names(VARIANTS))
RESID <- VARIANTS[[id]]

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))
msg("=== ", id, " starting ===")

run_dir <- file.path("network/runs", id)
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
dir.create("network/config", showWarnings = FALSE, recursive = TRUE)

SEED_BASE  <- 20261009L + 66000L
N_DRAWS    <- 20L
N_SUB      <- 37L            # largest n letting bp_nolith (74) split disjointly
POWER      <- 12L
THRESHOLDS <- c(0.02, 0.05, 0.10, 0.20)
N_SRC      <- 200L
REF_CACHE  <- "network/cache/ref_control_all_bicor_signed_p12_ds4_mch00_offsetrm.rds"

SLICE_ID <- suppressWarnings(as.integer(Sys.getenv("SLICE_ID", "0")))
if (is.na(SLICE_ID)) SLICE_ID <- 0L
SUBSET <- suppressWarnings(as.integer(Sys.getenv("NETMET_SUBSET", "0")))
if (is.na(SUBSET)) SUBSET <- 0L

if (!file.exists(REF_CACHE)) { msg("DEPENDENCY: missing ", REF_CACHE); quit(save = "no", status = 1L) }
suppressPackageStartupMessages({ library(WGCNA); library(igraph) })
enableWGCNAThreads(); options(stringsAsFactors = FALSE)
source("network/00_load.R")
d <- load_project(cache = TRUE, verbose = FALSE)
genes <- rownames(d$logtpm)

if (SUBSET > 0L) {
  msg("*** NETMET_SUBSET=", SUBSET, " -- TEST RUN, results are not valid ***")
  set.seed(1L); genes <- sort(sample(genes, min(SUBSET, length(genes))))
  run_dir <- paste0(run_dir, "_test"); dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
}

sid <- colnames(d$logtpm)
pool <- list(ctrl = sid[!is.na(d$groups) & d$groups == "control"],
             nol  = sid[!is.na(d$groups) & d$groups == "bp_nolith"],
             lit  = sid[!is.na(d$groups) & d$groups == "bp_lith"])
msg("pools: ", paste(sprintf("%s=%d", names(pool), lengths(pool)), collapse = " "))
stopifnot(all(lengths(pool) >= 2L * N_SUB))

if (is.null(RESID)) {
  Emat <- function(ids) expr_samples(d, ids)[, genes, drop = FALSE]
} else {
  msg("residualising on ", deparse(RESID), " ...")
  ERES <- expr_for(d, group = NULL, residualise = RESID)
  Emat <- function(ids) ERES[ids, genes, drop = FALSE]
}

if (SUBSET == 0L) {
  cp <- file.path("network/config", paste0(id, ".json"))
  if (!file.exists(cp)) {
    sha <- tryCatch(trimws(system("git rev-parse --short HEAD", intern = TRUE)[1]),
                    error = function(e) "unknown", warning = function(w) "unknown")
    writeLines(c("{", paste0('  "id": "', id, '", "family": "netmetrics",'),
      '  "spec": "graph shape with EVERY group subsampled, so every group has its own within-group null",',
      paste0('  "n_draws": ', N_DRAWS, ', "n_sub": ', N_SUB, ','),
      '  "why_n37": "largest n at which bp_nolith (74) supplies two disjoint subsamples",',
      '  "primary": "threshold-free: mean adjacency and correlation quantiles",',
      '  "secondary": "edge count at 4 thresholds, to expose threshold sensitivity",',
      paste0('  "residualise": ', if (is.null(RESID)) '"none"' else paste0('"', deparse(RESID), '"'), ','),
      paste0('  "seed_base": ', SEED_BASE, ', "git_sha": "', sha, '"'), "}"), cp)
  }
}

## ---- one subsample ----------------------------------------------------------
shape <- function(ids, label) {
  Ex <- Emat(ids)
  msg("  [", label, "] bicor on ", nrow(Ex), " x ", ncol(Ex), " ...")
  # The correlation matrix, kept once. Everything threshold-free is read off
  # it directly, so no cutoff can move those numbers.
  R <- WGCNA::bicor(Ex, use = "pairwise.complete.obs", nThreads = 1)
  ut <- upper.tri(R)
  r  <- R[ut]
  r  <- r[is.finite(r)]
  A  <- ((1 + R) / 2)^POWER
  a  <- A[ut]; a <- a[is.finite(a)]
  qs <- quantile(r, c(0.5, 0.9, 0.95, 0.99), names = FALSE)

  out <- data.frame(group = label, n = length(ids),
    mean_adjacency = mean(a), median_adjacency = median(a),
    mean_abs_bicor = mean(abs(r)), median_bicor = qs[1],
    q90_bicor = qs[2], q95_bicor = qs[3], q99_bicor = qs[4])
  for (th in THRESHOLDS) out[[paste0("edges_", sub("[.]", "", sprintf("%.2f", th)))]] <- sum(a > th)

  # path length on the 0.05 graph, for continuity with netmet-064
  A[A <= 0.05] <- 0
  diag(A) <- 0
  g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
  rm(A, R); gc(verbose = FALSE)
  set.seed(SEED_BASE + nchar(label) + ecount(g) %% 997L)
  D <- igraph::distances(g, v = sample.int(vcount(g), min(N_SRC, vcount(g))),
                         mode = "all", weights = NA)
  fin <- is.finite(D) & D > 0
  ecc <- apply(D, 1, function(x) { x <- x[is.finite(x) & x > 0]; if (length(x)) max(x) else NA_real_ })
  out$apl <- mean(D[fin]); out$efficiency <- mean(1 / D[fin])
  out$diameter_lb <- suppressWarnings(max(ecc, na.rm = TRUE))
  out$transitivity <- transitivity(g, type = "global")
  rm(g, D); gc(verbose = FALSE)
  out
}

do_draw <- function(k) {
  f <- file.path(run_dir, sprintf("draw_%02d.csv", k))
  if (file.exists(f)) { msg("draw ", k, " exists; skipping"); return(invisible(NULL)) }
  msg("--- draw ", k, " of ", N_DRAWS, " ---")
  set.seed(SEED_BASE + k)
  # Every group split into two DISJOINT subsamples of the same size. The
  # within-group pair is that group's own null.
  sel <- list()
  for (p in names(pool)) {
    x <- sample(pool[[p]], 2L * N_SUB, replace = FALSE)
    sel[[paste0(p, "A")]] <- x[seq_len(N_SUB)]
    sel[[paste0(p, "B")]] <- x[N_SUB + seq_len(N_SUB)]
  }
  stopifnot(all(vapply(sel, length, 1L) == N_SUB),
            length(intersect(sel$nolA, sel$nolB)) == 0L,
            length(intersect(sel$litA, sel$litB)) == 0L)
  R <- do.call(rbind, lapply(names(sel), function(nm) shape(sel[[nm]], nm)))
  R$draw <- k
  R$variant <- if (is.null(RESID)) "raw" else "resid"
  write.csv(R, f, row.names = FALSE)
  msg("draw ", k, " written")
}

if (SLICE_ID > 0L) {
  # 20 draws over 10 slices: two draws each.
  for (k in c(SLICE_ID, SLICE_ID + 10L)) if (k <= N_DRAWS) do_draw(k)
  msg("=== ", id, " slice ", SLICE_ID, " COMPLETE ===")
  quit(save = "no", status = 0L)
}

for (k in seq_len(N_DRAWS)) do_draw(k)
msg("=== ", id, " COMPLETE ===")
