#!/usr/bin/env Rscript
# ===========================================================================
# netmet-060  Seven network-metric families, four groups, one matched design
#
# THE QUESTION
#
#   Build bipolar and non-bipolar co-expression networks. Compute module
#   eigengene correlation, intramodular connectivity, gene significance,
#   topological overlap, differential connectivity, network-propagation
#   proximity, and degree/betweenness/maximal-clique centrality. Ask whether
#   each differs. Then do the same for lithium.
#
# WHY THIS IS NOT cent-040 AGAIN
#
# cent-040 computed three of these metrics on three groups and reported that
# 52 of 71 modules "differ" in degree between control and lithium. cent-056
# then split the control group in half and compared the two halves -- same
# people, same biology, nothing to find -- and got 43 of 71. Almost the whole
# cent-040 result is what this measurement does when nothing is happening.
#
# Two causes, both structural rather than statistical:
#
#   1. UNMATCHED n. control has 234 samples, bp_nolith 74. bicor variance
#      falls with n, so the 234-sample adjacency is systematically tighter
#      than any 74-sample one. Every metric derived from it inherits that.
#      A Wilcoxon test across the genes of a module then detects the sample
#      size, with a p-value that gets smaller the bigger the module is.
#   2. NO SAME-GROUP REFERENCE. Without knowing what the contrast returns on
#      two samples that cannot differ, "52 modules differ" has no scale.
#
# So this row fixes n = 74 everywhere and carries the null inside the design
# rather than bolting it on afterwards. Each draw builds FOUR networks:
#
#      ctrlA   74 controls
#      ctrlB   74 different controls, disjoint from ctrlA
#      nolith  all 74 bipolar-I subjects off lithium
#      lith    74 of the 152 bipolar-I subjects on lithium
#
# and scores three contrasts through the identical code path:
#
#      NULL     ctrlA  vs ctrlB    <- same biology; whatever this returns is noise
#      BIPOLAR  ctrlA  vs nolith
#      LITHIUM  nolith vs lith
#
# A metric is reported as differing only where the BIPOLAR or LITHIUM
# contrast exceeds what NULL returns on the same draw. That comparison is
# paired -- the same draw supplies all three -- so the draw-to-draw variance
# that swamped cent-040 cancels.
#
# THE SEVEN FAMILIES
#
#   1 module eigengene correlation   ME matrix per group; mean |r| off-diagonal
#                                    (the "modules move in lockstep" statistic)
#   2 intramodular connectivity      kIN (within-module weighted degree), kME
#   3 gene significance              GS = |cor(gene, trait)|, paired with kME
#   4 topological overlap            per-module mean within-module TOM, plus a
#                                    fixed pair subsample for cross-group r
#   5 differential connectivity      delta kIN, scaled (computed in netmet-061)
#   6 propagation proximity          personalised PageRank from DE seed sets
#   7 centrality                     degree (strength), betweenness, eigenvector,
#                                    and log maximal-clique centrality (MCC)
#
# SLICING. SLICE_ID picks one draw. Draws are exchangeable, so a partial set
# is still a valid design -- it only widens the confidence interval.
# Aggregation is netmet-061, which reads whatever draw files exist.
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "netmet-060"
stopifnot(id == "netmet-060")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))
msg("=== ", id, " starting ===")

run_dir <- file.path("network/runs", id)
cfg_dir <- "network/config"
for (d_ in c(run_dir, cfg_dir)) dir.create(d_, showWarnings = FALSE, recursive = TRUE)

## ---- pre-specified parameters ----------------------------------------------
SEED_BASE   <- 20261009L + 60000L
N_DRAWS     <- 10L
N_MATCH     <- 74L        # every group is cut to this; see "UNMATCHED n" above
POWER       <- 12L
NET_TYPE    <- "signed"
COR_FNC     <- "bicor"
BTW_THRESH  <- 0.05       # as cent-040, so the two rows remain comparable
N_SEED      <- 100L       # top-ranked DE genes used as propagation seeds
PPR_DAMP    <- 0.85
TOM_PAIRS   <- 2e6L       # fixed random gene pairs for the cross-group TOM r
MCC_DEG     <- 10L        # per-module edge budget: MCC_DEG * n_nodes, i.e. mean
                          # degree 2*MCC_DEG. Bounds clique enumeration, which an
                          # edge cap alone does not -- see mcc_log().
MIN_MOD     <- 15L
REF_CACHE   <- "network/cache/ref_control_all_bicor_signed_p12_ds4_mch00_offsetrm.rds"

SLICE_ID <- suppressWarnings(as.integer(Sys.getenv("SLICE_ID", "0")))
if (is.na(SLICE_ID)) SLICE_ID <- 0L
msg("SLICE_ID=", SLICE_ID, " (0 = run every draw locally)")

if (!file.exists(REF_CACHE)) { msg("DEPENDENCY: missing ", REF_CACHE); quit(save = "no", status = 1L) }

suppressPackageStartupMessages({ library(WGCNA); library(igraph) })
enableWGCNAThreads(); options(stringsAsFactors = FALSE)
source("network/00_load.R")
d <- load_project(cache = TRUE, verbose = FALSE)

ref_net <- readRDS(REF_CACHE)
genes   <- rownames(d$logtpm)
# base-001 caches moduleColors positionally and the preservation rows consume
# them that way, so the cache is left alone and this copy is named instead.
mc <- ref_net$moduleColors
names(mc) <- genes
# TEST HOOK, not a production path. A full draw is four 12368-gene networks
# and takes hours, so there is no way to exercise this file end to end before
# committing a 10-pod job without one. NETMET_SUBSET=N cuts the gene set to N
# and makes every number meaningless -- it exists to prove the code runs, not
# to produce a result. A subsetted run writes draw files under a _test suffix
# so it can never be mistaken for, or aggregated with, a real draw.
SUBSET <- suppressWarnings(as.integer(Sys.getenv("NETMET_SUBSET", "0")))
if (is.na(SUBSET)) SUBSET <- 0L
if (SUBSET > 0L) {
  msg("*** NETMET_SUBSET=", SUBSET, " -- TEST RUN, results are not valid ***")
  set.seed(1L)
  genes <- sort(sample(genes, min(SUBSET, length(genes))))
  mc <- mc[genes]
  run_dir <- paste0(run_dir, "_test")
  dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
}

mods <- sort(setdiff(unique(mc), "grey"))
mods <- mods[vapply(mods, function(m) sum(mc == m) >= MIN_MOD, TRUE)]
msg(length(genes), " genes | ", length(mods), " modules with >= ", MIN_MOD, " genes")

## ---- group membership -------------------------------------------------------
sid     <- colnames(d$logtpm)
ctrl_id <- sid[!is.na(d$groups) & d$groups == "control"]
nol_id  <- sid[!is.na(d$groups) & d$groups == "bp_nolith"]
lit_id  <- sid[!is.na(d$groups) & d$groups == "bp_lith"]
msg("pool: control ", length(ctrl_id), " | bp_nolith ", length(nol_id),
    " | bp_lith ", length(lit_id))
stopifnot(length(ctrl_id) >= 2L * N_MATCH, length(nol_id) == N_MATCH,
          length(lit_id) >= N_MATCH)

# Every matrix in this row goes through here, so the gene set and its ORDER are
# the same everywhere by construction rather than by assumption. group_metrics
# then asserts it, which would otherwise be an assertion about luck.
Emat <- function(ids) expr_samples(d, ids)[, genes, drop = FALSE]

## ---- propagation seeds ------------------------------------------------------
# Top N_SEED by the published DE rank. Fixed before any network is built, and
# identical for all four groups, so a difference in proximity is a difference
# in the graph rather than in the seed set.
read_seed <- function(tag) {
  p <- file.path("de_final/results", paste0("DE_", tag, ".csv.gz"))
  if (!file.exists(p)) { msg("WARNING: no ", p, "; seeds for ", tag, " skipped"); return(character(0)) }
  x <- read.csv(gzfile(p))
  g <- intersect(x$gene[order(x$rank)], genes)
  g[seq_len(min(N_SEED, length(g)))]
}
SEEDS <- list(BPD = read_seed("WB_BPD"), LI = read_seed("WB_LI"))
msg("seeds: BPD ", length(SEEDS$BPD), " | LI ", length(SEEDS$LI))

## ---- config -----------------------------------------------------------------
cfg_path <- file.path(cfg_dir, paste0(id, ".json"))
# A test run must not leave a config behind: config/netmet-060.json is the
# pre-specification record, and one written from a 600-gene smoke run would
# be indistinguishable from the real one.
if (SUBSET == 0L && !file.exists(cfg_path)) {
  sha <- tryCatch(trimws(system("git rev-parse --short HEAD", intern = TRUE)[1]),
                  error = function(e) "unknown", warning = function(w) "unknown")
  writeLines(c("{",
    paste0('  "id": "', id, '", "family": "netmetrics",'),
    '  "spec": "seven network-metric families over four n-matched groups with a control-split null",',
    paste0('  "n_draws": ', N_DRAWS, ', "n_matched": ', N_MATCH, ','),
    '  "groups": ["ctrlA", "ctrlB", "nolith", "lith"],',
    '  "contrasts": {"NULL": "ctrlA vs ctrlB", "BIPOLAR": "ctrlA vs nolith", "LITHIUM": "nolith vs lith"},',
    paste0('  "power": ', POWER, ', "type": "', NET_TYPE, '", "corFnc": "', COR_FNC, '",'),
    paste0('  "betweenness_threshold": ', BTW_THRESH, ', "ppr_damping": ', PPR_DAMP, ','),
    paste0('  "n_seed": ', N_SEED, ', "tom_pairs": ', TOM_PAIRS, ', "mcc_degree_budget": ', MCC_DEG, ','),
    '  "mcc": "log scale: logsumexp(lgamma(clique_size)) over maximal cliques, per module subgraph",',
    '  "decision_rule": "a contrast counts only where it exceeds the paired NULL contrast on the same draw",',
    paste0('  "seed_base": ', SEED_BASE, ', "git_sha": "', sha, '",'),
    paste0('  "r_version": "', paste(R.version$major, R.version$minor, sep = "."), '"'),
    "}"), cfg_path)
  msg("config written: ", cfg_path)
}

## ===========================================================================
## per-group metric computation
## ===========================================================================

# log maximal-clique centrality. cytoHubba defines MCC(v) as the sum over
# maximal cliques C containing v of (|C| - 1)!. On 12368 nodes that is both
# intractable and numerically hopeless -- (|C| - 1)! overflows a double at
# |C| = 171, and a dense module has cliques far larger. So: enumerate per
# module subgraph, and keep the log, logsumexp(lgamma(|C|)), which is exact
# and bounded. Ranking is unchanged by the monotone transform.
#
# Sparsification is not an optimisation, it is what makes the row finish.
# The number of maximal cliques is exponential in the worst case and grows
# with DENSITY, not with edge count, so a flat edge cap does not bound it: a
# 60k-edge graph on 600 nodes is far worse than one on 6000. Keeping the
# MCC_DEG * n highest-weighted edges bounds mean degree at 2 * MCC_DEG
# whatever the module's size, which is the quantity that actually controls
# the enumeration. The rule is fixed in config and applied identically to
# every group, so it cannot favour one.
# igraph renamed subgraph.edges() to subgraph_from_edges() in 2.1.0. The
# cluster image and this laptop are not on the same igraph, and the old name
# warns once per module per group, which is 284 warnings of pure noise.
.sub_edges <- if (exists("subgraph_from_edges", where = asNamespace("igraph"),
                         inherits = FALSE)) {
  igraph::subgraph_from_edges
} else {
  getExportedValue("igraph", "subgraph.edges")
}

# Eigenvector centrality by power iteration. The adjacency is non-negative and
# connected here, so Perron-Frobenius gives a unique positive leading
# eigenvector and the iteration converges to it. Scaled so the maximum is 1,
# which is igraph's convention, so the values stay comparable with cent-040.
eig_leading <- function(A, iter = 1000L, tol = 1e-12) {
  n <- nrow(A)
  v <- rep(1 / sqrt(n), n)
  for (i in seq_len(iter)) {
    w  <- as.numeric(A %*% v)
    nw <- sqrt(sum(w * w))
    if (!is.finite(nw) || nw == 0) break
    w <- w / nw
    if (max(abs(w - v)) < tol) { v <- w; break }
    v <- w
  }
  v <- abs(v)
  mx <- max(v)
  setNames(if (mx > 0) v / mx else v, rownames(A))
}

sparsify <- function(g_sub, deg_budget) {
  keep <- as.integer(deg_budget) * vcount(g_sub)
  if (ecount(g_sub) <= keep) return(g_sub)
  w  <- igraph::E(g_sub)$weight
  th <- sort(w, decreasing = TRUE)[keep]
  .sub_edges(g_sub, igraph::E(g_sub)[w >= th], delete.vertices = FALSE)
}

mcc_log <- function(g_sub) {
  nv <- vcount(g_sub)
  nm <- V(g_sub)$name
  if (nv < 3L) return(setNames(rep(-Inf, nv), nm))
  cl <- tryCatch(max_cliques(g_sub, min = 2L), error = function(e) NULL)
  if (is.null(cl)) return(setNames(rep(NA_real_, nv), nm))
  acc <- rep(-Inf, nv)
  for (C in cl) {
    w   <- lgamma(length(C))              # log((|C| - 1)!)
    ix  <- as.integer(C)
    cur <- acc[ix]
    mx  <- pmax(cur, w)
    acc[ix] <- mx + log(exp(cur - mx) + exp(w - mx))
  }
  setNames(acc, nm)
}

# Everything needing the full adjacency happens inside here, so only one
# 1.2 GB matrix is alive at a time.
group_metrics <- function(Ex, label, pair_i, pair_j) {
  msg("  [", label, "] adjacency on ", nrow(Ex), " x ", ncol(Ex), " ...")
  A <- adjacency(Ex, power = POWER, type = NET_TYPE, corFnc = COR_FNC)
  diag(A) <- 0
  stopifnot(identical(colnames(A), genes))

  ## -- 7. centrality ----------------------------------------------------------
  # Computed from the matrix, not from an igraph object. The full weighted
  # graph on 12368 nodes has 76.5M edges and costs about 1.8 GB on top of the
  # 1.2 GB adjacency, which is what put earlier rows near the pod's memory
  # ceiling. Neither measure needs it: weighted degree IS the row sum, and
  # eigenvector centrality IS the leading eigenvector. Checked against
  # igraph::strength and igraph::eigen_centrality on a 1500-gene control
  # network: max absolute difference 1.6e-13 and 1.1e-11, Spearman 1.000,
  # and 13x faster.
  msg("  [", label, "] degree + eigenvector ...")
  deg <- rowSums(A)
  eig <- tryCatch(eig_leading(A),
                  error = function(e) setNames(rep(NA_real_, nrow(A)), rownames(A)))

  ## -- 2. intramodular connectivity -------------------------------------------
  msg("  [", label, "] intramodular connectivity ...")
  kIN <- setNames(rep(NA_real_, length(genes)), genes)
  for (m in mods) {
    gs <- names(mc)[mc == m]
    kIN[gs] <- rowSums(A[gs, gs, drop = FALSE])
  }

  ## -- 1. module eigengene correlation ----------------------------------------
  msg("  [", label, "] module eigengenes ...")
  ME  <- moduleEigengenes(Ex, colors = mc, excludeGrey = TRUE)$eigengenes
  Rme <- cor(ME, use = "pairwise.complete.obs")
  me_mean_abs <- mean(abs(Rme[upper.tri(Rme)]))
  kme <- as.data.frame(signedKME(Ex, ME))
  kme_own <- vapply(genes, function(gn) {
    m  <- mc[[gn]]
    cn <- paste0("kME", m)
    if (m == "grey" || !cn %in% colnames(kme)) NA_real_ else kme[gn, cn]
  }, numeric(1))

  ## -- 4. topological overlap --------------------------------------------------
  msg("  [", label, "] TOM ...")
  TOM <- TOMsimilarity(A, TOMType = NET_TYPE, verbose = 0)
  dimnames(TOM) <- list(genes, genes)
  tom_mod <- vapply(mods, function(m) {
    gs <- names(mc)[mc == m]
    S  <- TOM[gs, gs]
    mean(S[upper.tri(S)])
  }, numeric(1))
  tom_vec <- TOM[cbind(pair_i, pair_j)]   # the same entries in every group
  rm(TOM); gc(verbose = FALSE)

  ## -- 7b. betweenness and MCC on the thresholded graph ------------------------
  msg("  [", label, "] thresholding at ", BTW_THRESH, " ...")
  A[A <= BTW_THRESH] <- 0
  gT <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
  rm(A); gc(verbose = FALSE)
  msg("  [", label, "] sparse graph: ", vcount(gT), " nodes, ", ecount(gT), " edges")
  btw <- tryCatch(betweenness(gT, normalized = TRUE),
                  error = function(e) setNames(rep(NA_real_, vcount(gT)), V(gT)$name))

  msg("  [", label, "] log-MCC per module ...")
  mccv <- setNames(rep(NA_real_, length(genes)), genes)
  mcc_capped <- character(0)
  for (m in mods) {
    gs <- intersect(names(mc)[mc == m], V(gT)$name)
    if (length(gs) < 3L) next
    sg <- induced_subgraph(gT, gs)
    e0 <- ecount(sg)
    sg <- sparsify(sg, MCC_DEG)
    if (ecount(sg) < e0) mcc_capped <- c(mcc_capped, m)
    v <- mcc_log(sg)
    mccv[names(v)] <- v
  }
  if (length(mcc_capped))
    msg("  [", label, "] MCC sparsified in ", length(mcc_capped), " of ",
        length(mods), " modules (mean degree capped at ", 2L * MCC_DEG, ")")

  ## -- 6. propagation proximity ------------------------------------------------
  msg("  [", label, "] personalised PageRank from ", length(SEEDS), " seed sets ...")
  ppr <- lapply(SEEDS, function(sgen) {
    s <- intersect(sgen, V(gT)$name)
    if (!length(s)) return(setNames(rep(NA_real_, vcount(gT)), V(gT)$name))
    p <- rep(0, vcount(gT)); names(p) <- V(gT)$name; p[s] <- 1 / length(s)
    page_rank(gT, damping = PPR_DAMP, personalized = p,
              weights = igraph::E(gT)$weight)$vector
  })
  rm(gT); gc(verbose = FALSE)

  list(
    gene = data.frame(
      gene        = genes,
      module      = unname(mc[genes]),
      degree      = unname(deg[genes]),
      eigenvector = unname(eig[genes]),
      betweenness = unname(btw[genes]),
      mcc_log     = unname(mccv[genes]),
      kIN         = unname(kIN[genes]),
      kME_own     = unname(kme_own[genes]),
      ppr_BPD     = unname(ppr$BPD[genes]),
      ppr_LI      = unname(ppr$LI[genes]),
      row.names   = NULL),
    me_cor       = Rme,
    me_mean_abs  = me_mean_abs,
    tom_mod      = tom_mod,
    tom_vec      = tom_vec,
    n_mcc_capped = length(mcc_capped))
}

## ===========================================================================
## one draw
## ===========================================================================

do_draw <- function(k) {
  out <- file.path(run_dir, sprintf("draw_%02d.rds", k))
  if (file.exists(out)) { msg("draw ", k, " checkpoint exists; skipping"); return(invisible(NULL)) }
  msg("--- draw ", k, " of ", N_DRAWS, " ---")
  set.seed(SEED_BASE + k)

  cs  <- sample(ctrl_id, 2L * N_MATCH, replace = FALSE)
  sel <- list(ctrlA  = cs[seq_len(N_MATCH)],
              ctrlB  = cs[N_MATCH + seq_len(N_MATCH)],
              nolith = nol_id,
              lith   = sample(lit_id, N_MATCH, replace = FALSE))
  stopifnot(length(intersect(sel$ctrlA, sel$ctrlB)) == 0L,
            all(vapply(sel, length, 1L) == N_MATCH))

  # One pair set for the whole draw, so the cross-group TOM correlation
  # compares the same entries in every group.
  ng <- length(genes)
  pair_i <- sample.int(ng, TOM_PAIRS, replace = TRUE)
  pair_j <- sample.int(ng, TOM_PAIRS, replace = TRUE)
  ok <- pair_i != pair_j
  pair_i <- pair_i[ok]; pair_j <- pair_j[ok]

  M <- lapply(names(sel), function(lab)
    group_metrics(Emat(sel[[lab]]), lab, pair_i, pair_j))
  names(M) <- names(sel)

  ## -- 3. gene significance ----------------------------------------------------
  # GS needs the trait to vary, so it is computed on the union of the two
  # groups a contrast compares, not within a single group.
  gs_for <- function(a, b) {
    Ei <- Emat(c(sel[[a]], sel[[b]]))
    tr <- c(rep(0, length(sel[[a]])), rep(1, length(sel[[b]])))
    abs(as.numeric(cor(Ei, tr, use = "pairwise.complete.obs")))
  }
  GS <- data.frame(gene    = genes,
                   gs_bpd  = gs_for("ctrlA", "nolith"),
                   gs_lith = gs_for("nolith", "lith"),
                   gs_null = gs_for("ctrlA", "ctrlB"))

  saveRDS(list(draw = k, sel = sel, metrics = M, gs = GS), out)
  msg("draw ", k, " written: ", out)
}

if (SLICE_ID > 0L) {
  if (SLICE_ID > N_DRAWS) {
    msg("SLICE_ID ", SLICE_ID, " > N_DRAWS ", N_DRAWS, "; nothing to do")
  } else {
    do_draw(SLICE_ID)
  }
  msg("=== ", id, " slice ", SLICE_ID, " COMPLETE ===")
  quit(save = "no", status = 0L)
}

for (k in seq_len(N_DRAWS)) do_draw(k)
msg("=== ", id, " all draws COMPLETE (aggregate with netmet-061) ===")
