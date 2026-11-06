#!/usr/bin/env Rscript
# ===========================================================================
# SUPERSEDED by grp-051.R, which runs 8 draws on the cluster instead of 1,
# adds GO enrichment, and fixes the 1-D table() bug below. Kept because its
# four networks reproduce grp-051's first draw exactly (69/64/83/78 modules),
# which is the cross-check that the cluster run is deterministic.
#
# grp-050 — group-specific networks: which modules and hub genes does each
#           group actually have, and how much do they overlap?
#
# Everything in the arm so far builds ONE reference and measures preservation
# of those modules elsewhere. That cannot answer "which modules do lithium
# users have", because the module set is never rebuilt per group. This does.
#
# THE SAMPLE-SIZE TRAP. The existing references are n = 234 (control), 74
# (bp_nolith) and 152 (bp_lith), and module count moves with n for reasons that
# have nothing to do with biology. Every network here is built on exactly
# n = 74 with identical parameters, so counts are comparable.
#
# THE CEILING. Two networks built on different 74-person samples of the SAME
# population also disagree, because module detection is unstable. So
# control_A vs control_B (disjoint halves of the controls) is computed first:
# that is how much disagreement means nothing. Any cross-group comparison is
# read against it, exactly as base-002-ceil does for preservation.
#
# Parameters copied from ref-003.R so results sit on the same footing as the
# rest of the arm: bicor, signed, power 12, deepSplit 4, minModuleSize 30,
# mergeCutHeight 0, offset-removed expression.
#
# Outputs (network/runs/grp-050/):
#   genes_<net>.csv        gene, module, kME        <- the file no run saved
#   hubs_<net>.csv         top 10 hub genes per module
#   module_sizes_<net>.csv module, n_genes
#   pair_<A>_vs_<B>.csv    per-module best match and Jaccard across two nets
#   summary.csv            one row per network pair: ARI, matched fraction
# ===========================================================================
suppressPackageStartupMessages(library(WGCNA))  # flashClust not installed; blockwiseModules clusters internally
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

source("network/00_load.R")
d <- load_project(cache = TRUE, verbose = FALSE)

OUT <- "network/runs/grp-050"; CK <- file.path(OUT, "ck")
dir.create(CK, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                   paste0(..., collapse = ""))); flush.console() }

POWER <- 12L; NET <- "signed"; DEEP <- 4L; MINMOD <- 30L; MERGE <- 0.0
N_PER_GROUP <- 74L
SEED <- 20261105L

# -- the four sample sets ---------------------------------------------------
set.seed(SEED)
ids   <- colnames(d$logtpm)
ctrl  <- ids[which(d$groups == "control")]     # 234
nolit <- ids[which(d$groups == "bp_nolith")]   # 74  — never subsampled
lith  <- ids[which(d$groups == "bp_lith")]     # 152
stopifnot(length(nolit) == N_PER_GROUP)

cs <- sample(ctrl, 2 * N_PER_GROUP)            # 148 controls, split disjointly
SETS <- list(control_A = cs[1:N_PER_GROUP],
             control_B = cs[(N_PER_GROUP + 1):(2 * N_PER_GROUP)],
             bp_nolith = nolit,
             bp_lith   = sample(lith, N_PER_GROUP))
stopifnot(length(intersect(SETS$control_A, SETS$control_B)) == 0)

# -- build one network, cached ----------------------------------------------
build <- function(nm) {
  f <- file.path(CK, paste0(nm, ".rds"))
  if (file.exists(f)) { msg(nm, ": checkpoint hit"); return(readRDS(f)) }
  msg(nm, ": building on ", length(SETS[[nm]]), " samples ...")
  E <- expr_for(d, group = NULL)[SETS[[nm]], , drop = FALSE]   # samples x genes
  stopifnot(nrow(E) == N_PER_GROUP, ncol(E) == nrow(d$logtpm))
  t0 <- Sys.time()
  bwm <- blockwiseModules(E, power = POWER, networkType = NET, TOMType = NET,
                          corType = "bicor", deepSplit = DEEP,
                          minModuleSize = MINMOD, mergeCutHeight = MERGE,
                          maxBlockSize = 15000L, numericLabels = FALSE,
                          saveTOMs = FALSE, verbose = 0)
  # kME: each gene's correlation with its own module eigengene = how central
  # it is. The top few by kME are what WGCNA calls the module's hub genes.
  kme <- signedKME(E, bwm$MEs, corFnc = "bicor")
  own <- paste0("kME", bwm$colors)
  kme_own <- vapply(seq_along(bwm$colors), function(i)
    if (own[i] %in% colnames(kme)) kme[i, own[i]] else NA_real_, numeric(1))
  res <- list(colors = bwm$colors, MEs = bwm$MEs, genes = colnames(E),
              kME = kme_own, secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  saveRDS(res, f)
  msg(nm, ": ", length(unique(res$colors)) - ("grey" %in% res$colors),
      " modules in ", round(res$secs), "s")
  rm(E, bwm, kme); gc(verbose = FALSE)
  res
}

NETS <- lapply(names(SETS), build); names(NETS) <- names(SETS)

# -- per-network tables -----------------------------------------------------
sym <- d$sym
for (nm in names(NETS)) {
  n <- NETS[[nm]]
  g <- data.frame(gene = n$genes,
                  symbol = ifelse(n$genes %in% names(sym), sym[n$genes], n$genes),
                  module = n$colors, kME = round(n$kME, 4))
  write.csv(g, file.path(OUT, paste0("genes_", nm, ".csv")), row.names = FALSE)
  sz <- as.data.frame(table(g$module[g$module != "grey"]))
  names(sz) <- c("module", "n_genes"); sz <- sz[order(-sz$n_genes), ]
  write.csv(sz, file.path(OUT, paste0("module_sizes_", nm, ".csv")), row.names = FALSE)
  hub <- do.call(rbind, lapply(split(g[g$module != "grey", ], g$module[g$module != "grey"]),
    function(x) { x <- x[order(-x$kME), ]; head(x, 10) }))
  write.csv(hub, file.path(OUT, paste0("hubs_", nm, ".csv")), row.names = FALSE)
}

# -- compare two partitions -------------------------------------------------
# Adjusted Rand index on the non-grey genes of both: 1 = identical grouping,
# 0 = no better than chance. Plus, per module of A, the module of B it overlaps
# most and the Jaccard of that match.
ari <- function(a, b) {
  tb <- table(a, b); n <- sum(tb)
  s <- function(x) sum(choose(x, 2))
  i <- s(tb); ea <- s(rowSums(tb)); eb <- s(colSums(tb)); tot <- choose(n, 2)
  exp <- ea * eb / tot
  (i - exp) / ((ea + eb) / 2 - exp)
}

pair <- function(A, B) {
  a <- NETS[[A]]; b <- NETS[[B]]
  stopifnot(identical(a$genes, b$genes))
  keep <- a$colors != "grey" & b$colors != "grey"
  rand <- ari(a$colors[keep], b$colors[keep])
  tb <- table(a$colors, b$colors)
  tb <- tb[rownames(tb) != "grey", colnames(tb) != "grey", drop = FALSE]
  sa <- table(a$colors)[rownames(tb)]   # table() of one factor is 1-D:
  sb <- table(b$colors)[colnames(tb)]   # rowSums/colSums error on it
  rows <- lapply(rownames(tb), function(m) {
    ov <- tb[m, ]; j <- ov / (sa[m] + sb - ov)
    k <- which.max(j)
    data.frame(module_A = m, n_A = as.integer(sa[m]),
               best_match_B = colnames(tb)[k], n_B = as.integer(sb[k]),
               n_shared = as.integer(ov[k]), jaccard = round(max(j), 3))
  })
  df <- do.call(rbind, rows); df <- df[order(-df$jaccard), ]
  write.csv(df, file.path(OUT, sprintf("pair_%s_vs_%s.csv", A, B)), row.names = FALSE)
  data.frame(comparison = sprintf("%s vs %s", A, B),
             n_modules_A = nrow(tb), n_modules_B = ncol(tb),
             ARI = round(rand, 4),
             median_jaccard = round(median(df$jaccard), 3),
             frac_matched_j50 = round(mean(df$jaccard >= 0.5), 3),
             frac_matched_j25 = round(mean(df$jaccard >= 0.25), 3))
}

PAIRS <- list(c("control_A", "control_B"),   # the ceiling
              c("control_A", "bp_nolith"),   # bipolar, unmedicated
              c("bp_nolith", "bp_lith"),     # lithium
              c("control_A", "bp_lith"))
S <- do.call(rbind, lapply(PAIRS, function(p) pair(p[1], p[2])))
write.csv(S, file.path(OUT, "summary.csv"), row.names = FALSE)

# -- hub stability: are the same genes central in each group? ---------------
hub_top <- function(nm, k = 200) {
  n <- NETS[[nm]]; ok <- n$colors != "grey"
  g <- n$genes[ok]; v <- n$kME[ok]
  g[order(-v)][1:min(k, length(g))]
}
H <- lapply(names(NETS), hub_top); names(H) <- names(NETS)
hs <- do.call(rbind, lapply(PAIRS, function(p) data.frame(
  comparison = sprintf("%s vs %s", p[1], p[2]),
  shared_top200_hubs = length(intersect(H[[p[1]]], H[[p[2]]])),
  pct = round(100 * length(intersect(H[[p[1]]], H[[p[2]]])) / 200, 1))))
write.csv(hs, file.path(OUT, "hub_overlap.csv"), row.names = FALSE)

msg("=== module counts (all at n = 74) ===")
for (nm in names(NETS))
  cat(sprintf("  %-10s %3d modules, %5d genes assigned, %5d grey\n", nm,
              length(unique(NETS[[nm]]$colors)) - 1,
              sum(NETS[[nm]]$colors != "grey"), sum(NETS[[nm]]$colors == "grey")))
msg("=== partition agreement ==="); print(S, row.names = FALSE)
msg("=== top-200 hub overlap ==="); print(hs, row.names = FALSE)
msg("done")
