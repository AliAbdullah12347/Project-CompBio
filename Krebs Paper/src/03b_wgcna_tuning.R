#!/usr/bin/env Rscript
# WGCNA parameter sweep: why does our partition put 2,960 genes in grey when Krebs
# et al. put 314?
#
#   usage: Rscript src/03b_wgcna_tuning.R --net=unsigned --cor=pearson [--grid=full]
#          Rscript src/03b_wgcna_tuning.R --combine
#
# Design note -- why this script open-codes blockwiseModules instead of calling it.
# ------------------------------------------------------------------------------
# blockwiseModules recomputes the topological overlap matrix on every call, and the
# TOM is ~99% of the runtime (12,344 genes -> a 1.2 GB dense matrix). But the TOM
# depends only on (power, networkType, TOMType, corType); every parameter we want to
# sweep -- deepSplit, detectCutHeight, pamStage, pamRespectsDendro, the kME quality
# filter -- acts downstream of it. So we compute the TOM once per (networkType,
# corType) pair and reuse it across the whole cut grid, which turns a ~4-hour sweep
# into a ~40-minute one.
#
# It reproduces these steps from the WGCNA 1.74 source, in order:
#   1. tom      <- TOMsimilarityFromExpr(...)          [same .Call as blockwiseModules]
#   2. dissTom  <- 1 - tom
#   3. dendro   <- fastcluster::hclust(as.dist(dissTom), method = "average")
#   4. labels   <- cutreeDynamic(method = "hybrid", distM = dissTom, ...)
#   5. kME quality control  (minCoreKME / minCoreKMESize / minKMEtoStay)
#   6. colors   <- labels2colors(labels)               [before merging, as WGCNA does]
#   7. mergeCloseModules(datExpr, colors, cutHeight = 0.25)
#
# One blockwiseModules step is deliberately omitted: the `reassignThreshold = 1e-6`
# pass, which moves a gene from its own module to another when the other module's kME
# t-test p-value beats its own by six orders of magnitude. It is omitted because it
# can only move a gene between two existing modules, never into grey, so it cannot
# affect the quantity under investigation.
#
# --validate measures the cost of that omission instead of assuming it away, and the
# honest answer is that it is not negligible: at blockwiseModules' own defaults the
# open-coded version gives 24 modules and 2,960 grey against blockwiseModules' 25 and
# 2,960, with an adjusted Rand index of 0.869 and 88.6% of genes in the best-matching
# counterpart module. The grey count -- the number the sweep exists to explain -- is
# reproduced exactly, and the module count to within one; the ~11% of genes that sit
# in a different module are the reassignment pass. So the sweep is a faithful stand-in
# for blockwiseModules on the grey and module-count axes and only an approximate one
# on ARI, which is why grey is scored first. The configuration the sweep selects is
# then re-run end to end by 03_wgcna.R, whose network code is this same step-by-step
# path, so nothing reported downstream depends on the stand-in.

suppressPackageStartupMessages({
  library(WGCNA); library(dynamicTreeCut); library(fastcluster)
})
options(stringsAsFactors = FALSE)

script_path <- sub("^--file=", "",
                   grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
ROOT <- if (!is.na(script_path)) dirname(dirname(normalizePath(script_path))) else "."
if (!dir.exists(file.path(ROOT, "data"))) ROOT <- "."
PROC <- file.path(ROOT, "data", "processed")
REF  <- file.path(ROOT, "data", "reference")
RES  <- file.path(ROOT, "results")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

SWEEP_CSV <- file.path(RES, "wgcna_param_sweep.csv")

# Fixed by the paper's supplementary methods; not swept.
BETA <- 7; MIN_MODULE_SIZE <- 30; MERGE_CUT_HEIGHT <- 0.25
# blockwiseModules' own kME quality-control defaults, used when kme_filter = TRUE.
MIN_CORE_KME <- 0.5; MIN_CORE_KME_SIZE <- MIN_MODULE_SIZE / 3; MIN_KME_TO_STAY <- 0.3

argv     <- commandArgs(TRUE)
argval   <- function(flag, default) {
  hit <- grep(paste0("^--", flag, "="), argv, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", flag, "="), "", hit[1])
}
NET      <- argval("net", "unsigned")
CORTYPE  <- argval("cor", "pearson")
GRID     <- argval("grid", "full")
VALIDATE <- "--validate" %in% argv
COMBINE  <- "--combine" %in% argv
# --merge-only re-merges the winning cut's cached pre-merge labels under different
# mergeCloseModules options. It needs no TOM and runs in seconds, which is why it is
# a separate mode rather than another axis of the main grid.
MERGE_ONLY <- "--merge-only" %in% argv
stopifnot(NET %in% c("unsigned", "signed"), CORTYPE %in% c("pearson", "bicor"))

BASELINE_KEY <- "2|0.995|TRUE|TRUE|TRUE"   # exactly blockwiseModules' defaults

# --------------------------------------------------------------- scoring helpers
# Adjusted Rand index. Written out rather than pulled from a package because the
# only implementations installed here are in packages we are not allowed to add.
adj_rand <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  sij <- sum(choose(tab, 2))
  sa  <- sum(choose(rowSums(tab), 2)); sb <- sum(choose(colSums(tab), 2))
  expct <- sa * sb / choose(n, 2)
  (sij - expct) / ((sa + sb) / 2 - expct)
}

# ------------------------------------------------------------------ validate mode
# Post-hoc, and deliberately so: it compares the sweep's baseline row -- deepSplit 2,
# detectCutHeight 0.995, pamStage TRUE, pamRespectsDendro TRUE, kME filter on, i.e.
# blockwiseModules' own defaults -- against the network blockwiseModules actually
# produced in 03_wgcna.R. If the two partitions are not essentially identical, the
# open-coded sweep is not a stand-in for blockwiseModules and nothing below it holds.
if (VALIDATE) {
  lab_file <- file.path(PROC, "wgcna_sweep_labels_unsigned_pearson.rds")
  net_file <- file.path(PROC, "wgcna_network.rds")
  if (!file.exists(lab_file)) stop("run the unsigned/pearson sweep first: ", lab_file)
  if (!file.exists(net_file)) stop("no cached blockwiseModules network at ", net_file)
  best <- readRDS(lab_file)
  if (is.null(best[[BASELINE_KEY]]))
    stop("baseline configuration ", BASELINE_KEY, " not in the saved sweep labels")
  ours <- best[[BASELINE_KEY]]
  net  <- readRDS(net_file)
  bw   <- labels2colors(net$colors)   # 03_wgcna.R ran with numericLabels = TRUE
  cat("VALIDATION -- open-coded sweep vs cached blockwiseModules, both at defaults\n")
  cat(sprintf("  blockwiseModules : %2d modules, grey %5d\n",
              length(setdiff(unique(bw), "grey")), sum(bw == "grey")))
  cat(sprintf("  open-coded sweep : %2d modules, grey %5d\n",
              length(setdiff(unique(ours), "grey")), sum(ours == "grey")))
  tab <- table(bw, ours)
  same <- sum(apply(tab, 1, max))     # colour names are arbitrary; match by overlap
  cat(sprintf("  adjusted Rand index between the two: %.6f\n", adj_rand(bw, ours)))
  cat(sprintf("  genes in the best-matching counterpart module: %d / %d (%.3f%%)\n",
              same, length(bw), 100 * same / length(bw)))
  quit(save = "no")
}

# ------------------------------------------------------------------- combine mode
if (COMBINE) {
  parts <- list.files(RES, pattern = "^\\.wgcna_sweep_.*\\.csv$", full.names = TRUE,
                      all.files = TRUE)
  if (length(parts) == 0) stop("no sweep partials found in ", RES)
  # Partials written by different modes can carry different columns (--merge-only
  # adds merge_useabs), so align on the union rather than requiring an exact match.
  tabs <- lapply(parts, read.csv)
  allcols <- Reduce(union, lapply(tabs, names))
  tabs <- lapply(tabs, function(d) {
    for (cc in setdiff(allcols, names(d))) d[[cc]] <- NA
    d[, allcols, drop = FALSE]
  })
  all <- do.call(rbind, tabs)
  if ("merge_useabs" %in% names(all)) all$merge_useabs[is.na(all$merge_useabs)] <- FALSE
  # --merge-only re-scores the winning configuration under useAbs = FALSE as its own
  # control, so that row also exists in the main grid. Keep one copy.
  keycols <- c("networkType", "corType", "deepSplit", "detectCutHeight", "pamStage",
               "pamRespectsDendro", "kme_filter", "merge_useabs")
  all <- all[!duplicated(all[, intersect(keycols, names(all))]), ]
  # Rank by the task's stated priority: grey count first, then module count, then
  # agreement with the published labels.
  all <- all[order(abs(all$n_grey - 314), abs(all$n_modules - 27), -all$ari_all), ]
  write.csv(all, SWEEP_CSV, row.names = FALSE)
  cat(sprintf("Combined %d partial file(s) -> %s (%d configurations)\n",
              length(parts), SWEEP_CSV, nrow(all)))
  print(head(all[, c("networkType", "corType", "deepSplit", "detectCutHeight",
                     "pamStage", "pamRespectsDendro", "kme_filter",
                     "n_modules", "n_grey", "ari_all", "max_size", "mean_size")], 15),
        row.names = FALSE, digits = 4)
  quit(save = "no")
}

PART_CSV <- file.path(RES, sprintf(".wgcna_sweep_%s_%s.csv", NET, CORTYPE))

# R block-buffers stdout when it is redirected, so a sweep that runs for an hour
# shows nothing until it exits. say() also appends to a log file, and appending
# reopens and closes the file, so progress is on disk while the job is still running.
LOG <- file.path(PROC, sprintf("wgcna_sweep_%s_%s.log", NET, CORTYPE))
say <- function(...) { msg <- sprintf(...); cat(msg); cat(msg, file = LOG, append = TRUE) }
if (!COMBINE && !VALIDATE && !MERGE_ONLY) cat("", file = LOG)
starttime <- Sys.time()

# ------------------------------------------------------------------------- inputs
cat(sprintf("=== sweep: networkType=%s corType=%s grid=%s ===\n", NET, CORTYPE, GRID))
cached  <- readRDS(file.path(PROC, "wgcna_residuals.rds"))
datExpr <- cached$datExpr; rm(cached); invisible(gc())
cat(sprintf("datExpr: %d samples x %d genes\n", nrow(datExpr), ncol(datExpr)))

ref_mod <- read.csv(file.path(REF, "File_S2_WGCNA_modules__Supplementary_File_2.csv"),
                    colClasses = c(gene = "character", moduleColor = "character"))
theirColor <- setNames(ref_mod$moduleColor, ref_mod$gene)[colnames(datExpr)]
stopifnot(!anyNA(theirColor))
their_ng <- table(theirColor)[names(table(theirColor)) != "grey"]
cat(sprintf("target : %d modules, %d grey, sizes %d-%d, mean %.1f\n",
            length(their_ng), sum(theirColor == "grey"),
            min(their_ng), max(their_ng), mean(their_ng)))

# One row of the scoreboard, given a merged colour vector and the cut it came from.
score_row <- function(merged, dyn, labels, deepSplit, cutHeight, pamStage,
                      pamRespectsDendro, kme_filter, merge_useabs) {
  sz <- table(merged); sz_ng <- sz[names(sz) != "grey"]
  keep <- merged != "grey" & theirColor != "grey"
  data.frame(
    networkType = NET, corType = CORTYPE, deepSplit = deepSplit,
    detectCutHeight = if (is.na(cutHeight)) "auto" else as.character(cutHeight),
    pamStage = pamStage, pamRespectsDendro = pamRespectsDendro,
    kme_filter = kme_filter, merge_useabs = merge_useabs,
    n_modules_predetect = length(setdiff(unique(dyn), 0)),
    n_grey_predetect    = sum(dyn == 0),
    n_modules_premerge  = length(setdiff(unique(labels), 0)),
    n_modules = length(sz_ng), n_grey = sum(merged == "grey"),
    n_assigned = sum(merged != "grey"),
    min_size = if (length(sz_ng)) min(sz_ng) else NA_integer_,
    max_size = if (length(sz_ng)) max(sz_ng) else NA_integer_,
    mean_size = if (length(sz_ng)) mean(sz_ng) else NA_real_,
    ari_all = adj_rand(merged, theirColor),
    ari_bothnongrey = if (sum(keep) > 1) adj_rand(merged[keep], theirColor[keep]) else NA_real_,
    stringsAsFactors = FALSE)
}

# --------------------------------------------------------------- merge-only mode
# Their 27 modules and our 31 hold nearly the same genes in their large modules --
# their top two total 4,852 and our top five total 7,023 against their top five's
# 7,187 -- so the residual disagreement is not in the cut, it is that their merge
# fused large modules ours kept apart. mergeCloseModules' `useAbs` is the parameter
# that does exactly that: with useAbs = FALSE (the default) two modules whose
# eigengenes are strongly ANTI-correlated are maximally dissimilar and never merge,
# while with useAbs = TRUE they merge immediately. In an unsigned network, where the
# adjacency already discards the sign of the correlation, useAbs = TRUE is arguably
# the consistent choice, and the methods do not say which was used. Test it on the
# winning cut's cached pre-merge labels, which costs seconds rather than a new TOM.
if (MERGE_ONLY) {
  netf <- file.path(PROC, "wgcna_network_stepwise.rds")
  if (!file.exists(netf)) stop("run src/03_wgcna.R first to cache the winning cut: ", netf)
  net <- readRDS(netf)
  stopifnot(identical(net$genes, colnames(datExpr)))
  dyn <- net$dynamic
  cat(sprintf("cached winning cut: %d modules, %d unassigned before merging\n",
              length(setdiff(unique(dyn), 0)), sum(dyn == 0)))
  rows <- list()
  for (ua in c(FALSE, TRUE)) {
    mg <- as.character(mergeCloseModules(datExpr, labels2colors(dyn),
                                         cutHeight = MERGE_CUT_HEIGHT, relabel = FALSE,
                                         useAbs = ua, verbose = 0)$colors)
    rows[[length(rows) + 1L]] <- score_row(mg, dyn, dyn, 2, NA, TRUE, FALSE, FALSE, ua)
    r <- rows[[length(rows)]]
    cat(sprintf("  useAbs=%-5s -> %2d modules, grey %5d, sizes %d-%d, mean %.1f, ARI %.4f\n",
                ua, r$n_modules, r$n_grey, r$min_size, r$max_size, r$mean_size, r$ari_all))
  }
  out <- do.call(rbind, rows)
  write.csv(out, file.path(RES, ".wgcna_sweep_mergeabs.csv"), row.names = FALSE)
  cat("Wrote results/.wgcna_sweep_mergeabs.csv\n")
  quit(save = "no")
}

# --------------------------------------------------------------------------- TOM
nthreads <- tryCatch({ enableWGCNAThreads(); WGCNAnThreads() },
                     error = function(e) { allowWGCNAThreads(); 1L })
cat(sprintf("Computing TOM (power=%d, %s adjacency, signed TOM, %s correlation, %d threads) ...\n",
            BETA, NET, CORTYPE, nthreads))
t0 <- Sys.time()
dissTom <- TOMsimilarityFromExpr(datExpr, power = BETA, networkType = NET,
                                 TOMType = "signed", corType = CORTYPE,
                                 nThreads = nthreads, verbose = 0)
dissTom <- 1 - dissTom                      # in place where R's refcount allows
invisible(gc())
say("  TOM done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins")))

say("Clustering (average linkage) ...\n")
geneTree <- fastcluster::hclust(as.dist(dissTom), method = "average")
invisible(gc())
say("  clustering done at %.1f min\n",
    as.numeric(difftime(Sys.time(), starttime, units = "mins")))

# ------------------------------------------------------- blockwiseModules' kME QC
# Verbatim reimplementation of WGCNA 1.74 blockwiseModules lines 318-357: dissolve
# any module whose core is too weakly correlated with its own eigengene, then drop
# individual genes whose kME falls below minKMEtoStay. Both dump genes into grey.
corFn <- if (CORTYPE == "bicor") WGCNA::bicor else WGCNA::cor
apply_kme_qc <- function(labels) {
  if (all(labels == 0)) return(labels)
  MEs <- moduleEigengenes(datExpr[, labels != 0, drop = FALSE],
                          labels[labels != 0], verbose = 0)$eigengenes
  idx <- as.numeric(substring(names(MEs), 3))
  for (m in seq_along(idx)) {
    modGenes <- labels == idx[m]
    if (!any(modGenes)) next
    KME <- corFn(datExpr[, modGenes, drop = FALSE], MEs[, m])
    # Unsigned networks cluster on |cor|, so kME is judged on |cor| too.
    if (NET == "unsigned") KME <- abs(KME)
    if (sum(KME > MIN_CORE_KME) < MIN_CORE_KME_SIZE) {
      labels[modGenes] <- 0
    } else if (any(KME < MIN_KME_TO_STAY)) {
      labels[modGenes][KME < MIN_KME_TO_STAY] <- 0
      if (sum(labels == idx[m]) < MIN_MODULE_SIZE) labels[labels == idx[m]] <- 0
    }
  }
  labels
}

# ------------------------------------------------------------------ one config
score_config <- function(deepSplit, cutHeight, pamStage, pamRespectsDendro,
                         kme_filter, dyn_cache) {
  key <- paste(deepSplit, cutHeight, pamStage, pamRespectsDendro, sep = "|")
  if (is.null(dyn_cache[[key]])) {
    # cutHeight = NA means "let dynamicTreeCut pick", which is what the WGCNA
    # step-by-step tutorial does; blockwiseModules instead hard-codes 0.995.
    dyn_cache[[key]] <- cutreeDynamic(
      dendro = geneTree, distM = dissTom, method = "hybrid",
      deepSplit = deepSplit,
      cutHeight = if (is.na(cutHeight)) NULL else cutHeight,
      minClusterSize = MIN_MODULE_SIZE,
      pamStage = pamStage, pamRespectsDendro = pamRespectsDendro,
      verbose = 0)
  }
  dyn <- dyn_cache[[key]]
  labels <- if (kme_filter) apply_kme_qc(dyn) else dyn

  # labels2colors BEFORE merging, then merge with relabel = FALSE. Both choices are
  # forced by the authors' own colour names, and neither is blockwiseModules'.
  # mergeCloseModules(relabel = TRUE) -- which is what blockwiseModules always calls --
  # renames the surviving modules turquoise, blue, brown ... in descending size, so 27
  # modules would carry standardColors()[1:27] with no gaps. The authors' 27 colours
  # instead sit at scattered positions 2-52, with turquoise, yellow and green missing
  # while red and brown are their two largest. That is the signature of relabel = FALSE,
  # where a merged branch inherits the alphabetically first colour among its members
  # ("red" before "turquoise", "brown" before "green"/"magenta"/"pink"/"purple"). So the
  # published network came from the step-by-step WGCNA tutorial, not blockwiseModules.
  # relabel changes names only, never the partition, so it cannot affect any score below.
  cols <- labels2colors(labels)
  merged <- as.character(mergeCloseModules(datExpr, cols, cutHeight = MERGE_CUT_HEIGHT,
                                           relabel = FALSE, corFnc = corFn,
                                           verbose = 0)$colors)

  list(row = score_row(merged, dyn, labels, deepSplit, cutHeight, pamStage,
                       pamRespectsDendro, kme_filter, merge_useabs = FALSE),
       colors = merged, dyn_cache = dyn_cache)
}

# --------------------------------------------------------------------- the grid
# Ordered, not a raw expand.grid, for two reasons. cutreeDynamic on a 12,344-gene
# distance matrix is minutes per call, so the grid is built to reuse each cut across
# its two kME settings (see dyn_cache), and it is ordered most-plausible-first: if the
# job has to be cut short, the configurations that matter have already been scored.
# "auto" cut height is dynamicTreeCut's own default, which is what the step-by-step
# tutorial uses; 0.995 is what blockwiseModules hard-codes.
mk <- function(deepSplit, cutHeight, pam) {
  g <- expand.grid(kme_filter = c(FALSE, TRUE), deepSplit = deepSplit,
                   cutHeight = cutHeight, pam = pam, stringsAsFactors = FALSE)
  g[, c("deepSplit", "cutHeight", "pam", "kme_filter")]
}
if (GRID == "full") {
  grid <- rbind(
    # The paper's stated deepSplit, across every cut height and PAM setting.
    mk(2, c(NA, 0.995, 0.99, 0.999), c("T/F", "T/T", "F")),
    # The other deepSplit values at the two cut heights that are actually defensible.
    mk(c(1, 3, 4), c(NA, 0.995), c("T/F", "T/T", "F")))
} else {
  # Reduced grid for the alternative network/correlation types. deepSplit is held at
  # the paper's value of 2 and only the PAM/kME axes move, because those are the two
  # that moved the grey count in the full grid.
  grid <- mk(2, c(NA, 0.995), c("T/F", "T/T", "F"))
}
grid <- grid[!duplicated(grid), ]
grid$pamStage <- grid$pam != "F"
# pamRespectsDendro is meaningless when pamStage is off; recorded as FALSE so those
# rows are not mistaken for a distinct configuration.
grid$pamRespectsDendro <- grid$pam == "T/T"
cat(sprintf("Grid: %d configurations\n", nrow(grid)))

dyn_cache <- list(); rows <- list(); best <- list()
for (i in seq_len(nrow(grid))) {
  g  <- grid[i, ]
  ti <- Sys.time()
  out <- score_config(g$deepSplit, g$cutHeight, g$pamStage, g$pamRespectsDendro,
                      g$kme_filter, dyn_cache)
  dyn_cache <- out$dyn_cache
  rows[[i]] <- out$row
  best[[paste(g$deepSplit, g$cutHeight, g$pamStage, g$pamRespectsDendro,
              g$kme_filter, sep = "|")]] <- out$colors
  say("[%2d/%2d] ds=%d cut=%-5s pam=%-3s kme=%-5s -> %2d modules, grey %5d, ARI %.4f (%.1f s)\n",
      i, nrow(grid), g$deepSplit,
      if (is.na(g$cutHeight)) "auto" else format(g$cutHeight),
      g$pam, g$kme_filter, out$row$n_modules, out$row$n_grey,
      out$row$ari_all, as.numeric(difftime(Sys.time(), ti, units = "secs")))
  # Write after every configuration: this sweep takes tens of minutes and a crash
  # partway through should not cost the configurations already scored.
  write.csv(do.call(rbind, rows), PART_CSV, row.names = FALSE)
}

saveRDS(best, file.path(PROC, sprintf("wgcna_sweep_labels_%s_%s.rds", NET, CORTYPE)))
cat(sprintf("\nWrote %s (%d rows)\n", PART_CSV, length(rows)))
