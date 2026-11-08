#!/usr/bin/env Rscript
# ===========================================================================
# Rebuild the base-001 reference network cache that seven experiment scripts
# depend on and that is NOT in the repository.
#
# network/cache/ is gitignored (*.rds), so a clone has no
#   ref_control_all_bicor_signed_p12_ds4_mch00_offsetrm.rds
# and perm-031, hub-036, hub-037, comm-038, comm-039, cent-040 and tob-045 all
# quit(status = 1) on sight of its absence. On the iMac it exists because
# base-001 built it there in October.
#
# IT MUST BE REBUILT THE WAY base-001 BUILT IT, WHICH IS NOT blockwiseModules.
# base-001.R takes a fast path: it reads pre-computed cutreeDynamic labels from
# mod-000-resolution/labels_p12_ds4.rds and only falls back to blockwiseModules
# when they are absent. The two are NOT interchangeable — blockwiseModules adds
# a kME reassignment pass and leaves far more genes unassigned. A first attempt
# here used it and produced the right module COUNT (71) with the wrong
# membership: 568 grey genes against base-001's 8, and turquoise 822 against
# 1145. The verification below caught it and refused to write.
#
# So this reproduces mod-000-resolution's path exactly:
#   adjacency(power 12, signed, bicor) -> TOMsimilarity -> 1 - TOM
#   -> fastcluster::hclust(average) -> cutreeDynamic(deepSplit 4,
#      minClusterSize 30, method "hybrid") -> labels2colors
# with seed 20261009 throughout, and mergeCutHeight 0 means no merge follows.
#
# THE REBUILD IS VERIFIED, NOT ASSUMED. base-001's own
# reference_module_sizes.csv IS tracked, so every module size is compared name
# by name. A mismatch stops the script: a reference that silently disagreed
# would score every downstream row against a decomposition the arm's existing
# results never used, which is worse than having no cache at all.
# ===========================================================================
suppressPackageStartupMessages({ library(WGCNA); library(fastcluster) })
options(stringsAsFactors = FALSE)
enableWGCNAThreads(as.integer(Sys.getenv("PREP_THREADS", "4")))

source("network/00_load.R")
d <- load_project(cache = TRUE, verbose = TRUE)

CACHE_KEY <- "control_all_bicor_signed_p12_ds4_mch00_offsetrm"
OUT    <- file.path("network/cache", paste0("ref_", CACHE_KEY, ".rds"))
LABELS <- "network/runs/mod-000-resolution/labels_p12_ds4.rds"
REF    <- "network/runs/base-001/reference_module_sizes.csv"
POWER <- 12L; NET_TYPE <- "signed"; COR_FNC <- "bicor"
DEEP_SPLIT <- 4L; MIN_MODULE <- 30L; MERGE_THRESH <- 0.0
REF_GROUP <- "control"; GENE_SET <- "all"; EXPR_MODE <- "offsetrm"
SEED <- 20261009L
dir.create(dirname(LABELS), recursive = TRUE, showWarnings = FALSE)
dir.create("network/cache", recursive = TRUE, showWarnings = FALSE)
msg <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                   paste0(..., collapse = ""))); flush.console() }

stopifnot(file.exists(REF))
expect <- read.csv(REF); names(expect)[1:2] <- c("module", "n")
msg("target: ", nrow(expect), " entries from base-001/reference_module_sizes.csv")

tmp <- expr_for(d, group = REF_GROUP)
msg("expression: ", nrow(tmp), " samples x ", ncol(tmp), " genes")
stopifnot(nrow(tmp) == 234L, ncol(tmp) == 12368L)

if (file.exists(LABELS)) {
  msg("labels checkpoint present, reusing")
  raw_labels <- readRDS(LABELS)
} else {
  set.seed(SEED)
  msg("adjacency (power ", POWER, ", ", COR_FNC, ", ", NET_TYPE, ") ...")
  t0 <- Sys.time()
  adj <- adjacency(tmp, power = POWER, type = NET_TYPE, corFnc = COR_FNC)
  diag(adj) <- 0
  msg("  mean k = ", round(mean(colSums(adj)), 2),
      "  [", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min]")

  msg("TOM ...")
  t0 <- Sys.time()
  TOM <- TOMsimilarity(adj, TOMType = NET_TYPE); rm(adj); invisible(gc())
  dissTOM <- 1 - TOM; rm(TOM); invisible(gc())
  msg("  [", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min]")

  msg("dendrogram (average linkage) ...")
  t0 <- Sys.time()
  gene_tree <- fastcluster::hclust(as.dist(dissTOM), method = "average")
  invisible(gc())
  msg("  [", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min]")

  msg("cutreeDynamic deepSplit=", DEEP_SPLIT, " minClusterSize=", MIN_MODULE, " ...")
  set.seed(SEED)
  raw_labels <- cutreeDynamic(dendro = gene_tree, distM = dissTOM,
                              deepSplit = DEEP_SPLIT, minClusterSize = MIN_MODULE,
                              method = "hybrid", verbose = 2)
  saveRDS(raw_labels, LABELS)
  rm(dissTOM, gene_tree); invisible(gc())
  msg("  ", max(raw_labels), " modules pre-merge")
}
stopifnot(length(raw_labels) == 12368L)
moduleColors <- labels2colors(raw_labels)   # mergeCutHeight 0 => no merge step

got <- as.data.frame(table(moduleColors)); names(got) <- c("module", "n")
got$module <- as.character(got$module); expect$module <- as.character(expect$module)
m <- merge(expect, got, by = "module", all = TRUE, suffixes = c("_base001", "_rebuilt"))
m$n_base001[is.na(m$n_base001)] <- 0L; m$n_rebuilt[is.na(m$n_rebuilt)] <- 0L
m$delta <- m$n_rebuilt - m$n_base001
write.csv(m[order(-abs(m$delta)), ], "network/cache/rebuild_verification.csv", row.names = FALSE)

n_mods <- sum(names(table(moduleColors)) != "grey")
agree  <- sum(m$delta == 0)
msg("modules: ", n_mods, " non-grey, ", sum(moduleColors == "grey"), " grey")
msg("verification: ", agree, " of ", nrow(m), " module sizes match base-001 exactly")

if (agree != nrow(m)) {
  bad <- m[m$delta != 0, ]
  msg("MISMATCHES (worst 10):"); print(utils::head(bad[order(-abs(bad$delta)), ], 10))
  stop("rebuilt reference does not reproduce base-001. Refusing to write the cache.")
}

ref_net <- list(expr_mode = EXPR_MODE, moduleColors = moduleColors,
                n_modules = n_mods, module_sizes = table(moduleColors),
                power = POWER, deep_split = DEEP_SPLIT,
                merge_cut_height = MERGE_THRESH, cache_key = CACHE_KEY,
                corFnc = COR_FNC, networkType = NET_TYPE,
                ref_group = REF_GROUP, gene_set = GENE_SET,
                source = "prep_cache.R via mod-000 cutreeDynamic path, verified vs base-001")
saveRDS(ref_net, OUT)
msg("wrote ", OUT, " (", round(file.size(OUT) / 1024), " KB)")
msg("done")
