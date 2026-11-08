#!/usr/bin/env Rscript
# ===========================================================================
# Rebuild the base-001 reference network cache that seven experiment scripts
# depend on and that is NOT in the repository.
#
# network/cache/ is gitignored (*.rds), so a clone has no
#   ref_control_all_bicor_signed_p12_ds4_mch00_offsetrm.rds
# and perm-031, hub-036, hub-037, comm-038, comm-039, cent-040, tob-045 all
# quit(status = 1) on sight of its absence. On the iMac it exists because
# base-001 built it there in October.
#
# Parameters are copied from base-001.R verbatim: control group (all 234),
# bicor, signed, power 12, deepSplit 4, minModuleSize 30, mergeCutHeight 0,
# maxBlockSize 15000, seed 20261009, offset-removed expression.
#
# THE REBUILD IS VERIFIED, NOT ASSUMED. base-001's own
# reference_module_sizes.csv IS tracked, so the rebuilt module sizes are
# compared against it name by name and count by count. A mismatch stops the
# script, because every downstream row would otherwise be scored against a
# different decomposition than the one the arm's existing results used.
# ===========================================================================
suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors = FALSE)
enableWGCNAThreads(as.integer(Sys.getenv("PREP_THREADS", "4")))

source("network/00_load.R")
d <- load_project(cache = TRUE, verbose = TRUE)

CACHE_KEY <- "control_all_bicor_signed_p12_ds4_mch00_offsetrm"
OUT  <- file.path("network/cache", paste0("ref_", CACHE_KEY, ".rds"))
REF  <- "network/runs/base-001/reference_module_sizes.csv"
POWER <- 12L; NET_TYPE <- "signed"; COR_FNC <- "bicor"
DEEP_SPLIT <- 4L; MIN_MODULE <- 30L; MERGE_THRESH <- 0.0
REF_GROUP <- "control"; GENE_SET <- "all"; EXPR_MODE <- "offsetrm"
SEED_BASE <- 20261009L
msg <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                   paste0(..., collapse = ""))); flush.console() }

stopifnot(file.exists(REF))
expect <- read.csv(REF)
names(expect)[1:2] <- c("module", "n")
msg("target: ", nrow(expect), " entries from base-001/reference_module_sizes.csv")

msg("building reference on the ", REF_GROUP, " group ...")
tmp <- expr_for(d, group = REF_GROUP)
msg("  expression matrix: ", nrow(tmp), " samples x ", ncol(tmp), " genes")
stopifnot(nrow(tmp) == 234L, ncol(tmp) == 12368L)
set.seed(SEED_BASE)
t0 <- Sys.time()
bwm <- blockwiseModules(tmp, power = POWER, networkType = NET_TYPE,
                        TOMType = NET_TYPE, corType = COR_FNC,
                        deepSplit = DEEP_SPLIT, minModuleSize = MIN_MODULE,
                        mergeCutHeight = MERGE_THRESH, maxBlockSize = 15000L,
                        numericLabels = FALSE, saveTOMs = FALSE, verbose = 3)
moduleColors <- bwm$colors
msg("built in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")

got <- as.data.frame(table(moduleColors)); names(got) <- c("module", "n")
got$module <- as.character(got$module); expect$module <- as.character(expect$module)
m <- merge(expect, got, by = "module", all = TRUE, suffixes = c("_base001", "_rebuilt"))
m$n_base001[is.na(m$n_base001)] <- 0L; m$n_rebuilt[is.na(m$n_rebuilt)] <- 0L
m$delta <- m$n_rebuilt - m$n_base001
write.csv(m[order(-abs(m$delta)), ], "network/cache/rebuild_verification.csv", row.names = FALSE)

n_mods <- sum(names(table(moduleColors)) != "grey")
agree  <- sum(m$delta == 0); exact <- identical(nrow(m), agree)
msg("modules: ", n_mods, " non-grey, ", sum(moduleColors == "grey"), " grey")
msg("verification: ", agree, " of ", nrow(m), " module sizes match base-001 exactly")

if (!exact) {
  bad <- m[m$delta != 0, ]
  msg("MISMATCHES (worst 10):"); print(utils::head(bad[order(-abs(bad$delta)), ], 10))
  # A rebuild that disagrees is worse than no rebuild: every downstream row
  # would be scored against a decomposition the arm's results never used.
  stop("rebuilt reference does not reproduce base-001. Refusing to write the cache.")
}

ref_net <- list(expr_mode = EXPR_MODE, moduleColors = moduleColors,
                n_modules = n_mods, module_sizes = table(moduleColors),
                power = POWER, deep_split = DEEP_SPLIT,
                merge_cut_height = MERGE_THRESH, cache_key = CACHE_KEY,
                corFnc = COR_FNC, networkType = NET_TYPE,
                ref_group = REF_GROUP, gene_set = GENE_SET,
                source = "prep_cache.R rebuild, verified against base-001")
saveRDS(ref_net, OUT)
msg("wrote ", OUT, " (", round(file.size(OUT) / 1024), " KB)")
msg("done")
