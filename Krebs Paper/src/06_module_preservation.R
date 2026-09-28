#!/usr/bin/env Rscript
# Module preservation analysis, reproducing Krebs et al. 2020 Supplementary Methods
# ("Weighted gene co-expression network module preservation analysis", Figure S6).
#
# Their specification, verbatim in substance:
#   - build WGCNA networks in four groups separately, each with its own soft power
#       BD cases beta = 4.2 | controls beta = 5.5 | lithium users beta = 7 |
#       non-lithium users beta = 7
#   - WGCNA::modulePreservation with 200 permutations, on four ref -> test pairs
#       cases -> controls, controls -> cases, lithium -> non-lithium, non-lithium -> lithium
#   - report the composite Zsummary and medianRank per module
#   - Z < 2 not preserved; 2 <= Z < 10 moderately preserved; Z >= 10 well preserved
#   - their reported result: "full preservation", i.e. essentially nothing below Z = 2
#
# USAGE
#   Rscript src/06_module_preservation.R [--perms=200] [--pair=casecontrol|lithium|both]
#
#   The four group networks are the expensive part (27-42 min each here) and are cached
#   to data/processed/wgcna_net_<group>.rds. The preservation result for each pair is
#   cached to data/processed/mp_rows_<pair>_<perms>.rds, so the two pairs can be run as
#   two concurrent processes (--pair=casecontrol and --pair=lithium) and the CSV is then
#   assembled from whichever pair caches exist for that permutation count. Delete the
#   .rds files to force a recompute.
#
# PERMUTATION BUDGET -- read before quoting any number from the output.
#   The count actually used is recorded in the nPermutations column of
#   results/module_preservation.csv and printed in the summary. The default is the
#   paper's 200. Zsummary is a Z score against the permutation mean and sd of each
#   component statistic, so its sampling error falls as 1/sqrt(nPerm) and it is normally
#   stable by ~50. The verdict thresholds (2 and 10) are coarse relative to that error,
#   so a module rarely changes verdict between 50 and 200 permutations -- but a module
#   within ~1 Z unit of a threshold should be treated as undecided if fewer than 200
#   permutations were used.
#
# THREE DEVIATIONS from a literal reading of their methods, all forced and all reported:
#   (1) Covariate residualisation is the cached all-444 fit from 03_wgcna.R, not a
#       within-group fit. Their pipeline residualises on assessment group, which is
#       constant within cases (every case is group A), so a within-case fit is singular
#       and would have to drop a covariate the other three groups keep -- which would
#       make the four networks non-comparable. One joint fit keeps all four on one scale.
#   (2) Gene universe is their published 12,344 (wgcna_residuals.rds is already
#       restricted to it), so this is like-for-like with their Supplementary File 2 even
#       though our own DE filter gives 12,353.
#   (3) "non-lithium users" is non-lithium CASES (n = 87 here), per their Table 1, not
#       the 292 people in the cohort who are not on lithium.
#
# Signed/unsigned is inherited from 03_wgcna.R (unsigned adjacency, signed TOM); the
# methods do not state it and beta in the 4-7 range is the unsigned convention.

suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors = FALSE)

script_path <- sub("^--file=", "",
                   grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
ROOT <- if (!is.na(script_path)) dirname(dirname(normalizePath(script_path))) else "."
if (!dir.exists(file.path(ROOT, "data"))) ROOT <- "."
PROC <- file.path(ROOT, "data", "processed")
RES  <- file.path(ROOT, "results")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

args   <- commandArgs(TRUE)
getarg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit[1])
}
NPERM <- as.integer(getarg("perms", "200"))
WHICH <- getarg("pair", "both")

# Betas are theirs; the rest is the main-pipeline setting from 03_wgcna.R, which their
# methods say the group networks also used ("the same pipeline as described in Methods").
BETA <- c(cases = 4.2, controls = 5.5, lithium = 7, nonlithium = 7)
MIN_MODULE_SIZE <- 30; DEEP_SPLIT <- 2; MERGE_CUT_HEIGHT <- 0.25
SEED <- 12345
PAIRS <- list(casecontrol = c("cases", "controls"),
              lithium     = c("lithium", "nonlithium"))
stopifnot(WHICH %in% c(names(PAIRS), "both"))
todo <- if (WHICH == "both") names(PAIRS) else WHICH

cat(sprintf("=== Module preservation | nPermutations = %d | pair(s): %s ===\n",
            NPERM, paste(todo, collapse = ", ")))

# ------------------------------------------------------------------ inputs
meta <- read.csv(file.path(PROC, "sample_metadata.csv"), check.names = FALSE)
meta <- meta[meta$qc_pass %in% c(TRUE, "True", "TRUE"), ]
stopifnot(nrow(meta) == 444, !anyDuplicated(meta$title))

datExpr <- readRDS(file.path(PROC, "wgcna_residuals.rds"))$datExpr  # samples x genes
# The cache was written in metadata order; assert rather than assume, because every
# group below is a metadata column indexed into this matrix.
stopifnot(identical(rownames(datExpr), as.character(meta$title)))
# quickCor = 1 below skips the NA-tolerant correlation path, which is exact only if
# there is nothing for that path to do.
stopifnot(!anyNA(datExpr), min(apply(datExpr, 2, sd)) > 0)
cat(sprintf("residual matrix: %d samples x %d genes\n", nrow(datExpr), ncol(datExpr)))

meta$lithium <- as.numeric(meta$lithium)
is_case <- meta$diagnosis != "Control"
grp_idx <- list(
  cases      = which(is_case),
  controls   = which(!is_case),
  lithium    = which(is_case & meta$lithium == 1),
  nonlithium = which(is_case & meta$lithium == 0))   # non-lithium CASES, not all non-users
stopifnot(all(meta$lithium[!is_case] == 0))          # controls are never on lithium

paper_n <- c(cases = 240, controls = 204, lithium = 152, nonlithium = 88)
cat("group sizes (paper's Table 1 in brackets):\n")
for (g in names(grp_idx))
  cat(sprintf("  %-11s n = %3d   [%3d]   beta = %.1f\n",
              g, length(grp_idx[[g]]), paper_n[[g]], BETA[[g]]))

nthreads <- tryCatch({ enableWGCNAThreads(); WGCNAnThreads() },
                     error = function(e) { allowWGCNAThreads(); 1L })
cat(sprintf("WGCNA threads: %s\n", nthreads))

# ------------------------------------------------------- per-group networks
# One block per network: a split network would cluster genes that never met, and the
# resulting module set would not be comparable across groups or with 03_wgcna.R.
build_net <- function(g) {
  f <- file.path(PROC, sprintf("wgcna_net_%s.rds", g))
  if (file.exists(f)) { cat(sprintf("  [%s] cached\n", g)); return(readRDS(f)) }
  cat(sprintf("  [%s] blockwiseModules beta = %.1f, n = %d -- slow ...\n",
              g, BETA[[g]], length(grp_idx[[g]])))
  t0 <- Sys.time()
  net <- blockwiseModules(datExpr[grp_idx[[g]], , drop = FALSE],
                          power             = BETA[[g]],
                          networkType       = "unsigned",
                          TOMType           = "signed",
                          maxBlockSize      = 13000,
                          minModuleSize     = MIN_MODULE_SIZE,
                          deepSplit         = DEEP_SPLIT,
                          mergeCutHeight    = MERGE_CUT_HEIGHT,
                          numericLabels     = TRUE,
                          pamRespectsDendro = TRUE,
                          saveTOMs          = FALSE,
                          randomSeed        = SEED,
                          verbose           = 0)
  cat(sprintf("  [%s] done in %.1f min, %d block(s)\n", g,
              as.numeric(difftime(Sys.time(), t0, units = "mins")),
              length(net$blockGenes)))
  # Keep only what is needed downstream; the TOM-adjacent slots are large and disk on
  # this machine is constrained.
  out <- list(colors = labels2colors(net$colors), nBlocks = length(net$blockGenes))
  names(out$colors) <- colnames(datExpr)
  saveRDS(out, f); rm(net); invisible(gc()); out
}

need <- unique(unlist(PAIRS[todo]))
cat("\n--- group networks ---\n")
nets <- lapply(need, build_net); names(nets) <- need
invisible(gc())

cat("\nmodule structure per group (grey = unassigned, not a module):\n")
for (g in need) {
  tb <- table(nets[[g]]$colors); nm <- setdiff(names(tb), "grey")
  cat(sprintf("  %-11s blocks %d | %2d modules | sizes %d-%d | grey %d\n",
              g, nets[[g]]$nBlocks, length(nm), min(tb[nm]), max(tb[nm]), tb[["grey"]]))
  if (nets[[g]]$nBlocks != 1)
    cat(sprintf("    WARNING: %s split into %d blocks -- its modules are not comparable\n",
                g, nets[[g]]$nBlocks), "    to a single-block network.\n")
}

# ---------------------------------------------------------- preservation
verdict_of <- function(z) ifelse(is.na(z), NA_character_,
                          ifelse(z < 2,  "not preserved",
                          ifelse(z < 10, "moderately preserved", "well preserved")))

# modulePreservation nests its output as [[position in referenceNetworks]][[set index of
# the test network]] and names those levels "ref.<set>" / "inColumnsAlsoPresentIn.<set>".
# Index numerically -- referenceNetworks = c(1,2) pins the outer position to the set
# number -- but assert the names so a WGCNA change cannot silently transpose a reference
# and a test network.
extract <- function(mp, sets, iref, itest) {
  ref <- sets[iref]; test <- sets[itest]
  Z   <- mp$preservation$Z[[iref]][[itest]]
  obs <- mp$preservation$observed[[iref]][[itest]]
  stopifnot(names(mp$preservation$Z)[iref] == paste0("ref.", ref),
            names(mp$preservation$Z[[iref]])[itest] ==
              paste0("inColumnsAlsoPresentIn.", test),
            is.data.frame(Z), is.data.frame(obs),
            identical(rownames(Z), rownames(obs)))
  data.frame(comparison     = sprintf("%s(ref) -> %s(test)", ref, test),
             reference      = ref,
             test           = test,
             beta_reference = BETA[[ref]],
             beta_test      = BETA[[test]],
             n_reference    = length(grp_idx[[ref]]),
             n_test         = length(grp_idx[[test]]),
             module         = rownames(Z),
             moduleSize     = Z$moduleSize,
             Zsummary       = Z$Zsummary.pres,
             Zdensity       = Z$Zdensity.pres,
             Zconnectivity  = Z$Zconnectivity.pres,
             medianRank     = obs$medianRank.pres,
             verdict        = verdict_of(Z$Zsummary.pres),
             nPermutations  = NPERM,
             row.names      = NULL)
}

# Two calls, each with referenceNetworks = c(1, 2), give exactly the four ref -> test
# pairs they specify while computing each set's correlation structure once, not twice.
run_pair <- function(pn) {
  cachef <- file.path(PROC, sprintf("mp_rows_%s_%d.rds", pn, NPERM))
  if (file.exists(cachef)) {
    cat(sprintf("\n--- %s: cached result for %d permutations ---\n", pn, NPERM))
    return(readRDS(cachef))
  }
  sets <- PAIRS[[pn]]
  multiData  <- lapply(sets, function(g) list(data = datExpr[grp_idx[[g]], , drop = FALSE]))
  multiColor <- lapply(sets, function(g) nets[[g]]$colors)
  names(multiData) <- names(multiColor) <- sets

  cat(sprintf("\n--- modulePreservation: %s <-> %s, %d permutations ---\n",
              sets[1], sets[2], NPERM))
  t0 <- Sys.time()
  mp <- modulePreservation(
          multiData, multiColor,
          dataIsExpr             = TRUE,
          networkType            = "unsigned",
          referenceNetworks      = c(1, 2),
          nPermutations          = NPERM,
          randomSeed             = SEED,
          quickCor               = 1,      # exact here: no NAs, no zero-variance genes
          parallelCalculation    = FALSE,  # serial keeps peak memory at one copy of the
                                           # data; this machine has < 1 GB free
          savePermutedStatistics = TRUE,
          permutedStatisticsFile = file.path(PROC, sprintf("mp_perm_%s_%d.RData",
                                                           pn, NPERM)),
          plotInterpolation      = FALSE,
          verbose                = 3, indent = 2)   # verbose 3 prints per-permutation
                                                    # progress, so the ETA is visible
  cat(sprintf("  elapsed %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  rows <- rbind(extract(mp, sets, 1, 2), extract(mp, sets, 2, 1))
  saveRDS(rows, cachef); rm(mp); invisible(gc())
  rows
}

for (pn in todo) invisible(run_pair(pn))

# --------------------------------------------- assemble the CSV from pair caches
have <- names(PAIRS)[file.exists(file.path(PROC, sprintf("mp_rows_%s_%d.rds",
                                                         names(PAIRS), NPERM)))]
if (!setequal(have, names(PAIRS)))
  cat(sprintf("\n*** INCOMPLETE: pair(s) %s have no result at %d permutations. The CSV\n",
              paste(setdiff(names(PAIRS), have), collapse = ", "), NPERM),
      "    below covers only the pairs listed above.\n")
res <- do.call(rbind, lapply(have, function(pn)
         readRDS(file.path(PROC, sprintf("mp_rows_%s_%d.rds", pn, NPERM)))))

# "grey" is the unassigned bin and "gold" is modulePreservation's own random-gene
# control; neither is a module, so neither counts toward the preservation verdict.
res$is_module <- !res$module %in% c("grey", "gold")
res <- res[order(res$comparison, !res$is_module, res$Zsummary), ]
write.csv(res, file.path(RES, "module_preservation.csv"), row.names = FALSE)
cat(sprintf("\nWrote results/module_preservation.csv (%d rows)\n", nrow(res)))

# ---------------------------------------------------------------- report
cat("\n=== Zsummary by comparison (real modules only; grey and gold excluded) ===\n")
real <- res[res$is_module, ]
for (cmp in unique(real$comparison)) {
  s <- real[real$comparison == cmp, ]
  cat(sprintf("\n%s   [%d modules, sizes %d-%d]\n", cmp, nrow(s),
              min(s$moduleSize), max(s$moduleSize)))
  cat(sprintf("  Zsummary %.2f - %.2f (median %.2f) | medianRank %.1f - %.1f\n",
              min(s$Zsummary), max(s$Zsummary), median(s$Zsummary),
              min(s$medianRank), max(s$medianRank)))
  cat(sprintf("  not preserved (Z<2) %d | moderate (2-10) %d | well (Z>=10) %d\n",
              sum(s$verdict == "not preserved"),
              sum(s$verdict == "moderately preserved"),
              sum(s$verdict == "well preserved")))
  # The weakest modules are the only ones where the verdict is in any doubt.
  low <- head(s[order(s$Zsummary),
                c("module", "moduleSize", "Zsummary", "Zdensity",
                  "Zconnectivity", "medianRank", "verdict")], 3)
  cat("  three weakest modules:\n"); print(low, row.names = FALSE, digits = 3)
}

nfail <- sum(real$verdict == "not preserved")
nmod  <- sum(real$verdict == "moderately preserved")
cat("\n=== Verdict vs the paper's claim of full preservation ===\n")
cat(sprintf("  %d of %d reference modules have Zsummary < 2; %d are between 2 and 10\n",
            nfail, nrow(real), nmod))
cat(if (nfail == 0)
      "  -> REPRODUCES: no reference module falls below Z = 2 in any comparison.\n"
    else
      "  -> DOES NOT fully reproduce: the modules above fall below Z = 2.\n")
cat(sprintf("  based on %d permutations (paper: 200)\n", NPERM))
# Worth saying out loud: a Zsummary in the tens or hundreds is not a stronger claim than
# one of 15, it just means the permutation null is easy to beat. All four networks here
# are built from the same covariate-residualised matrix on overlapping subjects, so high
# preservation is the expected result and is weak evidence about biology.
cat("  NOTE: Zsummary scales with module size and sample size. Values this far above\n")
cat("  the threshold mean the permuted null is easy to beat, not that the effect is\n")
cat("  large. Use medianRank to compare modules with each other.\n")
