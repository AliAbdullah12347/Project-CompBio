#!/usr/bin/env Rscript
# ==============================================================================
# 06_bmind_profiles.R -- cell-type expression estimates, the usable way.
#
# WHY THIS REPLACES bmind_de()
#
# Two defects make bmind_de() unable to answer our hypotheses, both verified
# directly rather than suspected:
#
#   1. ITS P-VALUES HAVE A HARD FLOOR AT 2e-4. Only ~600 distinct p-values
#      appear across 61,840 tests. Benjamini-Hochberg at FDR 0.05 over 12,368
#      genes requires the smallest p-value to fall below 4.0e-6. The floor is
#      50x too large, so NO gene can reach significance however strong the true
#      effect is. Raising the MCMC iterations tenfold moved the floor from
#      0.086 to 0.073 on a test set -- it is structural, not a resolution
#      problem that more sampling fixes.
#
#   2. IT SILENTLY IGNORES COVARIATES. bmind_de(covariate = NULL),
#      bmind_de(covariate = CV) and bmind_de(covariate_bulk = CV) return
#      byte-identical p-values. The argument is accepted and discarded. Every
#      result from it is therefore unadjusted for age, sex, tobacco, RIN,
#      plate, sequencing PCs and the cell-composition balances.
#
# Boltz et al. used bmind_de and reported zero cell-type lithium DEGs. Defect 1
# means that null was very likely unobtainable rather than observed.
#
# THE REPLACEMENT
#
# Run bMIND() WITHOUT the phenotype to get per-sample, per-lineage expression
# estimates, then test them with limma in a design we control. This gives:
#   * continuous p-values with no floor,
#   * covariates that actually enter the model,
#   * the same estimator family as the whole-blood arm, so the two levels are
#     comparable rather than measured on different instruments.
#
# Running bMIND without y is also LESS circular: the phenotype never informs
# the deconvolution, so the cell-type estimates are not shaped by the
# comparison they will later be used for. One run serves both contrasts.
#
# bMIND returns posterior SEs. Those are carried forward and used as limma
# precision weights in 07, so the deconvolution's uncertainty is propagated
# rather than discarded.
# ==============================================================================

suppressPackageStartupMessages({library(MIND)})
HERE <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/de_analysis"
setwd(HERE)
p <- readRDS("data/prep.rds")
OUT <- file.path("results", "bmind_profiles"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
CHUNK <- 400; NCORE <- 6

bulk <- p$logtpm                      # all 474 samples, all genes
frac <- p$lineage
cat(sprintf("== 06_bmind_profiles ==\nbulk %d x %d | frac %d x %d | chunk %d | ncore %d\n\n",
            nrow(bulk), ncol(bulk), nrow(frac), ncol(frac), CHUNK, NCORE))

genes <- rownames(bulk)
chunks <- split(genes, ceiling(seq_along(genes) / CHUNK))
t0 <- Sys.time()
for (i in seq_along(chunks)) {
  f <- file.path(OUT, sprintf("chunk_%03d.rds", i))
  if (file.exists(f)) next
  g <- chunks[[i]]
  r <- tryCatch(bMIND(bulk[g, , drop = FALSE], frac = frac, ncore = NCORE),
                error = function(e) { cat(sprintf("chunk %d FAILED: %s\n", i, conditionMessage(e))); NULL })
  if (is.null(r)) next
  # A is genes x cell types x samples; SE likewise. Both kept.
  saveRDS(list(genes = g, A = r$A, SE = r$SE), f)
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("  chunk %3d/%d  elapsed %5.1f min  ETA %5.1f min\n",
              i, length(chunks), el, el / i * (length(chunks) - i)))
  flush.console()
}

## ---- assemble ---------------------------------------------------------------
fs <- sort(list.files(OUT, "^chunk_.*\\.rds$", full.names = TRUE))
cat(sprintf("\nassembling %d chunks\n", length(fs)))
parts <- lapply(fs, readRDS)
gn <- unlist(lapply(parts, `[[`, "genes"), use.names = FALSE)
cts <- dimnames(parts[[1]]$A)[[2]]; sm <- dimnames(parts[[1]]$A)[[3]]
A  <- array(NA_real_, c(length(gn), length(cts), length(sm)), list(gn, cts, sm))
SE <- array(NA_real_, c(length(gn), length(cts), length(sm)), list(gn, cts, sm))
k <- 0
for (q in parts) {
  idx <- (k + 1):(k + length(q$genes)); k <- k + length(q$genes)
  A[idx, , ] <- q$A; SE[idx, , ] <- q$SE
}
cat(sprintf("  A %s | SE %s\n", paste(dim(A), collapse = " x "), paste(dim(SE), collapse = " x ")))
cat(sprintf("  any NA in A: %s | SE range %.4g to %.4g\n",
            any(is.na(A)), min(SE, na.rm = TRUE), max(SE, na.rm = TRUE)))
saveRDS(list(A = A, SE = SE, lineages = cts, samples = sm),
        file.path("results", "bmind_profiles.rds"))
cat(sprintf("wrote results/bmind_profiles.rds [%.1f min total]\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
