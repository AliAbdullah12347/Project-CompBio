#!/usr/bin/env Rscript
# ==============================================================================
# 08_celltype_mde.R -- what size of within-lineage effect could we have seen?
#
# WHY THIS IS THE MOST IMPORTANT SCRIPT FOR THE CELL-TYPE NULLS
#
# Bulk expression is the abundance-weighted average of its cell types:
#
#     bulk_i  =  sum_c  p_ci * x_ci
#
# so a within-lineage effect of size L in lineage c moves the bulk by only
# L * p_c. Granulocytes are 38.7% of this mixture and B cells are 1.45%. The
# SAME biological effect inside a B cell is therefore ~27x harder to see than
# inside a granulocyte, before any statistical consideration at all.
#
# That makes "no B-cell genes" and "no granulocyte genes" completely different
# statements, and reporting them as one number would be misleading. This script
# measures the detection floor for each lineage separately.
#
# CONSTRUCTION. A within-lineage effect is simulated where it actually lives --
# in the BULK, scaled by that lineage's per-sample fraction:
#
#     bulk_i  <-  bulk_i + L * p_ci * y_i          for spiked genes
#
# then the whole deconvolution-and-test pipeline runs on the spiked bulk. This
# is the honest construction: spiking bMIND's OUTPUT instead would skip the
# deconvolution entirely and measure nothing except limma.
#
# Labels are permuted first so no real signal contaminates the calibration.
#
# Effect sizes bracket the arithmetic prediction MDE_lineage ~= MDE_bulk / p_c,
# so the result either confirms or refutes that prediction rather than just
# reporting a number.
#
# ------------------------------------------------------------------------------
# WHY THREE EFFECT SIZES SHARE ONE bMIND RUN
#
# bMIND's model is per-gene: profile[i,] and covariance[i,,] are built from
# gene i alone and each gene gets its own MCMC. The ONLY cross-gene coupling is
# the final clip of posterior means to [min(X), max(X)] of the input matrix.
# 06b measures that clip; it affects a negligible fraction of values, so three
# effect sizes can be tested in ONE run by spiking three DISJOINT blocks of
# genes at three different levels, plus a fourth unspiked block as controls.
#
# This cuts the cost per lineage threefold and lets all FIVE lineages be
# covered -- which is what config.yaml asks for ("applied_to: each lineage").
# The previous version covered only three and would have left the mono and NK
# nulls without a detection floor, i.e. uninterpretable.
#
# If 06b had found clipping to be common, this sharing would be unsafe and the
# fallback is SHARE_RUN <- FALSE, which runs each (lineage, level) separately.
# ==============================================================================

suppressPackageStartupMessages({library(MIND); library(limma)})
HERE <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/de_analysis"
setwd(HERE); set.seed(481)
RES <- "results"
p <- readRDS("data/prep.rds")

# NCORE is 4, not 6. At 6 workers this machine (7.7 GB) ran out of physical
# memory and began paging: 1.5 GB of pagefile in use and only ~1.2 of the 6
# cores actually busy, which made the run roughly 5x SLOWER than fewer workers
# would have been. More parallelism is not free when each PSOCK worker carries
# its own copy of the data.
NBLOCK   <- 200        # genes per spiked block, and per control block
NCORE    <- 4
MDE_BULK <- 0.203      # measured in 03, whole blood, LI contrast
SHARE_RUN <- TRUE      # see header; set FALSE to isolate every (lineage, level)
PARTS    <- file.path("results", "celltype_mde_parts")   # per-lineage checkpoints
dir.create(PARTS, showWarnings = FALSE, recursive = TRUE)

# Effect sizes bracketing MDE_BULK / p_c for each lineage.
GRID <- list(
  gran = c(0.25, 0.50, 1.00),    # predicted 0.52
  mono = c(0.40, 0.80, 1.60),    # predicted 0.75
  T    = c(0.40, 0.80, 1.60),    # predicted 0.82
  NK   = c(1.00, 2.50, 5.00),    # predicted 2.59
  B    = c(5.00, 15.0, 30.0))    # predicted 13.97

s  <- p$sel$LI
mm <- p$meta[s$keep, , drop = FALSE]; mm$grp <- s$grp
mm$plate <- droplevels(factor(mm$plate)); mm$sex <- droplevels(factor(mm$sex))
cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]])) > 1, logical(1))]
X <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse = " + "))), data = mm)
# Permute the group column: the calibration must measure the spike alone.
X[, "grpcase"] <- X[sample(nrow(X)), "grpcase"]
y <- X[, "grpcase"]

bulk0 <- p$logtpm[, mm$title, drop = FALSE]
frac  <- p$lineage[mm$title, , drop = FALSE]

cat(sprintf("== 08_celltype_mde ==\nn=%d | %d genes per block | bulk MDE %.3f | share run: %s\n",
            nrow(mm), NBLOCK, MDE_BULK, SHARE_RUN))
cat("\nlineage  mean_frac  predicted MDE = bulk_MDE / frac\n")
for (ct in names(GRID))
  cat(sprintf("  %-5s   %.4f     %.2f log2FC\n", ct, mean(frac[, ct]), MDE_BULK / mean(frac[, ct])))

# Disjoint gene blocks: three spiked, one control. Drawn once and reused for
# every lineage so differences between lineages are not gene-set differences.
allg   <- sample(rownames(bulk0))
blocks <- list(allg[1:NBLOCK],
               allg[(NBLOCK + 1):(2 * NBLOCK)],
               allg[(2 * NBLOCK + 1):(3 * NBLOCK)])
ctrl   <- allg[(3 * NBLOCK + 1):(4 * NBLOCK)]
keepg  <- c(unlist(blocks), ctrl)

test_one <- function(r, ct, genes) {
  Y <- r$A[, ct, ]; W <- 1 / (r$SE[, ct, ]^2)
  fin <- is.finite(W) & W > 0
  if (any(!fin)) W[!fin] <- min(W[fin])
  W <- W / mean(W)
  tt <- topTable(eBayes(lmFit(Y[, mm$title], X, weights = W[, mm$title]), trend = TRUE),
                 coef = "grpcase", number = Inf, sort.by = "none")
  q <- setNames(p.adjust(tt$P.Value, "BH"), rownames(tt))
  sum(q[genes] < 0.05, na.rm = TRUE)
}

OUT <- list()
for (ct in names(GRID)) {
  # Checkpoint per lineage: a lineage already on disk is not recomputed. Each
  # lineage is ~10 minutes of deconvolution, so losing the lot to an
  # interruption is avoidable and was avoided.
  ckpt <- file.path(PARTS, sprintf("%s.rds", ct))
  if (file.exists(ckpt)) {
    OUT[[length(OUT) + 1]] <- readRDS(ckpt)
    cat(sprintf("  %-5s loaded from checkpoint\n", ct)); next
  }
  L3 <- GRID[[ct]]
  t0 <- Sys.time()
  if (SHARE_RUN) {
    B <- bulk0
    for (j in 1:3) {
      add <- L3[j] * frac[, ct] * y
      B[blocks[[j]], ] <- B[blocks[[j]], ] + matrix(rep(add, each = NBLOCK), nrow = NBLOCK)
    }
    r <- bMIND(B[keepg, , drop = FALSE], frac = frac, ncore = NCORE)
    det <- vapply(1:3, function(j) test_one(r, ct, blocks[[j]]), numeric(1))
    fp  <- test_one(r, ct, ctrl)
    el  <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    rows <- data.frame(
      lineage = ct, mean_frac = round(mean(frac[, ct]), 4), log2FC = L3,
      bulk_equivalent = round(L3 * mean(frac[, ct]), 4),
      detected = det, n_spiked = NBLOCK, power = round(det / NBLOCK, 3),
      false_pos = fp, n_control = NBLOCK, stringsAsFactors = FALSE)
    saveRDS(rows, ckpt)
    OUT[[length(OUT) + 1]] <- rows
    cat(sprintf("  %-5s power %.3f / %.3f / %.3f at L = %.2f / %.2f / %.2f | FP %d/%d  [%.1f min]\n",
                ct, det[1] / NBLOCK, det[2] / NBLOCK, det[3] / NBLOCK, L3[1], L3[2], L3[3],
                fp, NBLOCK, el))
  } else {
    for (j in 1:3) {
      tj <- Sys.time(); B <- bulk0
      add <- L3[j] * frac[, ct] * y
      B[blocks[[1]], ] <- B[blocks[[1]], ] + matrix(rep(add, each = NBLOCK), nrow = NBLOCK)
      g <- c(blocks[[1]], ctrl)
      r <- bMIND(B[g, , drop = FALSE], frac = frac, ncore = NCORE)
      d <- test_one(r, ct, blocks[[1]]); f <- test_one(r, ct, ctrl)
      OUT[[length(OUT) + 1]] <- data.frame(
        lineage = ct, mean_frac = round(mean(frac[, ct]), 4), log2FC = L3[j],
        bulk_equivalent = round(L3[j] * mean(frac[, ct]), 4),
        detected = d, n_spiked = NBLOCK, power = round(d / NBLOCK, 3),
        false_pos = f, n_control = NBLOCK, stringsAsFactors = FALSE)
      cat(sprintf("  %-5s L=%6.2f power %.3f FP %d  [%.1f min]\n", ct, L3[j], d / NBLOCK, f,
                  as.numeric(difftime(Sys.time(), tj, units = "mins"))))
    }
  }
  flush.console()
}
MDE <- do.call(rbind, OUT)
cat("\n=== cell-type minimum detectable effect ===\n"); print(MDE, row.names = FALSE)

mde80 <- sapply(names(GRID), function(ct) {
  m <- MDE[MDE$lineage == ct, ]
  if (all(m$power < 0.8)) return(NA_real_)
  if (all(m$power >= 0.8)) return(min(m$log2FC))
  approx(m$power, m$log2FC, xout = 0.8, ties = "ordered")$y })

cat("\n  MDE at 80% power, by lineage:\n")
for (ct in names(mde80))
  cat(sprintf("    %-5s %-18s (predicted %.2f from abundance alone)\n", ct,
              if (is.na(mde80[ct])) "> largest tested" else sprintf("%.2f log2FC", mde80[ct]),
              MDE_BULK / mean(frac[, ct])))

cat("\n  config meaningfulness floor = 0.5 log2FC\n")
cat("  A lineage whose MDE exceeds that floor cannot support a null claim:\n")
for (ct in names(mde80)) {
  v <- mde80[ct]
  cat(sprintf("    %-5s %s\n", ct,
              if (is.na(v) || v > 0.5) "NULL NOT INTERPRETABLE -- underpowered by construction"
              else "null interpretable"))
}
write.csv(MDE, file.path(RES, "celltype_mde.csv"), row.names = FALSE)
write.csv(data.frame(lineage = names(mde80), mde80 = as.numeric(mde80),
                     predicted = as.numeric(MDE_BULK / colMeans(frac)[names(mde80)]),
                     interpretable = !is.na(mde80) & mde80 <= 0.5),
          file.path(RES, "celltype_mde80.csv"), row.names = FALSE)
cat("\nwrote results/celltype_mde*.csv\n")
