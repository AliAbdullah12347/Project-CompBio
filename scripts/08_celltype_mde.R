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
# L * p_c. Granulocytes are 39% of this mixture and B cells are 1.5%. The SAME
# biological effect inside a B cell is therefore ~27x harder to see than inside
# a granulocyte, before any statistical consideration at all.
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
# Effect sizes are chosen per lineage to bracket the arithmetic prediction
# MDE_lineage ~= MDE_bulk / p_c, so the result either confirms or refutes that
# prediction rather than just reporting a number.
# ==============================================================================

suppressPackageStartupMessages({library(MIND); library(limma)})
HERE <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/de_analysis"
setwd(HERE); set.seed(481)
RES <- "results"
p <- readRDS("data/prep.rds")
NSPIKE <- 600; NCORE <- 6
MDE_BULK <- 0.203                      # measured in 03, whole blood, LI contrast

s  <- p$sel$LI
mm <- p$meta[s$keep, , drop = FALSE]; mm$grp <- s$grp
mm$plate <- droplevels(mm$plate); mm$sex <- droplevels(mm$sex)
cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]])) > 1, logical(1))]
X <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse = " + "))), data = mm)
# Permute the group column: the calibration must measure the spike alone.
X[, "grpcase"] <- X[sample(nrow(X)), "grpcase"]
y <- X[, "grpcase"]

bulk0 <- p$logtpm[, mm$title, drop = FALSE]
frac  <- p$lineage[mm$title, , drop = FALSE]
spk   <- sample(rownames(bulk0), NSPIKE)
bg    <- setdiff(rownames(bulk0), spk)
keepg <- c(spk, sample(bg, NSPIKE))          # spiked + unspiked controls

# Per lineage, effect sizes bracketing MDE_BULK / p_c.
GRID <- list(
  gran = c(0.25, 0.5, 1.0),
  T    = c(0.4, 0.8, 1.6),
  B    = c(5, 15, 30))

cat(sprintf("== 08_celltype_mde ==\nn=%d | spiked %d genes + %d controls | bulk MDE %.3f\n",
            nrow(mm), NSPIKE, NSPIKE, MDE_BULK))
cat("\nlineage  mean_frac  predicted MDE = bulk_MDE / frac\n")
for (ct in names(GRID))
  cat(sprintf("  %-5s   %.4f     %.2f log2FC\n", ct, mean(frac[, ct]), MDE_BULK / mean(frac[, ct])))

run_one <- function(ct, L) {
  B <- bulk0
  add <- L * frac[, ct] * y                     # the within-lineage effect, as it reaches bulk
  B[spk, ] <- B[spk, ] + matrix(rep(add, each = NSPIKE), nrow = NSPIKE)
  r <- bMIND(B[keepg, , drop = FALSE], frac = frac, ncore = NCORE)
  Y <- r$A[, ct, ]; W <- 1 / (r$SE[, ct, ]^2); W[!is.finite(W)] <- NA
  W <- W / mean(W, na.rm = TRUE)
  tt <- topTable(eBayes(lmFit(Y[, mm$title], X, weights = W[, mm$title]), trend = TRUE),
                 coef = "grpcase", number = Inf, sort.by = "none")
  q <- p.adjust(tt$P.Value, "BH")
  det <- sum(q[match(spk, rownames(tt))] < 0.05, na.rm = TRUE)
  fp  <- sum(q[match(setdiff(keepg, spk), rownames(tt))] < 0.05, na.rm = TRUE)
  data.frame(lineage = ct, mean_frac = round(mean(frac[, ct]), 4), log2FC = L,
             bulk_equivalent = round(L * mean(frac[, ct]), 4),
             detected = det, power = round(det / NSPIKE, 3),
             false_pos = fp, stringsAsFactors = FALSE)
}

OUT <- list()
for (ct in names(GRID)) for (L in GRID[[ct]]) {
  t0 <- Sys.time()
  r <- run_one(ct, L)
  OUT[[length(OUT) + 1]] <- r
  cat(sprintf("  %-5s L=%5.2f  bulk-equiv %.3f  power %.3f  FP %d   [%.1f min]\n",
              ct, L, r$bulk_equivalent, r$power, r$false_pos,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
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
  cat(sprintf("    %-5s %s log2FC   (predicted %.2f from abundance alone)\n", ct,
              if (is.na(mde80[ct])) "> largest tested" else sprintf("%.2f", mde80[ct]),
              MDE_BULK / mean(frac[, ct])))
cat(sprintf("\n  meaningfulness floor from config = 0.5 log2FC\n"))
cat("  A lineage whose MDE exceeds that floor cannot support a null claim:\n")
for (ct in names(mde80)) {
  v <- mde80[ct]
  cat(sprintf("    %-5s %s\n", ct,
              if (is.na(v) || v > 0.5) "NULL NOT INTERPRETABLE -- underpowered by construction"
              else "null interpretable"))
}
write.csv(MDE, file.path(RES, "celltype_mde.csv"), row.names = FALSE)
write.csv(data.frame(lineage = names(mde80), mde80 = as.numeric(mde80),
                     predicted = MDE_BULK / colMeans(frac)[names(mde80)]),
          file.path(RES, "celltype_mde80.csv"), row.names = FALSE)
cat("\nwrote results/celltype_mde*.csv\n")
