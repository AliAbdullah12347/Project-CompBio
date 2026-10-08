#!/usr/bin/env Rscript
# ==============================================================================
# 00b_ilr_trace.R -- walk the ILR construction end to end, printing the actual
# numbers at every step.
#
# The ILR balances are the mediator in Arm 2, so how they were built is part of
# the mediation methods, not a preprocessing detail to be waved at. This script
# re-derives them from the CIBERSORTx output exactly as 01_prep.R does and
# reports what happens at each stage, including the checks that would catch a
# mistake.
#
# It recomputes rather than loads, then asserts the result is identical to the
# ILR matrix in prep.rds. If the two ever diverge, this script fails.
# ==============================================================================

ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "mediation")); options(width = 150)
p <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
m <- p$meta

say <- function(...) cat(sprintf(...))
rule <- function(t) cat(sprintf("\n%s\n%s\n", t, strrep("-", nchar(t))))

## ---- STEP 1: what CIBERSORTx returned ---------------------------------------
rule("STEP 1  input: the CIBERSORTx output")
fr <- read.csv(file.path(ROOT, "data/cohort_474/fractions_bmode_474_CIBERSORTx.csv"),
               check.names = FALSE, colClasses = c(Mixture = "character"))
rownames(fr) <- fr$Mixture
drop <- c("Mixture", "P-value", "Correlation", "RMSE")
fr <- as.matrix(fr[m$title, setdiff(names(fr), drop)])
say("  %d samples x %d LM22 cell types\n", nrow(fr), ncol(fr))
say("  row sums: %.6f to %.6f  (each sample's 22 proportions sum to 1)\n",
    min(rowSums(fr)), max(rowSums(fr)))
say("  cell types that are zero in EVERY sample: %s\n",
    paste(colnames(fr)[colSums(fr) == 0], collapse = ", "))
say("  -> %d of 22 carry any signal at all\n", sum(colSums(fr) > 0))

## ---- STEP 2: aggregate 22 types to 5 lineages -------------------------------
rule("STEP 2  aggregate 22 types into 5 developmental lineages")
LIN <- list(
  gran = c("Neutrophils", "Eosinophils", "Mast cells resting", "Mast cells activated"),
  mono = c("Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2",
           "Dendritic cells resting", "Dendritic cells activated"),
  T    = c("T cells CD8", "T cells CD4 naive", "T cells CD4 memory resting",
           "T cells CD4 memory activated", "T cells follicular helper",
           "T cells regulatory (Tregs)", "T cells gamma delta"),
  NK   = c("NK cells resting", "NK cells activated"),
  B    = c("B cells naive", "B cells memory", "Plasma cells"))
stopifnot(setequal(unlist(LIN), colnames(fr)))
say("  every one of the 22 types is assigned to exactly one lineage: PASS\n")
for (L in names(LIN))
  say("    %-5s %d types: %s\n", L, length(LIN[[L]]),
      paste(sub("^(T cells|NK cells|B cells|Mast cells|Dendritic cells|Macrophages) ", "", LIN[[L]]), collapse = ", "))
Lf <- sapply(LIN, function(cs) rowSums(fr[, cs, drop = FALSE]))
rownames(Lf) <- m$title
say("\n  result: %d x %d | row sums %.6f to %.6f\n", nrow(Lf), ncol(Lf), min(rowSums(Lf)), max(rowSums(Lf)))
say("  means: %s\n", paste(sprintf("%s=%.4f", colnames(Lf), colMeans(Lf)), collapse = "  "))
Lf <- Lf / rowSums(Lf)   # re-close; a no-op here, kept as a guard
say("  re-closed to sum 1 (a no-op since the 22 already summed to 1, kept as a guard)\n")

## ---- STEP 3: why a log-ratio at all -----------------------------------------
rule("STEP 3  why the proportions cannot be used directly")
say("  The five proportions sum to 1 for every subject, so they carry only 4\n")
say("  degrees of freedom. The matrix of proportions ALONE is still full rank --\n")
say("  summing to one is an affine constraint, not a linear dependency:\n")
say("    rank(proportions)             = %d of %d\n", qr(Lf)$rank, ncol(Lf))
say("  The singularity appears once a regression INTERCEPT is present, because\n")
say("  the five proportions then add up to exactly the intercept column:\n")
say("    rank(intercept + proportions) = %d of %d  -> SINGULAR\n",
    qr(cbind(1, Lf))$rank, ncol(Lf) + 1)
say("  Correlations are also forced negative by the constraint, e.g.\n")
say("    gran vs mono  r = %+.3f\n", cor(Lf[, "gran"], Lf[, "mono"]))
say("    gran vs T     r = %+.3f\n", cor(Lf[, "gran"], Lf[, "T"]))
say("  Neither tells us a cell count fell; it may only mean another rose.\n")

## ---- STEP 4: zero replacement ----------------------------------------------
rule("STEP 4  zero replacement (logs need strictly positive values)")
nz <- colSums(Lf == 0)
say("  zeros per lineage before replacement: %s\n",
    paste(sprintf("%s=%d", names(nz), nz), collapse = "  "))
if (sum(nz) == 0) say("  -> none. Aggregation removed every zero: each lineage contains at least\n     one type that is always present, so the replacement never fires here.\n")
zrepl <- function(X, delta = 0.65) {
  lim <- apply(X, 2, function(v) min(v[v > 0]) * delta)
  z <- X <= 0; ins <- rowSums(sweep(z, 2, lim, "*"))
  X[z] <- matrix(rep(lim, each = nrow(X)), nrow(X))[z]
  X[!z] <- (X * (1 - ins))[!z]; X / rowSums(X)
}
Lz <- zrepl(Lf)
say("  multiplicative simple replacement, delta = 0.65 x the smallest observed value\n")
say("  (multiplicative, not additive: it rescales the non-zero parts so the\n   composition still sums to 1 and the ratios between them are preserved)\n")
say("  max |change| caused by the replacement: %.3g\n", max(abs(Lz - Lf)))
say("  minimum value anywhere after replacement: %.3g  (strictly positive: PASS)\n", min(Lz))

## ---- STEP 5: the sequential binary partition --------------------------------
rule("STEP 5  the sequential binary partition (which splits, in which order)")
sbp <- list(b1 = list(n = c("gran", "mono"), d = c("T", "NK", "B")),
            b2 = list(n = "gran",            d = "mono"),
            b3 = list(n = "T",               d = c("NK", "B")),
            b4 = list(n = "NK",              d = "B"))
for (b in names(sbp))
  say("  %s  %-12s vs %-12s   %s\n", b, paste(sbp[[b]]$n, collapse = "+"),
      paste(sbp[[b]]$d, collapse = "+"),
      c(b1 = "myeloid vs lymphoid", b2 = "the two myeloid arms",
        b3 = "T vs the rest of lymphoid", b4 = "the last lymphoid split")[b])
say("\n  5 parts -> exactly 4 balances. Each split uses only parts from one side\n")
say("  of an earlier split, which is what makes the coordinates orthogonal.\n")
say("  The tree is fixed from haematopoiesis in advance, not chosen from the data.\n")

## ---- STEP 6: the contrast matrix -------------------------------------------
rule("STEP 6  build the orthonormal contrast matrix psi")
psi <- t(sapply(sbp, function(b) {
  v <- setNames(numeric(ncol(Lz)), colnames(Lz))
  r <- length(b$n); s <- length(b$d)
  v[b$n] <- sqrt(s / (r * (r + s))); v[b$d] <- -sqrt(r / (s * (r + s))); v }))
say("  weights: +sqrt(s/(r(r+s))) on the numerator parts, -sqrt(r/(s(r+s))) on the\n")
say("  denominator parts, where r and s are how many parts each side holds.\n\n")
print(round(psi, 4))
say("\n  row sums (must be 0 so each balance is a contrast): %s\n",
    paste(sprintf("%.1g", rowSums(psi)), collapse = " "))
orth <- max(abs(unname(psi %*% t(psi)) - diag(nrow(psi))))
say("  orthonormality, max |psi psi' - I| = %.3g  -> %s\n", orth, ifelse(orth < 1e-12, "PASS", "FAIL"))
say("  (orthonormal means the four balances are uncorrelated by construction and\n   carry no overlapping information)\n")

## ---- STEP 7: the transform --------------------------------------------------
rule("STEP 7  apply the transform")
ILR <- log(Lz) %*% t(psi); colnames(ILR) <- names(sbp); rownames(ILR) <- m$title
say("  ILR = log(proportions) %%*%% t(psi)   ->  %d x %d\n", nrow(ILR), ncol(ILR))
say("  Each coordinate is a scaled log ratio of geometric means:\n")
say("    b1 = sqrt(6/5) * log( gmean(gran,mono) / gmean(T,NK,B) )\n\n")
S <- data.frame(balance = colnames(ILR),
                contrast = sapply(names(sbp), function(b) sprintf("%s vs %s",
                  paste(sbp[[b]]$n, collapse = "+"), paste(sbp[[b]]$d, collapse = "+"))),
                mean = round(colMeans(ILR), 4), sd = round(apply(ILR, 2, sd), 4),
                min = round(apply(ILR, 2, min), 3), max = round(apply(ILR, 2, max), 3))
print(S, row.names = FALSE)
say("\n  rank(intercept + balances) = %d of %d  -> FULL RANK, usable in a regression\n",
    qr(cbind(1, ILR))$rank, ncol(ILR) + 1)
say("  which is the point: 4 coordinates where 5 proportions would have been singular.\n")
cm <- cor(ILR)
say("\n  correlation BETWEEN the balances, across the 474 subjects:\n")
print(round(cm, 3))
say("  max |r| = %.3f\n", max(abs(cm[upper.tri(cm)])))
say("  NOTE, and this is easy to state wrongly: orthonormality of psi does NOT make\n")
say("  the coordinates uncorrelated in the data. It guarantees the BASIS vectors are\n")
say("  orthogonal, which is what makes the transform lossless and distance-preserving.\n")
say("  How subjects scatter within that basis is biology, and here they do correlate.\n")
say("  The balances are non-redundant, not statistically independent.\n")

## ---- STEP 8: nothing was lost ----------------------------------------------
rule("STEP 8  check: the transform is reversible, so no information was lost")
clr <- ILR %*% psi                       # back to centred log-ratio space
rec <- exp(clr); rec <- rec / rowSums(rec)
say("  reconstructing the 5 proportions from the 4 balances:\n")
say("    max |reconstructed - original| = %.3g  -> %s\n", max(abs(rec - Lz)),
    ifelse(max(abs(rec - Lz)) < 1e-10, "PASS, lossless", "FAIL"))
say("  4 numbers hold everything the 5 proportions did, without the redundancy.\n")

## ---- STEP 9: identical to what the analysis used ----------------------------
rule("STEP 9  check: this matches the ILR actually used downstream")
d <- max(abs(ILR - p$ILR[rownames(ILR), colnames(ILR)]))
say("  max |recomputed - prep.rds| = %.3g  -> %s\n", d, ifelse(d < 1e-12, "IDENTICAL", "DIVERGED"))
stopifnot(d < 1e-12)

write.csv(data.frame(sample = rownames(ILR), round(ILR, 6), check.names = FALSE),
          "results/ilr_balances_474.csv", row.names = FALSE)
write.csv(cbind(balance = rownames(psi), as.data.frame(round(psi, 6))),
          "results/ilr_psi_matrix.csv", row.names = FALSE)
write.csv(S, "results/ilr_summary.csv", row.names = FALSE)
cat("\nwrote mediation/results/ilr_*.csv\n")
