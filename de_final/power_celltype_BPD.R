#!/usr/bin/env Rscript
# ==============================================================================
# power_celltype_BPD.R -- cell-type detection floors for the BIPOLAR contrast.
#
# WHY THIS RUN IS NEEDED
#
# The per-lineage floors quoted so far (granulocytes 0.73, B cells 26.11) were
# measured on the LITHIUM contrast, n = 226 with 74 vs 152. They have been used
# as a proxy for the bipolar contrast, n = 308 with 74 vs 234. The two share a
# group of 74, which is the limiting factor, so the proxy is reasonable -- but
# it is a proxy, and the claim "CT_BPD's zero is a power problem" rests on it.
# This measures the bipolar floors directly.
#
# WHY A SPIKE-IN RATHER THAN A FORMULA
#
# At the whole-blood level an analytical calculation works: limma hands you the
# standard error of the group coefficient per gene. At the cell-type level it
# does not, because a within-lineage effect does not arrive at the bMIND
# estimate one-for-one. The deconvolution attenuates it by an amount that
# depends on the lineage's abundance and on bMIND's shrinkage, and neither is
# in a closed form. The only honest measurement runs the effect through the
# whole pipeline.
#
# CONSTRUCTION. A within-lineage effect of size L is injected where it actually
# lives -- in the bulk, scaled by that sample's lineage fraction:
#
#     bulk_i  <-  bulk_i + L * p_ci * y_i          for spiked genes
#
# then the deconvolution and the test are re-run on the spiked bulk. Group
# labels are permuted first, so the real bipolar contrast cannot contaminate
# the calibration and whatever is detected is attributable to the spike alone.
#
# Three effect sizes per lineage, spiked into DISJOINT gene blocks so one bMIND
# run covers all three, plus a fourth unspiked block that measures false
# positives. Checkpointed per lineage.
# ==============================================================================

suppressPackageStartupMessages({library(MIND); library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "de_final")); set.seed(481); options(width = 170)
RES <- "results"; dir.create(RES, showWarnings = FALSE, recursive = TRUE)
PARTS <- file.path(RES, "power_ct_parts"); dir.create(PARTS, showWarnings = FALSE, recursive = TRUE)

p <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
NBLOCK <- 200; NCORE <- 1      # 1, not 4: only ~0.6 GB free; each worker copies the data
MDE_BULK_BPD <- 0.208        # measured for whole blood, bipolar contrast

s  <- p$sel$BPD
mm <- p$meta[s$keep, , drop = FALSE]; mm$grp <- s$grp
mm$plate <- droplevels(factor(mm$plate)); mm$sex <- droplevels(factor(mm$sex))
mm$group <- droplevels(factor(mm$group))
cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]])) > 1, logical(1))]
X <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse = " + "))), data = mm)
stopifnot(qr(X)$rank == ncol(X))
X[, "grpcase"] <- X[sample(nrow(X)), "grpcase"]     # permute: measure the spike alone
y <- X[, "grpcase"]

bulk0 <- p$logtpm[, mm$title, drop = FALSE]
frac  <- p$lineage[mm$title, , drop = FALSE]
LINS  <- colnames(frac)

cat(sprintf("== power_celltype_BPD ==\nn = %d (%d controls vs %d cases) | %d genes per block | whole-blood floor %.3f\n",
            nrow(mm), sum(s$grp == "ctrl"), sum(s$grp == "case"), NBLOCK, MDE_BULK_BPD))
cat("\nlineage  mean_frac  predicted floor = bulk_floor / fraction\n")
for (ct in LINS)
  cat(sprintf("  %-5s   %.4f     %.2f log2FC\n", ct, mean(frac[, ct]), MDE_BULK_BPD / mean(frac[, ct])))

GRID <- list(gran = c(0.25, 0.50, 1.00), mono = c(0.40, 0.80, 1.60),
             T = c(0.40, 0.80, 1.60), NK = c(1.00, 2.50, 5.00), B = c(5.0, 15.0, 30.0))

allg   <- sample(rownames(bulk0))
blocks <- list(allg[1:NBLOCK], allg[(NBLOCK+1):(2*NBLOCK)], allg[(2*NBLOCK+1):(3*NBLOCK)])
ctrl   <- allg[(3*NBLOCK+1):(4*NBLOCK)]
keepg  <- c(unlist(blocks), ctrl)
# Subset ONCE, before the loop. bMIND only ever sees these 800 rows, so
# spiking into the full matrix and subsetting afterwards is identical
# arithmetic at ~70x the memory. This is what the earlier run did, and it is
# what the OS killed it for.
bulk_sub <- bulk0[keepg, , drop = FALSE]
rm(bulk0); invisible(gc())

test_one <- function(r, ct, genes) {
  Y <- r$A[, ct, ]; W <- 1 / (r$SE[, ct, ]^2)
  fin <- is.finite(W) & W > 0; if (any(!fin)) W[!fin] <- min(W[fin]); W <- W / mean(W)
  tt <- topTable(eBayes(lmFit(Y[, mm$title], X, weights = W[, mm$title]), trend = TRUE),
                 coef = "grpcase", number = Inf, sort.by = "none")
  q <- setNames(p.adjust(tt$P.Value, "BH"), rownames(tt))
  sum(q[genes] < 0.05, na.rm = TRUE) }

OUT <- list(); t0 <- Sys.time()
for (ct in LINS) {
  ck <- file.path(PARTS, sprintf("%s.rds", ct))
  if (file.exists(ck)) { OUT[[length(OUT)+1]] <- readRDS(ck)
    cat(sprintf("  %-5s from checkpoint\n", ct)); next }
  L3 <- GRID[[ct]]; tL <- Sys.time()
  B <- bulk_sub
  for (j in 1:3) {
    add <- L3[j] * frac[, ct] * y
    B[blocks[[j]], ] <- B[blocks[[j]], ] + matrix(rep(add, each = NBLOCK), nrow = NBLOCK) }
  r <- bMIND(B, frac = frac, ncore = NCORE)
  det <- vapply(1:3, function(j) test_one(r, ct, blocks[[j]]), numeric(1))
  fp  <- test_one(r, ct, ctrl)
  rows <- data.frame(lineage = ct, mean_frac = round(mean(frac[, ct]), 4), log2FC = L3,
                     bulk_equivalent = round(L3 * mean(frac[, ct]), 4),
                     detected = det, n_spiked = NBLOCK, power = round(det / NBLOCK, 3),
                     false_pos = fp, n_control = NBLOCK, stringsAsFactors = FALSE)
  saveRDS(rows, ck); OUT[[length(OUT)+1]] <- rows
  rm(r, B); invisible(gc())        # release the posterior arrays before the next lineage
  cat(sprintf("  %-5s power %.3f / %.3f / %.3f at L = %.2f / %.2f / %.2f | FP %d/%d  [%.1f min]\n",
              ct, det[1]/NBLOCK, det[2]/NBLOCK, det[3]/NBLOCK, L3[1], L3[2], L3[3],
              fp, NBLOCK, as.numeric(difftime(Sys.time(), tL, units = "mins"))))
  flush.console()
}
MDE <- do.call(rbind, OUT)
cat("\n=== cell-type power, bipolar contrast ===\n"); print(MDE, row.names = FALSE)

mde80 <- sapply(LINS, function(ct) {
  m <- MDE[MDE$lineage == ct, ]
  if (all(m$power < 0.8)) return(NA_real_)
  if (all(m$power >= 0.8)) return(min(m$log2FC))
  approx(m$power, m$log2FC, xout = 0.8, ties = "ordered")$y })

cat("\n=== detection floor at 80% power, bipolar contrast ===\n")
for (ct in LINS) {
  v <- mde80[ct]
  cat(sprintf("  %-5s %-20s predicted %6.2f | bulk-equivalent %.3f | %s\n", ct,
      if (is.na(v)) "> largest tested" else sprintf("%.2f log2FC", v),
      MDE_BULK_BPD / mean(frac[, ct]),
      if (is.na(v)) NA_real_ else v * mean(frac[, ct]),
      if (is.na(v) || v > 0.5) "NULL NOT INTERPRETABLE" else "null interpretable")) }
cat(sprintf("\n  whole blood, same contrast: %.3f log2FC. Meaningfulness floor set in advance: 0.5.\n",
            MDE_BULK_BPD))

write.csv(MDE, file.path(RES, "power_celltype_BPD.csv"), row.names = FALSE)
write.csv(data.frame(lineage = LINS, mde80 = as.numeric(mde80),
                     predicted = as.numeric(MDE_BULK_BPD / colMeans(frac)[LINS]),
                     bulk_equivalent = as.numeric(mde80 * colMeans(frac)[LINS]),
                     interpretable = !is.na(mde80) & mde80 <= 0.5),
          file.path(RES, "power_celltype_BPD_floors.csv"), row.names = FALSE)
cat(sprintf("\nwrote results/power_celltype_BPD*.csv [%.1f min]\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
