#!/usr/bin/env Rscript
# ==============================================================================
# 06b_calibration_proposal.R -- close a gap: S3 (regression calibration) had only
# ever been run on the 4-balance LINEAGE tree, because reliability had only been
# measured for those four balances. The reported corrected figure of 0.587 was
# therefore a lineage-tree number quoted as if it covered both trees.
#
# This script measures reliability for the proposal's FIVE balances the same way
# 00_reliability.R did for the lineage ones -- by putting Krebs et al.'s
# deposited CIBERSORT fractions through the identical 6-part aggregation and the
# identical sequential binary partition -- and then runs the identical
# calibration on the proposal tree.
# ==============================================================================

ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "mediation")); set.seed(481)
suppressPackageStartupMessages({library(limma); library(edgeR)})
source("scripts/med.R"); options(width = 160)
say <- function(...) cat(sprintf(...))
rule <- function(t) cat(sprintf("\n%s\n%s\n", t, strrep("=", nchar(t))))

prep <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
dp <- readRDS("data/med_input_proposal.rds")
md <- read.csv(file.path(ROOT, "data/cohort_474/metadata_474_imputed.csv"),
               check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(identical(md$title, prep$meta$title))

## ---- Krebs' deposited fractions, aggregated to the proposal's six parts -----
KPARTS <- list(
  neut = "neutrophils", mono = "monocytes", cd4n = "t.cells.cd4.naive",
  cd4m = c("t.cells.cd4.memory.resting", "t.cells.cd4.memory.activated"),
  NK   = c("nk.cells.resting", "nk.cells.activated"),
  B    = c("b.cells.naive", "b.cells.memory"))
stopifnot(all(unlist(KPARTS) %in% names(md)))
KP <- sapply(KPARTS, function(cs)
  rowSums(sapply(cs, function(g) as.numeric(md[[g]])), na.rm = FALSE))
rownames(KP) <- md$title
say("== 06b: reliability of the proposal's five balances ==\n")
say("  Krebs six-part mass retained: mean %.3f | zeros per part: %s\n",
    mean(rowSums(KP)), paste(sprintf("%s=%d", names(KPARTS), colSums(KP == 0)), collapse = " "))
KP <- KP / rowSums(KP)

zrepl <- function(X, delta = 0.65) {
  lim <- apply(X, 2, function(v) { v <- v[v > 0]; if (!length(v)) 1e-8 else min(v) }) * delta
  z <- X <= 0; ins <- rowSums(sweep(z, 2, lim, "*"))
  X[z] <- matrix(rep(lim, each = nrow(X)), nrow(X))[z]
  X[!z] <- (X * (1 - ins))[!z]; X / rowSums(X) }

psi <- dp$psi
I_kr <- log(zrepl(KP)) %*% t(psi)
I_cx <- log(zrepl(dp$parts)) %*% t(psi)        # ours, same construction
colnames(I_kr) <- colnames(I_cx) <- rownames(psi)

icc_A1 <- function(x, y) {
  M <- cbind(x, y); n <- nrow(M); k <- ncol(M)
  gm <- mean(M); rm <- rowMeans(M); cm <- colMeans(M)
  MSR <- k * sum((rm - gm)^2) / (n - 1); MSC <- n * sum((cm - gm)^2) / (k - 1)
  MSE <- (sum((M - gm)^2) - k * sum((rm - gm)^2) - n * sum((cm - gm)^2)) / ((n - 1) * (k - 1))
  (MSR - MSE) / (MSR + (k - 1) * MSE + k * (MSC - MSE) / n) }

R <- data.frame(balance = colnames(I_cx), contrast = dp$summary$contrast,
                icc_absolute = round(vapply(seq_len(ncol(I_cx)),
                  function(j) icc_A1(I_cx[, j], I_kr[, j]), numeric(1)), 4),
                stringsAsFactors = FALSE)
print(R, row.names = FALSE)
write.csv(R, "results/balance_reliability_proposal.csv", row.names = FALSE)

## ---- the identical calibration, on the proposal tree ------------------------
rule("S3 on the proposal tree")
d <- dp; K <- ncol(d$M)
f <- readRDS("results/primary/proposal_fit.rds"); tb <- f$tb; mf <- f$mf
cM <- centre(d$M); cC <- centre(d$C); x <- d$lithium
Xd <- med_design_full(x, cM$X, cC$X); ix <- med_idx(K)
of <- outcome_fit(d$v_main$E, d$v_main$weights, Xd, se = TRUE, moderate = TRUE)
T2 <- of$theta[, ix$iM, drop = FALSE]; T3 <- of$theta[, ix$iI, drop = FALSE]
lam <- R$icc_absolute
sig <- tb$TE_q_BH < 0.05

Sw <- mf$Sigma; Su0 <- diag((1 - lam) * apply(d$M, 2, var), K, K)
minev <- function(cc) min(eigen(Sw - cc * Su0, symmetric = TRUE)$values)
if (minev(1) > 0) { cshrink <- 1 } else {
  lo <- 0; hi <- 1
  for (it in 1:60) { mid <- (lo + hi) / 2; if (minev(mid) > 1e-10) lo <- mid else hi <- mid }
  cshrink <- lo * 0.99 }
say("  admissible multiple of the implied error covariance: c = %.4f\n", cshrink)

# the REPORTED correction: per-balance, stable (same choice as the lineage tree)
NIEc <- as.vector(((T2 + T3) %*% diag(1 / lam, K, K)) %*% mf$b1)
PMc  <- NIEc / tb$TE
say("\n  median proportion mediated, proposal tree:\n")
say("    uncorrected : %.4f\n", median(tb$PM[sig]))
say("    corrected   : %.4f   (change %+.4f)\n", median(PMc[sig]), median(PMc[sig]) - median(tb$PM[sig]))
say("    share with PM > 0.50: %.1f%% -> %.1f%%\n",
    100 * mean(tb$PM[sig] > 0.5), 100 * mean(PMc[sig] > 0.5))

L <- readRDS("results/sens/calibration.rds")
rule("BOTH TREES SIDE BY SIDE")
B <- data.frame(
  tree = c("lineage (4 balances)", "proposal (5 balances)"),
  n_balances = c(4, 5),
  median_icc = c(median(read.csv("results/balance_reliability.csv")$icc_absolute),
                 median(lam)),
  PM_uncorrected = c(L$cal$uncorrected[2], median(tb$PM[sig])),
  PM_corrected   = c(L$cal$corrected[2],   median(PMc[sig])),
  stringsAsFactors = FALSE)
B$raises <- ifelse(B$PM_corrected > B$PM_uncorrected, "YES", "NO")
print(transform(B, median_icc = round(median_icc, 3),
                PM_uncorrected = round(PM_uncorrected, 4),
                PM_corrected = round(PM_corrected, 4)), row.names = FALSE)
write.csv(B, "results/sens/S3_both_trees.csv", row.names = FALSE)
saveRDS(list(reliability = R, PM_corr = PMc, summary = B), "results/sens/calibration_proposal.rds")
cat("\n06b complete.\n")
