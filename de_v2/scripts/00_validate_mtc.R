#!/usr/bin/env Rscript
# ==============================================================================
# 00_validate_mtc.R -- prove the hand-written corrections are correct before any
# of them touches real data.
#
# WHY THIS SCRIPT EXISTS
#
# mtc.R implements Storey q-values, Benjamini-Bogomolov and permutation FDR by
# hand, because this machine has no qvalue/IHW/fdrtool and CLAUDE.md forbids
# adding a dependency without asking. A hand-written FDR method that is subtly
# wrong does not announce itself: it just returns a plausible number of genes.
#
# So each one is pinned to an identity that must hold EXACTLY or to a simulation
# with known ground truth. The script stops on the first failure. If it prints
# ALL CHECKS PASSED, the numbers in the rest of this folder rest on something.
#
# Part 1 is exact identities. Part 2 is an operating-characteristics simulation:
# 200 datasets with a known set of true effects, measuring whether each method
# actually delivers the error rate it promises. That second part is the one that
# would catch a method that is self-consistent but miscalibrated.
# ==============================================================================

HERE <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/de_v2"
setwd(HERE)
source("scripts/mtc.R")
set.seed(481)
RES <- "results"; dir.create(RES, showWarnings = FALSE, recursive = TRUE)

fails <- 0
ok <- function(label, pass, detail = "") {
  cat(sprintf("  [%s] %-58s %s\n", if (pass) "PASS" else "FAIL", label, detail))
  if (!pass) fails <<- fails + 1
}

cat("== 00_validate_mtc ==\n\nPART 1: exact identities\n")

## ---- 1. Storey with pi0 = 1 IS Benjamini-Hochberg ---------------------------
# This is the definitional anchor. Storey q = pi0 * BH, so setting pi0 = 1 must
# reproduce p.adjust exactly -- not approximately. Any disagreement means the
# ordering, the monotonicity step, or the m/i factor is wrong.
p <- runif(5000)^2
d <- max(abs(qvalue_storey(p, pi0 = 1) - p.adjust(p, "BH")))
ok("Storey q at pi0=1 equals BH exactly", d < 1e-12, sprintf("max diff %.3g", d))

## ---- 2. q-values are monotone in p and bounded ------------------------------
q <- qvalue_storey(p, pi0 = 0.7)
ok("q-values monotone in p", all(diff(q[order(p)]) >= -1e-12))
ok("q-values within [0,1]", all(q >= 0 & q <= 1))

## ---- 3. ordering of the FWER / FDR panel ------------------------------------
# These inequalities are theorems, so they hold for ANY p-vector. If one breaks,
# the panel is wired up wrong.
a <- adjust_panel(p)
ok("Bonferroni >= Holm",   all(a$bonferroni >= a$holm - 1e-12))
ok("Holm >= Hommel",       all(a$holm >= a$hommel - 1e-12))
ok("Holm >= BH (FWER >= FDR)", all(a$holm >= a$BH - 1e-12))
ok("BY >= BH",             all(a$BY >= a$BH - 1e-12))
ok("Storey q <= BH (pi0<=1 so adaptive is never worse)",
   all(qvalue_storey(p) <= a$BH + 1e-12))

## ---- 4. pi0 estimators on a vector that is entirely null ---------------------
# Uniform p-values mean every hypothesis is null, so both estimators must land
# near 1. An estimator that returns 0.6 here would silently inflate every
# downstream discovery count by 40%.
pu <- runif(20000)
s <- pi0_spline(pu); b <- pi0_bootstrap(pu)
ok("pi0 spline ~ 1 on pure null",    abs(s - 1) < 0.05, sprintf("got %.4f", s))
ok("pi0 bootstrap ~ 1 on pure null", abs(b - 1) < 0.05, sprintf("got %.4f", b))

## ---- 5. pi0 estimators recover a known mixture ------------------------------
# 70% null, 30% strong alternative. Both should land near 0.70. Being slightly
# HIGH is acceptable and conservative; being low is not.
mix <- c(runif(14000), rbeta(6000, 0.3, 8))
s <- pi0_spline(mix); b <- pi0_bootstrap(mix)
ok("pi0 spline recovers 0.70",    abs(s - 0.70) < 0.08, sprintf("got %.4f", s))
ok("pi0 bootstrap recovers 0.70", abs(b - 0.70) < 0.08, sprintf("got %.4f", b))

## ---- 6. Simes ----------------------------------------------------------------
# Under the null the Simes p-value is uniform, so it should exceed 0.05 about
# 95% of the time. Checked over 2000 null families.
sim <- replicate(2000, simes(runif(50)))
r <- mean(sim < 0.05)
ok("Simes p uniform under the null", abs(r - 0.05) < 0.02, sprintf("rejection rate %.3f", r))

## ---- 7. Benjamini-Bogomolov collapses to BH when K = 1 ----------------------
# With one family that gets selected, the selection penalty R/K = 1 and the
# procedure must be ordinary BH. This checks the level-2 rescaling arithmetic.
pl <- list(f1 = c(rbeta(200, 0.2, 8), runif(800)))
bb <- bb_hierarchical(pl, q = 0.05)
d <- max(abs(bb$adjusted$f1 - p.adjust(pl$f1, "BH")), na.rm = TRUE)
ok("BB with K=1 reduces to BH", bb$n_selected == 1 && d < 1e-12, sprintf("max diff %.3g", d))

## ---- 8. permutation FDR against a uniform null ------------------------------
# If the null really is uniform then E[#{p<=t}] = m*t, and the permutation
# estimate must agree with what BH implies at the same threshold. This checks
# the plug-in formula and the monotonicity step.
m <- 12368; grid <- c(1e-5, 1e-4, 1e-3, 1e-2, 0.05)
pobs <- c(rbeta(600, 0.15, 10), runif(m - 600))
nullc <- t(replicate(60, vapply(grid, function(t) sum(runif(m) <= t), numeric(1))))
pf <- perm_fdr(pobs, grid, nullc, pi0 = 1)
bh_implied <- vapply(grid, function(t) m * t / max(1, sum(pobs <= t)), numeric(1))
bh_implied <- pmin(1, rev(cummin(rev(bh_implied))))
d <- max(abs(pf$fdr - bh_implied))
ok("permutation FDR matches BH under a uniform null", d < 0.02, sprintf("max diff %.4f", d))

cat("\nPART 2: operating characteristics on simulated data with known truth\n")
cat("  150 datasets, 8,000 genes, 5% truly non-null, target FDR 0.05\n\n")

## ---- 9. does each method deliver the error rate it promises? ----------------
# SCOPE. Only methods written BY HAND need this: Storey q-values and the pi0
# estimator behind them. BH, BY, Holm, Bonferroni and Hommel come from base R's
# p.adjust and are not ours to validate -- Part 1 already confirms this file
# wires them up in the right order. BH is carried along purely as the reference
# Storey must not beat on error rate.
#
# Hommel is deliberately EXCLUDED from this loop. Its algorithm is quadratic in
# the number of tests and takes minutes at m = 12,368; running it 150 times
# would cost hours to re-validate a base-R function. It is still reported on the
# real data, where it runs once per contrast.
#
# Self-consistency is not calibration: here the truth is known, so the realised
# false discovery proportion is measured directly rather than assumed.
NSIM <- 150; M <- 8000; PI1 <- 0.05; ALPHA <- 0.05
n1 <- round(M * PI1)
meth <- c("BH", "Storey")
res <- matrix(0, NSIM, 2, dimnames = list(NULL, meth))
pwr <- matrix(0, NSIM, 2, dimnames = list(NULL, meth))
pi0hat <- numeric(NSIM)
truth <- c(rep(TRUE, n1), rep(FALSE, M - n1))
for (i in seq_len(NSIM)) {
  z <- c(stats::rnorm(n1, 3.2), stats::rnorm(M - n1, 0))
  pv <- 2 * stats::pnorm(-abs(z))
  pi0hat[i] <- pi0_spline(pv)
  adj <- list(BH = stats::p.adjust(pv, "BH"), Storey = qvalue_storey(pv, pi0 = pi0hat[i]))
  for (mth in meth) {
    rej <- adj[[mth]] < ALPHA
    res[i, mth] <- if (sum(rej) == 0) 0 else sum(rej & !truth) / sum(rej)
    pwr[i, mth] <- sum(rej & truth) / n1
  }
}
cat("  mean realised false discovery proportion (target <= 0.05):\n")
for (mth in meth) cat(sprintf("    %-8s %.4f\n", mth, mean(res[, mth])))
cat(sprintf("  mean power (of %d true effects):\n", n1))
for (mth in meth) cat(sprintf("    %-8s %.4f\n", mth, mean(pwr[, mth])))
cat(sprintf("  mean estimated pi0: %.4f   (true value %.4f)\n", mean(pi0hat), 1 - PI1))

ok("BH controls FDR at 0.05 (reference)", mean(res[, "BH"]) <= 0.055, sprintf("%.4f", mean(res[, "BH"])))
ok("Storey controls FDR at 0.05",         mean(res[, "Storey"]) <= 0.055, sprintf("%.4f", mean(res[, "Storey"])))
ok("Storey at least as powerful as BH",   mean(pwr[, "Storey"]) >= mean(pwr[, "BH"]) - 1e-9,
   sprintf("%.4f vs %.4f", mean(pwr[, "Storey"]), mean(pwr[, "BH"])))
ok("pi0 estimator not biased downward",   mean(pi0hat) >= (1 - PI1) - 0.03,
   sprintf("%.4f vs true %.4f", mean(pi0hat), 1 - PI1))

SUM <- data.frame(
  method = meth, type = "FDR",
  realised_FDP = round(colMeans(res), 4),
  power = round(colMeans(pwr), 4),
  n_sim = NSIM, m_tests = M, true_pi0 = 1 - PI1,
  mean_pi0_hat = round(mean(pi0hat), 4), stringsAsFactors = FALSE)
write.csv(SUM, file.path(RES, "mtc_validation.csv"), row.names = FALSE)

cat(sprintf("\n%s\n", strrep("=", 72)))
if (fails == 0) {
  cat("ALL CHECKS PASSED -- the corrections in mtc.R behave as advertised.\n")
} else {
  cat(sprintf("%d CHECK(S) FAILED -- do not use these results.\n", fails))
  quit(status = 1)
}
cat("wrote results/mtc_validation.csv\n")
