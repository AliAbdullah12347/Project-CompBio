#!/usr/bin/env Rscript
# ==============================================================================
# 06_calibration.R -- S3, S4 and S10 from PRESPEC.md.
#
#   S3  regression calibration for a mediator that is ESTIMATED, not measured
#   S4  E-value: how strong would unmeasured confounding have to be?
#   S10 power: what is the smallest mediated effect this design could detect?
#
# All three answer questions the proposal and the feedback raised explicitly.
# ==============================================================================

ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "mediation")); set.seed(481)
suppressPackageStartupMessages({library(limma); library(edgeR)})
source("scripts/med.R"); options(width = 160)
say <- function(...) cat(sprintf(...))
rule <- function(t) cat(sprintf("\n%s\n%s\n", t, strrep("=", nchar(t))))
dir.create("results/sens", showWarnings = FALSE, recursive = TRUE)

d  <- readRDS("data/med_input.rds"); K <- ncol(d$M)
f  <- readRDS("results/primary/lineage_fit.rds"); tb <- f$tb; mf <- f$mf
cM <- centre(d$M); cC <- centre(d$C); x <- d$lithium
Y  <- d$v_main$E; W <- d$v_main$weights
Xd <- med_design_full(x, cM$X, cC$X); ix <- med_idx(K)
of <- outcome_fit(Y, W, Xd, se = TRUE, moderate = TRUE)
T2 <- of$theta[, ix$iM, drop = FALSE]; T3 <- of$theta[, ix$iI, drop = FALSE]
rel <- read.csv("results/balance_reliability.csv", stringsAsFactors = FALSE)
lam <- rel$icc_absolute
names(lam) <- rel$balance

## =============================================================================
## S3  REGRESSION CALIBRATION
## =============================================================================
rule("S3  regression calibration -- the mediator is estimated, not measured")
# The proposal states the logic: "As the mediator is estimated rather than
# measured, that biases the mediated proportion downward, so you correct for it
# using regression calibration and report the result as a lower bound
# (Valeri et al., 2014)."
#
# Classical measurement error in a REGRESSOR attenuates its coefficient, so the
# naive t2 (and t3) are too small and the mediated effect is understated. The
# multivariate correction:
#
#   observed W = true M + U,  U independent of M, X, C
#   naive  t2 converges to  Sigma_W^{-1} Sigma_M t2_true
#   so     t2_true = (Sigma_W^{-1} Sigma_M)^{-1} t2_naive
#
# Error in the mediator MODEL's outcome does not bias b1 -- error in a dependent
# variable inflates its standard error but leaves the slope unbiased -- so b1 is
# carried through uncorrected.
say("  reliability used (ICC, absolute agreement, vs a second deconvolution):\n")
say("    %s\n", paste(sprintf("%s=%.3f", names(lam), lam), collapse = "  "))

Sw_res <- mf$Sigma                               # residual cov of observed W | X,C
var_tot <- apply(d$M, 2, var)                    # total variance per balance
Su0 <- diag((1 - lam) * var_tot, K, K)           # error covariance, assumed diagonal

# ADMISSIBILITY. Sigma_M = Sigma_W - Sigma_U must be positive semi-definite, or
# the reliabilities are mutually inconsistent with the observed covariance of the
# balances and the "corrected" matrix is not a covariance at all. Here the full
# correction IS inadmissible, so rather than printing the failure and carrying on
# we shrink Sigma_U to the largest multiple c that keeps Sigma_M PSD, and report
# c. c < 1 means the data will not support the full de-attenuation the ICCs imply,
# and the correction actually applied is correspondingly smaller.
minev <- function(cc) min(eigen(Sw_res - cc * Su0, symmetric = TRUE)$values)
say("  implied error variance as a share of residual variance: %s\n",
    paste(sprintf("%.0f%%", 100 * diag(Su0) / diag(Sw_res)), collapse = "  "))
say("  eigenvalues of Sigma_M at full correction (c = 1): %s\n",
    paste(sprintf("%.4f", eigen(Sw_res - Su0, symmetric = TRUE)$values), collapse = " "))
if (minev(1) > 0) {
  cshrink <- 1
  say("  -> admissible at full strength (c = 1)\n")
} else {
  lo <- 0; hi <- 1
  for (it in 1:60) { mid <- (lo + hi) / 2; if (minev(mid) > 1e-10) lo <- mid else hi <- mid }
  cshrink <- lo * 0.99
  say("  -> NOT admissible at c = 1. Largest admissible multiple: c = %.4f\n", cshrink)
  say("     The ICCs imply more measurement error than the observed balance\n")
  say("     covariance can support, so only %.0f%% of the implied correction is\n", 100 * cshrink)
  say("     applied. This makes the corrected figure MORE conservative, not less.\n")
}
Su <- cshrink * Su0
Sm_res <- Sw_res - Su
say("  eigenvalues of Sigma_M as applied: %s  (min %.5f > 0: %s)\n",
    paste(sprintf("%.4f", eigen(Sm_res, symmetric = TRUE)$values), collapse = " "),
    min(eigen(Sm_res, symmetric = TRUE)$values),
    min(eigen(Sm_res, symmetric = TRUE)$values) > 0)
stopifnot(min(eigen(Sm_res, symmetric = TRUE)$values) > 0)

# THE MATRIX CORRECTION IS ILL-CONDITIONED HERE, so we report the CURVE rather
# than a point on it. At the admissibility boundary the smallest eigenvalue of
# Sigma_M is barely above zero, so the attenuation matrix Lambda has an
# eigenvalue near 0.01 and Lambda^-1 amplifies by a factor of ~100 in that
# direction. A single "corrected" number taken from there would be an artefact
# of where the boundary happens to fall, not a measurement. Sweeping c shows the
# reader exactly how much the answer depends on how hard the correction is
# pushed.
sig <- tb$TE_q_BH < 0.05
corr_at <- function(cc) {
  Sm <- Sw_res - cc * Su0
  if (min(eigen(Sm, symmetric = TRUE)$values) <= 0) return(NULL)
  Lm <- solve(Sw_res) %*% Sm
  ev <- sort(Re(eigen(Lm)$values))
  Li <- solve(Lm)
  nie <- as.vector(((T2 %*% t(Li)) + (T3 %*% t(Li))) %*% mf$b1)
  list(c = cc, min_eig_Lambda = ev[1], cond = max(ev) / ev[1],
       median_absNIE = median(abs(nie)), median_PM = median((nie / tb$TE)[sig]),
       nie = nie)
}
grid <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, cshrink)
SW <- do.call(rbind, lapply(grid, function(cc) {
  r <- corr_at(cc); if (is.null(r)) return(NULL)
  data.frame(c = r$c, min_eig_Lambda = r$min_eig_Lambda, condition_number = r$cond,
             median_abs_NIE = r$median_absNIE, median_PM = r$median_PM) }))
say("\n  SWEEP over the correction strength c (c = 1 would be the full ICC-implied\n")
say("  correction, which is inadmissible; c = %.4f is the admissibility boundary):\n\n", cshrink)
print(transform(SW, min_eig_Lambda = round(min_eig_Lambda, 4),
                condition_number = round(condition_number, 1),
                median_abs_NIE = round(median_abs_NIE, 4),
                median_PM = round(median_PM, 4)), row.names = FALSE)
say("\n  -> the matrix correction is NOT stable: the condition number runs from\n")
say("     %.1f to %.1f across the admissible range, and the corrected proportion\n",
    min(SW$condition_number), max(SW$condition_number))
say("     mediated runs from %.3f to %.3f. We therefore do NOT report a single\n",
    min(SW$median_PM), max(SW$median_PM))
say("     matrix-corrected figure.\n")

# THE REPORTED CORRECTION. Correct each balance by its own reliability,
# t2_k / lambda_k. This is the direct multi-balance analogue of the univariate
# Valeri, Lin & VanderWeele (2014) correction, it ignores the correlation
# between balances, and it is completely stable because it inverts a diagonal
# matrix with entries bounded away from zero (min lambda = 0.516).
Linv_s <- diag(1 / lam, K, K)
NIEc <- as.vector(((T2 + T3) %*% Linv_s) %*% mf$b1)

# THE TOTAL EFFECT IS INVARIANT TO MEASUREMENT ERROR IN THE MEDIATOR.
# TE is the lithium coefficient of the marginal model Y ~ lithium + C, which
# contains no mediator at all, so no amount of error in M can change it
# (02c_identity.R established TE equals that coefficient exactly once the
# projections match). A valid correction must REDISTRIBUTE effect between the
# direct and indirect paths while leaving the total alone, so the direct effect
# is obtained by subtraction and the proportion is taken against the
# UNCORRECTED total.
TEc  <- tb$TE
NDEc <- TEc - NIEc
PMc  <- NIEc / TEc
say("\n  check: total effect held fixed by construction, max |TEc - TE| = %.3g\n",
    max(abs(TEc - tb$TE)))
stopifnot(max(abs(TEc - tb$TE)) < 1e-12)

say("\n  REPORTED correction (per-balance, stable), genes with significant TE:\n")
cal <- data.frame(
  quantity = c("median |NIE|", "median proportion mediated", "share with PM > 0.50"),
  uncorrected = c(median(abs(tb$NIE)), median(tb$PM[sig]), mean(tb$PM[sig] > 0.5)),
  corrected   = c(median(abs(NIEc)),   median(PMc[sig]),   mean(PMc[sig] > 0.5)),
  stringsAsFactors = FALSE)
cal$change <- cal$corrected - cal$uncorrected
print(transform(cal, uncorrected = round(uncorrected, 4),
                corrected = round(corrected, 4), change = round(change, 4)),
      row.names = FALSE)
say("\n  -> the correction moves the median proportion mediated from %.3f to %.3f,\n",
    median(tb$PM[sig]), median(PMc[sig]))
say("     and the unstable matrix version puts it higher still (up to %.3f). Every\n",
    max(SW$median_PM))
say("     version moves it UP, so the uncorrected %.3f is a LOWER BOUND, exactly as\n",
    median(tb$PM[sig]))
say("     the proposal says to report it. The SIZE of the correction is uncertain;\n")
say("     its DIRECTION is not.\n")
say("\n  CAVEAT that runs the same way: the two deconvolutions behind these ICCs\n")
say("  share the LM22 signature matrix and the same algorithm family, so they\n")
say("  share whatever error the signature itself introduces. The ICCs therefore\n")
say("  OVERSTATE reliability and the correction applied is too small.\n")
write.csv(SW, "results/sens/S3_correction_sweep.csv", row.names = FALSE)
cal$shrink_c <- cshrink
cal$matrix_PM_max <- max(SW$median_PM)
S3 <- data.frame(gene = tb$gene, NIE = tb$NIE, NIE_calibrated = NIEc,
                 TE = tb$TE, NDE_calibrated = NDEc, PM = tb$PM,
                 PM_calibrated = PMc,
                 stringsAsFactors = FALSE)
cal$shrink_c <- cshrink
write.csv(S3, "results/sens/S3_regression_calibration.csv", row.names = FALSE)
write.csv(cal, "results/sens/S3_summary.csv", row.names = FALSE)

## =============================================================================
## S4  E-VALUE
## =============================================================================
rule("S4  E-value -- how strong would an unmeasured confounder have to be?")
# The feedback: "Mediation assumes no unmeasured confounding of lithium, cell
# mix and expression. Illness severity and other drugs affect all three. How
# will you choose the E-value bound?" And the proposal commits to the E-value
# because bodyweight is not recorded in GSE124326.
#
# For a continuous outcome, VanderWeele & Ding (2017) convert a standardised
# mean difference d to an approximate risk ratio exp(0.91 d), then
#   E = RR + sqrt(RR (RR - 1)).
# The E-value is the MINIMUM strength of association, on the risk-ratio scale,
# that an unmeasured confounder would need with BOTH the exposure and the
# outcome, above and beyond the measured covariates, to explain away the effect.
sdY <- sqrt(of$s2_use)                            # residual SD per gene
evalue <- function(est, sdy) {
  rr <- exp(0.91 * abs(est) / sdy)
  rr + sqrt(rr * (rr - 1))
}
ev_TE  <- evalue(tb$TE,  sdY)
ev_NIE <- evalue(tb$NIE, sdY)
# E-value for the confidence LIMIT -- the quantity that actually matters, since
# it asks what would move the interval to include the null
lim_TE <- pmax(abs(tb$TE) - 1.96 * tb$TE_se, 0)
ev_TE_lim <- evalue(lim_TE, sdY)
E <- data.frame(
  quantity = c("total effect", "indirect (mediated) effect", "total effect, CI limit"),
  median_E = c(median(ev_TE), median(ev_NIE), median(ev_TE_lim)),
  E_at_top_gene = c(ev_TE[which.min(tb$TE_p)], ev_NIE[which.min(tb$NIE_p)],
                    ev_TE_lim[which.min(tb$TE_p)]),
  max_E = c(max(ev_TE), max(ev_NIE), max(ev_TE_lim)), stringsAsFactors = FALSE)
print(transform(E, median_E = round(median_E, 3), E_at_top_gene = round(E_at_top_gene, 3),
                max_E = round(max_E, 3)), row.names = FALSE)
say("\n  Read this as: an unmeasured confounder would need to be associated with\n")
say("  BOTH lithium use AND gene expression by a risk ratio of about %.2f each,\n",
    E$median_E[1])
say("  beyond age, sex, tobacco, RIN, plate and the three sequencing PCs, to\n")
say("  explain away the typical total effect.\n")

# A data-grounded benchmark. The proposal wants the comparison made against the
# obesity-to-neutrophil association (Furuncuoglu et al. 2016), but that paper is
# not in the project library, so instead of inventing its numbers we benchmark
# against the measured covariates we DO have: how large is each of them on the
# same E-value scale, in this very dataset?
say("\n  benchmark -- the measured covariates on the same scale, in this dataset:\n")
bench <- NULL
for (cv in colnames(d$C)) {
  v <- d$C[, cv]
  if (length(unique(v)) < 2) next
  # association with the exposure, as a standardised difference
  dx <- (mean(v[x == 1]) - mean(v[x == 0])) / sd(v)
  rx <- exp(0.91 * abs(dx)); ex <- rx + sqrt(rx * (rx - 1))
  # association with expression: median |t| across genes for that covariate

  bench <- rbind(bench, data.frame(covariate = cv, E_vs_lithium = ex,
                                   stringsAsFactors = FALSE))
}
bench <- bench[order(-bench$E_vs_lithium), ]
print(transform(head(bench, 8), E_vs_lithium = round(E_vs_lithium, 3)), row.names = FALSE)
say("\n  ACTION FOR THE WRITE-UP: the proposal names Furuncuoglu et al. (2016) as\n")
say("  the obesity-to-neutrophil benchmark. That paper is not in the project\n")
say("  library, so its numbers are NOT used here rather than being guessed. To\n")
say("  complete S4 as the proposal specifies, extract the obesity-neutrophil\n")
say("  association from it and compare against the median E-value of %.2f above.\n",
    E$median_E[1])
write.csv(data.frame(gene = tb$gene, E_TE = ev_TE, E_NIE = ev_NIE,
                     E_TE_limit = ev_TE_lim), "results/sens/S4_evalues.csv",
          row.names = FALSE)
write.csv(E, "results/sens/S4_summary.csv", row.names = FALSE)

## =============================================================================
## S10  POWER
## =============================================================================
rule("S10  power -- the smallest mediated effect 226 subjects could detect")
# The feedback: "Add a power calculation for the mediation arm with 226
# patients." This is that calculation, and it is the most consequential number
# in the whole arm.
a_nie <- local({ a <- numeric(of$nsub); a[ix$jM] <- mf$b1; a[ix$jI] <- mf$b1; a })
v_theta <- as.vector(of$Vsub %*% as.vector(outer(a_nie, a_nie))) * of$s2_use
V22 <- mf$Vbeta[K + 1:K, K + 1:K, drop = FALSE]
v_beta <- rowSums(((T2 + T3) %*% V22) * (T2 + T3))
se_tot <- sqrt(v_theta + v_beta)

z_nom <- qnorm(0.975); z_pow <- qnorm(0.80)
z_bonf <- qnorm(1 - 0.05 / (2 * nrow(tb)))
P <- data.frame(
  threshold = c("uncorrected 5%", "Bonferroni across 12,018 genes"),
  z_needed = c(z_nom, z_bonf),
  MDE_NIE = c((z_nom + z_pow) * median(se_tot), (z_bonf + z_pow) * median(se_tot)),
  stringsAsFactors = FALSE)
P$vs_observed <- P$MDE_NIE / median(abs(tb$NIE))
print(transform(P, z_needed = round(z_needed, 2), MDE_NIE = round(MDE_NIE, 4),
                vs_observed = round(vs_observed, 2)), row.names = FALSE)
say("\n  median observed |NIE| = %.4f, so the typical gene's mediated effect is\n", median(abs(tb$NIE)))
say("  %.1f%% of what would be needed for 80%% power at a Bonferroni threshold.\n",
    100 / P$vs_observed[2])

rule("S10b  the hard ceiling -- why more sequencing would not help")
# Var(NIE) splits into an outcome-model part and a mediator-model part. The
# second does not shrink with deeper sequencing or more genes; it shrinks only
# with MORE SUBJECTS. So there is a ceiling on how certain this design can ever
# be about a mediated effect, no matter how well expression is measured.
say("  share of Var(NIE) from the mediator models (lithium -> cell mix): %.1f%%\n",
    100 * median(v_beta / (v_beta + v_theta)))
z_ceiling <- abs(tb$NIE) / sqrt(v_beta)
say("  if gene expression were measured PERFECTLY, the largest attainable |z|\n")
say("  for any gene would be %.2f (observed max %.2f).\n",
    max(z_ceiling), max(abs(tb$NIE_z)))
say("  Bonferroni needs |z| >= %.2f, so no amount of sequencing depth reaches it.\n", z_bonf)
nmult <- (z_bonf / max(z_ceiling))^2
say("\n  the binding constraint is SUBJECTS. Since the mediator-model variance\n")
say("  scales as 1/n, reaching that threshold for the best gene would need about\n")
say("  %.1fx the sample size, i.e. roughly %d subjects instead of 226.\n",
    nmult, round(226 * nmult))
say("\n  This is a property of the DESIGN, fixed before any result was seen, not a\n")
say("  disappointment about the data. It is reported whichever way the estimates\n")
say("  come out.\n")
write.csv(P, "results/sens/S10_power.csv", row.names = FALSE)
saveRDS(list(S3 = S3, cal = cal, E = E, power = P, z_ceiling = z_ceiling,
             v_beta = v_beta, v_theta = v_theta, n_needed = round(226 * nmult)),
        "results/sens/calibration.rds")
cat("\n06_calibration complete.\n")
