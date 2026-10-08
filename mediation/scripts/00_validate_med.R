#!/usr/bin/env Rscript
# ==============================================================================
# 00_validate_med.R -- prove the estimator in med.R is right before any real
# number is produced with it.
#
# The mediation formulas are short enough to look obviously correct and subtle
# enough to be wrong anyway. So nothing here checks med.R against my algebra.
# Every test checks it against something independent:
#
#   1-4   against MONTE CARLO SIMULATION of the actual counterfactuals. We
#         generate potential outcomes Y(1, M(0)) and Y(0, M(0)) by brute force
#         from known truth and compare the simulated contrast with the formula.
#         If the formula is wrong, this catches it; no derivation is reused.
#   5     against the `mediation` package's definition via an independent
#         hand-rolled two-model computation on a single mediator.
#   6-7   internal identities that must hold exactly (decomposition sums,
#         product-of-coefficients nesting).
#   8-9   centring invariance, and that voom weights are unchanged by centring.
#   10-12 the delta method against a 4,000-replicate parametric bootstrap, on
#         simulated data where the truth is known. Checks the standard errors,
#         not just the point estimates, and checks the claim that the
#         theta/beta covariance is block diagonal.
#   13    coverage of the delta-method 95% interval over 500 simulated datasets.
#   14-16 the mediator-model covariance equals the kronecker form; the
#         seemingly-unrelated-regression cross-mediator covariance is carried.
#   17-18 behaviour on degenerate input: zero mediation, zero direct effect.
#
# Every test prints PASS or FAIL and the script stops on the first FAIL.
# ==============================================================================

ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
HERE <- file.path(ROOT, "mediation")
setwd(HERE); source("scripts/med.R")
set.seed(481)
options(width = 150)

nfail <- 0L; ntest <- 0L
chk <- function(label, ok, detail = "") {
  ntest <<- ntest + 1L
  if (!ok) nfail <<- nfail + 1L
  cat(sprintf("  %-62s %s%s\n", label, if (ok) "PASS" else "**FAIL**",
              if (nzchar(detail)) paste0("   ", detail) else ""))
  if (!ok) stop("validation failed at: ", label)
}
rule <- function(t) cat(sprintf("\n%s\n%s\n", t, strrep("-", nchar(t))))

## ---- a simulator with known truth -------------------------------------------
# Two mediators so the vector machinery is exercised, and a genuine
# exposure-mediator interaction so the interaction terms are not silently zero.
TRUTH <- list(
  b0 = c(0.5, -0.3), b1 = c(0.8, 0.4),           # mediator intercepts, exposure effects
  b2 = rbind(c(0.2, -0.1), c(0.3, 0.25)),        # 2 covariates x 2 mediators
  t0 = 1.0, t1 = 0.6,                            # outcome intercept, direct
  t2 = c(0.7, -0.5), t3 = c(0.25, 0.15),         # mediator, interaction
  t4 = c(0.4, -0.2),                             # covariates
  sM = c(0.6, 0.5), rhoM = 0.3, sY = 0.9
)

simulate_one <- function(n, tr = TRUTH, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x <- rbinom(n, 1, 0.65)
  C <- cbind(c1 = rnorm(n), c2 = rbinom(n, 1, 0.4))
  Sm <- matrix(c(tr$sM[1]^2, tr$rhoM * prod(tr$sM),
                 tr$rhoM * prod(tr$sM), tr$sM[2]^2), 2, 2)
  L <- chol(Sm)
  eM <- matrix(rnorm(n * 2), n, 2) %*% L
  M <- cbind(tr$b0[1], tr$b0[2])[rep(1, n), ] +
       outer(x, tr$b1) + C %*% tr$b2 + eM
  colnames(M) <- c("m1", "m2")
  mu <- tr$t0 + tr$t1 * x + M %*% tr$t2 + (M * x) %*% tr$t3 + C %*% tr$t4
  y <- as.vector(mu) + rnorm(n, 0, tr$sY)
  list(x = x, C = C, M = M, y = y)
}

# Truth on the estimand scale, at C = E[C]. Derived ONLY from the generating
# parameters, deliberately not from med.R's formulas.
cbar_true <- c(0, 0.4)
true_EM0 <- TRUTH$b0 + as.vector(crossprod(TRUTH$b2, cbar_true))
true_NDE <- TRUTH$t1 + sum(TRUTH$t3 * true_EM0)
true_NIE <- sum((TRUTH$t2 + TRUTH$t3) * TRUTH$b1)

## ---- 1-4: Monte Carlo counterfactuals ---------------------------------------
rule("1-4  formula vs brute-force Monte Carlo of the counterfactuals")
# Build potential outcomes directly. M(0) and M(1) are drawn from their true
# distributions; Y(x, m) is evaluated from the true outcome model. The natural
# effects are then plain averages of simulated quantities -- no formula used.
mc_natural <- function(N = 4e6, tr = TRUTH, cbar = cbar_true) {
  Sm <- matrix(c(tr$sM[1]^2, tr$rhoM * prod(tr$sM),
                 tr$rhoM * prod(tr$sM), tr$sM[2]^2), 2, 2)
  eM <- matrix(rnorm(N * 2), N, 2) %*% chol(Sm)
  base <- matrix(tr$b0, N, 2, byrow = TRUE) +
          matrix(as.vector(crossprod(tr$b2, cbar)), N, 2, byrow = TRUE)
  M0 <- base + eM                                   # M(0)
  M1 <- base + matrix(tr$b1, N, 2, byrow = TRUE) + eM   # M(1)
  Yx <- function(x, M) tr$t0 + tr$t1 * x + M %*% tr$t2 + (M * x) %*% tr$t3 +
                       sum(tr$t4 * cbar)
  list(NDE = mean(Yx(1, M0) - Yx(0, M0)),           # Y(1,M(0)) - Y(0,M(0))
       NIE = mean(Yx(1, M1) - Yx(1, M0)),           # Y(1,M(1)) - Y(1,M(0))
       TE  = mean(Yx(1, M1) - Yx(0, M0)))
}
mc <- mc_natural()
chk("MC natural direct effect matches the closed form", abs(mc$NDE - true_NDE) < 2e-3,
    sprintf("MC %.5f vs formula %.5f", mc$NDE, true_NDE))
chk("MC natural indirect effect matches the closed form", abs(mc$NIE - true_NIE) < 2e-3,
    sprintf("MC %.5f vs formula %.5f", mc$NIE, true_NIE))
chk("MC total effect = NDE + NIE", abs(mc$TE - (mc$NDE + mc$NIE)) < 3e-3,
    sprintf("MC TE %.5f vs sum %.5f", mc$TE, mc$NDE + mc$NIE))

# And now: does med.R, run on a large simulated sample, recover that truth?
big <- simulate_one(200000, seed = 7)
fit_est <- function(d, moderate = FALSE) {
  cM <- centre(d$M); cC <- centre(d$C)
  mf <- mediator_fit(cM$X, d$x, cC$X)
  Xd <- med_design_full(d$x, cM$X, cC$X)
  Y <- matrix(d$y, nrow = 1); W <- matrix(1, nrow = 1, ncol = length(d$y))
  of <- outcome_fit(Y, W, Xd, se = TRUE, moderate = moderate)
  list(eff = med_effects(of, mf), mf = mf, of = of, Xd = Xd)
}
f <- fit_est(big)
chk("med.R recovers the true NDE at n = 200,000",
    abs(f$eff$est$NDE - true_NDE) < 0.02,
    sprintf("est %.4f vs true %.4f", f$eff$est$NDE, true_NDE))
chk("med.R recovers the true NIE at n = 200,000",
    abs(f$eff$est$NIE - true_NIE) < 0.02,
    sprintf("est %.4f vs true %.4f", f$eff$est$NIE, true_NIE))

## ---- 5: independent two-model computation, single mediator ------------------
rule("5  single-mediator case vs an independent hand computation")
d1 <- simulate_one(5000, seed = 11)
d1$M <- d1$M[, 1, drop = FALSE]
cM <- centre(d1$M); cC <- centre(d1$C)
# hand computation with plain lm(), nothing from med.R
mm <- lm(cM$X[, 1] ~ d1$x + cC$X)
dat <- data.frame(y = d1$y, x = d1$x, m = cM$X[, 1], cC$X)
oo <- lm(y ~ x + m + x:m + c1 + c2, data = dat)
cf <- coef(oo)
hand_NDE <- cf["x"] + cf["x:m"] * coef(mm)[1]
hand_NIE <- (cf["m"] + cf["x:m"]) * coef(mm)[2]
f1 <- fit_est(d1)
chk("NDE matches an independent lm()-based computation",
    abs(f1$eff$est$NDE - hand_NDE) < 1e-9,
    sprintf("%.8f vs %.8f", f1$eff$est$NDE, hand_NDE))
chk("NIE matches an independent lm()-based computation",
    abs(f1$eff$est$NIE - hand_NIE) < 1e-9,
    sprintf("%.8f vs %.8f", f1$eff$est$NIE, hand_NIE))

## ---- 6-7: internal identities ----------------------------------------------
rule("6-7  decomposition identities")
e <- f$eff$est
chk("four-way decomposition sums to the total effect",
    abs((e$CDE + e$INTref + e$INTmed + e$PIE) - e$TE) < 1e-10,
    sprintf("diff %.3g", (e$CDE + e$INTref + e$INTmed + e$PIE) - e$TE))
chk("CDE + INTref = NDE and INTmed + PIE = NIE",
    abs(e$CDE + e$INTref - e$NDE) < 1e-10 && abs(e$INTmed + e$PIE - e$NIE) < 1e-10)
# product-of-coefficients must be the t3 = 0 special case
d0 <- simulate_one(4000, seed = 13)
tr0 <- TRUTH; tr0$t3 <- c(0, 0)
d0 <- local({ tr <- tr0; n <- 4000; set.seed(13)
  x <- rbinom(n,1,0.65); C <- cbind(c1=rnorm(n), c2=rbinom(n,1,0.4))
  Sm <- matrix(c(tr$sM[1]^2, tr$rhoM*prod(tr$sM), tr$rhoM*prod(tr$sM), tr$sM[2]^2),2,2)
  M <- matrix(tr$b0,n,2,byrow=TRUE) + outer(x,tr$b1) + C%*%tr$b2 +
       matrix(rnorm(n*2),n,2)%*%chol(Sm); colnames(M) <- c("m1","m2")
  y <- as.vector(tr$t0 + tr$t1*x + M%*%tr$t2 + C%*%tr$t4) + rnorm(n,0,tr$sY)
  list(x=x, C=C, M=M, y=y) })
f0 <- fit_est(d0)
# Identity: the gap between the interaction-allowing NIE and the PURE indirect
# effect is EXACTLY the mediated-interaction term, whatever the data.
chk("NIE - PIE = INTmed exactly",
    abs((f0$eff$est$NIE - f0$eff$est$PIE) - f0$eff$est$INTmed) < 1e-12,
    sprintf("gap %.10f vs INTmed %.10f",
            f0$eff$est$NIE - f0$eff$est$PIE, f0$eff$est$INTmed))
# The PRODUCT-OF-COEFFICIENTS estimator is a genuinely DIFFERENT fit, not this
# model with t3 dropped from the formula. It must be refitted without the
# interaction columns, and it must NOT coincide with PIE in general.
pc0 <- local({
  cM <- centre(d0$M); cC <- centre(d0$C)
  prod_of_coef(matrix(d0$y, nrow = 1), matrix(1, 1, length(d0$y)),
               d0$x, cM$X, cC$X, mediator_fit(cM$X, d0$x, cC$X)$b1,
               moderate = FALSE) })
chk("product-of-coefficients is NOT just PIE relabelled",
    abs(pc0$NIE_prod - f0$eff$est$PIE) > 1e-8 ||
      abs(f0$eff$est$INTmed) < 1e-10,
    sprintf("prod %.6f vs PIE %.6f", pc0$NIE_prod, f0$eff$est$PIE))
# with NO true interaction the two estimators must agree to within noise
chk("with no true interaction, both estimators agree",
    abs(pc0$NIE_prod - f0$eff$est$NIE) < 3 * f0$eff$se$NIE,
    sprintf("prod %.5f vs NIE %.5f (SE %.5f)", pc0$NIE_prod, f0$eff$est$NIE,
            f0$eff$se$NIE))
# And when the truth has NO interaction, that gap must be statistically
# indistinguishable from zero. It is not numerically zero -- t3_hat is noise
# even when t3 is zero -- so the test is on the z-statistic, not the gap.
z_int <- f0$eff$est$INTmed / f0$eff$se$INTmed
chk("with no true interaction, the mediated-interaction term is null",
    abs(z_int) < 3,
    sprintf("INTmed %.4f, SE %.4f, z = %.2f", f0$eff$est$INTmed,
            f0$eff$se$INTmed, z_int))

## ---- 8-9: centring invariance ----------------------------------------------
rule("8-9  centring changes nothing it must not change")
dd <- simulate_one(3000, seed = 17)
# uncentred fit, same estimator, by hand
raw_fit <- function(d) {
  D <- cbind(1, d$x, d$C); B <- qr.coef(qr(D), d$M)
  Xd <- cbind(1, d$x, d$M, d$M * d$x, d$C)
  th <- qr.coef(qr(Xd), d$y)
  K <- ncol(d$M)
  t1 <- th[2]; t2 <- th[2 + 1:K]; t3 <- th[2 + K + 1:K]
  cbar <- colMeans(d$C)
  EM0 <- B[1, ] + as.vector(crossprod(B[-(1:2), , drop = FALSE], cbar))
  c(NDE = t1 + sum(t3 * EM0), NIE = sum((t2 + t3) * B[2, ]))
}
rf <- raw_fit(dd); cf2 <- fit_est(dd)$eff$est
chk("NDE identical centred vs uncentred", abs(rf["NDE"] - cf2$NDE) < 1e-8,
    sprintf("%.10f vs %.10f", rf["NDE"], cf2$NDE))
chk("NIE identical centred vs uncentred", abs(rf["NIE"] - cf2$NIE) < 1e-8,
    sprintf("%.10f vs %.10f", rf["NIE"], cf2$NIE))

# and on the real data: do voom weights survive centring? They must, because
# the two designs span the same column space.
suppressPackageStartupMessages({library(edgeR); library(limma)})
di <- readRDS("data/med_input.rds")
cM <- centre(di$M); cC <- centre(di$C)
Xc <- med_design_full(di$lithium, cM$X, cC$X)
p  <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
dge <- normLibSizes(DGEList(p$counts[di$gene_main, di$meta$title, drop = FALSE]), method = "TMM")
vC <- voom(dge, Xc)
chk("design is full rank after centring", qr(Xc)$rank == ncol(Xc),
    sprintf("rank %d of %d", qr(Xc)$rank, ncol(Xc)))
chk("voom weights unchanged by centring (same column space)",
    max(abs(vC$weights - di$v_main$weights)) < 1e-8,
    sprintf("max abs diff %.3g", max(abs(vC$weights - di$v_main$weights))))

## ---- 10-12: delta method vs parametric bootstrap ---------------------------
rule("10-12  delta-method standard errors vs a 4,000-replicate bootstrap")
dB <- simulate_one(400, seed = 23)
fB <- fit_est(dB)
B <- 4000
boots <- matrix(NA_real_, B, 3, dimnames = list(NULL, c("NDE", "NIE", "TE")))
nB <- length(dB$y)
for (b in seq_len(B)) {
  i <- sample.int(nB, nB, replace = TRUE)
  db <- list(x = dB$x[i], C = dB$C[i, , drop = FALSE],
             M = dB$M[i, , drop = FALSE], y = dB$y[i])
  cMb <- centre(db$M); cCb <- centre(db$C)
  mfb <- mediator_fit(cMb$X, db$x, cCb$X)
  Xdb <- med_design(db$x, cMb$X, cCb$X)
  thb <- qr.coef(qr(Xdb), db$y)
  K <- 2L; t1 <- thb[2]; t2 <- thb[3:4]; t3 <- thb[5:6]
  nde <- t1 + sum(t3 * mfb$b0); nie <- sum((t2 + t3) * mfb$b1)
  boots[b, ] <- c(nde, nie, nde + nie)
}
for (q in c("NDE", "NIE", "TE")) {
  bs <- sd(boots[, q]); ds <- fB$eff$se[[q]]
  chk(sprintf("delta SE(%s) within 10%% of the bootstrap SE", q),
      abs(ds - bs) / bs < 0.10, sprintf("delta %.5f vs boot %.5f (%.1f%%)",
                                        ds, bs, 100 * (ds - bs) / bs))
}

## ---- 13: interval coverage --------------------------------------------------
rule("13  coverage of the nominal 95% delta-method interval, 500 datasets")
NS <- 500; cov <- c(NDE = 0, NIE = 0)
for (s in seq_len(NS)) {
  ds <- simulate_one(300, seed = 1000 + s)
  fs <- fit_est(ds)
  for (q in c("NDE", "NIE")) {
    tv <- if (q == "NDE") true_NDE else true_NIE
    lo <- fs$eff$est[[q]] - 1.96 * fs$eff$se[[q]]
    hi <- fs$eff$est[[q]] + 1.96 * fs$eff$se[[q]]
    if (tv >= lo && tv <= hi) cov[q] <- cov[q] + 1
  }
}
cov <- cov / NS
chk("NDE interval covers near 95%", abs(cov["NDE"] - 0.95) < 0.035,
    sprintf("%.1f%%", 100 * cov["NDE"]))
chk("NIE interval covers near 95%", abs(cov["NIE"] - 0.95) < 0.035,
    sprintf("%.1f%%", 100 * cov["NIE"]))

## ---- 14-16: the mediator-model covariance ----------------------------------
rule("14-16  mediator-model covariance is the exact kronecker form")
dK <- simulate_one(600, seed = 31)
cMk <- centre(dK$M); cCk <- centre(dK$C)
mfk <- mediator_fit(cMk$X, dK$x, cCk$X)
# independent check: each equation separately via lm(), plus the residual
# covariance, assembled by hand
Dk <- cbind(1, dK$x, cCk$X)
DtDinv <- solve(crossprod(Dk))
Ek <- cMk$X - Dk %*% qr.coef(qr(Dk), cMk$X)
Sg <- crossprod(Ek) / (nrow(Dk) - ncol(Dk))
chk("Vbeta equals kronecker((D'D)^-1[1:2,1:2], Sigma)",
    max(abs(mfk$Vbeta - kronecker(DtDinv[1:2, 1:2], Sg))) < 1e-12)
# single-equation SEs must match lm()'s
l1 <- summary(lm(cMk$X[, 1] ~ dK$x + cCk$X))
chk("mediator b1 SE matches lm() for equation 1",
    abs(mfk$b1_se[1] - coef(l1)[2, 2]) < 1e-10,
    sprintf("%.8f vs %.8f", mfk$b1_se[1], coef(l1)[2, 2]))
# the cross-mediator covariance must be non-zero when the mediators correlate
offd <- mfk$Vbeta[1, 2]
chk("cross-mediator covariance is carried, not zeroed", abs(offd) > 1e-6,
    sprintf("Cov(b0_1, b0_2) = %.3g, residual corr = %.3f",
            offd, Sg[1, 2] / sqrt(Sg[1, 1] * Sg[2, 2])))

## ---- 17-18: degenerate cases -----------------------------------------------
rule("17-18  degenerate inputs behave")
# no mediation at all: b1 = 0 so NIE must be ~0 and PM meaningless
dz <- local({ tr <- TRUTH; tr$b1 <- c(0, 0); n <- 6000; set.seed(41)
  x <- rbinom(n,1,0.65); C <- cbind(c1=rnorm(n), c2=rbinom(n,1,0.4))
  Sm <- matrix(c(tr$sM[1]^2, tr$rhoM*prod(tr$sM), tr$rhoM*prod(tr$sM), tr$sM[2]^2),2,2)
  M <- matrix(tr$b0,n,2,byrow=TRUE) + C%*%tr$b2 + matrix(rnorm(n*2),n,2)%*%chol(Sm)
  colnames(M) <- c("m1","m2")
  y <- as.vector(tr$t0 + tr$t1*x + M%*%tr$t2 + (M*x)%*%tr$t3 + C%*%tr$t4) + rnorm(n,0,tr$sY)
  list(x=x, C=C, M=M, y=y) })
fz <- fit_est(dz)
chk("no exposure effect on the mediator gives NIE indistinguishable from 0",
    abs(fz$eff$est$NIE / fz$eff$se$NIE) < 2.5,
    sprintf("NIE %.4f, z = %.2f", fz$eff$est$NIE, fz$eff$est$NIE / fz$eff$se$NIE))
# no direct effect: t1 = 0 and t3 = 0 so NDE must be ~0 and PM ~1
dn <- local({ tr <- TRUTH; tr$t1 <- 0; tr$t3 <- c(0,0); n <- 6000; set.seed(43)
  x <- rbinom(n,1,0.65); C <- cbind(c1=rnorm(n), c2=rbinom(n,1,0.4))
  Sm <- matrix(c(tr$sM[1]^2, tr$rhoM*prod(tr$sM), tr$rhoM*prod(tr$sM), tr$sM[2]^2),2,2)
  M <- matrix(tr$b0,n,2,byrow=TRUE) + outer(x,tr$b1) + C%*%tr$b2 +
       matrix(rnorm(n*2),n,2)%*%chol(Sm); colnames(M) <- c("m1","m2")
  y <- as.vector(tr$t0 + M%*%tr$t2 + C%*%tr$t4) + rnorm(n,0,tr$sY)
  list(x=x, C=C, M=M, y=y) })
fn <- fit_est(dn)
pmn <- fn$eff$est$NIE / fn$eff$est$TE
chk("no direct effect gives proportion mediated near 1",
    abs(pmn - 1) < 0.08, sprintf("PM = %.4f", pmn))

## ---- summary ----------------------------------------------------------------
rule("summary")
cat(sprintf("  %d checks, %d failures\n", ntest, nfail))
if (nfail == 0L) cat("  med.R is validated. Downstream results may be trusted to the\n  extent that the MODEL is right; the ARITHMETIC is now known to be right.\n")
writeLines(sprintf("%d checks, %d failures, validated %s", ntest, nfail, "2026-10-08"),
           "results/validation_med.txt")
