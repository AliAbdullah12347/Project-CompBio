#!/usr/bin/env Rscript
# ==============================================================================
# 02b_diagnose.R -- find out why the two mediator definitions give NIE estimates
# that correlate at 0.980 but wildly different numbers of significant genes.
#
# Either one of the standard errors is wrong, or something real and important is
# going on. This script decides which, using checks that do not reuse the
# estimator's own algebra:
#
#   A  EXTERNAL: the mediation total effect must equal the lithium coefficient
#      from the plain marginal model Y ~ lithium + C. If it does, the fitting
#      machinery is sound and the problem is confined to the NIE variance.
#   B  decompose Var(NIE) into its outcome-model part and its mediator-model
#      part, and see which dominates in each variant.
#   C  check the theoretical ceiling on the NIE z-statistic: you cannot be more
#      certain of a mediated effect than of the exposure-to-mediator link that
#      carries it.
#   D  compare the delta-method SE against a direct subject bootstrap on a
#      sample of genes, for both variants. This is the arbiter.
# ==============================================================================

ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "mediation")); set.seed(481)
suppressPackageStartupMessages({library(limma); library(edgeR)})
source("scripts/med.R"); options(width = 160)
say <- function(...) cat(sprintf(...))
rule <- function(t) cat(sprintf("\n%s\n%s\n", t, strrep("-", nchar(t))))

VAR <- list(lineage = "data/med_input.rds", proposal = "data/med_input_proposal.rds")

for (vn in names(VAR)) {
  rule(sprintf("VARIANT %s", vn))
  d  <- readRDS(VAR[[vn]]); K <- ncol(d$M)
  f  <- readRDS(sprintf("results/primary/%s_fit.rds", vn)); tb <- f$tb; mf <- f$mf
  cM <- centre(d$M); cC <- centre(d$C); x <- d$lithium
  Y  <- d$v_main$E; W <- d$v_main$weights
  ix <- med_idx(K)

  ## A -- EXTERNAL CHECK: TE vs the marginal model ----------------------------
  Xm <- cbind(`(Intercept)` = 1, X = x, cC$X)
  marg <- matrix(NA_real_, nrow(Y), 2)
  for (g in seq_len(nrow(Y))) {
    w <- W[g, ]; ch <- chol(crossprod(Xm, Xm * w))
    bb <- backsolve(ch, backsolve(ch, crossprod(Xm, Y[g, ] * w), transpose = TRUE))
    marg[g, 1] <- bb[2]
  }
  say("A  total effect vs the marginal model Y ~ lithium + C (no mediator):\n")
  say("     correlation %.6f | max abs diff %.3g | median abs diff %.3g\n",
      cor(tb$TE, marg[, 1]), max(abs(tb$TE - marg[, 1])), median(abs(tb$TE - marg[, 1])))
  say("     -> %s\n", if (cor(tb$TE, marg[, 1]) > 0.999 &&
                          max(abs(tb$TE - marg[, 1])) < 1e-8)
        "EXACT. The fitting machinery is sound." else
        "NOT exact; investigate (expected when the mediator model is misspecified)")

  ## B -- variance decomposition ----------------------------------------------
  th <- f$tb

  # recompute theta to get t2/t3 (not stored in the saved table)
  Xd <- med_design_full(x, cM$X, cC$X)
  of <- outcome_fit(Y, W, Xd, se = TRUE, moderate = TRUE)
  T2 <- of$theta[, ix$iM, drop = FALSE]; T3 <- of$theta[, ix$iI, drop = FALSE]
  S  <- T2 + T3
  Vb <- mf$Vbeta; V22 <- Vb[K + 1:K, K + 1:K, drop = FALSE]
  mkA <- function(m = numeric(K), i = numeric(K), xx = 0) {
    a <- numeric(of$nsub); a[ix$jX] <- xx; a[ix$jM] <- m; a[ix$jI] <- i; a }
  a_nie <- mkA(m = mf$b1, i = mf$b1)
  v_theta <- as.vector(of$Vsub %*% as.vector(outer(a_nie, a_nie))) * of$s2_use
  v_beta  <- rowSums((S %*% V22) * S)
  say("\nB  Var(NIE) decomposition, median across genes:\n")
  say("     outcome-model part  %.3e  (%.1f%%)\n", median(v_theta),
      100 * median(v_theta) / median(v_theta + v_beta))
  say("     mediator-model part %.3e  (%.1f%%)\n", median(v_beta),
      100 * median(v_beta) / median(v_theta + v_beta))
  say("     median SE(NIE) %.5f | median |NIE| %.5f | median |z| %.3f | max |z| %.3f\n",
      median(sqrt(v_theta + v_beta)), median(abs(tb$NIE)), median(abs(tb$NIE_z)),
      max(abs(tb$NIE_z)))
  say("     min NIE p %.3e | BH threshold for 1 discovery p <= %.3e\n",
      min(tb$NIE_p), 0.05 / nrow(tb))
  say("     genes at BH 5%%: %d\n", sum(tb$NIE_q_BH < 0.05))

  ## C -- the ceiling ----------------------------------------------------------
  tstat <- mf$b1 / mf$b1_se
  say("\nC  mediator-model t-statistics (lithium -> each balance):\n")
  say("     %s\n", paste(sprintf("%s t=%.2f", colnames(d$M), tstat), collapse = "  "))
  say("     root-sum-square of the t's (loose ceiling on |z(NIE)|): %.2f\n",
      sqrt(sum(tstat^2)))
  say("     observed max |z(NIE)| %.2f -> ceiling %s\n", max(abs(tb$NIE_z)),
      if (max(abs(tb$NIE_z)) <= sqrt(sum(tstat^2)) * 1.05) "RESPECTED" else "**VIOLATED**")

  ## D -- the arbiter: a direct subject bootstrap on a sample of genes ---------
  set.seed(481)
  gs <- sort(sample.int(nrow(Y), 40))
  B <- 400; n <- ncol(Y)
  bn <- matrix(NA_real_, B, length(gs))
  for (b in seq_len(B)) {
    i <- sample.int(n, n, replace = TRUE)
    cMb <- centre(d$M[i, , drop = FALSE]); cCb <- centre(d$C[i, , drop = FALSE])
    mfb <- mediator_fit(cMb$X, x[i], cCb$X)
    Xb  <- med_design(x[i], cMb$X, cCb$X)
    for (jj in seq_along(gs)) {
      g <- gs[jj]; w <- W[g, i]
      ch <- chol(crossprod(Xb, Xb * w))
      bb <- backsolve(ch, backsolve(ch, crossprod(Xb, Y[g, i] * w), transpose = TRUE))
      bn[b, jj] <- sum((bb[ix$iM] + bb[ix$iI]) * mfb$b1)
    }
  }
  bse <- apply(bn, 2, sd); dse <- tb$NIE_se[gs]
  say("\nD  delta-method SE vs a %d-replicate subject bootstrap, %d sampled genes:\n", B, length(gs))
  say("     median ratio delta/bootstrap  %.3f\n", median(dse / bse))
  say("     range of that ratio           %.3f to %.3f\n", min(dse / bse), max(dse / bse))
  say("     -> %s\n", if (abs(median(dse / bse) - 1) < 0.15)
        "delta method AGREES with the bootstrap" else
        "**delta method DISAGREES with the bootstrap -- the SE is wrong**")
}
