# ==============================================================================
# med.R -- the mediation estimator. One implementation, used by every script in
# this arm, so there is exactly one place for it to be wrong and exactly one
# place that gets tested (00_validate_med.R).
#
# WHAT IT COMPUTES
#
# Natural direct and natural indirect effects for a VECTOR-valued mediator,
# allowing exposure-mediator interaction. These are Pearl's (2001) natural
# effects, after Robins & Greenland (1992), in the regression form of
# VanderWeele & Vansteelandt (2009, 2010), extended to multiple mediators by
# VanderWeele & Vansteelandt (2014).
#
#   mediator models, k = 1..K, sharing one design D = [1, X, C]:
#       M_k = b0_k + b1_k X + b2_k' C + e_k
#
#   outcome model, per gene, design Xd = [1, X, M, X*M, C]:
#       Y = t0 + t1 X + sum_k t2_k M_k + sum_k t3_k X M_k + t4' C + e
#
#   NDE = t1 + sum_k t3_k b0_k
#   NIE = sum_k (t2_k + t3_k) b1_k
#   TE  = NDE + NIE
#   PM  = NIE / TE
#
# and the four-way decomposition of VanderWeele (2014), TE = CDE + INTref +
# INTmed + PIE, which separates the part of the total effect that needs neither
# mediation nor interaction (CDE), interaction alone (INTref), both (INTmed),
# and mediation alone (PIE).
#
# M AND C ARE CENTRED. With M and C centred at their sample means, the mediator
# model intercept b0_k IS E[M_k | X=0, C=Cbar], which is the quantity the NDE
# formula needs. Centring provably leaves NDE, NIE and TE unchanged -- the
# design's column space is identical -- and it removes the near-collinearity
# between X and X*M that uncentred ILR values (means 1.88, 0.24, 1.76, 1.37)
# would otherwise create. 00_validate_med.R checks the invariance numerically
# rather than trusting this paragraph.
#
# WHY THE DELTA-METHOD VARIANCE IS BLOCK-DIAGONAL, AND WHY THAT IS EXACT
#
# The estimates come from two separate regressions, so one might expect a cross
# term. There is none, asymptotically. The outcome model's estimating equation
# is driven by e_Y = Y - E[Y | X, M, C], the mediator models' by
# e_M = M - E[M | X, C]. Since E[e_Y | X, M, C] = 0 and e_M is a function of
# (M, X, C), iterated expectation gives E[e_Y e_M f(X,M,C)] = 0. The two score
# functions are uncorrelated, so the joint covariance is block diagonal:
#
#       Var(g) = a' V_theta a  +  b' V_beta b
#
# with no cross term. This is not an approximation we are making for
# convenience; it follows from the outcome model being correctly specified.
#
# THE MEDIATOR-MODEL COVARIANCE IS EXACT TOO
#
# The K mediator models share one design matrix, so this is a seemingly
# unrelated regression in which OLS equation-by-equation is fully efficient, and
#
#       Cov(beta_k, beta_l) = Sigma[k,l] * (D'D)^{-1}
#
# with Sigma the residual covariance across mediators. Stacking the (b0, b1)
# pairs in the order (b0_1..b0_K, b1_1..b1_K) makes the 2K x 2K block exactly
#
#       V_beta = kronecker( (D'D)^{-1}[1:2, 1:2] , Sigma )
#
# which is what mediator_fit() returns. Nothing is assumed independent that is
# not independent -- the across-mediator correlation is carried in full.
# ==============================================================================

suppressPackageStartupMessages(library(limma))

## ---- small helpers ----------------------------------------------------------

# Centre columns and remember the centres, so a test set or a bootstrap
# replicate can be put on the same scale if ever needed.
centre <- function(X) {
  mu <- colMeans(X)
  list(X = sweep(X, 2, mu, "-"), mu = mu)
}

# Column layout of the outcome design. Fixed here so no script can disagree
# with another about which column is which.
#   1          intercept
#   2          exposure X
#   2+(1..K)   mediators
#   2+K+(1..K) exposure x mediator
#   rest       covariates
med_design <- function(x, Mc, Cc) {
  K <- ncol(Mc)
  Xd <- cbind(`(Intercept)` = 1, X = x, Mc, Mc * x, Cc)
  colnames(Xd)[(2 + K + 1):(2 + 2 * K)] <- paste0("X_x_", colnames(Mc))
  Xd
}

med_idx <- function(K) list(
  iX = 2L, iM = 2L + seq_len(K), iI = 2L + K + seq_len(K),
  # positions inside the extracted (1 + 2K) sub-block, in the same order
  jX = 1L, jM = 1L + seq_len(K), jI = 1L + K + seq_len(K)
)

## ---- stage 1: the mediator models -------------------------------------------
# M is passed ALREADY CENTRED; C likewise. Returns b0, b1 and the exact 2K x 2K
# covariance of (b0, b1).
mediator_fit <- function(Mc, x, Cc) {
  K <- ncol(Mc); n <- nrow(Mc)
  D <- cbind(`(Intercept)` = 1, X = x, Cc)
  qd <- qr(D)
  if (qd$rank < ncol(D)) stop("mediator design is rank deficient")
  B <- qr.coef(qd, Mc)                       # (2 + pC) x K
  E <- Mc - D %*% B
  dfD <- n - ncol(D)
  Sigma <- crossprod(E) / dfD                # K x K residual covariance
  DtDinv <- chol2inv(chol(crossprod(D)))
  Vbeta <- kronecker(DtDinv[1:2, 1:2, drop = FALSE], Sigma)
  dimnames(Vbeta) <- list(
    c(paste0("b0_", colnames(Mc)), paste0("b1_", colnames(Mc))),
    c(paste0("b0_", colnames(Mc)), paste0("b1_", colnames(Mc))))
  list(b0 = unname(B[1, ]), b1 = unname(B[2, ]), B = B, Sigma = Sigma,
       Vbeta = Vbeta, df = dfD, resid = E, D = D,
       # per-balance Wald test of the exposure effect on each mediator
       b1_se = sqrt(DtDinv[2, 2] * diag(Sigma)),
       names = colnames(Mc))
}

## ---- stage 2: the outcome models --------------------------------------------
# Weighted least squares, one fit per gene, with voom precision weights.
#
# `se = FALSE` is the fast path used by the bootstrap and the permutation null,
# where only point estimates are needed: it skips the 9x9 inverse entirely.
#
# `moderate = TRUE` applies limma's empirical-Bayes variance shrinkage
# (squeezeVar, the engine inside eBayes) to the per-gene residual variance, for
# consistency with the differential-expression arm, which used limma-voom with
# eBayes throughout. Mediation on an unmoderated variance scale would not be
# comparable with the DE results it is meant to explain.
outcome_fit <- function(Y, W, Xd, se = TRUE, moderate = TRUE) {
  G <- nrow(Y); n <- ncol(Y); p <- ncol(Xd)
  stopifnot(ncol(W) == n, nrow(W) == G, nrow(Xd) == n)
  theta <- matrix(NA_real_, G, p, dimnames = list(rownames(Y), colnames(Xd)))
  s2 <- numeric(G)
  nsub <- NULL; Vsub <- NULL
  if (se) {
    nsub <- length(attr(Xd, "sub"))
    Vsub <- matrix(NA_real_, G, nsub * nsub)
  }
  sub <- attr(Xd, "sub")
  for (g in seq_len(G)) {
    w <- W[g, ]
    XtWX <- crossprod(Xd, Xd * w)
    ch <- chol(XtWX)
    b <- backsolve(ch, backsolve(ch, crossprod(Xd, Y[g, ] * w), transpose = TRUE))
    theta[g, ] <- b
    r <- Y[g, ] - Xd %*% b
    s2[g] <- sum(w * r^2) / (n - p)
    if (se) {
      inv <- chol2inv(ch)
      Vsub[g, ] <- as.vector(inv[sub, sub, drop = FALSE])
    }
  }
  out <- list(theta = theta, s2 = s2, df = n - p, Vsub = Vsub, nsub = nsub,
              sub = sub, p = p, n = n)
  if (moderate) {
    sq <- limma::squeezeVar(s2, df = n - p)
    out$s2_use <- sq$var.post
    out$df_prior <- sq$df.prior
    out$df_total <- (n - p) + ifelse(is.finite(sq$df.prior), sq$df.prior, 0)
  } else {
    out$s2_use <- s2; out$df_prior <- 0; out$df_total <- n - p
  }
  out
}

# Attach the sub-block index to a design so outcome_fit knows what to keep.
med_design_full <- function(x, Mc, Cc) {
  Xd <- med_design(x, Mc, Cc)
  K <- ncol(Mc); ix <- med_idx(K)
  attr(Xd, "sub") <- c(ix$iX, ix$iM, ix$iI)
  Xd
}

## ---- the estimands, with delta-method standard errors ------------------------
# Every gradient with respect to theta is CONSTANT across genes (it depends only
# on b0 and b1, which come from the mediator models). Every gradient with
# respect to beta varies by gene (it depends on that gene's t2, t3). Both are
# therefore vectorised: one 81-length matrix-vector product per estimand for the
# theta part, one quadratic form per estimand for the beta part.
med_effects <- function(ofit, mfit, want_se = TRUE) {
  K <- length(mfit$b1); ix <- med_idx(K)
  th <- ofit$theta
  t1 <- th[, ix$iX]
  T2 <- th[, ix$iM, drop = FALSE]
  T3 <- th[, ix$iI, drop = FALSE]
  b0 <- mfit$b0; b1 <- mfit$b1
  S  <- T2 + T3                                   # (t2_k + t3_k), per gene

  est <- list(
    NDE    = as.vector(t1 + T3 %*% b0),
    NIE    = as.vector(S  %*% b1),
    CDE    = as.vector(t1),
    INTref = as.vector(T3 %*% b0),
    INTmed = as.vector(T3 %*% b1),
    PIE    = as.vector(T2 %*% b1)
  )
  est$TE <- est$NDE + est$NIE
  # NOTE: the product-of-coefficients estimator is NOT computed here. Forcing t3
  # to zero in the FORMULA is not the same as removing the interaction columns
  # from the DESIGN -- t2 itself changes when they are dropped, because t2 in
  # this model is the mediator slope at X = 0 whereas the Baron-Kenny / Sobel
  # estimator uses the pooled slope. Taking T2 %*% b1 here would simply
  # reproduce PIE under a second name and make any "cost of assuming no
  # interaction" tautological. prod_of_coef() below refits properly.

  if (!want_se) return(est)

  nsub <- ofit$nsub
  # a' V_theta a  =  s2 * (a' Sinv a), with a constant across genes
  qform_theta <- function(a) {
    as.vector(ofit$Vsub %*% as.vector(outer(a, a))) * ofit$s2_use
  }
  mkA <- function(x = 0, m = numeric(K), i = numeric(K)) {
    a <- numeric(nsub); a[ix$jX] <- x; a[ix$jM] <- m; a[ix$jI] <- i; a
  }
  Vb <- mfit$Vbeta
  V11 <- Vb[1:K, 1:K, drop = FALSE]                 # b0 block
  V22 <- Vb[K + 1:K, K + 1:K, drop = FALSE]         # b1 block
  qform_beta <- function(Bmat, V) rowSums((Bmat %*% V) * Bmat)

  se <- list(
    NDE    = sqrt(qform_theta(mkA(x = 1, i = b0))      + qform_beta(T3, V11)),
    NIE    = sqrt(qform_theta(mkA(m = b1, i = b1))     + qform_beta(S,  V22)),
    CDE    = sqrt(qform_theta(mkA(x = 1))),
    INTref = sqrt(qform_theta(mkA(i = b0))             + qform_beta(T3, V11)),
    INTmed = sqrt(qform_theta(mkA(i = b1))             + qform_beta(T3, V22)),
    PIE    = sqrt(qform_theta(mkA(m = b1))             + qform_beta(T2, V22))
  )
  # TE needs the COMBINED gradient, not the sum of the two variances, because
  # NDE and NIE are correlated through theta.
  se$TE <- sqrt(qform_theta(mkA(x = 1, m = b1, i = b0 + b1)) +
                qform_beta(cbind(T3, S), Vb))

  list(est = est, se = se)
}

# Wrap estimates and SEs into a tidy per-gene table with z, p and CI.
med_table <- function(eff, genes) {
  nm <- names(eff$est)
  out <- data.frame(gene = genes, stringsAsFactors = FALSE)
  for (k in nm) {
    e <- eff$est[[k]]; s <- eff$se[[k]]
    out[[paste0(k)]]        <- e
    out[[paste0(k, "_se")]] <- s
    z <- e / s
    out[[paste0(k, "_z")]]  <- z
    out[[paste0(k, "_p")]]  <- 2 * pnorm(-abs(z))
  }
  out$PM <- out$NIE / out$TE
  out
}

## ---- joint test of the exposure-mediator interaction ------------------------
# H0: t3_1 = ... = t3_K = 0, per gene. A moderated F on K df, computed from the
# limma fit so the variance shrinkage matches everything else in this arm.
#
# This exists because the feedback on the pre-proposal said, of the proposal's
# original product-of-coefficients plan: "product of coefficients assumes no
# lithium-by-cell-mix interaction. Test that interaction before reporting a
# proportion mediated." This is that test.
interaction_test <- function(Y, W, Xd, K) {
  ix <- med_idx(K)
  fit <- limma::lmFit(Y, Xd, weights = W)
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = ix$iI, number = nrow(Y), sort.by = "none")
  data.frame(gene = rownames(Y), F = tt$F, p = tt$P.Value,
             stringsAsFactors = FALSE)
}

## ---- fast path for bootstrap / permutation ----------------------------------
# Point estimates only. Returns NDE, NIE, TE per gene, nothing else, because
# the resampling supplies the uncertainty and the delta method is not needed
# inside the loop.
med_point <- function(Y, W, x, Mc, Cc) {
  Xd <- med_design(x, Mc, Cc)
  G <- nrow(Y); n <- ncol(Y); p <- ncol(Xd); K <- ncol(Mc)
  ix <- med_idx(K)
  mfit <- mediator_fit(Mc, x, Cc)
  th <- matrix(NA_real_, G, p)
  for (g in seq_len(G)) {
    w <- W[g, ]
    ch <- chol(crossprod(Xd, Xd * w))
    th[g, ] <- backsolve(ch, backsolve(ch, crossprod(Xd, Y[g, ] * w),
                                       transpose = TRUE))
  }
  t1 <- th[, ix$iX]; T2 <- th[, ix$iM, drop = FALSE]; T3 <- th[, ix$iI, drop = FALSE]
  NDE <- as.vector(t1 + T3 %*% mfit$b0)
  NIE <- as.vector((T2 + T3) %*% mfit$b1)
  list(NDE = NDE, NIE = NIE, TE = NDE + NIE,
       b1 = mfit$b1, b0 = mfit$b0)
}

## ---- the product-of-coefficients estimator, properly refitted ---------------
# The proposal's original flowchart specified product-of-coefficients
# (Baron & Kenny 1986; Sobel 1982), which assumes no exposure-mediator
# interaction. Prof. Ay's feedback required testing that assumption before
# reporting any proportion mediated, so the primary estimator allows the
# interaction. To report "both estimators side by side", as PRESPEC.md S8
# promises, the no-interaction estimator must be fitted from its OWN design:
#
#     Y = t0 + t1 X + sum_k t2_k M_k + t4' C + e        (no X*M columns)
#     NIE_prod = sum_k t2_k b1_k        NDE_prod = t1
#
# For a binary exposure the no-interaction slope is a within-group weighted
# average of the two group-specific slopes, so t2_ni = t2 + w1 t3 with w1 the
# share of exposed subjects -- which is why dropping t3 from the formula of the
# INTERACTION model gives a different (and wrong) answer.
prod_of_coef <- function(Y, W, x, Mc, Cc, b1, moderate = TRUE) {
  K <- ncol(Mc)
  Xd0 <- cbind(`(Intercept)` = 1, X = x, Mc, Cc)
  attr(Xd0, "sub") <- c(2L, 2L + seq_len(K))
  of0 <- outcome_fit(Y, W, Xd0, se = TRUE, moderate = moderate)
  T2 <- of0$theta[, 2L + seq_len(K), drop = FALSE]
  nsub <- of0$nsub
  qt <- function(a) as.vector(of0$Vsub %*% as.vector(outer(a, a))) * of0$s2_use
  aN <- numeric(nsub); aN[1L + seq_len(K)] <- b1
  aD <- numeric(nsub); aD[1L] <- 1
  list(NIE_prod = as.vector(T2 %*% b1), NDE_prod = as.vector(of0$theta[, 2L]),
       NIE_prod_se_theta = sqrt(qt(aN)), NDE_prod_se = sqrt(qt(aD)), T2_ni = T2)
}
