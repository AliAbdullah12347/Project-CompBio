# ==============================================================================
# mtc.R -- multiple-testing correction library.
#
# WHY THESE ARE HAND-WRITTEN
#
# The machine has limma, edgeR, statmod and MIND. It does NOT have qvalue, IHW,
# fdrtool, ashr, locfdr, sva or swfdr, and CLAUDE.md section 7 says to ask before
# adding a dependency. So every method beyond base R's p.adjust() is implemented
# here from its published definition.
#
# That is only acceptable if the implementations are CHECKED rather than
# trusted. 00_validate_mtc.R exercises each one against an identity that must
# hold exactly, and refuses to proceed if any fails. The identities are stated
# beside each function below.
#
# WHAT IS IMPLEMENTED, AND WHY EACH IS HERE
#
#   Bonferroni   FWER, any dependence. The strictest thing anyone asks for.
#   Holm         FWER, any dependence, uniformly more powerful than Bonferroni.
#                There is never a reason to prefer Bonferroni to Holm; both are
#                reported only because reviewers ask for Bonferroni by name.
#   Hommel       FWER, more powerful than Holm, valid under the same positive
#                dependence BH needs.
#   BH           FDR under independence or positive regression dependence
#                (PRDS). The field standard and our primary.
#   BY           FDR under ARBITRARY dependence. Costs a factor of
#                sum(1/i) ~ 9.9 at m = 12,368. This is the honest bound if one
#                refuses to assume PRDS for co-expressed genes.
#   Storey q     Adaptive BH: multiplies by the estimated share of true nulls,
#                so it is uniformly at least as powerful as BH. Two pi0
#                estimators are provided because they disagree when the p-value
#                histogram is irregular, and reporting one alone would hide that.
#   Permutation  Empirical FDR built from label permutations. The only method
#     FDR        here that respects the ACTUAL correlation between genes rather
#                than assuming something about it.
#   Benjamini-   Hierarchical correction across families (lineages, contrasts).
#     Bogomolov  Answers "I looked at four comparisons, what does that cost me".
#
# NOTHING here chooses a threshold or returns a verdict. Thresholds are a
# reporting choice made in 05_tables.R, and this file deliberately has no
# opinion about 0.05.
# ==============================================================================

## ---- pi0 estimators ---------------------------------------------------------
# Share of hypotheses that are truly null. Both estimators exploit the same
# fact: p-values from true nulls are uniform, so the flat right-hand part of
# the histogram is nearly all null. pi0(lambda) = #{p > lambda} / (m(1-lambda)).
# The two differ only in how they pick lambda, which is the hard part -- small
# lambda is biased upward, large lambda is unbiased but noisy.

pi0_lambda <- function(p, lambda) sum(p > lambda) / (length(p) * (1 - lambda))

# Storey & Tibshirani (2003) PNAS: fit a cubic spline through pi0(lambda) over a
# grid and extrapolate to lambda = 1, where the estimator would be unbiased but
# has infinite variance. The spline trades the two off.
# IDENTITY: with a perfectly uniform p-value vector this must return ~1.
pi0_spline <- function(p, lambda = seq(0.05, 0.95, 0.05)) {
  m <- length(p)
  pl <- vapply(lambda, function(l) pi0_lambda(p, l), numeric(1))
  if (all(!is.finite(pl))) return(1)
  fit <- stats::smooth.spline(lambda, pl, df = 3)
  min(1, max(0, stats::predict(fit, x = 1)$y))
}

# Storey, Taylor & Siegmund (2004) JRSS-B: choose lambda by bootstrap, minimising
# estimated mean squared error against the smallest (least biased) pi0 on the
# grid. More stable than the spline when the histogram is irregular, which is
# exactly the situation in our rare lineages.
pi0_bootstrap <- function(p, lambda = seq(0.05, 0.95, 0.05), B = 100, seed = 481) {
  set.seed(seed)
  m <- length(p)
  pl <- vapply(lambda, function(l) pi0_lambda(p, l), numeric(1))
  minpi0 <- stats::quantile(pl, 0.1, names = FALSE)
  mse <- numeric(length(lambda))
  for (b in seq_len(B)) {
    pb <- sample(p, m, replace = TRUE)
    plb <- vapply(lambda, function(l) pi0_lambda(pb, l), numeric(1))
    mse <- mse + (plb - minpi0)^2
  }
  min(1, max(0, pl[which.min(mse)]))
}

## ---- Storey q-values --------------------------------------------------------
# q(i) = min over j >= i of ( pi0 * m * p(j) / j ), i.e. BH scaled by pi0 and
# made monotone from the top down.
# IDENTITY: qvalue_storey(p, pi0 = 1) must equal p.adjust(p, "BH") EXACTLY.
qvalue_storey <- function(p, pi0 = NULL) {
  if (is.null(pi0)) pi0 <- pi0_spline(p)
  m <- length(p)
  o <- order(p); ro <- order(o)
  ps <- p[o]
  q <- pi0 * m * ps / seq_len(m)
  q <- rev(cummin(rev(q)))          # enforce monotonicity from the largest p down
  pmin(1, q[ro])
}

## ---- the full base-R panel --------------------------------------------------
# Every method base R offers that is defensible here, in one table, so no result
# depends on which correction a reader happens to prefer.
adjust_panel <- function(p) {
  data.frame(
    p_raw       = p,
    bonferroni  = stats::p.adjust(p, "bonferroni"),
    holm        = stats::p.adjust(p, "holm"),
    hommel      = stats::p.adjust(p, "hommel"),
    BH          = stats::p.adjust(p, "BH"),
    BY          = stats::p.adjust(p, "BY"),
    stringsAsFactors = FALSE)
}

## ---- Simes global p-value ---------------------------------------------------
# One p-value summarising a whole family. Valid under PRDS. Used as the
# family-level statistic in the hierarchical procedure below.
simes <- function(p) {
  p <- sort(p[is.finite(p)])
  if (!length(p)) return(NA_real_)
  min(length(p) * p / seq_along(p), 1)
}

## ---- Benjamini-Bogomolov hierarchical correction ---------------------------
# Benjamini & Bogomolov (2014) JRSS-B. The question it answers: we did not test
# one family, we tested K of them, and we only report findings from the families
# that looked interesting. What does that selection cost?
#
#   Stage 1  one Simes p-value per family; BH across the K families at level q.
#            R families are selected.
#   Stage 2  inside each selected family, BH at the REDUCED level q * R / K.
#
# This controls the expected average FDR over the selected families. A family
# not selected in stage 1 reports nothing, which is the point: it is what stops
# four comparisons turning into four free chances.
# IDENTITY: with K = 1 and that family selected, stage 2 reduces to plain BH.
bb_hierarchical <- function(plist, q = 0.05) {
  K <- length(plist)
  fam_p <- vapply(plist, simes, numeric(1))
  fam_adj <- stats::p.adjust(fam_p, "BH")
  selected <- fam_adj <= q
  R <- sum(selected)
  level2 <- if (R > 0) q * R / K else NA_real_
  out <- lapply(seq_len(K), function(i) {
    if (!selected[i]) return(rep(NA_real_, length(plist[[i]])))
    # BH at level2 expressed as an adjusted p-value on the ORIGINAL scale, so a
    # reader can compare it against q directly instead of against level2.
    stats::p.adjust(plist[[i]], "BH") * (q / level2)
  })
  names(out) <- names(plist)
  list(family_p = fam_p, family_adj = fam_adj, selected = selected,
       n_selected = R, K = K, level2 = level2, adjusted = out)
}

## ---- permutation (empirical) FDR -------------------------------------------
# Storey & Tibshirani's permutation plug-in, the same idea SAM uses.
#
#   FDR(t) = pi0 * E_perm[ #{null p <= t} ] / max(1, #{observed p <= t})
#
# The numerator is measured, not assumed. Every permutation leaves the
# expression matrix and therefore every gene-gene correlation intact, so unlike
# BH this makes no independence or PRDS assumption at all. That is the whole
# reason for paying for it.
#
# null_counts: B x length(grid) matrix, counts of null p <= grid[j].
# IDENTITY: if the null behaves exactly uniformly, E[count] = m*t and the result
# must agree with BH at the same threshold to within Monte Carlo error.
perm_fdr <- function(p_obs, grid, null_counts, pi0 = 1) {
  obs <- vapply(grid, function(t) sum(p_obs <= t), numeric(1))
  expected <- colMeans(null_counts)
  fdr <- pi0 * expected / pmax(1, obs)
  fdr <- pmin(1, rev(cummin(rev(fdr))))      # monotone in t
  data.frame(threshold = grid, n_observed = obs,
             n_null_expected = round(expected, 2), fdr = fdr,
             stringsAsFactors = FALSE)
}

# Map each gene's p-value to the empirical FDR of the grid threshold that
# captures it, so a per-gene q-like value can be reported alongside BH.
#
# BOTH ENDS OF THE GRID NEED CARE, and the first version of this function got
# both wrong:
#   * p BELOW the smallest grid point returned 1, i.e. the most significant
#     genes in the study were assigned the worst possible q-value. Since FDR is
#     monotone increasing in t, the right answer is fdr[1] -- conservative and
#     correct.
#   * p ABOVE the largest grid point returned fdr[last], which is a number well
#     under 1 and so credited genes the grid never calibrated. The right answer
#     there is 1: outside the calibrated range we know nothing.
perm_qvalue <- function(p_obs, pf) {
  idx <- findInterval(p_obs, pf$threshold, left.open = FALSE)
  out <- rep(1, length(p_obs))
  below <- p_obs <= pf$threshold[1]
  inside <- idx >= 1 & p_obs <= pf$threshold[nrow(pf)]
  out[inside] <- pf$fdr[idx[inside]]
  out[below] <- pf$fdr[1]
  pmin(1, out)
}

## ---- independent filtering --------------------------------------------------
# Bourgon, Gentleman & Huber (2010) PNAS. Discard the genes least likely to be
# detectable using a statistic that is independent of the p-value under the
# null -- mean expression qualifies -- which reduces m and so raises power for
# everything that survives. DESeq2 does this by default.
#
# The filter must be chosen WITHOUT looking at the direction of the effect; mean
# expression is computed across all samples ignoring group, which satisfies that.
independent_filter <- function(p, filter_stat, alpha = 0.05,
                               quantiles = seq(0, 0.5, 0.05)) {
  res <- do.call(rbind, lapply(quantiles, function(qq) {
    thr <- stats::quantile(filter_stat, qq, names = FALSE)
    keep <- filter_stat >= thr
    n <- sum(stats::p.adjust(p[keep], "BH") < alpha, na.rm = TRUE)
    data.frame(filter_quantile = qq, n_kept = sum(keep), n_rejected = n,
               stringsAsFactors = FALSE)
  }))
  res$is_best <- res$n_rejected == max(res$n_rejected)
  res
}
