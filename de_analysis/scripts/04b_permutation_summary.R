#!/usr/bin/env Rscript
# ==============================================================================
# 04b_permutation_summary.R -- correct summary of the permutation nulls from 04.
#
# WHY THIS SUPERSEDES THE FIRST SUMMARY
#
# 04 judged the observed pi0 by whether it fell inside the null's central 95%
# interval [q2.5, q97.5]. That is a TWO-SIDED criterion, and it is the wrong
# one. Signal can only push pi0 DOWN -- a contrast with real differential
# expression has a p-value distribution with excess mass near zero, so pi0 is
# smaller than 1. There is no alternative hypothesis under which pi0 is too
# LARGE. The correct test is therefore one-sided, and the correct reference
# point is the null's 5th percentile, not its 2.5th.
#
# Using the two-sided interval made the test more conservative than intended
# and, for the BPD contrast, flipped the reading: pi0 = 0.7091 sits ABOVE the
# 2.5th percentile (0.6901) and so looked "inside the null", but only 3.5% of
# permutations reach a pi0 that low, so the one-sided p is 0.0398.
#
# THE NULL IS BIMODAL, AND THE MEAN IS A USELESS SUMMARY
#
# 195 of 200 permutations produce exactly ZERO differentially expressed genes.
# The remaining five produce as many as 1,683. The mean DEG count (17.1 for
# LI) describes neither mode. Medians and exact tail counts are reported
# instead. The five heavy permutations are those that happened to land close to
# a real axis of variation in the expression matrix; they are legitimate draws
# from the null and are kept, but they are why a mean must not be quoted.
# ==============================================================================

HERE <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/de_analysis"
setwd(HERE)
RES <- "results"
wb <- read.csv(file.path(RES, "wholeblood_summary.csv"), row.names = 1)

OUT <- do.call(rbind, lapply(c("LI", "BPD"), function(cn) {
  d <- read.csv(file.path(RES, sprintf("permnull_%s.csv", cn)))
  B <- nrow(d)
  lab <- sprintf("WB_%s_raw", cn)
  o_deg <- wb[lab, "deg_fdr05"]; o_pi0 <- wb[lab, "pi0"]
  data.frame(
    contrast = cn, n_perm = B,
    obs_deg = o_deg,
    null_deg_median = median(d$deg),
    null_deg_q95 = quantile(d$deg, 0.95, names = FALSE),
    null_deg_max = max(d$deg),
    null_deg_nonzero = sprintf("%d/%d", sum(d$deg > 0), B),
    n_null_ge_obs_deg = sum(d$deg >= o_deg),
    p_deg = round((sum(d$deg >= o_deg) + 1) / (B + 1), 4),
    obs_pi0 = o_pi0,
    null_pi0_median = round(median(d$pi0), 4),
    null_pi0_q05 = round(quantile(d$pi0, 0.05, names = FALSE), 4),
    null_pi0_min = round(min(d$pi0), 4),
    n_null_le_obs_pi0 = sum(d$pi0 <= o_pi0),
    p_pi0_onesided = round((sum(d$pi0 <= o_pi0) + 1) / (B + 1), 4),
    boot_pi0_lo = wb[lab, "pi0_lo"], boot_pi0_hi = wb[lab, "pi0_hi"],
    stringsAsFactors = FALSE) }))

cat("== 04b_permutation_summary ==\n\n")
cat("--- DEG count ---\n")
print(OUT[, c("contrast", "obs_deg", "null_deg_median", "null_deg_q95", "null_deg_max",
              "null_deg_nonzero", "n_null_ge_obs_deg", "p_deg")], row.names = FALSE)
cat("\n--- pi0 (one-sided: signal can only lower pi0) ---\n")
print(OUT[, c("contrast", "obs_pi0", "null_pi0_median", "null_pi0_q05", "null_pi0_min",
              "n_null_le_obs_pi0", "p_pi0_onesided", "boot_pi0_lo", "boot_pi0_hi")], row.names = FALSE)

cat("\n--- how far apart are the two uncertainty estimates for pi0? ---\n")
for (i in seq_len(nrow(OUT))) {
  w <- OUT$boot_pi0_hi[i] - OUT$boot_pi0_lo[i]
  nw <- 1 - OUT$null_pi0_q05[i]
  cat(sprintf("  %-4s bootstrap CI width %.4f | distance from null median(=1) to its 5th pct %.4f | ratio %.1fx\n",
              OUT$contrast[i], w, nw, nw / w))
}

cat("\n--- verdicts ---\n")
for (i in seq_len(nrow(OUT))) {
  cn <- OUT$contrast[i]
  sig <- OUT$p_deg[i] < 0.05 | OUT$p_pi0_onesided[i] < 0.05
  cat(sprintf("  %-4s DEG p=%.4f, pi0 p=%.4f -> signal distinguishable from chance: %s\n",
              cn, OUT$p_deg[i], OUT$p_pi0_onesided[i], ifelse(sig, "YES", "no")))
}
cat("\n  NOTE: for BPD this says the signal is REAL but does not say it is LARGE.\n")
cat("  4 genes is 0.032% of the transcriptome and nothing survives TREAT at |logFC|>=0.1.\n")

write.csv(OUT, file.path(RES, "permutation_vs_observed.csv"), row.names = FALSE)
cat("\nrewrote results/permutation_vs_observed.csv\n")
