#!/usr/bin/env Rscript
# ==============================================================================
# norm_sweep_perm.R -- permutation calibration for the normalisation sweep.
#
# The sweep found that skipping normalisation nearly triples the DEG count
# (1382 -> 4013). Two explanations are on the table and they have opposite
# consequences:
#
#   (a) NO NORMALISATION IS GENERICALLY ANTI-CONSERVATIVE. The unmodelled
#       per-sample scale variation inflates the test statistic for everyone, so
#       any exposure would look significant.
#   (b) NO NORMALISATION IS CALIBRATED, but on THIS dataset the scaling factors
#       happen to be correlated with lithium, so omitting them leaves a real
#       global offset that is then read off as thousands of gene effects.
#
# A permutation of the exposure column separates them. Permuting breaks the
# lithium/factor association while leaving the factor variation itself fully
# intact. Under (a) the null DEG counts stay inflated for norm="none"; under (b)
# every method collapses to roughly zero and the real-data gap is attributable
# entirely to the confound.
#
# Covariates stay attached to their own samples -- only the exposure column
# moves -- so the null tested is "lithium carries no signal", not "these
# samples are exchangeable".
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})
source("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation/scripts/common.R")

OUTDIR <- file.path(EXP_HOME, "runs", "norm_sweep")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

NPERM <- 50
set.seed(20261218)
p <- load_prep()
n <- build_design(p, "lithium")$n
stopifnot(n == 226)

# One shared set of permutations across all configurations. Paired, not
# independent: any difference between methods is then a difference in method
# and not a difference in which shuffles they happened to draw.
PERMS <- replicate(NPERM, sample.int(n), simplify = FALSE)

CFG <- list(voom_TMM           = list(norm = "TMM",           method = "voom"),
            voom_none_CPM      = list(norm = "none",          method = "voom"),
            voom_upperquartile = list(norm = "upperquartile", method = "voom"))

out <- list()
for (id in names(CFG)) {
  cf <- CFG[[id]]
  for (k in seq_len(NPERM)) {
    f <- de_fit(p, "lithium", norm = cf$norm, method = cf$method,
                permute_exposure = PERMS[[k]])
    out[[length(out) + 1]] <- data.frame(
      config = id, perm = k,
      deg_fdr05 = n_deg(f$tt, 0.05), deg_fdr10 = n_deg(f$tt, 0.10),
      pi0 = pi0_storey(f$tt$P.Value),
      p_lt_05_frac = mean(f$tt$P.Value < 0.05),
      p_lt_01_frac = mean(f$tt$P.Value < 0.01),
      min_p = min(f$tt$P.Value),
      median_abs_t = median(abs(f$tt$t)),
      # a global offset shows up as a location shift in logFC; under the null
      # there should be none
      median_lfc = median(f$tt$logFC),
      stringsAsFactors = FALSE)
    cat(sprintf("  %-20s perm %2d/%d  DEG=%4d  pi0=%.3f  p<.05=%.4f\n",
                id, k, NPERM, n_deg(f$tt), pi0_storey(f$tt$P.Value),
                mean(f$tt$P.Value < 0.05)))
    flush(stdout())
  }
}
perm <- do.call(rbind, out)
write.csv(perm, file.path(OUTDIR, "permutation_null.csv"), row.names = FALSE)

agg <- do.call(rbind, lapply(split(perm, perm$config), function(d) data.frame(
  config = d$config[1], nperm = nrow(d),
  deg_mean = mean(d$deg_fdr05), deg_median = median(d$deg_fdr05),
  deg_max = max(d$deg_fdr05), deg_p95 = quantile(d$deg_fdr05, .95),
  frac_perms_with_any_deg = mean(d$deg_fdr05 > 0),
  pi0_mean = mean(d$pi0),
  # under a correct test this is 0.05 by construction; above it is inflation
  p_lt_05_mean = mean(d$p_lt_05_frac), p_lt_01_mean = mean(d$p_lt_01_frac),
  median_abs_t_mean = mean(d$median_abs_t),
  abs_median_lfc_mean = mean(abs(d$median_lfc)),
  stringsAsFactors = FALSE)))
write.csv(agg, file.path(OUTDIR, "permutation_null_summary.csv"), row.names = FALSE)
cat("\n-- permutation null summary --\n"); print(format(agg, digits = 4), row.names = FALSE)
cat("\ndone.\n")
