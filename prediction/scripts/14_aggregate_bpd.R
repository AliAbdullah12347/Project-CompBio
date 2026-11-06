#!/usr/bin/env Rscript
# Pool the out-of-fold predictions and report the contest.
#
# AUC is computed per repeat, over that repeat's pooled out-of-fold
# predictions -- every subject scored exactly once by a model that never saw
# them. Pooling across all 25 repeats first would average predictions for the
# same subject, which flatters the estimate; averaging the 25 per-repeat AUCs
# keeps each one honest.
#
# The confidence interval bootstraps over SUBJECTS, not over predictions: a
# subject appears 25 times (once per repeat) and those appearances are not
# independent, so resampling rows would understate the interval badly.
#
# Usage: Rscript 14_aggregate_bpd.R [results_dir]

args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args)) args[1] else "results/bpd_cluster"
set.seed(481)
NBOOT <- 2000

fs <- list.files(dir, pattern = "\\.csv$", full.names = TRUE)
d <- do.call(rbind, lapply(fs, read.csv, stringsAsFactors = FALSE))
d <- unique(d)
cat(sprintf("%d files, %d rows\n", length(fs), nrow(d)))

auc <- function(score, lab) {
  if (length(unique(lab)) < 2) return(NA_real_)
  r <- rank(score); n1 <- sum(lab == 1); n0 <- sum(lab == 0)
  (sum(r[lab == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# One AUC per repeat, then the mean. Each repeat is a complete out-of-fold
# pass over the cohort.
per_repeat <- function(x) {
  a <- tapply(seq_len(nrow(x)), x$repeat_no,
              function(i) auc(x$pred[i], x$y[i]))
  as.numeric(a)
}

boot_ci <- function(x) {
  subj <- unique(x$sample_id)
  idx <- split(seq_len(nrow(x)), x$sample_id)
  est <- numeric(NBOOT)
  for (b in seq_len(NBOOT)) {
    pick <- sample(subj, length(subj), replace = TRUE)
    rows <- unlist(idx[pick], use.names = FALSE)
    xb <- x[rows, ]
    est[b] <- mean(per_repeat(xb), na.rm = TRUE)
  }
  quantile(est, c(0.025, 0.975), na.rm = TRUE)
}

real <- d[d$perm == 0, ]
cat("\n")
for (coh in c("groupA", "groupA_plate", "full")) {
  s <- real[real$cohort == coh, ]
  if (!nrow(s)) next
  n_case <- length(unique(s$sample_id[s$y == 1]))
  n_ctrl <- length(unique(s$sample_id[s$y == 0]))
  cat(sprintf("=== %s : %d cases / %d controls ===\n", coh, n_case, n_ctrl))
  cat(sprintf("  %-12s %-5s %8s %8s   %s\n", "features", "model", "AUC", "sd", "95% CI (subject bootstrap)"))
  for (fset in c("covariates", "genes", "genes_cov", "genes_resid")) {
    for (md in c("enet", "rf")) {
      x <- s[s$feature_set == fset & s$model == md, ]
      if (!nrow(x)) next
      a <- per_repeat(x)
      ci <- boot_ci(x)
      cat(sprintf("  %-12s %-5s %8.3f %8.3f   [%.3f, %.3f]\n",
                  fset, md, mean(a, na.rm = TRUE), sd(a, na.rm = TRUE), ci[1], ci[2]))
    }
  }
  cat("\n")
}

## ---- permutation null --------------------------------------------------------
pm <- d[d$perm > 0, ]
if (nrow(pm)) {
  cat("=== permutation null (primary cohort, labels shuffled) ===\n")
  for (fset in unique(pm$feature_set)) {
    x <- pm[pm$feature_set == fset, ]
    a <- tapply(seq_len(nrow(x)), x$perm, function(i) auc(x$pred[i], x$y[i]))
    a <- as.numeric(a)
    obs <- mean(per_repeat(real[real$cohort == "groupA" &
                                real$feature_set == fset & real$model == "enet", ]), na.rm = TRUE)
    p <- (sum(a >= obs, na.rm = TRUE) + 1) / (sum(!is.na(a)) + 1)
    cat(sprintf("  %-12s null mean %.3f  95th pct %.3f  (n=%d)   observed %.3f   p = %.3f\n",
                fset, mean(a, na.rm = TRUE), quantile(a, 0.95, na.rm = TRUE),
                sum(!is.na(a)), obs, p))
  }
  cat("\n  A null that does not centre on 0.5 is the point of running it:\n")
  cat("  small-n CV with thousands of features does not behave like a coin.\n")
}
