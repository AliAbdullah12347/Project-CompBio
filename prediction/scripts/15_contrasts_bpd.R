#!/usr/bin/env Rscript
# The pre-declared test: does adding expression beat demographics alone?
#
# Paired bootstrap over SUBJECTS -- the same resampled subjects are used for
# both feature sets, so the difference is measured on the same people and the
# interval reflects uncertainty in the comparison rather than in the two AUCs
# separately.
args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args)) args[1] else "results/bpd_cluster"
set.seed(481); NBOOT <- 2000
d <- unique(do.call(rbind, lapply(list.files(dir, "[.]csv$", full.names = TRUE),
                                  read.csv, stringsAsFactors = FALSE)))
d <- d[d$perm == 0, ]
auc <- function(s, l) { if (length(unique(l)) < 2) return(NA_real_)
  r <- rank(s); n1 <- sum(l == 1); n0 <- sum(l == 0)
  (sum(r[l == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
mean_auc <- function(x) mean(as.numeric(tapply(seq_len(nrow(x)), x$repeat_no,
                              function(i) auc(x$pred[i], x$y[i]))), na.rm = TRUE)

for (coh in c("groupA", "groupA_plate", "full")) {
  cat(sprintf("=== %s ===\n", coh))
  for (md in c("enet", "rf")) {
    for (pair in list(c("genes_cov", "covariates"), c("genes", "covariates"),
                      c("genes_resid", "covariates"))) {
      A <- d[d$cohort == coh & d$model == md & d$feature_set == pair[1], ]
      B <- d[d$cohort == coh & d$model == md & d$feature_set == pair[2], ]
      if (!nrow(A) || !nrow(B)) next
      subj <- unique(A$sample_id)
      ia <- split(seq_len(nrow(A)), A$sample_id); ib <- split(seq_len(nrow(B)), B$sample_id)
      obs <- mean_auc(A) - mean_auc(B)
      bs <- numeric(NBOOT)
      for (b in seq_len(NBOOT)) {
        pick <- sample(subj, length(subj), replace = TRUE)
        bs[b] <- mean_auc(A[unlist(ia[pick], use.names = FALSE), ]) -
                 mean_auc(B[unlist(ib[pick], use.names = FALSE), ])
      }
      ci <- quantile(bs, c(0.025, 0.975), na.rm = TRUE)
      verdict <- if (ci[1] > 0) "expression ADDS" else if (ci[2] < 0) "expression HURTS" else "no added value"
      cat(sprintf("  %-5s %-12s vs covariates  %+0.3f  [%+.3f, %+.3f]  %s\n",
                  md, pair[1], obs, ci[1], ci[2], verdict))
    }
  }
  cat("\n")
}
