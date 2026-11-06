#!/usr/bin/env Rscript
# Tidy summary of every cohort x feature set x model, written once so the
# figures do not each re-run a 2,000-sample bootstrap.
args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args)) args[1] else "results/bpd_cluster"
out <- if (length(args) > 1) args[2] else "results/summary_auc.csv"
set.seed(481); NBOOT <- 1000
d <- unique(do.call(rbind, lapply(list.files(dir, "[.]csv$", full.names = TRUE),
                                  read.csv, stringsAsFactors = FALSE)))
auc <- function(s, l) { if (length(unique(l)) < 2) return(NA_real_)
  r <- rank(s); n1 <- sum(l == 1); n0 <- sum(l == 0)
  (sum(r[l == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
mean_auc <- function(x) mean(as.numeric(tapply(seq_len(nrow(x)), x$repeat_no,
                              function(i) auc(x$pred[i], x$y[i]))), na.rm = TRUE)
real <- d[d$perm == 0, ]; pm <- d[d$perm > 0, ]
rows <- list(); n <- 0L
for (coh in unique(real$cohort)) for (fs in unique(real$feature_set)) for (md in unique(real$model)) {
  x <- real[real$cohort == coh & real$feature_set == fs & real$model == md, ]
  if (!nrow(x)) next
  subj <- unique(x$sample_id); idx <- split(seq_len(nrow(x)), x$sample_id)
  bs <- replicate(NBOOT, { p <- sample(subj, length(subj), TRUE)
                           mean_auc(x[unlist(idx[p], use.names = FALSE), ]) })
  ci <- quantile(bs, c(0.025, 0.975), na.rm = TRUE)
  # Null from the permutation runs, where they exist for this cohort/feature set
  q <- pm[pm$cohort == coh & pm$feature_set == fs, ]
  nullmean <- if (nrow(q)) mean(as.numeric(tapply(seq_len(nrow(q)), q$perm,
                      function(i) auc(q$pred[i], q$y[i]))), na.rm = TRUE) else NA_real_
  null95 <- if (nrow(q)) quantile(as.numeric(tapply(seq_len(nrow(q)), q$perm,
                      function(i) auc(q$pred[i], q$y[i]))), 0.95, na.rm = TRUE) else NA_real_
  n <- n + 1L
  rows[[n]] <- data.frame(cohort = coh, feature_set = fs, model = md,
    n_case = length(unique(x$sample_id[x$y == 1])),
    n_ctrl = length(unique(x$sample_id[x$y == 0])),
    auc = mean_auc(x), lo = ci[1], hi = ci[2],
    null_mean = nullmean, null_p95 = null95, row.names = NULL)
}
res <- do.call(rbind, rows)
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
write.csv(res, out, row.names = FALSE)
cat(sprintf("wrote %s (%d rows)\n", out, nrow(res)))
print(res[res$model == "enet", c("cohort","feature_set","auc","lo","hi","null_mean")], digits = 3)
