#!/usr/bin/env Rscript
# Interim analysis for perm-031 — works with however many draws are available.
# Run during or after the run to monitor progress.

run_dir  <- "network/runs/perm-031"
base_dir <- "network/runs/base-001"

done_k <- vapply(1:20, function(i)
  file.exists(file.path(run_dir, sprintf("perm_draw_%03d.rds", i))), logical(1))
n_done <- sum(done_k)
cat(sprintf("Permuted draws available: %d / 20\n", n_done))
if (n_done == 0) { cat("No draws yet.\n"); quit(save="no", status=0) }

# Load null draws
null_all <- do.call(rbind, lapply(which(done_k), function(i)
  readRDS(file.path(run_dir, sprintf("perm_draw_%03d.rds", i)))))

# Load observed (base-001 bp_lith)
obs_csv <- file.path(base_dir, "results_lith_alldraws.csv")
if (!file.exists(obs_csv)) { cat("base-001 lith CSV missing\n"); quit(save="no", status=1) }
obs_lith <- read.csv(obs_csv, stringsAsFactors=FALSE)
obs_mean <- tapply(obs_lith$Zsummary, obs_lith$module, mean, na.rm=TRUE)

cat("\n=== Null distribution (gold module) ===\n")
gold_null <- null_all[null_all$module == "gold", ]
cat(sprintf("null gold: mean=%.2f  sd=%.2f  [%.2f, %.2f]  (n=%d draws)\n",
            mean(gold_null$Zsummary), sd(gold_null$Zsummary),
            min(gold_null$Zsummary), max(gold_null$Zsummary), nrow(gold_null)))
cat(sprintf("observed gold (bp_lith mean 20 draws): %.2f\n", obs_mean["gold"]))
cat(sprintf("gold observed < null mean: %s (specific reduction in bp_lith)\n",
            obs_mean["gold"] < mean(gold_null$Zsummary)))
gold_emp_p <- mean(gold_null$Zsummary <= obs_mean["gold"])
cat(sprintf("gold empirical p (lower-tail, K=%d): %.3f\n", n_done, gold_emp_p))

cat("\n=== Observed vs null per module (preliminary, K=", n_done, ") ===\n", sep="")
modules <- sort(unique(null_all$module)); modules <- modules[modules != "gold"]
direction <- sapply(modules, function(m) {
  nz <- null_all[null_all$module == m, "Zsummary"]
  oz <- unname(obs_mean[m])
  if (is.na(oz)) return(NA)
  ifelse(oz < mean(nz, na.rm=TRUE), "less_preserved", "more_preserved")
})
n_less <- sum(direction == "less_preserved", na.rm=TRUE)
n_more <- sum(direction == "more_preserved", na.rm=TRUE)
cat(sprintf("Less preserved in bp_lith than null: %d / %d (%.0f%%)\n",
            n_less, length(modules), 100*n_less/length(modules)))
cat(sprintf("More preserved in bp_lith than null: %d / %d (%.0f%%)\n",
            n_more, length(modules), 100*n_more/length(modules)))

emp_p <- sapply(modules, function(m) {
  nz <- null_all[null_all$module == m, "Zsummary"]
  oz <- unname(obs_mean[m])
  if (is.na(oz)||!length(nz)) return(NA_real_)
  mean(nz <= oz)
})
cat(sprintf("\nModules with empirical p = 0.05 (= 1/%d, the minimum): %d\n",
            n_done, sum(emp_p <= 0.05, na.rm=TRUE)))
padj <- p.adjust(emp_p, method="BH")
cat(sprintf("BH-adjusted significant at 5%%: %d\n", sum(padj <= 0.05, na.rm=TRUE)))
cat(sprintf("Note: minimum achievable p = 1/%d = %.3f; BH 5%% over %d modules requires p_(1) <= %.4f — impossible\n",
            n_done, 1/n_done, length(modules)+1, 0.05*1/(length(modules)+1)))

cat("\n=== Top 5 most reduced below null ===\n")
delta <- sapply(modules, function(m) {
  unname(obs_mean[m]) - mean(null_all[null_all$module==m,"Zsummary"],na.rm=TRUE)
})
top5 <- sort(delta)[1:5]
for (i in seq_along(top5)) cat(sprintf("  %s: observed-null=%.2f  p=%.3f\n",
  names(top5)[i], top5[i], emp_p[names(top5)[i]]))

cat("\n=== Top 5 most elevated above null ===\n")
top5u <- sort(delta, decreasing=TRUE)[1:5]
for (i in seq_along(top5u)) cat(sprintf("  %s: observed-null=%.2f  p=%.3f\n",
  names(top5u)[i], top5u[i], emp_p[names(top5u)[i]]))

cat("\n=== DONE ===\n")
