#!/usr/bin/env Rscript
# Post-hoc analysis for perm-031 (group-label permutation null)
# Run after all K=20 permuted draws complete.
# Prints a summary table ready to paste into the finding file.

run_dir <- "network/runs/perm-031"
run_dir_base <- "network/runs/base-001"

cat_section <- function(title) cat("\n===", title, "===\n")

done_k <- vapply(1:20, function(i)
  file.exists(file.path(run_dir, sprintf("perm_draw_%03d.rds", i))), logical(1))
cat("Permuted draws complete:", sum(done_k), "/ 20\n")

results_file <- file.path(run_dir, "results_permutation_test.csv")
null_file    <- file.path(run_dir, "null_zsummary_alldraws.csv")

if (!file.exists(results_file)) {
  if (sum(done_k) < 20) {
    cat("Results not yet aggregated (", sum(done_k), "/20 draws done)\n")
    quit(save="no", status=0)
  }
  cat("results_permutation_test.csv missing despite all draws done — rerun perm-031\n")
  quit(save="no", status=1)
}

cmp <- read.csv(results_file, stringsAsFactors=FALSE)
null_all <- read.csv(null_file, stringsAsFactors=FALSE)

cat_section("perm-031 overall")
n_mods <- nrow(cmp)
n_less <- sum(cmp$direction == "less_preserved_than_random", na.rm=TRUE)
n_more <- sum(cmp$direction == "more_preserved_than_random", na.rm=TRUE)
n_p05  <- sum(cmp$empirical_pval_lower <= 0.05, na.rm=TRUE)
n_sig  <- sum(cmp$padj <= 0.05, na.rm=TRUE)
cat(sprintf("Total modules: %d\n", n_mods))
cat(sprintf("Less preserved than null: %d of %d (%.0f%%)\n", n_less, n_mods, 100*n_less/n_mods))
cat(sprintf("More preserved than null: %d of %d (%.0f%%)\n", n_more, n_mods, 100*n_more/n_mods))
cat(sprintf("Empirical p=0.05 (min possible): %d\n", n_p05))
cat(sprintf("BH significant at 5%%: %d (min p with K=20 is 1/20=0.05, BH over %d modules)\n", n_sig, n_mods))

cat_section("Null Zsummary distribution")
cat(sprintf("null mean across all modules and draws: %.2f\n",
            mean(null_all$Zsummary, na.rm=TRUE)))
cat(sprintf("null sd: %.2f\n", sd(null_all$Zsummary, na.rm=TRUE)))
cat(sprintf("null range: [%.2f, %.2f]\n",
            min(null_all$Zsummary, na.rm=TRUE), max(null_all$Zsummary, na.rm=TRUE)))

cat_section("Gold module in permutation null")
gold_null <- null_all[null_all$module == "gold", ]
if (nrow(gold_null) > 0) {
  cat(sprintf("null gold Zsummary: mean=%.2f  sd=%.2f  [%.2f, %.2f]\n",
              mean(gold_null$Zsummary), sd(gold_null$Zsummary),
              min(gold_null$Zsummary), max(gold_null$Zsummary)))
}
gold_cmp <- cmp[cmp$module == "gold", ]
if (nrow(gold_cmp) > 0) {
  cat(sprintf("observed gold (base-001 lith mean): %.2f\n",
              gold_cmp$obs_Zsummary_mean))
  cat(sprintf("gold empirical p (lower-tail): %.3f  padj=%.3f\n",
              gold_cmp$empirical_pval_lower, gold_cmp$padj))
}

cat_section("Observed vs null by module (top 10 most below null)")
cmp_sorted <- cmp[order(cmp$obs_Zsummary_mean - cmp$null_Zsummary_mean), ]
print(head(cmp_sorted[, c("module","obs_Zsummary_mean","null_Zsummary_mean",
                           "null_Zsummary_sd","empirical_pval_lower","padj","direction")], 10))

cat_section("Observed vs null by module (top 10 most above null)")
cmp_sorted2 <- cmp[order(-(cmp$obs_Zsummary_mean - cmp$null_Zsummary_mean)), ]
print(head(cmp_sorted2[, c("module","obs_Zsummary_mean","null_Zsummary_mean",
                            "null_Zsummary_sd","empirical_pval_lower","padj","direction")], 10))

cat("\n=== DONE ===\n")
