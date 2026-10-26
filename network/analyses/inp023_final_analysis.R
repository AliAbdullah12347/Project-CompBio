#!/usr/bin/env Rscript
# Run this script after all 20 lith draws of inp-023 are complete.
# Prints statistics for the finding file.
# Usage: Rscript network/analyses/inp023_final_analysis.R

run_dir <- "network/runs/inp-023"

# ---- nolith ----------------------------------------------------------------
nolith_ckpt <- file.path(run_dir, "nolith_preservation.rds")
stopifnot(file.exists(nolith_ckpt))
nr <- readRDS(nolith_ckpt)

cat("\n--- bp_nolith ---\n")
gold_nr <- nr[nr$module == "gold", ]
real_nr <- nr[nr$module != "gold", ]
cat(sprintf("gold Zsummary: %.2f\n", gold_nr$Zsummary))
cat(sprintf("gold medianRank: %.0f / %.0f\n", gold_nr$medianRank, nrow(nr)))
cat(sprintf("modules above gold (Zsummary): %d of %d (%.0f%%)\n",
    sum(real_nr$Zsummary > gold_nr$Zsummary, na.rm=TRUE), nrow(real_nr),
    100 * mean(real_nr$Zsummary > gold_nr$Zsummary, na.rm=TRUE)))
cat(sprintf("top modules: %s\n",
    paste(head(nr$module[order(-nr$Zsummary)], 5), collapse=", ")))

# ---- lith draws ------------------------------------------------------------
DRAWS <- 20L
lith_files <- file.path(run_dir, sprintf("lith_draw_%03d.rds", seq_len(DRAWS)))
missing_draws <- which(!file.exists(lith_files))
if (length(missing_draws) > 0)
  cat(sprintf("\nWARNING: missing draws %s\n", paste(missing_draws, collapse=",")))

existing <- which(file.exists(lith_files))
lith_list <- lapply(existing, function(i) readRDS(lith_files[i]))
lith_all  <- do.call(rbind, lith_list)

cat(sprintf("\n--- bp_lith (%d of %d draws) ---\n", length(existing), DRAWS))

# per-draw statistics
gold_per_draw <- lapply(lith_list, function(dr) dr[dr$module == "gold", ])
gold_Z <- sapply(gold_per_draw, function(x) x$Zsummary)
gold_rank <- sapply(gold_per_draw, function(x) x$medianRank)
n_total  <- sapply(lith_list, nrow)  # includes gold

cat(sprintf("gold mean Z: %.2f (sd %.2f, range %.2f–%.2f)\n",
    mean(gold_Z), sd(gold_Z), min(gold_Z), max(gold_Z)))
cat(sprintf("gold medianRank: last (%d/%d) in all %d draws? %s\n",
    round(mean(n_total)), round(mean(n_total)),
    length(existing),
    if (all(gold_rank == n_total)) "YES" else sprintf("NO (%d draws not last)", sum(gold_rank != n_total))))

# modules above gold per draw
above_gold_per_draw <- sapply(seq_along(lith_list), function(i) {
  dr <- lith_list[[i]]
  gold_Z_i <- dr$Zsummary[dr$module == "gold"]
  real_Z   <- dr$Zsummary[dr$module != "gold"]
  sum(real_Z > gold_Z_i, na.rm=TRUE)
})
n_real_modules <- n_total - 1L  # exclude gold

cat(sprintf("modules above gold per draw: mean %.1f/%.0f (%.1f%%)\n",
    mean(above_gold_per_draw), mean(n_real_modules),
    100 * mean(above_gold_per_draw / n_real_modules)))

# modules above gold in ALL draws
modules_real <- unique(lith_all$module[lith_all$module != "gold"])
above_in_all <- sapply(modules_real, function(m) {
  rows <- lith_all[lith_all$module == m, ]
  # get gold Z for each draw this module appears in
  all_above <- sapply(existing, function(i) {
    dr <- lith_list[[which(existing == i)]]
    mZ <- dr$Zsummary[dr$module == m]
    gZ <- dr$Zsummary[dr$module == "gold"]
    if (length(mZ)==0 || length(gZ)==0) return(NA)
    mZ > gZ
  })
  all(all_above, na.rm=TRUE)
})
cat(sprintf("modules above gold in ALL %d draws: %d of %d (%.0f%%)\n",
    length(existing), sum(above_in_all, na.rm=TRUE), length(modules_real),
    100 * mean(above_in_all, na.rm=TRUE)))

# per-draw gold Z breakdown
cat("\nPer-draw gold Z:\n")
for (i in seq_along(existing)) {
  cat(sprintf("  draw %2d: Z=%.2f, rank=%.0f/%.0f\n",
      existing[i], gold_Z[i], gold_rank[i], n_total[i]))
}

# summary for above-gold modules
cat(sprintf("\nModule mean Z (top 5):\n"))
mod_meanZ <- tapply(lith_all$Zsummary, lith_all$module, mean, na.rm=TRUE)
top5 <- head(sort(mod_meanZ, decreasing=TRUE), 5)
for (nm in names(top5)) cat(sprintf("  %s: %.1f\n", nm, top5[nm]))

cat(sprintf("\nModule mean Z (bottom 5 excluding gold):\n"))
real_meanZ <- mod_meanZ[names(mod_meanZ) != "gold"]
bot5 <- head(sort(real_meanZ), 5)
for (nm in names(bot5)) cat(sprintf("  %s: %.1f\n", nm, bot5[nm]))

cat(sprintf("\ngold mean Z across draws: %.2f (for comparison)\n", mod_meanZ["gold"]))
