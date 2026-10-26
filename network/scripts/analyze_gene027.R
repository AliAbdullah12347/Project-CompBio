#!/usr/bin/env Rscript
# Post-hoc analysis for gene-027 (LM22 circularity check)
# Run after all 20 lith draws + nolith complete.
# Prints a summary table ready to paste into the finding file.

run_dir_gene027  <- "network/runs/gene-027"
run_dir_base001  <- "network/runs/base-001"

# ---- helpers ---------------------------------------------------------------
cat_section <- function(title) cat("\n===", title, "===\n")

# ---- load gene-027 checkpoints ----------------------------------------------
nolith_ckpt <- file.path(run_dir_gene027, "nolith_preservation.rds")
stopifnot(file.exists(nolith_ckpt))
nr <- readRDS(nolith_ckpt)

done_lith <- vapply(1:20, function(i)
  file.exists(file.path(run_dir_gene027, sprintf("lith_draw_%03d.rds", i))), logical(1))
cat("Lith draws complete:", sum(done_lith), "/ 20\n")

lith_all <- do.call(rbind, lapply(which(done_lith), function(i)
  readRDS(file.path(run_dir_gene027, sprintf("lith_draw_%03d.rds", i)))))

# ---- base-001 nolith for comparison ----------------------------------------
b1_nolith <- read.csv(file.path(run_dir_base001, "results_nolith.csv"),
                      stringsAsFactors = FALSE)
b1_lith   <- read.csv(file.path(run_dir_base001, "results_lith_alldraws.csv"),
                      stringsAsFactors = FALSE)

# ---- gene-027 nolith -------------------------------------------------------
cat_section("gene-027 nolith (50 perms)")
gold_n <- nr[nr$module == "gold", ]
cat(sprintf("gold: Zsummary=%.2f  medianRank=%.0f of %.0f\n",
            gold_n$Zsummary, gold_n$medianRank, max(nr$medianRank)))
cat(sprintf("n modules (incl gold): %d\n", nrow(nr)))
cat(sprintf("gold ranks last (medianRank == max): %s\n",
            gold_n$medianRank == max(nr$medianRank)))

# ---- base-001 nolith for comparison ----------------------------------------
cat_section("base-001 nolith (500 perms) [for comparison]")
b1_gold_n <- b1_nolith[b1_nolith$module == "gold", ]
cat(sprintf("gold: Zsummary=%.2f  medianRank=%.0f of %.0f\n",
            b1_gold_n$Zsummary, b1_gold_n$medianRank, max(b1_nolith$medianRank)))

# ---- gene-027 lith ---------------------------------------------------------
cat_section("gene-027 lith (20 draws, 50 perms each)")
n_draws_done <- sum(done_lith)
if (n_draws_done > 0) {
  gold_L <- lith_all[lith_all$module == "gold", ]
  gold_Z <- gold_L$Zsummary
  cat(sprintf("gold lith: mean Z=%.2f  sd=%.2f  [%.2f, %.2f]\n",
              mean(gold_Z), sd(gold_Z), min(gold_Z), max(gold_Z)))
  gold_rank_vals <- gold_L$medianRank
  n_mods_per_draw <- nrow(lith_all) / n_draws_done  # includes gold
  cat(sprintf("gold medianRank: always last (== %d): %s  (range %.0f–%.0f)\n",
              round(mean(gold_rank_vals)), all(gold_rank_vals == max(gold_rank_vals)),
              min(gold_rank_vals), max(gold_rank_vals)))

  # Beat-gold count per module
  mods <- unique(lith_all$module)
  mods <- mods[mods != "gold"]
  beat_counts <- sapply(mods, function(m) {
    mod_Z <- lith_all[lith_all$module == m, "Zsummary"]
    sum(mod_Z > gold_Z, na.rm = TRUE)
  })
  cat(sprintf("\nModules beating gold in all %d draws: %d of %d\n",
              n_draws_done, sum(beat_counts == n_draws_done), length(beat_counts)))
  cat(sprintf("Modules beating gold in 0 draws: %d\n", sum(beat_counts == 0)))
  cat(sprintf("Modules beating gold in < 10 draws: %d\n", sum(beat_counts < 10)))

  cat("\nBeat-gold count distribution:\n")
  print(table(beat_counts))
} else {
  cat("No lith draws yet.\n")
}

# ---- base-001 lith for comparison -------------------------------------------
cat_section("base-001 lith (20 draws, 50 perms) [for comparison]")
b1_gold_L  <- b1_lith[b1_lith$module == "gold", ]
b1_gold_Z  <- b1_gold_L$Zsummary
cat(sprintf("gold lith: mean Z=%.2f  sd=%.2f  [%.2f, %.2f]\n",
            mean(b1_gold_Z), sd(b1_gold_Z), min(b1_gold_Z), max(b1_gold_Z)))
b1_mods <- unique(b1_lith$module); b1_mods <- b1_mods[b1_mods != "gold"]
b1_beats <- sapply(b1_mods, function(m) {
  mZ <- b1_lith[b1_lith$module == m, "Zsummary"]
  sum(mZ > b1_gold_Z, na.rm = TRUE)
})
cat(sprintf("Modules beating gold in all 20: %d of %d\n",
            sum(b1_beats == 20), length(b1_beats)))
cat(sprintf("Modules beating gold in 0: %d\n", sum(b1_beats == 0)))

# ---- module-level comparison (by colour): are rankings similar? ------------
cat_section("Module colour overlap: gene-027 vs base-001")
g27_mods <- unique(nr$module); g27_mods <- g27_mods[g27_mods != "gold"]
b1_mods  <- unique(b1_nolith$module); b1_mods <- b1_mods[b1_mods != "gold"]
shared   <- intersect(g27_mods, b1_mods)
cat(sprintf("gene-027 modules (non-grey, non-gold): %d\n", length(g27_mods)))
cat(sprintf("base-001 modules: %d\n", length(b1_mods)))
cat(sprintf("Shared colour names: %d\n", length(shared)))
cat(sprintf("In base-001 only: %d\n", length(setdiff(b1_mods, g27_mods))))
cat(sprintf("In gene-027 only: %d\n", length(setdiff(g27_mods, b1_mods))))

if (length(shared) > 0) {
  g27_Z <- setNames(nr$Zsummary[match(shared, nr$module)], shared)
  b1_Z  <- setNames(b1_nolith$Zsummary[match(shared, b1_nolith$module)], shared)
  cor_Z <- cor(g27_Z, b1_Z, use = "complete.obs")
  cat(sprintf("Correlation of Zsummary (shared modules, nolith): %.3f\n", cor_Z))
}

cat("\n=== DONE ===\n")
