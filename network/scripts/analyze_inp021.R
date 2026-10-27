#!/usr/bin/env Rscript
# Post-hoc analysis for inp-021 (full technical residualisation network)
# Run after nolith + all lith draws complete.
# Prints a summary table ready to paste into the finding file.

run_dir_inp021 <- "network/runs/inp-021"
run_dir_base001 <- "network/runs/base-001"
run_dir_inp023  <- "network/runs/inp-023"

cat_section <- function(title) cat("\n===", title, "===\n")

# ---- nolith ---------------------------------------------------------------
nolith_ckpt <- file.path(run_dir_inp021, "nolith_preservation.rds")
if (!file.exists(nolith_ckpt)) { cat("nolith checkpoint not yet done\n"); quit() }
nr <- readRDS(nolith_ckpt)

cat("Lith draws complete:")
done_lith <- vapply(1:20, function(i)
  file.exists(file.path(run_dir_inp021, sprintf("lith_draw_%03d.rds", i))), logical(1))
cat(sum(done_lith), "/ 20\n")

cat_section("inp-021 nolith (50 perms)")
gold_n <- nr[nr$module == "gold", ]
cat(sprintf("gold: Zsummary=%.2f  medianRank=%.0f of %.0f\n",
            gold_n$Zsummary, gold_n$medianRank, max(nr$medianRank)))
cat(sprintf("n modules (incl gold): %d\n", nrow(nr)))
cat(sprintf("gold ranks last: %s\n", gold_n$medianRank == max(nr$medianRank)))
cat(sprintf("modules above gold by Zsummary: %d of %d\n",
            sum(nr$Zsummary[nr$module != "gold"] > gold_n$Zsummary, na.rm=TRUE),
            nrow(nr) - 1))

# ---- base-001 nolith for comparison ----------------------------------------
b1_nolith <- read.csv(file.path(run_dir_base001, "results_nolith.csv"), stringsAsFactors=FALSE)
cat_section("base-001 nolith (500 perms) [for comparison]")
b1_gold_n <- b1_nolith[b1_nolith$module == "gold", ]
cat(sprintf("gold: Zsummary=%.2f  medianRank=%.0f of %.0f\n",
            b1_gold_n$Zsummary, b1_gold_n$medianRank, max(b1_nolith$medianRank)))

# ---- inp-021 lith -----------------------------------------------------------
if (sum(done_lith) > 0) {
  lith_all <- do.call(rbind, lapply(which(done_lith), function(i)
    readRDS(file.path(run_dir_inp021, sprintf("lith_draw_%03d.rds", i)))))

  cat_section(sprintf("inp-021 lith (%d draws, 50 perms each)", sum(done_lith)))
  gold_L <- lith_all[lith_all$module == "gold", ]
  gold_Z <- gold_L$Zsummary
  cat(sprintf("gold lith: mean Z=%.2f  sd=%.2f  [%.2f, %.2f]\n",
              mean(gold_Z), sd(gold_Z), min(gold_Z), max(gold_Z)))
  n_mods_total <- nrow(lith_all) / sum(done_lith)
  cat(sprintf("gold medianRank: always last (%.0f): %s  (range %.0f–%.0f)\n",
              max(gold_L$medianRank), all(gold_L$medianRank == max(gold_L$medianRank)),
              min(gold_L$medianRank), max(gold_L$medianRank)))

  mods <- unique(lith_all$module); mods <- mods[mods != "gold"]
  beat_counts <- sapply(mods, function(m) {
    mod_Z <- lith_all[lith_all$module == m, "Zsummary"]
    sum(mod_Z > gold_Z, na.rm=TRUE)
  })
  n_draws <- sum(done_lith)
  cat(sprintf("\nModules beating gold in all %d draws: %d of %d\n",
              n_draws, sum(beat_counts == n_draws), length(beat_counts)))
  cat(sprintf("Modules beating gold in 0 draws: %d\n", sum(beat_counts == 0)))
  cat(sprintf("Modules beating gold in < 10 draws: %d\n", sum(beat_counts < 10)))

  cat("\nBeat-gold count distribution:\n")
  print(table(beat_counts))
}

# ---- base-001 lith for comparison -------------------------------------------
b1_lith <- read.csv(file.path(run_dir_base001, "results_lith_alldraws.csv"), stringsAsFactors=FALSE)
cat_section("base-001 lith (20 draws, 50 perms) [for comparison]")
b1_gold_L <- b1_lith[b1_lith$module == "gold", ]
b1_gold_Z <- b1_gold_L$Zsummary
cat(sprintf("gold lith: mean Z=%.2f  sd=%.2f  [%.2f, %.2f]\n",
            mean(b1_gold_Z), sd(b1_gold_Z), min(b1_gold_Z), max(b1_gold_Z)))
b1_mods <- unique(b1_lith$module); b1_mods <- b1_mods[b1_mods != "gold"]
b1_beats <- sapply(b1_mods, function(m) {
  mZ <- b1_lith[b1_lith$module == m, "Zsummary"]
  sum(mZ > b1_gold_Z, na.rm=TRUE)
})
cat(sprintf("Modules beating gold in all 20: %d of %d\n",
            sum(b1_beats == 20), length(b1_beats)))

# ---- inp-023 for comparison -------------------------------------------------
if (file.exists(file.path(run_dir_inp023, "results_nolith.csv"))) {
  i23_nolith <- read.csv(file.path(run_dir_inp023, "results_nolith.csv"), stringsAsFactors=FALSE)
  cat_section("inp-023 nolith (50 perms, tech+ILR) [for comparison]")
  i23_gold_n <- i23_nolith[i23_nolith$module == "gold", ]
  cat(sprintf("gold: Zsummary=%.2f  medianRank=%.0f of %.0f\n",
              i23_gold_n$Zsummary, i23_gold_n$medianRank, max(i23_nolith$medianRank)))
}
if (file.exists(file.path(run_dir_inp023, "results_lith_alldraws.csv"))) {
  i23_lith <- read.csv(file.path(run_dir_inp023, "results_lith_alldraws.csv"), stringsAsFactors=FALSE)
  cat_section("inp-023 lith (20 draws, tech+ILR) [for comparison]")
  i23_gold <- i23_lith[i23_lith$module == "gold", ]
  cat(sprintf("gold lith: mean Z=%.2f  sd=%.2f  [%.2f, %.2f]\n",
              mean(i23_gold$Zsummary), sd(i23_gold$Zsummary),
              min(i23_gold$Zsummary), max(i23_gold$Zsummary)))
  i23_mods <- unique(i23_lith$module); i23_mods <- i23_mods[i23_mods != "gold"]
  i23_beats <- sapply(i23_mods, function(m) {
    mZ <- i23_lith[i23_lith$module == m, "Zsummary"]
    sum(mZ > i23_gold$Zsummary, na.rm=TRUE)
  })
  cat(sprintf("Modules beating gold in all 20: %d of %d\n",
              sum(i23_beats == 20), length(i23_mods)))
}

# ---- SFT comparison ---------------------------------------------------------
cat_section("SFT R2 comparison across arms")
sft021 <- read.csv(file.path(run_dir_inp021, "sft_r2_curve.csv"), stringsAsFactors=FALSE)
cat(sprintf("inp-021 max R2: %.3f at power %d\n",
            max(sft021$SFT.R.sq), sft021$Power[which.max(sft021$SFT.R.sq)]))

cat("\n=== DONE ===\n")
