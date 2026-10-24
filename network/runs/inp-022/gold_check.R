#!/usr/bin/env Rscript
run_dir <- 'network/runs/inp-022'
draws <- list()
for (i in 1:20) {
  fn <- file.path(run_dir, sprintf('lith_draw_%03d.rds', i))
  if (file.exists(fn)) {
    d <- readRDS(fn)
    gold_row <- d[d$module == 'gold', ]
    draws[[i]] <- data.frame(draw=i, Zsummary=gold_row$Zsummary,
                              medianRank=gold_row$medianRank)
  }
}
df <- do.call(rbind, draws)
cat('inp-022 gold Z across 20 draws:\n')
cat('  Mean:', mean(df$Zsummary), '\n')
cat('  SD:  ', sd(df$Zsummary), '\n')
cat('  Min: ', min(df$Zsummary), '\n')
cat('  Max: ', max(df$Zsummary), '\n')
cat('  Per draw:\n')
print(df, row.names=FALSE)

# compare to base-001 gold
b1 <- read.csv('network/runs/base-001/results_lith_summary.csv')
b1_gold <- b1[b1$module == 'gold', ]
cat('\nbase-001 gold Z stats:\n')
cat('  Mean:', b1_gold$Zsummary_mean, '\n')
cat('  SD:  ', b1_gold$Zsummary_sd, '\n')
cat('  Min: ', b1_gold$Zsummary_min, '\n')
cat('  Max: ', b1_gold$Zsummary_max, '\n')
