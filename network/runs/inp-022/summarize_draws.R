#!/usr/bin/env Rscript
run_dir <- 'network/runs/inp-022'

# Read all lith draws (1-19)
draws <- list()
for (i in 1:19) {
  fn <- file.path(run_dir, sprintf('lith_draw_%03d.rds', i))
  if (file.exists(fn)) draws[[i]] <- readRDS(fn)
}
all_draws <- do.call(rbind, draws)

agg <- aggregate(cbind(Zsummary, medianRank) ~ module, data=all_draws,
  FUN=function(x) c(mean=mean(x), sd=sd(x), min=min(x), max=max(x)))
df <- data.frame(module=agg$module,
  Zsummary_mean=agg$Zsummary[,'mean'], Zsummary_sd=agg$Zsummary[,'sd'],
  Zsummary_min=agg$Zsummary[,'min'], Zsummary_max=agg$Zsummary[,'max'],
  medianRank_mean=agg$medianRank[,'mean'],
  stringsAsFactors=FALSE)
df <- df[order(-df$Zsummary_mean),]
write.csv(df, file.path(run_dir, 'lith_draws1to19_summary.csv'), row.names=FALSE)

# nolith
nolith <- readRDS(file.path(run_dir, 'nolith_preservation.rds'))
ndf <- nolith[order(-nolith$Zsummary),]
write.csv(ndf, file.path(run_dir, 'nolith_summary.csv'), row.names=FALSE)

# Print key results
cat('=== inp-022 lith draws 1-19 summary ===\n')
print(head(df, 15))
gold_z <- df$Zsummary_mean[df$module=='gold']
cat('\nGold mean Z:', gold_z, '\n')
cat('Modules beating gold (lith, draws 1-19):', sum(df$Zsummary_mean > gold_z, na.rm=TRUE), '\n')
cat('All modules with mean Z > 10:\n')
print(df[df$Zsummary_mean > 10, c('module','Zsummary_mean','Zsummary_sd')])

cat('\n=== inp-022 nolith summary (top 15) ===\n')
print(head(ndf, 15))
gold_z_no <- nolith$Zsummary[nolith$module=='gold']
cat('\nnolith gold Z:', gold_z_no, '\n')
cat('Modules beating gold (nolith):', sum(nolith$Zsummary > gold_z_no, na.rm=TRUE), '\n')
