#!/usr/bin/env Rscript
x1 <- readRDS('network/runs/inp-022/lith_draw_001.rds')
cat('Class:', class(x1), '\n')
cat('Names:', names(x1), '\n')
if (is.list(x1)) {
  for (nm in names(x1)) {
    cat('\n$', nm, ': class=', class(x1[[nm]]), '\n', sep='')
    if (is.list(x1[[nm]])) cat('  Names:', names(x1[[nm]]), '\n')
  }
}
