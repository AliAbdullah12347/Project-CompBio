#!/usr/bin/env Rscript
# ===========================================================================
# 004-decomposition-robustness
#
# Finding 002 (lithium-associated PC2, 73% composition) rests on ONE
# decomposition: deepSplit=4, mergeCutHeight=0. Its robustness score is held
# down for exactly that reason.
#
# Test: rebuild all ten cells of the mod-000-resolution grid - deepSplit {2,4} x
# mergeCutHeight {0, 0.05, 0.10, 0.15, 0.25} - and ask the same question of each.
# Modules are re-merged on the CONTROL reference (IMAC_README 6.2: defined once
# on the reference), eigengenes computed on BP1, lithium tested against the Li &
# Ji effective dimensions with BH across them. Same a-priori rule every time; no
# post-hoc PC picking.
#
# If a lithium-associated dimension appears across decompositions - including the
# degenerate ones - the finding is a property of the data. If it appears only at
# mergeCutHeight=0, it is a property of that parameter choice.
#
# Writes only into this folder.
# ===========================================================================
setwd("/Users/Aaylab/Project-CompBio")
suppressPackageStartupMessages(library(WGCNA)); options(stringsAsFactors=FALSE); enableWGCNAThreads()
OUT <- "network/analyses/004-decomposition-robustness/results"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), paste0(...,collapse="")))
set.seed(20261014L)

source("network/00_load.R"); d <- load_project(cache=TRUE, verbose=FALSE)
ctrl <- !is.na(d$groups) & d$groups=="control"          # reference: modules defined here
bp   <- !is.na(d$groups) & d$groups %in% c("bp_nolith","bp_lith")
stopifnot(sum(ctrl)==234L, sum(bp)==226L)
Eref <- t(d$logtpm[, ctrl]); Ebp <- t(d$logtpm[, bp])
md <- d$meta[bp,]; lith <- as.numeric(md$lithium)
il <- d$ilr[match(md$title,d$ilr$sample), c("b1","b2","b3","b4")]

meff_li_ji <- function(cm) { lam <- abs(eigen(cm, only.values=TRUE)$values)
                             sum(ifelse(lam>=1,1,0) + (lam-floor(lam))) }

GRID <- expand.grid(deepSplit=c(2,4), mch=c(0, 0.05, 0.10, 0.15, 0.25))
rows <- list()
for (i in seq_len(nrow(GRID))) {
  ds <- GRID$deepSplit[i]; mch <- GRID$mch[i]
  lab <- readRDS(sprintf("network/runs/mod-000-resolution/labels_p12_ds%d.rds", ds))
  cols <- labels2colors(lab)
  if (mch > 0) cols <- mergeCloseModules(Eref, cols, cutHeight=mch, verbose=0)$colors
  nmod <- length(setdiff(unique(cols),"grey"))
  if (nmod < 2) { msg("ds=",ds," mch=",mch," -> only ",nmod," module(s), skipping"); next }

  me <- moduleEigengenes(Ebp, colors=cols, excludeGrey=TRUE)$eigengenes
  cm <- cor(me); npc <- max(2, ceiling(meff_li_ji(cm)))
  npc <- min(npc, ncol(me))
  pc <- prcomp(me, scale.=TRUE)

  r <- do.call(rbind, lapply(seq_len(npc), function(j) {
    s <- pc$x[,j]
    c1 <- summary(lm(s ~ lith))$coefficients["lith", c("Estimate","Pr(>|t|)")]
    c2 <- summary(lm(s ~ lith + il$b1+il$b2+il$b3+il$b4))$coefficients["lith", c("Estimate","Pr(>|t|)")]
    data.frame(pc=j, var=100*pc$sdev[j]^2/sum(pc$sdev^2), beta=c1[1], p=c1[2],
               beta_adj=c2[1], p_adj=c2[2],
               r2c=100*summary(lm(s ~ il$b1+il$b2+il$b3+il$b4))$r.squared)
  }))
  r$p_BH <- p.adjust(r$p, "BH")
  best <- r[which.min(r$p_BH), ]
  rows[[length(rows)+1]] <- data.frame(
    deepSplit=ds, mergeCutHeight=mch, n_modules=nmod,
    median_abs_r=round(median(abs(cm[upper.tri(cm)])),3),
    Meff=round(meff_li_ji(cm),1), tested_PCs=npc,
    lith_PC=best$pc, lith_var_pct=round(best$var,1),
    beta=round(best$beta,4), p=signif(best$p,3), p_BH=signif(best$p_BH,3),
    sig=ifelse(best$p_BH<0.05,"YES","no"),
    beta_adj=round(best$beta_adj,4), p_adj=signif(best$p_adj,3),
    atten_pct=round(100*(1-abs(best$beta_adj)/abs(best$beta)),1),
    pc_r2_comp=round(best$r2c,1))
  msg(sprintf("ds=%d mch=%.2f -> %2d modules | lithium on PC%d BH=%.4g | comp r2=%.0f%% | atten=%.0f%%",
              ds, mch, nmod, best$pc, best$p_BH, best$r2c, 100*(1-abs(best$beta_adj)/abs(best$beta))))
}
R <- do.call(rbind, rows)
write.csv(R, file.path(OUT,"robustness_grid.csv"), row.names=FALSE)

cat("\n================ LITHIUM DIMENSION ACROSS 10 DECOMPOSITIONS ================\n")
print(R[, c("deepSplit","mergeCutHeight","n_modules","Meff","lith_PC","lith_var_pct",
            "beta","p_BH","sig","atten_pct","pc_r2_comp")], row.names=FALSE)

cat("\n================ SUMMARY ================\n")
cat(sprintf("  decompositions where a lithium dimension reaches BH<0.05 : %d / %d\n", sum(R$sig=="YES"), nrow(R)))
cat(sprintf("  of those, median composition share of that dimension     : %.1f%%\n", median(R$pc_r2_comp[R$sig=="YES"])))
cat(sprintf("  of those, median attenuation after composition adjustment: %.1f%%\n", median(R$atten_pct[R$sig=="YES"])))
cat(sprintf("  of those, still BH<0.05 after composition adjustment      : %d\n",
            sum(R$sig=="YES" & p.adjust(R$p_adj[R$sig=="YES"],"BH") < 0.05)))
cat(sprintf("\n  usable cells (>=10 modules, largest<25%%) that replicate   : %d\n",
            sum(R$sig=="YES" & R$n_modules>=10)))
cat(sprintf("  degenerate cells (<10 modules) that also replicate        : %d\n",
            sum(R$sig=="YES" & R$n_modules<10)))
saveRDS(R, file.path(OUT,"summary.rds"))
cat("\nANALYSIS_004_DONE\n")
