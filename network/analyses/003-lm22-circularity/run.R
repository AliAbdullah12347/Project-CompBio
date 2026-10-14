#!/usr/bin/env Rscript
# ===========================================================================
# 003-lm22-circularity
#
# Analysis 002 found PC2 of the module eigengenes is lithium-associated and 73%
# explained by the ILR composition balances. Those balances come from CIBERSORTx
# fractions, which are estimated FROM the LM22 marker genes IN THIS SAME MATRIX.
# So the 73% may be partly circular: the composition estimate and the expression
# PC share genes. IMAC_README catalogue item 7 names exactly this check.
#
# Test: drop the 350 LM22 markers present among the 12,368 genes, recompute the
# eigengenes from the SAME module assignments, redo the PCA, and re-ask both
# questions. If the lithium association and the composition share survive, the
# finding is not an artifact of the marker genes.
#
# Writes only into this folder.
# ===========================================================================
setwd("/Users/Aaylab/Project-CompBio")
suppressPackageStartupMessages(library(WGCNA)); options(stringsAsFactors=FALSE); enableWGCNAThreads()
OUT <- "network/analyses/003-lm22-circularity/results"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), paste0(...,collapse="")))
set.seed(20261014L)

source("network/00_load.R"); d <- load_project(cache=TRUE, verbose=FALSE)
bp <- !is.na(d$groups) & d$groups %in% c("bp_nolith","bp_lith"); stopifnot(sum(bp)==226L)
md <- d$meta[bp,]; lith <- as.numeric(md$lithium)
il <- d$ilr[match(md$title,d$ilr$sample), c("b1","b2","b3","b4")]
lab <- readRDS("network/runs/mod-000-resolution/labels_p12_ds4.rds")
colors <- labels2colors(lab)

sym <- d$sym[rownames(d$logtpm)]
is_lm22 <- !is.na(sym) & sym %in% d$lm22
msg("LM22 markers among the 12,368 genes: ", sum(is_lm22))

# ---- helper: eigengenes -> PCA -> lithium + composition on a gene subset ----
assess <- function(keep, tag) {
  E <- t(d$logtpm[keep, bp])
  me <- moduleEigengenes(E, colors = colors[keep], excludeGrey = TRUE)$eigengenes
  pc <- prcomp(me, scale.=TRUE)
  # Identify the lithium-associated dimension by the SAME rule as 002: test the
  # effective dimensions, BH across them. Do not pick the best PC post hoc.
  meff <- { cm <- cor(me); lam <- abs(eigen(cm, only.values=TRUE)$values)
            ceiling(sum(ifelse(lam>=1,1,0) + (lam-floor(lam)))) }
  npc <- max(2, meff)
  rows <- do.call(rbind, lapply(seq_len(npc), function(i) {
    s <- pc$x[,i]
    m1 <- lm(s ~ lith); m2 <- lm(s ~ lith + il$b1+il$b2+il$b3+il$b4)
    g <- function(f) summary(f)$coefficients["lith", c("Estimate","Pr(>|t|)")]
    c1 <- g(m1); c2 <- g(m2)
    data.frame(set=tag, PC=paste0("PC",i), var_pct=round(100*pc$sdev[i]^2/sum(pc$sdev^2),1),
               beta=round(c1[1],4), p=signif(c1[2],3),
               beta_adj=round(c2[1],4), p_adj=signif(c2[2],3),
               atten=round(100*(1-abs(c2[1])/abs(c1[1])),1),
               r2_comp=round(100*summary(lm(s ~ il$b1+il$b2+il$b3+il$b4))$r.squared,1))
  }))
  rows$p_BH <- signif(p.adjust(rows$p,"BH"),3)
  list(tag=tag, n_genes=sum(keep), meff=meff, npc=npc, rows=rows, pc=pc, me=me)
}

full <- assess(rep(TRUE, nrow(d$logtpm)), "all 12,368 genes")
excl <- assess(!is_lm22,                  "LM22 markers EXCLUDED (12,018)")
msg("both assessments done")

cat("\n================ FULL GENE SET (replicates analysis 002) ================\n")
print(full$rows[, c("PC","var_pct","beta","p","p_BH","beta_adj","p_adj","atten","r2_comp")], row.names=FALSE)
cat("\n================ LM22 MARKERS EXCLUDED ================\n")
print(excl$rows[, c("PC","var_pct","beta","p","p_BH","beta_adj","p_adj","atten","r2_comp")], row.names=FALSE)

# which PC carries lithium in each set, by smallest BH
pick <- function(a) a$rows[which.min(a$rows$p_BH), ]
pf <- pick(full); pe <- pick(excl)
cat("\n================ SIDE BY SIDE: the lithium-associated dimension ================\n")
cmp <- data.frame(
  set        = c(pf$set, pe$set),
  genes      = c(full$n_genes, excl$n_genes),
  Meff       = c(full$meff, excl$meff),
  PC         = c(pf$PC, pe$PC),
  var_pct    = c(pf$var_pct, pe$var_pct),
  beta       = c(pf$beta, pe$beta),
  p          = c(pf$p, pe$p),
  p_BH       = c(pf$p_BH, pe$p_BH),
  atten_pct  = c(pf$atten, pe$atten),
  pc_r2_comp = c(pf$r2_comp, pe$r2_comp))
print(cmp, row.names=FALSE)
write.csv(rbind(full$rows, excl$rows), file.path(OUT,"pc_lithium_by_geneset.csv"), row.names=FALSE)
write.csv(cmp, file.path(OUT,"side_by_side.csv"), row.names=FALSE)

# how similar are the two PC2s?
k <- min(ncol(full$pc$x), ncol(excl$pc$x), 4)
cat("\n  |correlation| between the two gene sets' leading PCs (same samples):\n")
for (i in seq_len(k)) cat(sprintf("    PC%d: %.3f\n", i, abs(cor(full$pc$x[,i], excl$pc$x[,i]))))

cat("\n================ VERDICT ================\n")
circ <- (pe$p_BH >= 0.05) || (pe$r2_comp < 0.5*pf$r2_comp)
cat(sprintf("  lithium association without LM22 markers: beta %+.4f, BH %.4g  -> %s\n",
            pe$beta, pe$p_BH, ifelse(pe$p_BH<0.05,"SURVIVES","LOST")))
cat(sprintf("  composition share without LM22 markers:   %.1f%% (was %.1f%%) -> %s\n",
            pe$r2_comp, pf$r2_comp, ifelse(pe$r2_comp >= 0.5*pf$r2_comp,"HOLDS","COLLAPSES")))
cat(sprintf("  => finding 002 is %s by the circularity check\n",
            ifelse(circ, "WEAKENED", "NOT explained by the LM22 marker genes")))
saveRDS(list(full=full$rows, excl=excl$rows, cmp=cmp), file.path(OUT,"summary.rds"))
cat("\nANALYSIS_003_DONE\n")
