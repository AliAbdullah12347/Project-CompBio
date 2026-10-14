#!/usr/bin/env Rscript
# ===========================================================================
# 002-effective-dimensionality
#
# Analysis 001 found the 33 module eigengenes correlate at median |r| = 0.822.
# If they are not independent, two things follow:
#   (a) BH over 33 modules is mis-specified (conservative), and
#   (b) testing 33 correlated modules is a WEAKER test than testing the few
#       dimensions they actually span.
#
# So: how many effective dimensions are there, and does lithium associate with
# any of THEM? That is a properly powered version of the question 001 asked.
#
# Also compares the usable decomposition (33 modules, deepSplit=4, mch=0) with
# the degenerate one base-001 used (mergeCutHeight=0.25), to see whether the
# dimensionality is a property of the DATA or of the parameter choice.
#
# Writes only into this folder.
# ===========================================================================
setwd("/Users/Aaylab/Project-CompBio")
suppressPackageStartupMessages({ library(WGCNA) })
options(stringsAsFactors = FALSE); enableWGCNAThreads()
OUT <- "network/analyses/002-effective-dimensionality/results"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), paste0(...,collapse="")))
set.seed(20261014L)

source("network/00_load.R"); d <- load_project(cache = TRUE, verbose = FALSE)
bp <- !is.na(d$groups) & d$groups %in% c("bp_nolith","bp_lith")
stopifnot(sum(bp) == 226L)
E  <- t(d$logtpm[, bp]); md <- d$meta[bp, ]
il <- d$ilr[match(md$title, d$ilr$sample), c("b1","b2","b3","b4")]
lith <- as.numeric(md$lithium)

# Li & Ji (2005) effective number of independent tests from an eigenvalue spectrum.
meff_li_ji <- function(cm) {
  lam <- abs(eigen(cm, only.values = TRUE)$values)
  sum(ifelse(lam >= 1, 1, 0) + (lam - floor(lam)))
}

analyse <- function(labels_file, tag) {
  lab <- readRDS(labels_file); colors <- labels2colors(lab)
  me  <- moduleEigengenes(E, colors = colors, excludeGrey = TRUE)$eigengenes
  k   <- ncol(me)
  cm  <- cor(me); offd <- cm[upper.tri(cm)]
  ev  <- eigen(cm, only.values = TRUE)$values
  vexp <- cumsum(ev)/sum(ev)
  list(tag = tag, k = k,
       median_abs_r = median(abs(offd)), mean_r = mean(offd),
       meff = meff_li_ji(cm),
       pc_for_80 = which(vexp >= 0.80)[1], pc_for_90 = which(vexp >= 0.90)[1],
       pc1_var = ev[1]/sum(ev), ev = ev, me = me)
}

a33 <- analyse("network/runs/mod-000-resolution/labels_p12_ds4.rds", "deepSplit=4, mergeCutHeight=0 (33 modules)")
msg("usable decomposition analysed: ", a33$k, " modules")

# degenerate comparison: re-cut the cached dendrogram at mergeCutHeight 0.25 is
# already stored as a grid cell; use base-001's own reference colours instead.
ref <- readRDS("network/cache/ref_control_all_bicor_signed_p12.rds")
me4 <- moduleEigengenes(E, colors = ref$moduleColors, excludeGrey = TRUE)$eigengenes
cm4 <- cor(me4); ev4 <- eigen(cm4, only.values=TRUE)$values
a4 <- list(tag="mergeCutHeight=0.25 (base-001, 4 modules)", k=ncol(me4),
           median_abs_r=median(abs(cm4[upper.tri(cm4)])), mean_r=mean(cm4[upper.tri(cm4)]),
           meff=meff_li_ji(cm4), pc_for_80=which(cumsum(ev4)/sum(ev4)>=0.80)[1],
           pc_for_90=which(cumsum(ev4)/sum(ev4)>=0.90)[1], pc1_var=ev4[1]/sum(ev4))

cat("\n================ EFFECTIVE DIMENSIONALITY ================\n")
tab <- do.call(rbind, lapply(list(a33, a4), function(a) data.frame(
  decomposition = a$tag, n_modules = a$k,
  median_abs_r = round(a$median_abs_r,3), mean_r = round(a$mean_r,3),
  Meff_LiJi = round(a$meff,2), PCs_for_80pct = a$pc_for_80, PCs_for_90pct = a$pc_for_90,
  PC1_var_pct = round(100*a$pc1_var,1))))
print(tab, row.names = FALSE)
write.csv(tab, file.path(OUT, "effective_dimensionality.csv"), row.names = FALSE)

cat("\n  eigenvalue spectrum of the 33 eigengenes (first 10):\n  ")
cat(paste(sprintf("%.2f", a33$ev[1:10]), collapse="  "), "\n")

# ---- the properly powered test: lithium against the PCs -------------------
pc <- prcomp(a33$me, scale. = TRUE)
npc <- max(2, ceiling(a33$meff))
cat(sprintf("\n================ LITHIUM vs THE %d EFFECTIVE DIMENSIONS ================\n", npc))
rows <- do.call(rbind, lapply(seq_len(npc), function(i) {
  s <- pc$x[, i]
  m1 <- lm(s ~ lith); m2 <- lm(s ~ lith + il$b1 + il$b2 + il$b3 + il$b4)
  g <- function(f) summary(f)$coefficients["lith", c("Estimate","Pr(>|t|)")]
  c1 <- g(m1); c2 <- g(m2)
  data.frame(PC = paste0("PC", i),
             var_pct = round(100*pc$sdev[i]^2/sum(pc$sdev^2), 1),
             beta_marginal = round(c1[1],4), p_marginal = signif(c1[2],3),
             beta_compadj  = round(c2[1],4), p_compadj  = signif(c2[2],3),
             pct_attenuation = round(100*(1-abs(c2[1])/abs(c1[1])),1))
}))
rows$p_marginal_BH <- signif(p.adjust(rows$p_marginal, "BH"), 3)
rows$p_compadj_BH  <- signif(p.adjust(rows$p_compadj,  "BH"), 3)
print(rows, row.names = FALSE)
write.csv(rows, file.path(OUT, "lithium_vs_effective_dimensions.csv"), row.names = FALSE)

cat(sprintf("\n  PCs with BH<0.05 for lithium: marginal %d / %d | composition-adjusted %d\n",
            sum(rows$p_marginal_BH<0.05), nrow(rows), sum(rows$p_compadj_BH<0.05)))

# what IS PC1, if not lithium?
cat("\n================ WHAT DOES PC1 TRACK? ================\n")
cand <- data.frame(age=as.numeric(md$age), rin=as.numeric(md$rin),
                   seqpc1=as.numeric(md$seqpc1), gran=d$frac$gran[match(md$title,d$frac$sample)],
                   mono=d$frac$mono[match(md$title,d$frac$sample)],
                   Tcell=d$frac$T[match(md$title,d$frac$sample)],
                   b1=il$b1, b2=il$b2, lithium=lith)
cc <- sapply(cand, function(x) cor(pc$x[,1], x, use="pairwise.complete.obs"))
cc <- sort(cc[order(-abs(cc))])
for (n in names(cc)) cat(sprintf("  PC1 vs %-8s r = %+.3f\n", n, cc[n]))
saveRDS(list(tab=tab, rows=rows, pc1_cor=cc, ev=a33$ev), file.path(OUT,"summary.rds"))
cat("\nANALYSIS_002_DONE\n")
