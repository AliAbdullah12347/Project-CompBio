#!/usr/bin/env Rscript
# ==============================================================================
# filter_sweep_deep.R -- decompose the filter effect. No new model fits; this
# reads the topTables filter_sweep.R cached, so every number here is the same
# fit that produced the sweep table.
#
# Four questions the sweep table raises but cannot settle:
#
#  A. WHERE ON THE EXPRESSION AXIS IS THE LITHIUM SIGNAL?
#     Bin the loosest universe by average expression and by detection rate --
#     the quantity the filter actually thresholds -- and read pi0 per bin. A
#     bin with pi0 near 1 contains no detectable signal, whatever its size.
#
#  B. WHY DOES LOOSENING THE FILTER *LOSE* DEGs?
#     Three candidate mechanisms: extra genes add Benjamini-Hochberg burden;
#     TMM factors move; the voom mean-variance trend is refitted. Separate them
#     by splicing p-value vectors together, which isolates one mechanism at a
#     time without refitting anything.
#
#  C. DOES THE FILTER CHANGE THE RANKING, OR ONLY THE CUTOFF?
#     FDR-set overlap confounds the two. A rank-based concordance at the top
#     (CAT) curve does not.
#
#  D. IS THE BURDEN STORY ENOUGH ON ITS OWN?
#     Fixed nominal p-value counts on a fixed gene set remove multiple testing
#     from the picture entirely.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})
source("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation/scripts/common.R")

RUN   <- file.path(EXP_HOME, "runs", "filter_sweep")
CACHE <- file.path(RUN, "cache")
FDR   <- 0.05
p     <- load_prep()

get <- function(tag) readRDS(file.path(CACHE, paste0(tag, ".rds")))
CANON <- "mc10_mp090"; LOOSE <- "mc0_mp010"
ref <- get(CANON); loo <- get(LOOSE)
refU <- ref$tt$gene; refDEG <- deg_ids(ref$tt, FDR)
cat(sprintf("canonical %s: %d genes, %d DEG\nloosest   %s: %d genes, %d DEG\n",
            CANON, ref$n_genes, length(refDEG), LOOSE, loo$n_genes, n_deg(loo$tt, FDR)))

d0  <- build_design(p, "lithium")
cnt <- p$counts[, d0$samples, drop = FALSE]

# ===========================================================================
# A. Where the signal lives on the expression axis
# ===========================================================================
tt <- loo$tt
tt$detect10 <- rowMeans(cnt[tt$gene, ] > 10)      # the canonical filter's own statistic
tt$detect0  <- rowMeans(cnt[tt$gene, ] > 0)
tt$in_canon <- tt$gene %in% refU

# Zero-padded labels so the deciles sort in numeric order when split() turns
# them into names.
bins <- function(x, k = 10) sprintf("D%02d", cut(x, quantile(x, seq(0, 1, length.out = k + 1)),
                                include.lowest = TRUE, labels = FALSE))
tt$expr_bin   <- bins(tt$AveExpr)
tt$detect_bin <- cut(tt$detect10, c(-.01, .01, .1, .25, .5, .75, .9, 1.0001),
                     labels = c("0", "0-10%", "10-25%", "25-50%", "50-75%", "75-90%", ">=90%"))

strat <- function(g, lab) {
  s <- split(seq_len(nrow(tt)), g)
  s <- s[lengths(s) > 0]   # cut() keeps empty levels; they would give NaN pi0
  do.call(rbind, lapply(names(s), function(k) {
    i <- s[[k]]
    data.frame(axis = lab, bin = k, n = length(i),
               med_aveexpr = median(tt$AveExpr[i]),
               med_detect10 = median(tt$detect10[i]),
               pi0 = pi0_storey(tt$P.Value[i]),
               pct_p05 = 100 * mean(tt$P.Value[i] < 0.05),
               pct_p001 = 100 * mean(tt$P.Value[i] < 0.001),
               n_deg_globalFDR = sum(tt$adj.P.Val[i] < FDR),
               n_deg_binFDR = sum(p.adjust(tt$P.Value[i], "BH") < FDR),
               frac_in_canon = mean(tt$in_canon[i]),
               stringsAsFactors = FALSE)
  }))
}
A <- rbind(strat(tt$expr_bin, "AveExpr decile"), strat(tt$detect_bin, "detection rate (>10 counts)"))
cat("\n== A. signal by expression stratum (loosest universe, 37k genes) ==\n")
print(A, row.names = FALSE, digits = 4)
write.csv(A, file.path(RUN, "strata.csv"), row.names = FALSE)

# ===========================================================================
# B. Why loosening loses DEGs -- splice the p-value vectors
# ===========================================================================
extra   <- setdiff(loo$tt$gene, refU)
p_extra <- loo$tt$P.Value[match(extra, loo$tt$gene)]
p_canon <- ref$tt$P.Value                           # canonical genes, canonical fit
p_loose_on_canon <- loo$tt$P.Value[match(refU, loo$tt$gene)]  # same genes, loose fit

nd <- function(pv) sum(p.adjust(pv, "BH") < FDR)
set.seed(20261218)
sim <- replicate(500, nd(c(p_canon, runif(length(extra)))))

B <- data.frame(
  scenario = c("canonical as published",
               "canonical p-values + uniform noise genes (median of 500 draws)",
               "canonical p-values + the real extra genes' p-values",
               "loose fit restricted to the canonical genes, BH inside them",
               "loose fit as run (all 37k)"),
  n_genes  = c(length(p_canon), length(p_canon) + length(extra),
               length(p_canon) + length(extra), length(p_canon),
               nrow(loo$tt)),
  n_deg    = c(nd(p_canon), median(sim), nd(c(p_canon, p_extra)),
               nd(p_loose_on_canon), nd(loo$tt$P.Value)),
  stringsAsFactors = FALSE)
cat("\n== B. decomposition of the DEG loss ==\n")
print(B, row.names = FALSE)
cat(sprintf("uniform-padding simulation: median %.0f  (2.5%%-97.5%%: %.0f-%.0f)\n",
            median(sim), quantile(sim, .025), quantile(sim, .975)))
cat(sprintf("extra genes: %d tested, %d reach FDR<0.05 in the loose fit, pi0=%.3f\n",
            length(extra), sum(loo$tt$adj.P.Val[match(extra, loo$tt$gene)] < FDR),
            pi0_storey(p_extra)))
write.csv(B, file.path(RUN, "decomposition.csv"), row.names = FALSE)

# how far did the canonical genes' own p-values move?
cat(sprintf("canonical genes re-fitted in the loose matrix: Spearman(p) = %.4f, Pearson(t) = %.4f, Pearson(logFC) = %.4f\n",
            cor(p_canon, p_loose_on_canon, method = "spearman"),
            cor(ref$tt$t, loo$tt$t[match(refU, loo$tt$gene)]),
            cor(ref$tt$logFC, loo$tt$logFC[match(refU, loo$tt$gene)])))

# ===========================================================================
# C. Rank concordance -- threshold-free
# ===========================================================================
cat_curve <- function(tag, ks = c(100, 250, 500, 1000, 2000, 5000)) {
  f <- get(tag); if (!isTRUE(f$ok)) return(NULL)
  a <- f$tt$gene[order(f$tt$P.Value)]
  b <- refU[order(ref$tt$P.Value)]
  data.frame(tag = tag, n_genes = f$n_genes, k = ks,
             cat = vapply(ks, function(k)
               if (min(length(a), length(b)) < k) NA_real_
               else length(intersect(a[1:k], b[1:k])) / k, numeric(1)))
}
tags <- sub("\\.rds$", "", list.files(CACHE, "\\.rds$"))
C <- do.call(rbind, lapply(tags, cat_curve))
Cw <- reshape(C, idvar = c("tag", "n_genes"), timevar = "k", direction = "wide")
Cw <- Cw[order(-Cw$n_genes), ]
cat("\n== C. concordance at the top: fraction of the top-k shared with the canonical ranking ==\n")
print(Cw, row.names = FALSE, digits = 3)
write.csv(C, file.path(RUN, "cat_curves.csv"), row.names = FALSE)

# ===========================================================================
# D. Fixed gene set, fixed nominal threshold -- no multiple testing at all
# ===========================================================================
D <- do.call(rbind, lapply(tags, function(tg) {
  f <- get(tg); if (!isTRUE(f$ok)) return(NULL)
  cm <- intersect(f$tt$gene, refU)
  a  <- f$tt$P.Value[match(cm, f$tt$gene)]
  b  <- ref$tt$P.Value[match(cm, refU)]
  data.frame(tag = tg, n_genes = f$n_genes, n_common = length(cm),
             cell_p001 = sum(a < 0.001), canon_p001 = sum(b < 0.001),
             cell_p01 = sum(a < 0.01),   canon_p01 = sum(b < 0.01),
             ratio_p001 = sum(a < 0.001) / max(1, sum(b < 0.001)),
             stringsAsFactors = FALSE)
}))
D <- D[order(-D$n_genes), ]
cat("\n== D. same genes, same nominal cutoff: how much does the ESTIMATE move? ==\n")
print(D, row.names = FALSE, digits = 4)
write.csv(D, file.path(RUN, "fixed_nominal.csv"), row.names = FALSE)

# ===========================================================================
# E. Who are the genes the filter decides?
# ===========================================================================
loose_only <- setdiff(deg_ids(loo$tt, FDR), refDEG)
canon_only <- setdiff(refDEG, deg_ids(loo$tt, FDR))
shared     <- intersect(deg_ids(loo$tt, FDR), refDEG)
desc <- function(g, lab) data.frame(
  set = lab, n = length(g),
  med_aveexpr = if (length(g)) median(loo$tt$AveExpr[match(g, loo$tt$gene)], na.rm = TRUE) else NA,
  med_detect10 = if (length(g)) median(rowMeans(cnt[g, , drop = FALSE] > 10)) else NA,
  med_abs_lfc = if (length(g)) median(abs(loo$tt$logFC[match(g, loo$tt$gene)]), na.rm = TRUE) else NA,
  in_canon_universe = if (length(g)) mean(g %in% refU) else NA)
E <- rbind(desc(shared, "DEG in both"), desc(loose_only, "DEG only in loose fit"),
           desc(canon_only, "DEG only in canonical fit"))
cat("\n== E. character of the DEGs the filter decides ==\n")
print(E, row.names = FALSE, digits = 4)
write.csv(E, file.path(RUN, "deg_character.csv"), row.names = FALSE)
cat("\nwrote strata.csv decomposition.csv cat_curves.csv fixed_nominal.csv deg_character.csv\n")

# ===========================================================================
# F. One figure. The panel that matters is the middle one: if the DEG count
#    were bookkeeping, the fixed-universe curve would be flat and the raw
#    curve would climb with the universe. Neither happens.
# ===========================================================================
S <- read.csv(file.path(RUN, "summary.csv"), stringsAsFactors = FALSE)
G <- S[S$ok & !is.na(S$min_count), ]
G <- G[order(G$n_genes), ]
png(file.path(RUN, "filter_sweep.png"), width = 1500, height = 1150, res = 130)
op <- par(mfrow = c(2, 2), mar = c(4.2, 4.4, 2.6, 1.2), cex = 0.85)
cols <- c("#1f78b4", "#e31a1c")

plot(G$n_genes, G$n_deg_05, pch = 19, col = cols[1], log = "x",
     xlab = "genes tested (log scale)", ylab = "DEG at FDR 0.05",
     main = "DEG count vs size of testing universe", ylim = range(0, G$n_deg_05, G$n_deg_fixedU))
points(G$n_genes, G$n_deg_fixedU, pch = 1, col = cols[2])
abline(h = S$n_deg_05[S$is_canon], lty = 3)
# what a pure-bookkeeping result would look like: DEG proportional to genes
k <- S$n_deg_05[S$is_canon] / S$n_genes[S$is_canon]
lines(G$n_genes, k * G$n_genes, lty = 2, col = "grey40")
legend("topleft", bty = "n", cex = 0.8,
       legend = c("as run", "rescored on the canonical 12,173 genes",
                  "canonical result", "if DEG count were proportional to genes tested"),
       pch = c(19, 1, NA, NA), lty = c(NA, NA, 3, 2),
       col = c(cols[1], cols[2], "black", "grey40"))

plot(G$n_genes, G$pct_deg_05, pch = 19, col = cols[1], log = "x",
     xlab = "genes tested (log scale)", ylab = "DEG as % of genes tested",
     main = "The same numbers as a percentage")
abline(h = S$pct_deg_05[S$is_canon], lty = 3)

plot(G$n_genes, G$pi0, pch = 19, col = cols[1], log = "x", ylim = c(0.4, 1),
     xlab = "genes tested (log scale)", ylab = expression(pi[0]),
     main = "Estimated null fraction")
abline(h = 1, lty = 2, col = "grey60"); abline(h = S$pi0[S$is_canon], lty = 3)

plot(G$n_genes, G$jaccard_vs_canon, pch = 19, col = cols[1], log = "x", ylim = c(0, 1),
     xlab = "genes tested (log scale)", ylab = "Jaccard vs canonical DEG set",
     main = "Stability of the DEG list")
points(G$n_genes, G$frac_canon_deg_recovered, pch = 1, col = cols[2])
legend("bottomleft", bty = "n", cex = 0.8, pch = c(19, 1), col = cols,
       legend = c("Jaccard", "fraction of canonical DEGs recovered"))
par(op); dev.off()
cat("wrote runs/filter_sweep/filter_sweep.png\n")
