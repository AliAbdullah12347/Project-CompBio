#!/usr/bin/env Rscript
# ==============================================================================
# meeting_followup.R -- the three checks raised at the 18 Oct check-in, run on
# the cell-type results for BOTH contrasts, each reported separately.
#
#   A. Recentre the bMIND fold changes so the median is zero, and ask whether
#      gene-specific signal appears once the lineage-wide offset is removed.
#   B. Rank genes instead of thresholding them, and ask whether the same genes
#      rank highly in more than one lineage.
#   C. Filter genes on across-sample standard deviation, and choose the
#      threshold by a stated rule rather than by eye.
#
# WHY A IS WORTH DOING -- the diagnosis that motivates it
#
# Under a null where most genes do not change, roughly half should move up.
# Whole blood behaves: 44.7% up for lithium, 49.1% for bipolar. The cell-type
# estimates do not:
#
#     CT_LI  granulocytes   median logFC -0.118    7.6% of genes up
#     CT_BPD T cells        median logFC -0.066   14.8% up
#     CT_BPD B cells        median logFC +0.012   94.9% up
#
# In granulocytes the lineage-wide shift (-0.118) is LARGER than the entire
# gene-to-gene spread (SD 0.090). That is not biology: a real biological effect
# does not move every gene in a lineage by the same amount. It is an additive
# per-sample offset in the deconvolution that happens to correlate with group
# -- unsurprising, since bMIND's estimate for a lineage depends on that
# sample's fraction for that lineage, and both contrasts move those fractions.
#
# WHAT RECENTRING DOES TO THE HYPOTHESIS. Removing the offset changes the
# question from "does this gene differ between groups?" to "does this gene
# differ between groups MORE THAN THE TYPICAL GENE IN ITS LINEAGE DOES?" That
# is a different and, here, better-posed question -- but it must be stated,
# not slipped in.
# ==============================================================================

suppressPackageStartupMessages({library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "de_final")); set.seed(481); options(width = 185)
RES <- "results"
FDR <- 0.05

p  <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
bp <- readRDS(file.path(ROOT, "de_analysis/results/bmind_profiles.rds"))
A <- bp$A; SE <- bp$SE; LINS <- bp$lineages
sym <- local({ m <- read.delim(file.path(ROOT, "de_analysis/data/ens2sym.tsv"),
                               header = FALSE, col.names = c("ens", "sym"))
               setNames(m$sym, m$ens) })
GENES <- dimnames(A)[[1]]

design_for <- function(contrast) {
  s <- p$sel[[contrast]]
  mm <- p$meta[s$keep, , drop = FALSE]; mm$grp <- s$grp
  mm$plate <- droplevels(factor(mm$plate)); mm$sex <- droplevels(factor(mm$sex))
  if ("group" %in% s$cov) mm$group <- droplevels(factor(mm$group))
  cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]])) > 1, logical(1))]
  X <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse = " + "))), data = mm)
  stopifnot(qr(X)$rank == ncol(X))
  list(X = X, samples = mm$title) }

mkw <- function(S) { W <- 1/(S^2); f <- is.finite(W) & W > 0
  if (any(!f)) W[!f] <- min(W[f]); W/mean(W) }

# One lineage fit. `centre` chooses the correction:
#   "none"   fit the bMIND estimates as they are
#   "sample" subtract each SAMPLE's median across genes before fitting (A1)
fit_lin <- function(ct, d, centre = "none") {
  Y <- A[, ct, d$samples]
  if (centre == "sample") Y <- sweep(Y, 2, apply(Y, 2, median), "-")
  f <- eBayes(lmFit(Y, d$X, weights = mkw(SE[, ct, d$samples])), trend = TRUE)
  tt <- topTable(f, coef = "grpcase", number = Inf, sort.by = "none")
  data.frame(gene = rownames(tt), logFC = tt$logFC, t = tt$t, P = tt$P.Value,
             SE = f$stdev.unscaled[, "grpcase"] * sqrt(f$s2.post),
             df = f$df.total, stringsAsFactors = FALSE)[match(GENES, rownames(tt)), ] }

sumline <- function(x, q) sprintf("median %+7.4f | %4.1f%% up | %5d at BH %.2f (%d up, %d down)",
  median(x$logFC), 100*mean(x$logFC > 0), sum(q < FDR), FDR,
  sum(q < FDR & x$logFC > 0), sum(q < FDR & x$logFC < 0))

ALL <- list()
for (CN in c("LI", "BPD")) {
  d <- design_for(CN)
  cat(sprintf("\n%s\n== CONTRAST: %s  (n = %d) ==\n%s\n", strrep("=", 110), CN,
              length(d$samples), strrep("=", 110)))

  ## ======================= PART A: recentring ==============================
  cat("\n---- PART A: recentre the fold changes ----\n")
  cat("Three versions of the same fit, differing only in how the lineage-wide offset is handled.\n")
  cat("  raw      : bMIND estimates as they are\n")
  cat("  post-hoc : t recomputed as (logFC - median logFC) / SE, same SEs  [diagnostic]\n")
  cat("  per-sample: each sample's median across genes removed BEFORE fitting  [principled]\n\n")
  AOUT <- list()
  for (ct in LINS) {
    r0 <- fit_lin(ct, d, "none")
    q0 <- p.adjust(r0$P, "BH")
    # post-hoc recentring: shift the coefficient, keep SE and df, recompute
    sh <- median(r0$logFC)
    t1 <- (r0$logFC - sh) / r0$SE
    p1 <- 2 * pt(-abs(t1), df = r0$df)
    q1 <- p.adjust(p1, "BH")
    # principled: centre the data, refit
    r2 <- fit_lin(ct, d, "sample")
    q2 <- p.adjust(r2$P, "BH")
    cat(sprintf("  %-5s raw       %s\n", ct, sumline(r0, q0)))
    cat(sprintf("  %-5s post-hoc  %s   (shift removed %+.4f)\n", "",
        sumline(data.frame(logFC = r0$logFC - sh), q1), sh))
    cat(sprintf("  %-5s per-sample%s\n\n", "", sumline(r2, q2)))
    AOUT[[ct]] <- data.frame(contrast = CN, lineage = ct,
      median_raw = round(median(r0$logFC), 5), pct_up_raw = round(100*mean(r0$logFC > 0), 1),
      deg_raw = sum(q0 < FDR),
      shift_removed = round(sh, 5),
      deg_posthoc = sum(q1 < FDR),
      deg_posthoc_up = sum(q1 < FDR & (r0$logFC - sh) > 0),
      deg_posthoc_down = sum(q1 < FDR & (r0$logFC - sh) < 0),
      median_centred = round(median(r2$logFC), 5), pct_up_centred = round(100*mean(r2$logFC > 0), 1),
      deg_centred = sum(q2 < FDR),
      deg_centred_up = sum(q2 < FDR & r2$logFC > 0),
      deg_centred_down = sum(q2 < FDR & r2$logFC < 0),
      stringsAsFactors = FALSE)
    ALL[[paste(CN, ct, "raw", sep = "_")]] <- r0
    ALL[[paste(CN, ct, "cen", sep = "_")]] <- r2
  }
  AT <- do.call(rbind, AOUT)
  cat("  SUMMARY -- does removing the offset reveal genes?\n")
  print(AT[, c("lineage", "median_raw", "pct_up_raw", "deg_raw", "shift_removed",
               "deg_posthoc", "deg_centred", "pct_up_centred")], row.names = FALSE)
  write.csv(AT, file.path(RES, sprintf("followup_A_recentre_%s.csv", CN)), row.names = FALSE)

  ## ======================= PART B: ranking =================================
  cat("\n---- PART B: rank genes instead of thresholding them ----\n")
  cat("Ranking by the recentred t-statistic, which is both offset-free and precision-weighted.\n")
  cat("A constant offset does not change the ORDER of logFC, but it does change the order of t,\n")
  cat("because t divides by a gene-specific SE. So the ranking is taken after recentring.\n\n")
  TT <- sapply(LINS, function(ct) ALL[[paste(CN, ct, "cen", sep = "_")]]$t)
  rownames(TT) <- GENES
  cat("  Spearman correlation of t between lineages (observed):\n")
  print(round(cor(TT, method = "spearman"), 3))

  # Top-N overlap, observed
  cat("\n  genes in the top N of two or more lineages (|t|):\n")
  topn <- function(M, N) lapply(LINS, function(ct) rownames(M)[order(-abs(M[, ct]))][1:N])
  OV <- do.call(rbind, lapply(c(50, 100, 500), function(N) {
    tl <- topn(TT, N); u <- table(unlist(tl))
    data.frame(topN = N, in_2plus = sum(u >= 2), in_3plus = sum(u >= 3),
               in_all5 = sum(u == 5), stringsAsFactors = FALSE) }))
  print(OV, row.names = FALSE)

  # Rank product across lineages: a gene's combined evidence without any threshold
  RK <- apply(-abs(TT), 2, rank)
  rp <- exp(rowMeans(log(RK)))                       # geometric mean of ranks
  ord <- order(rp)
  cat("\n  top 15 genes by rank product across all five lineages:\n")
  tab <- data.frame(symbol = unname(ifelse(is.na(sym[GENES[ord]]), GENES[ord], sym[GENES[ord]])),
                    rank_product = round(rp[ord], 1),
                    round(RK[ord, , drop = FALSE]), stringsAsFactors = FALSE)[1:15, ]
  print(tab, row.names = FALSE)
  write.csv(data.frame(gene = GENES, symbol = unname(ifelse(is.na(sym[GENES]), GENES, sym[GENES])),
                       rank_product = round(rp, 2), RK)[order(rp), ],
            gzfile(file.path(RES, sprintf("followup_B_rankproduct_%s.csv.gz", CN))), row.names = FALSE)
  saveRDS(list(t = TT, rankprod = rp), file.path(RES, sprintf("followup_B_%s.rds", CN)))

  ## ======================= PART C: SD filter ===============================
  cat("\n---- PART C: filter on across-sample standard deviation ----\n")
  cat("The filter statistic is each gene's SD across all samples within a lineage, computed\n")
  cat("IGNORING GROUP. That independence is what makes filtering legitimate (Bourgon 2010):\n")
  cat("it cannot bias the group test under the null. The check is printed below.\n\n")
  SDm <- sapply(LINS, function(ct) apply(A[, ct, d$samples], 1, sd))
  rownames(SDm) <- GENES
  # independence check: among genes with no signal, SD must not predict the p-value
  cat("  independence check -- Spearman(SD, p-value) among genes with p > 0.5:\n")
  for (ct in LINS) {
    pv <- ALL[[paste(CN, ct, "cen", sep = "_")]]$P; k <- pv > 0.5
    cat(sprintf("    %-5s rho = %+.4f  (n = %d)   %s\n", ct,
        cor(SDm[k, ct], pv[k], method = "spearman"), sum(k),
        if (abs(cor(SDm[k, ct], pv[k], method = "spearman")) < 0.1) "independent" else "NOT independent -- filter invalid")) }

  cat("\n  discoveries retained at each SD threshold (BH recomputed on survivors):\n")
  QS <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  CT <- do.call(rbind, lapply(LINS, function(ct) {
    pv <- ALL[[paste(CN, ct, "cen", sep = "_")]]$P
    do.call(rbind, lapply(QS, function(qq) {
      thr <- quantile(SDm[, ct], qq, names = FALSE); keep <- SDm[, ct] >= thr
      data.frame(contrast = CN, lineage = ct, quantile_cut = qq,
                 sd_threshold = signif(thr, 4), n_kept = sum(keep),
                 deg = sum(p.adjust(pv[keep], "BH") < FDR), stringsAsFactors = FALSE) })) }))
  W <- reshape(CT[, c("lineage", "quantile_cut", "deg")], idvar = "quantile_cut",
               timevar = "lineage", direction = "wide")
  names(W) <- sub("deg\\.", "", names(W))
  W$pooled <- vapply(QS, function(qq) {
    pv <- unlist(lapply(LINS, function(ct) {
      k <- SDm[, ct] >= quantile(SDm[, ct], qq, names = FALSE)
      setNames(ALL[[paste(CN, ct, "cen", sep = "_")]]$P[k], GENES[k]) }))
    q <- p.adjust(pv, "BH"); length(unique(names(pv)[q < FDR])) }, integer(1))
  print(W, row.names = FALSE)
  best <- QS[which.max(W$pooled)]
  cat(sprintf("\n  threshold maximising pooled unique genes: SD quantile %.1f (%d genes, vs %d unfiltered)\n",
              best, max(W$pooled), W$pooled[QS == 0]))
  write.csv(CT, file.path(RES, sprintf("followup_C_sdfilter_%s.csv", CN)), row.names = FALSE)
}
cat("\nwrote results/followup_{A,B,C}_*.csv\n")
