#!/usr/bin/env Rscript
# ==============================================================================
# sdfilter_200_compare.R -- is standard deviation the right filter statistic?
#
# sdfilter_200.R answers the question as posed: filter on SD, keep ~200, test.
# This script asks the follow-up that determines whether the idea is worth
# keeping -- given that we are going to spend an independent filter, which
# statistic should we spend it on?
#
# THE ARGUMENT FOR FILTERING AT ALL IS SOUND, AND IS NOT IN DISPUTE
#
# Fewer tests means a less punishing Benjamini-Hochberg threshold. The gain is
# real and can be large. But a filter only helps if the genes it DISCARDS are
# overwhelmingly null. If it discards true positives, the relaxed threshold is
# bought with signal, and the trade can easily be a net loss.
#
# WHY SD MAY BE THE WRONG STATISTIC
#
# Total variance = between-group variance + within-group variance. For the
# effects in this data the second term dominates by an order of magnitude: a
# gene with the largest lithium effect in the cohort (log2FC 0.79, groups split
# 152/74) gains only about 0.37 in SD from that effect, while the top-200 cut
# sits at 0.93. A gene therefore cannot reach the top 200 on group difference
# alone -- it gets there by being NOISY. Ranking on SD ranks on within-group
# noise, which is the denominator of the t-statistic, not the numerator.
#
# THE ALTERNATIVE
#
# Bourgon, Gentleman & Huber (PNAS 2010) show that overall MEAN expression is,
# for a two-group comparison, independent of the t-statistic under the null --
# the precondition SD does not cleanly satisfy -- while still being strongly
# related to power, because low-count genes carry most of the measurement
# noise. It is the standard independent filter for exactly this reason.
#
# Both are run here at identical cuts, with identical permutation validation,
# so the comparison is like-for-like and the choice is made on evidence.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "de_final")); set.seed(481); options(width = 200)
RES <- "results"

FDR   <- 0.05
CUTS  <- c(200, 500, 1000, 2000, 4000, 6000, 8000, 10000, 12368)
NPERM <- 1000
LINS  <- c("gran", "mono", "T", "NK", "B")

p  <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
bp <- readRDS(file.path(ROOT, "de_analysis/results/bmind_profiles.rds"))
sym <- local({ m <- read.delim(file.path(ROOT, "de_analysis/data/ens2sym.tsv"),
                               header = FALSE, col.names = c("ens", "sym"))
               setNames(m$sym, m$ens) })
A <- bp$A; SE <- bp$SE

design_for <- function(contrast) {
  s <- p$sel[[contrast]]
  mm <- p$meta[s$keep, , drop = FALSE]; mm$grp <- s$grp
  mm$plate <- droplevels(factor(mm$plate)); mm$sex <- droplevels(factor(mm$sex))
  if ("group" %in% s$cov) mm$group <- droplevels(factor(mm$group))
  cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]])) > 1, logical(1))]
  X <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse = " + "))), data = mm)
  ne <- nonEstimable(X); if (!is.null(ne)) X <- X[, setdiff(colnames(X), ne), drop = FALSE]
  stopifnot(qr(X)$rank == ncol(X), "grpcase" %in% colnames(X))
  list(X = X, samples = mm$title) }

mkw <- function(S) { W <- 1/(S^2); f <- is.finite(W) & W > 0
  if (any(!f)) W[!f] <- min(W[f]); W/mean(W) }

nsig <- function(Y, X, W = NULL, trend = FALSE) {
  f <- eBayes(if (is.null(W)) lmFit(Y, X) else lmFit(Y, X, weights = W), trend = trend)
  tt <- topTable(f, coef = "grpcase", number = Inf, sort.by = "none")
  q <- p.adjust(tt$P.Value, "BH")
  list(n = sum(q < FDR), q = q, p = tt$P.Value, lfc = tt$logFC, genes = rownames(tt)) }

cat(sprintf("== sdfilter_200_compare ==\nSD vs mean-expression filtering | BH %.2f | %d permutations at the 200 cut\n",
            FDR, NPERM))

CURVE <- list(); VAL <- list(); GEN <- list()

## =============================================================== WHOLE BLOOD
for (CN in c("LI", "BPD")) {
  d <- design_for(CN)
  dge <- calcNormFactors(DGEList(p$counts[, d$samples, drop = FALSE]), method = "TMM")
  v <- voom(dge, d$X); E <- v$E; VW <- v$weights; dimnames(VW) <- dimnames(E)
  base <- nsig(E, d$X, VW)
  base_sig <- base$genes[base$q < FDR]

  stat <- list(SD = apply(E, 1, sd), MEAN = rowMeans(E))
  cat(sprintf("\n%s\n== WHOLE BLOOD %s (n = %d) -- unfiltered: %d significant ==\n%s\n",
              strrep("=", 160), CN, length(d$samples), base$n, strrep("=", 160)))
  cat("  filter | genes kept | significant | true-positive recall (of the unfiltered hits kept) | enrichment vs chance\n")

  for (FS in names(stat)) {
    ord <- order(stat[[FS]], decreasing = TRUE)
    for (m in CUTS) {
      keep <- ord[seq_len(min(m, nrow(E)))]
      r <- nsig(E[keep, , drop = FALSE], d$X, VW[keep, , drop = FALSE])
      kept_tp <- sum(base_sig %in% rownames(E)[keep])
      # If the filter were blind to the signal, it would keep m/12368 of the
      # true positives. Above 1 means it concentrates signal; below 1 means it
      # is throwing signal away faster than it throws away genes.
      enr <- (kept_tp / max(1, length(base_sig))) / (length(keep) / nrow(E))
      cat(sprintf("  %-6s | %10d | %11d | %5d of %-5d | %5.2fx\n",
                  FS, length(keep), r$n, kept_tp, length(base_sig), enr))
      CURVE[[length(CURVE)+1]] <- data.frame(
        level = "WB", contrast = CN, lineage = NA_character_, filter = FS,
        n_kept = length(keep), n_sig = r$n, base_sig = length(base_sig),
        kept_true_pos = kept_tp, enrichment = enr, stringsAsFactors = FALSE) }
    cat("\n") }

  ## permutation validation of BOTH filters at the 200 cut
  for (FS in names(stat)) {
    keep <- order(stat[[FS]], decreasing = TRUE)[seq_len(200)]
    Ys <- E[keep, , drop = FALSE]; Ws <- VW[keep, , drop = FALSE]
    obs <- nsig(Ys, d$X, Ws)$n
    NU <- integer(NPERM)
    for (b in seq_len(NPERM)) {
      Xp <- d$X; Xp[, "grpcase"] <- Xp[sample(nrow(Xp)), "grpcase"]
      NU[b] <- nsig(Ys, Xp, Ws)$n }
    pv <- (sum(NU >= obs) + 1) / (NPERM + 1)
    cat(sprintf("  %-6s 200-gene cut: observed %3d | null median %2.0f, q95 %4.1f, max %4d | p = %.4f\n",
                FS, obs, median(NU), quantile(NU, .95, names = FALSE), max(NU), pv))
    VAL[[length(VAL)+1]] <- data.frame(
      level = "WB", contrast = CN, lineage = NA_character_, filter = FS, n_kept = 200,
      observed = obs, null_median = median(NU), null_q95 = quantile(NU, .95, names = FALSE),
      null_max = max(NU), perm_p = pv, stringsAsFactors = FALSE) }

  ## the best cut each filter can reach, and the genes it buys
  for (FS in names(stat)) {
    ord <- order(stat[[FS]], decreasing = TRUE)
    best <- CURVE[vapply(CURVE, function(z) z$level=="WB" && z$contrast==CN && z$filter==FS, logical(1))]
    bb <- do.call(rbind, best); bb <- bb[which.max(bb$n_sig), ]
    keep <- ord[seq_len(bb$n_kept)]
    r <- nsig(E[keep, , drop = FALSE], d$X, VW[keep, , drop = FALSE])
    gg <- data.frame(level = "WB", contrast = CN, lineage = NA_character_, filter = FS,
                     n_kept = bb$n_kept, gene = r$genes,
                     symbol = unname(ifelse(is.na(sym[r$genes]), r$genes, sym[r$genes])),
                     logFC = r$lfc, P = r$p, BH = r$q, significant = r$q < FDR,
                     novel = !(r$genes %in% base_sig), stringsAsFactors = FALSE)
    GEN[[paste("WB", CN, FS, sep = "_")]] <- gg[order(gg$P), ]
    s <- gg[gg$significant, ]
    cat(sprintf("  %-6s best cut = %5d genes -> %4d significant (%d of them NOT found unfiltered)\n",
                FS, bb$n_kept, nrow(s), sum(s$novel))) }
}

## ================================================================ CELL TYPE
for (CN in c("LI", "BPD")) {
  d <- design_for(CN)
  cat(sprintf("\n%s\n== CELL TYPE %s (n = %d) ==\n%s\n",
              strrep("=", 160), CN, length(d$samples), strrep("=", 160)))
  cat("  lineage | filter | unfiltered | top-200 | null q95 | null max | perm p | recall of unfiltered hits\n")
  for (ct in LINS) {
    Y <- A[, ct, d$samples]; W <- mkw(SE[, ct, d$samples])
    base <- nsig(Y, d$X, W, trend = TRUE); base_sig <- base$genes[base$q < FDR]
    stat <- list(SD = apply(Y, 1, sd), MEAN = rowMeans(Y))
    for (FS in names(stat)) {
      keep <- order(stat[[FS]], decreasing = TRUE)[seq_len(200)]
      Ys <- Y[keep, , drop = FALSE]; Ws <- W[keep, , drop = FALSE]
      obs <- nsig(Ys, d$X, Ws, trend = TRUE)$n
      NU <- integer(NPERM)
      for (b in seq_len(NPERM)) {
        Xp <- d$X; Xp[, "grpcase"] <- Xp[sample(nrow(Xp)), "grpcase"]
        NU[b] <- nsig(Ys, Xp, Ws, trend = TRUE)$n }
      pv <- (sum(NU >= obs) + 1) / (NPERM + 1)
      kept_tp <- sum(base_sig %in% rownames(Y)[keep])
      cat(sprintf("  %-7s | %-6s | %10d | %7d | %8.1f | %8d | %6.4f | %d of %d\n",
                  ct, FS, base$n, obs, quantile(NU, .95, names = FALSE), max(NU), pv,
                  kept_tp, length(base_sig)))
      VAL[[length(VAL)+1]] <- data.frame(
        level = "CT", contrast = CN, lineage = ct, filter = FS, n_kept = 200,
        observed = obs, null_median = median(NU), null_q95 = quantile(NU, .95, names = FALSE),
        null_max = max(NU), perm_p = pv, stringsAsFactors = FALSE)
      CURVE[[length(CURVE)+1]] <- data.frame(
        level = "CT", contrast = CN, lineage = ct, filter = FS, n_kept = 200,
        n_sig = obs, base_sig = length(base_sig), kept_true_pos = kept_tp,
        enrichment = (kept_tp / max(1, length(base_sig))) / (200 / nrow(Y)),
        stringsAsFactors = FALSE) } }
}

CU <- do.call(rbind, CURVE); VV <- do.call(rbind, VAL); GG <- do.call(rbind, GEN)
write.csv(CU, file.path(RES, "sdfilter200_compare_curve.csv"), row.names = FALSE)
write.csv(VV, file.path(RES, "sdfilter200_compare_validation.csv"), row.names = FALSE)
write.csv(GG, gzfile(file.path(RES, "sdfilter200_compare_genes.csv.gz")), row.names = FALSE)
cat("\nwrote results/sdfilter200_compare_*.csv\n")
