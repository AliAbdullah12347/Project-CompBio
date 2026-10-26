#!/usr/bin/env Rscript
# ==============================================================================
# meanfilter.R -- mean-expression filtering, the principled alternative to SD.
#
# WHERE THIS COMES FROM
#
# sdfilter_200.R established that filtering on standard deviation is valid but
# harmful: it cut the lithium whole-blood result from 1,426 genes to 10 and
# discarded the three strongest genes in the study, because total variance is
# dominated by WITHIN-group noise rather than between-group difference, so
# ranking on it ranks on the denominator of the t-statistic.
#
# The underlying idea -- spend an independent filter to buy a less punishing
# multiple-testing correction -- is still sound. It just needs a filter
# statistic that is related to POWER rather than to noise.
#
# WHY MEAN EXPRESSION IS THE TEXTBOOK CHOICE
#
# Bourgon, Gentleman & Huber (PNAS 2010) show that for a two-group comparison
# the overall mean is independent of the t-statistic under the null -- the
# precondition SD only approximately satisfies -- while still being strongly
# related to power, because low-count genes carry most of the measurement
# noise. Discarding them removes tests that had little chance of succeeding
# without removing candidates that did.
#
# A PREDICTION, STATED BEFORE LOOKING
#
# We expect the gain here to be SMALL, for a reason specific to this pipeline:
#   1. The baseline filter (>10 counts in >=90% of samples) has already removed
#      the low-expression genes. The 12,368 that remain are the well-expressed
#      ones, so the pool mean filtering normally targets is largely gone.
#   2. voom already down-weights low-count observations through the
#      mean-variance trend, so part of the benefit is being collected twice.
# If mean filtering still helps, it helps on top of both of those.
#
# TWO ORDERS, BOTH REPORTED, BECAUSE THEY ISOLATE DIFFERENT MECHANISMS
#
#   filter BEFORE fitting -- subset, then eBayes. Matches what sdfilter_200.R
#       did, so the two filters are comparable. Changes the variance moderation
#       as well as the correction, because eBayes borrows across whichever
#       genes are present.
#   filter AFTER fitting  -- eBayes on all 12,368, then keep the subset's
#       p-values and apply BH within it. This is what DESeq2's independent
#       filtering does, and it isolates the PURE multiple-testing benefit,
#       which is the actual claim being tested.
#
# This script runs the cheap diagnostics only. Permutation nulls are run
# afterwards, and only where a filter appears to gain something -- there is
# nothing to defend where it does not.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "de_final")); set.seed(481); options(width = 200)
RES <- "results"
FDR  <- 0.05
CUTS <- c(200, 500, 1000, 2000, 4000, 6000, 8000, 10000, 12368)
LINS <- c("gran", "mono", "T", "NK", "B")

p  <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
bp <- readRDS(file.path(ROOT, "de_analysis/results/bmind_profiles.rds"))
A <- bp$A; SE <- bp$SE
sym <- local({ m <- read.delim(file.path(ROOT, "de_analysis/data/ens2sym.tsv"),
                 header = FALSE, col.names = c("ens","sym")); setNames(m$sym, m$ens) })
LM22G <- toupper(trimws(read.delim(file.path(ROOT, "data/LM22.txt"), check.names = FALSE)[[1]]))
ALLSYM <- toupper(ifelse(is.na(sym[dimnames(A)[[1]]]), dimnames(A)[[1]], sym[dimnames(A)[[1]]]))
BG <- mean(ALLSYM %in% LM22G)

design_for <- function(cn) {
  s <- p$sel[[cn]]; mm <- p$meta[s$keep, , drop = FALSE]; mm$grp <- s$grp
  mm$plate <- droplevels(factor(mm$plate)); mm$sex <- droplevels(factor(mm$sex))
  if ("group" %in% s$cov) mm$group <- droplevels(factor(mm$group))
  cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]])) > 1, logical(1))]
  X <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse = " + "))), data = mm)
  ne <- nonEstimable(X); if (!is.null(ne)) X <- X[, setdiff(colnames(X), ne), drop = FALSE]
  list(X = X, samples = mm$title) }
mkw <- function(S) { W <- 1/(S^2); f <- is.finite(W) & W > 0
  if (any(!f)) W[!f] <- min(W[f]); W/mean(W) }
pvals <- function(Y, X, W, trend) {
  tt <- topTable(eBayes(lmFit(Y, X, weights = W), trend = trend),
                 coef = "grpcase", number = Inf, sort.by = "none")
  setNames(tt$P.Value, rownames(tt)) }

cat(sprintf("== meanfilter ==\nSD vs mean-expression filtering | BH %.2f\nLM22 background %.2f%% of the %d tested genes\n",
            FDR, 100*BG, length(ALLSYM)))

CURVE <- list(); DIAG <- list()

run_one <- function(level, CN, lineage, Y, X, W, trend) {
  # per-gene p-values on the FULL gene set, computed once
  pfull <- pvals(Y, X, W, trend)
  base_n <- sum(p.adjust(pfull, "BH") < FDR)
  base_sig <- names(pfull)[p.adjust(pfull, "BH") < FDR]
  stat <- list(SD = apply(Y, 1, sd), MEAN = rowMeans(Y))

  # independence diagnostic: under PERMUTED labels, does the filter statistic
  # predict significance? It must not.
  Xp <- X; Xp[, "grpcase"] <- Xp[sample(nrow(Xp)), "grpcase"]
  pperm <- pvals(Y, Xp, W, trend)
  for (FS in names(stat))
    DIAG[[length(DIAG)+1]] <<- data.frame(level = level, contrast = CN,
      lineage = lineage, filter = FS,
      rho_null = cor(stat[[FS]], -log10(pperm), method = "spearman"),
      rho_obs  = cor(stat[[FS]], -log10(pfull), method = "spearman"),
      stringsAsFactors = FALSE)

  for (FS in names(stat)) {
    ord <- order(stat[[FS]], decreasing = TRUE)
    for (m in CUTS) {
      k <- ord[seq_len(min(m, nrow(Y)))]
      g <- rownames(Y)[k]
      n_after  <- sum(p.adjust(pfull[g], "BH") < FDR)                 # pure MTC effect
      n_before <- sum(p.adjust(pvals(Y[k,,drop=FALSE], X,
                      W[k,,drop=FALSE], trend), "BH") < FDR)          # refit on the subset
      ks <- toupper(ifelse(is.na(sym[g]), g, sym[g]))
      CURVE[[length(CURVE)+1]] <<- data.frame(level = level, contrast = CN,
        lineage = lineage, filter = FS, n_kept = length(k),
        n_sig_filter_after_fit = n_after, n_sig_filter_before_fit = n_before,
        base_n = base_n, kept_true_pos = sum(base_sig %in% g),
        tp_enrichment = (sum(base_sig %in% g)/max(1,length(base_sig)))/(length(k)/nrow(Y)),
        lm22_frac = mean(ks %in% LM22G), lm22_enrichment = mean(ks %in% LM22G)/BG,
        stringsAsFactors = FALSE) } }
  base_n }

## ---------------------------------------------------------------- whole blood
for (CN in c("LI", "BPD")) {
  d <- design_for(CN)
  v <- voom(calcNormFactors(DGEList(p$counts[, d$samples]), method = "TMM"), d$X)
  VW <- v$weights; dimnames(VW) <- dimnames(v$E)
  bn <- run_one("WB", CN, NA_character_, v$E, d$X, VW, FALSE)
  cat(sprintf("\n%s\n== WHOLE BLOOD %s -- unfiltered %d significant ==\n%s\n",
              strrep("=",168), CN, bn, strrep("=",168)))
  cat("  filter | kept  | sig (filter after fit) | sig (filter before fit) | true pos kept | TP enrich | LM22 enrich\n")
  for (FS in c("MEAN","SD")) {
    for (r in CURVE) if (r$level=="WB" && r$contrast==CN && r$filter==FS)
      cat(sprintf("  %-6s | %5d | %22d | %23d | %5d of %-5d | %8.2fx | %9.2fx\n",
          r$filter, r$n_kept, r$n_sig_filter_after_fit, r$n_sig_filter_before_fit,
          r$kept_true_pos, r$base_n, r$tp_enrichment, r$lm22_enrichment))
    cat("\n") } }

## ----------------------------------------------------------------- cell type
for (CN in c("LI", "BPD")) {
  d <- design_for(CN)
  cat(sprintf("\n%s\n== CELL TYPE %s ==\n%s\n", strrep("=",168), CN, strrep("=",168)))
  cat("  lineage | filter | kept  | sig (after fit) | sig (before fit) | TP enrich | LM22 enrich\n")
  for (ct in LINS) {
    Y <- A[, ct, d$samples]; W <- mkw(SE[, ct, d$samples])
    bn <- run_one("CT", CN, ct, Y, d$X, W, TRUE)
    for (FS in c("MEAN","SD"))
      for (r in CURVE) if (r$level=="CT" && r$contrast==CN && r$lineage==ct && r$filter==FS &&
                           r$n_kept %in% c(200, 2000, 6000, 12368))
        cat(sprintf("  %-7s | %-6s | %5d | %15d | %16d | %8.2fx | %9.2fx\n",
            ct, r$filter, r$n_kept, r$n_sig_filter_after_fit, r$n_sig_filter_before_fit,
            r$tp_enrichment, r$lm22_enrichment)) } }

CU <- do.call(rbind, CURVE); DG <- do.call(rbind, DIAG)
write.csv(CU, file.path(RES, "meanfilter_curve.csv"), row.names = FALSE)
write.csv(DG, file.path(RES, "meanfilter_independence.csv"), row.names = FALSE)

cat(sprintf("\n%s\nINDEPENDENCE CHECK -- Spearman(filter statistic, -log10 p) under PERMUTED labels\n%s\n",
            strrep("=",168), strrep("=",168)))
agg <- aggregate(cbind(rho_null, rho_obs) ~ level + filter, data = DG,
                 FUN = function(z) round(mean(z), 4))
print(agg, row.names = FALSE)
cat("  (near zero = the precondition for independent filtering is satisfied)\n")

cat(sprintf("\n%s\nDOES EITHER FILTER EVER BEAT NO FILTER?\n%s\n", strrep("=",168), strrep("=",168)))
for (lv in c("WB","CT")) for (FS in c("MEAN","SD")) {
  z <- CU[CU$level==lv & CU$filter==FS, ]
  unf <- z[z$n_kept==12368, ]
  beat_a <- sum(z$n_sig_filter_after_fit  > rep(unf$n_sig_filter_after_fit,  each = length(CUTS)))
  beat_b <- sum(z$n_sig_filter_before_fit > rep(unf$n_sig_filter_before_fit, each = length(CUTS)))
  cat(sprintf("  %-3s %-5s : %2d of %2d cutoffs beat unfiltered (filter after fit) | %2d of %2d (before fit)\n",
              lv, FS, beat_a, nrow(z), beat_b, nrow(z))) }
cat("\nwrote results/meanfilter_curve.csv and meanfilter_independence.csv\n")
