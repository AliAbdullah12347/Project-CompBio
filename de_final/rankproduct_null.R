#!/usr/bin/env Rscript
# ==============================================================================
# rankproduct_null.R -- the permutation null for the cross-lineage ranking.
#
# THE GAP THIS CLOSES
#
# Section 7 ranks genes by rank product across the five lineages and reports
# that the lineages agree strongly -- cross-lineage Spearman 0.556 to 0.871.
# Section 7.3 then states the obvious objection and leaves it open: a rank
# product has no built-in null, so agreement that looks impressive may be
# exactly what scrambled labels produce. Until that is measured the ranking is
# an ordering, not evidence.
#
# WHY THE ANSWER IS NOT OBVIOUS IN EITHER DIRECTION
#
# Section 10.2 established that one sample-level axis explains 71% of the
# variance in the bMIND estimates and any two genes correlate at 0.674. The
# five lineage estimates are therefore partly the SAME measurement, and two
# measurements of the same thing agree whether or not the labels mean
# anything. That predicts a null with high cross-lineage agreement -- in which
# case the observed agreement is not evidence of replication.
#
# But it could go the other way: shared structure raises the null's agreement
# without necessarily raising it to the observed level, and a genuine effect
# present in several lineages would push the observed above it.
#
# WHAT IS TESTED
#
# Three statistics, each computed identically on the real labels and on 100
# permutations of them:
#
#   1. mean Spearman correlation of the t-statistics over the 10 lineage pairs
#   2. how many genes land in the top 100 of two or more lineages
#   3. the best (smallest) rank product achieved by any gene
#
# BOTH RANKINGS ARE TESTED. Section 7 built its ranking on the RECENTRED
# t-statistics, because a constant offset changes the order of t even though
# it does not change the order of logFC. Section 10.2 then found centring to
# be poorly calibrated and recommended the raw analysis. So the raw ranking is
# tested too, and the reported conclusion will follow whichever the paper
# uses. Centring is label-blind -- each sample's own median across genes -- so
# it is legitimate to apply inside a permutation.
# ==============================================================================

suppressPackageStartupMessages({library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "de_final")); set.seed(481); options(width = 185)
RES <- "results"
NPERM <- 100
TOPN  <- 100

p  <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
bp <- readRDS(file.path(ROOT, "de_analysis/results/bmind_profiles.rds"))
A <- bp$A; SE <- bp$SE; LINS <- bp$lineages
sym <- local({ m <- read.delim(file.path(ROOT, "de_analysis/data/ens2sym.tsv"),
                 header = FALSE, col.names = c("ens","sym")); setNames(m$sym, m$ens) })

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

# the three statistics, from a genes x lineages matrix of t-statistics
stats_of <- function(TT) {
  C  <- cor(TT, method = "spearman")
  mc <- mean(C[upper.tri(C)])
  tl <- unlist(lapply(LINS, function(ct) rownames(TT)[order(-abs(TT[, ct]))][seq_len(TOPN)]))
  u  <- table(tl)
  RK <- apply(-abs(TT), 2, rank)
  c(mean_rho = mc, top_2plus = sum(u >= 2), top_3plus = sum(u >= 3),
    best_rp = min(exp(rowMeans(log(RK))))) }

cat(sprintf("== rankproduct_null ==\n%d permutations | top-%d overlap | raw and recentred rankings\n",
            NPERM, TOPN))

OUT <- list()
for (CN in c("LI", "BPD")) {
  d <- design_for(CN)
  # precompute per lineage: raw matrix, centred matrix, weights. None depends
  # on the labels, so all are fixed across permutations -- which is also what
  # makes permuting the labels a valid null here.
  PRE <- lapply(LINS, function(ct) {
    Y <- A[, ct, d$samples]
    list(raw = Y, cen = sweep(Y, 2, apply(Y, 2, median), "-"),
         W = mkw(SE[, ct, d$samples])) })
  names(PRE) <- LINS

  tmat <- function(X, which) {
    TT <- sapply(LINS, function(ct) {
      f <- eBayes(lmFit(PRE[[ct]][[which]], X, weights = PRE[[ct]]$W), trend = TRUE)
      f$t[, "grpcase"] })
    rownames(TT) <- dimnames(A)[[1]]; TT }

  obs <- list(raw = stats_of(tmat(d$X, "raw")), cen = stats_of(tmat(d$X, "cen")))
  cat(sprintf("\n%s\n== %s (n = %d) ==\n%s\n", strrep("=",185), CN, length(d$samples), strrep("=",185)))
  cat(sprintf("  observed  raw: mean rho %+.3f | %4d genes in top-%d of 2+ lineages | best rank product %.2f\n",
              obs$raw["mean_rho"], obs$raw["top_2plus"], TOPN, obs$raw["best_rp"]))
  cat(sprintf("  observed  cen: mean rho %+.3f | %4d genes in top-%d of 2+ lineages | best rank product %.2f\n",
              obs$cen["mean_rho"], obs$cen["top_2plus"], TOPN, obs$cen["best_rp"]))

  NUL <- list(raw = matrix(NA_real_, NPERM, 4), cen = matrix(NA_real_, NPERM, 4))
  t0 <- Sys.time()
  for (b in seq_len(NPERM)) {
    Xp <- d$X; Xp[, "grpcase"] <- Xp[sample(nrow(Xp)), "grpcase"]
    NUL$raw[b, ] <- stats_of(tmat(Xp, "raw"))
    NUL$cen[b, ] <- stats_of(tmat(Xp, "cen"))
    if (b %% 20 == 0) { el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      cat(sprintf("    %3d/%d  %.1f min, ETA %.1f\n", b, NPERM, el, el/b*(NPERM-b)))
      flush.console() } }

  for (w in c("raw","cen")) {
    colnames(NUL[[w]]) <- names(obs[[w]])
    cat(sprintf("\n  --- %s ranking, against %d permutations ---\n",
                if (w == "raw") "RAW" else "RECENTRED", NPERM))
    cat("   statistic                         observed | null median | null 95th pct |      p\n")
    for (k in c("mean_rho","top_2plus","top_3plus")) {
      o <- obs[[w]][k]; nn <- NUL[[w]][, k]
      pv <- (sum(nn >= o) + 1)/(NPERM + 1)
      cat(sprintf("   %-32s %8.3f | %11.3f | %13.3f | %6.4f%s\n", k, o, median(nn),
                  quantile(nn,.95,names=FALSE), pv, if (pv < 0.05) "  <- above the null" else ""))
      OUT[[length(OUT)+1]] <- data.frame(contrast = CN, ranking = w, statistic = k,
        observed = o, null_median = median(nn), null_q95 = quantile(nn,.95,names=FALSE),
        null_min = min(nn), null_max = max(nn), p = pv,
        exceeds_null = pv < 0.05, stringsAsFactors = FALSE) }
    # best rank product: SMALLER is better, so the tail is the other way
    o <- obs[[w]]["best_rp"]; nn <- NUL[[w]][, "best_rp"]
    pv <- (sum(nn <= o) + 1)/(NPERM + 1)
    cat(sprintf("   %-32s %8.3f | %11.3f | %13.3f | %6.4f%s   (smaller is better)\n",
                "best_rp", o, median(nn), quantile(nn,.05,names=FALSE), pv,
                if (pv < 0.05) "  <- better than the null" else ""))
    OUT[[length(OUT)+1]] <- data.frame(contrast = CN, ranking = w, statistic = "best_rp",
      observed = o, null_median = median(nn), null_q95 = quantile(nn,.05,names=FALSE),
      null_min = min(nn), null_max = max(nn), p = pv, exceeds_null = pv < 0.05,
      stringsAsFactors = FALSE) }
}
R <- do.call(rbind, OUT)
write.csv(R, file.path(RES, "rankproduct_null.csv"), row.names = FALSE)

cat(sprintf("\n%s\nVERDICT\n%s\n", strrep("=",185), strrep("=",185)))
for (CN in c("LI","BPD")) for (w in c("raw","cen")) {
  z <- R[R$contrast == CN & R$ranking == w, ]
  k <- sum(z$exceeds_null)
  cat(sprintf("  %-4s %-4s : %d of %d statistics exceed their null -> %s\n", CN, w, k, nrow(z),
      if (k == 0) "the cross-lineage agreement is NOT above chance"
      else if (k == nrow(z)) "agreement exceeds chance on every statistic"
      else "mixed; see the table")) }
cat("\nwrote results/rankproduct_null.csv\n")
