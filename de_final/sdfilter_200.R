#!/usr/bin/env Rscript
# ==============================================================================
# sdfilter_200.R -- variance-based gene filtering, cut to ~200 genes.
#
# THE IDEA (as proposed)
#
# A gene whose expression barely varies across the whole cohort cannot differ
# much between two subgroups of that cohort either, so testing it spends a
# multiple-testing penalty on a gene that cannot win. Drop those genes and the
# surviving ones face a far smaller correction: with Benjamini-Hochberg the
# threshold for the top-ranked gene is alpha/m, so going from m = 12,368 to
# m = 200 relaxes it by a factor of 62. Genuine but modest effects that were
# buried under the correction should surface.
#
# The logic is sound and it has a name -- independent filtering (Bourgon,
# Gentleman & Huber, PNAS 2010). It also has a precise precondition, which is
# the reason this script is as long as it is.
#
# THE PRECONDITION
#
# The filter statistic must be independent of the test statistic UNDER THE NULL.
# Overall variance is not obviously so: total variance = within-group variance
# + between-group variance, and the between-group part is the numerator of the
# very t-statistic being filtered on. Selecting high-variance genes therefore
# pre-selects, weakly, for genes with group differences. With n = 226 and small
# effects that contamination is slight -- but we are discarding 98.4% of genes,
# which is exactly the regime where a slight bias is amplified.
#
# So this script does not assume the filter is valid. It measures it, at the
# 200-gene threshold itself, against 1,000 label permutations.
#
# THE SECOND HAZARD, SPECIFIC TO CELL TYPE
#
# diag_null_mechanism.R established that in the bMIND estimates a single
# sample-level axis explains 71.4% of all variance and any two genes correlate
# at 0.674. The highest-variance genes in that matrix are plausibly just the
# genes most loaded on that nuisance axis. Ranking by SD could therefore select
# *for* the artifact. That is measured too.
#
# WHAT IS FILTERED, AND WHEN
#
# The SD is computed AFTER normalisation and BEFORE testing, across all samples
# in the contrast, never using the group labels -- that is what makes it
# label-blind and legitimate to hold fixed inside a permutation.
#
# Normalisation (TMM, voom) is computed on the full 12,368 genes and only then
# subset. TMM scaling factors estimated from 200 genes would be unstable, and
# re-running voom on a high-variance subset would mis-fit the mean-variance
# trend, since that subset is not representative of it.
#
# One real cost is accepted and reported: eBayes borrows variance information
# across genes, and with 200 genes there is far less to borrow. The moderation
# is weaker than in the 12,368-gene analysis.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "de_final")); set.seed(481); options(width = 200)
RES <- "results"

FDR    <- 0.05
NTOP   <- 200                                  # the requested target
SWEEP  <- c(50, 100, 200, 500, 1000, 2000, 5000, 12368)
NPERM  <- 1000
LINS   <- c("gran", "mono", "T", "NK", "B")

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

# One fit. Y and W are already subset to the genes being tested.
nsig <- function(Y, X, W = NULL, trend = FALSE) {
  f <- eBayes(if (is.null(W)) lmFit(Y, X) else lmFit(Y, X, weights = W), trend = trend)
  tt <- topTable(f, coef = "grpcase", number = Inf, sort.by = "none")
  q <- p.adjust(tt$P.Value, "BH")
  list(n = sum(q < FDR), q = q, p = tt$P.Value, lfc = tt$logFC, genes = rownames(tt)) }

cat(sprintf("== sdfilter_200 ==\ntarget %d genes | BH %.2f | %d label permutations\n",
            NTOP, FDR, NPERM))

SWEEPOUT <- list(); VALOUT <- list(); GENEOUT <- list(); DIAG <- list()

## =============================================================== WHOLE BLOOD
for (CN in c("LI", "BPD")) {
  d <- design_for(CN)
  dge <- calcNormFactors(DGEList(p$counts[, d$samples, drop = FALSE]), method = "TMM")
  v   <- voom(dge, d$X)                     # trend fitted on ALL 12,368 genes
  E   <- v$E
  VW  <- v$weights                          # voom's per-observation precision
  dimnames(VW) <- dimnames(E)               # lmFit(v$E, X) would silently drop these
  sdg <- apply(E, 1, sd)                    # label-blind filter statistic
  ord <- order(sdg, decreasing = TRUE)

  base <- nsig(E, d$X, VW)                  # the unfiltered reference
  cat(sprintf("\n%s\n== WHOLE BLOOD %s (n = %d) ==\n%s\n",
              strrep("=", 150), CN, length(d$samples), strrep("=", 150)))
  cat(sprintf("  unfiltered: %d of %d genes significant\n", base$n, nrow(E)))

  cat("\n  --- how the count moves as the filter tightens ---\n")
  cat("   genes kept | SD cutoff | significant | of the unfiltered hits, how many survive the filter\n")
  for (m in SWEEP) {
    keep <- ord[seq_len(min(m, nrow(E)))]
    r <- nsig(E[keep, , drop = FALSE], d$X, VW[keep, , drop = FALSE])
    base_sig <- base$genes[base$q < FDR]
    retained <- sum(base_sig %in% rownames(E)[keep])
    cat(sprintf("   %10d | %9.4f | %11d | %d of %d\n", length(keep),
                min(sdg[keep]), r$n, retained, length(base_sig)))
    SWEEPOUT[[length(SWEEPOUT) + 1]] <- data.frame(
      level = "WB", contrast = CN, lineage = NA_character_, n_kept = length(keep),
      sd_cutoff = min(sdg[keep]), n_sig = r$n,
      base_sig = length(base_sig), base_sig_retained = retained, stringsAsFactors = FALSE) }

  ## ---- the requested cut, in detail ---------------------------------------
  keep <- ord[seq_len(NTOP)]
  r <- nsig(E[keep, , drop = FALSE], d$X, VW[keep, , drop = FALSE])
  g <- data.frame(level = "WB", contrast = CN, lineage = NA_character_,
                  gene = r$genes, symbol = unname(ifelse(is.na(sym[r$genes]), r$genes, sym[r$genes])),
                  sd = sdg[r$genes], logFC = r$lfc, P = r$p, BH = r$q,
                  significant = r$q < FDR, stringsAsFactors = FALSE)
  GENEOUT[[paste0("WB_", CN)]] <- g[order(g$P), ]

  ## ---- is the filter independent of the test statistic under the null? -----
  cat(sprintf("\n  --- validating the %d-gene cut against %d label permutations ---\n", NTOP, NPERM))
  Ysub <- E[keep, , drop = FALSE]; Wsub <- VW[keep, , drop = FALSE]
  NU <- integer(NPERM)
  for (b in seq_len(NPERM)) {
    Xp <- d$X; Xp[, "grpcase"] <- Xp[sample(nrow(Xp)), "grpcase"]
    NU[b] <- nsig(Ysub, Xp, Wsub)$n }
  pv <- (sum(NU >= r$n) + 1) / (NPERM + 1)
  cat(sprintf("   observed %d | null median %.0f, 95th pct %.1f, max %d | permutation p = %.4f\n",
              r$n, median(NU), quantile(NU, .95, names = FALSE), max(NU), pv))

  # Direct test of the precondition: under permuted labels, does a gene's SD
  # predict its p-value? If the filter were truly independent it should not.
  Xp <- d$X; Xp[, "grpcase"] <- Xp[sample(nrow(Xp)), "grpcase"]
  rp <- nsig(E, Xp, VW)
  rho_null <- cor(sdg, -log10(rp$p), method = "spearman")
  rho_obs  <- cor(sdg, -log10(base$p), method = "spearman")
  cat(sprintf("   Spearman(SD, -log10 p): under permuted labels %+.3f (should be ~0)  |  observed %+.3f\n",
              rho_null, rho_obs))

  VALOUT[[length(VALOUT) + 1]] <- data.frame(
    level = "WB", contrast = CN, lineage = NA_character_, n_kept = NTOP,
    observed = r$n, null_median = median(NU), null_q95 = quantile(NU, .95, names = FALSE),
    null_max = max(NU), perm_p = pv, rho_sd_p_null = rho_null, rho_sd_p_obs = rho_obs,
    stringsAsFactors = FALSE)
}

## ================================================================ CELL TYPE
for (CN in c("LI", "BPD")) {
  d <- design_for(CN)
  cat(sprintf("\n%s\n== CELL TYPE %s (n = %d) ==\n%s\n",
              strrep("=", 150), CN, length(d$samples), strrep("=", 150)))

  for (ct in LINS) {
    Y <- A[, ct, d$samples]; W <- mkw(SE[, ct, d$samples])
    sdg <- apply(Y, 1, sd); ord <- order(sdg, decreasing = TRUE)
    base <- nsig(Y, d$X, W, trend = TRUE)

    # Does ranking by SD select for the dominant nuisance axis?
    Yc <- Y - rowMeans(Y)
    sv <- svd(Yc, nu = 1, nv = 1)
    load1 <- abs(sv$u[, 1])                       # each gene's loading on axis 1
    vexp  <- sv$d[1]^2 / sum(sv$d^2)
    keep  <- ord[seq_len(NTOP)]
    enrich <- mean(load1[keep]) / mean(load1)

    r <- nsig(Y[keep, , drop = FALSE], d$X, W[keep, , drop = FALSE], trend = TRUE)
    base_sig <- base$genes[base$q < FDR]
    retained <- sum(base_sig %in% rownames(Y)[keep])

    # permutation null at the same cut
    Ysub <- Y[keep, , drop = FALSE]; Wsub <- W[keep, , drop = FALSE]
    NU <- integer(NPERM)
    for (b in seq_len(NPERM)) {
      Xp <- d$X; Xp[, "grpcase"] <- Xp[sample(nrow(Xp)), "grpcase"]
      NU[b] <- nsig(Ysub, Xp, Wsub, trend = TRUE)$n }
    pv <- (sum(NU >= r$n) + 1) / (NPERM + 1)

    cat(sprintf("  %-5s unfiltered %4d -> top-%d %4d | null med %3.0f q95 %5.1f max %4d | p %.4f | axis1 %.0f%% var, loading enrichment %.2fx | kept %d of %d unfiltered hits\n",
                ct, base$n, NTOP, r$n, median(NU), quantile(NU, .95, names = FALSE),
                max(NU), pv, 100*vexp, enrich, retained, length(base_sig)))

    g <- data.frame(level = "CT", contrast = CN, lineage = ct,
                    gene = r$genes, symbol = unname(ifelse(is.na(sym[r$genes]), r$genes, sym[r$genes])),
                    sd = sdg[r$genes], logFC = r$lfc, P = r$p, BH = r$q,
                    significant = r$q < FDR, stringsAsFactors = FALSE)
    GENEOUT[[paste0("CT_", CN, "_", ct)]] <- g[order(g$P), ]

    SWEEPOUT[[length(SWEEPOUT) + 1]] <- data.frame(
      level = "CT", contrast = CN, lineage = ct, n_kept = NTOP,
      sd_cutoff = min(sdg[keep]), n_sig = r$n,
      base_sig = length(base_sig), base_sig_retained = retained, stringsAsFactors = FALSE)
    VALOUT[[length(VALOUT) + 1]] <- data.frame(
      level = "CT", contrast = CN, lineage = ct, n_kept = NTOP,
      observed = r$n, null_median = median(NU), null_q95 = quantile(NU, .95, names = FALSE),
      null_max = max(NU), perm_p = pv, rho_sd_p_null = NA_real_, rho_sd_p_obs = NA_real_,
      stringsAsFactors = FALSE)
    DIAG[[length(DIAG) + 1]] <- data.frame(
      contrast = CN, lineage = ct, var_explained_axis1 = vexp,
      loading_enrichment_top200 = enrich, stringsAsFactors = FALSE) }
}

SW <- do.call(rbind, SWEEPOUT); VA <- do.call(rbind, VALOUT); DG <- do.call(rbind, DIAG)
GE <- do.call(rbind, GENEOUT)
write.csv(SW, file.path(RES, "sdfilter200_sweep.csv"), row.names = FALSE)
write.csv(VA, file.path(RES, "sdfilter200_validation.csv"), row.names = FALSE)
write.csv(DG, file.path(RES, "sdfilter200_axis_diagnostic.csv"), row.names = FALSE)
write.csv(GE, gzfile(file.path(RES, "sdfilter200_genes.csv.gz")), row.names = FALSE)

cat(sprintf("\n%s\nVALIDATION SUMMARY -- is the 200-gene cut calibrated?\n%s\n",
            strrep("=", 150), strrep("=", 150)))
VA$verdict <- ifelse(VA$null_q95 > 10, "NOT CALIBRATED",
              ifelse(VA$observed > 0 & VA$perm_p < 0.05, "calibrated, beats its null",
              ifelse(VA$observed > 0, "calibrated, does NOT beat its null", "calibrated, nothing found")))
print(VA[, c("level","contrast","lineage","observed","null_median","null_q95","null_max","perm_p","verdict")],
      row.names = FALSE)

cat("\nsignificant genes found at the 200-gene cut:\n")
for (k in names(GENEOUT)) {
  z <- GENEOUT[[k]]; s <- z[z$significant, ]
  cat(sprintf("  %-16s %3d significant%s\n", k, nrow(s),
              if (nrow(s)) paste0(": ", paste(head(s$symbol, 15), collapse = ", "),
                                  if (nrow(s) > 15) ", ..." else "") else "")) }
cat("\nwrote results/sdfilter200_*.csv\n")
