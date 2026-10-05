#!/usr/bin/env Rscript
# ==============================================================================
# 10_celltype_null.R -- permutation null for the cell-type level. Two jobs.
#
# ------------------------------------------------------------------------------
# JOB 1 (the reason this script got bigger): IS pi0 EVEN CALIBRATED HERE?
#
# 07 produced per-lineage pi0 values that do not make sense at face value:
#
#     BPD, B cells : pi0 = 0.4395  but ZERO genes at FDR 0.05 and min p = 8e-5
#     LI,  B cells : pi0 = 0.4166  but ZERO genes at FDR 0.05
#     BPD, T cells : pi0 = 0.5268  but ZERO genes at FDR 0.05
#
# pi0 = 0.44 asserts that 56% of genes carry signal. A gene set with 56% true
# effects and n=308 should produce far more than zero FDR-significant genes.
# The two statements are not compatible, and the most likely explanation is
# not biology but MISCALIBRATION: bMIND shrinks every sample toward a per-gene
# prior, which makes the posterior means more similar to each other than the
# underlying data are. limma then estimates a residual variance that is too
# small, p-values drift low, and pi0 drops -- without any gene becoming
# individually significant.
#
# Taking those pi0 values at face value would mean reporting "widespread
# hidden signal in B cells", which would be an artifact of the estimator.
#
# The test: permute the labels and recompute pi0. If the NULL also produces
# pi0 around 0.44, then 0.44 means nothing. This is the single most important
# check on the cell-type arm.
#
# ------------------------------------------------------------------------------
# JOB 2: the pre-specified MECHANISM test for CT_BPD.
#
# config.yaml, CT_BPD: mechanism_test:
#     what:          sign discordance across lineages among genes null in whole blood
#     supported_if:  >=10% of bulk-null genes show opposite-sign lineage effects
#                    at nominal p<0.05
#     compared_against: permutation null
#
# CT_BPD is otherwise decided on a COUNT comparison, and a count cannot
# demonstrate a mechanism. Masking makes a sharper prediction: a gene hidden in
# bulk BECAUSE two cell types move it in opposite directions should show
# opposite-sign effects in two lineages. That is a claim about SIGNS.
#
# ------------------------------------------------------------------------------
# WHY PERMUTING LABELS HERE IS A COMPLETE NULL
#
# 06 ran bMIND WITHOUT the phenotype, so the cell-type profiles do not depend
# on the group labels in any way. Permuting labels afterwards and re-running
# only the limma step is therefore a full and valid null -- the deconvolution
# would be bit-identical and does not need repeating. (Had bMIND been run with
# y, as bmind_de does, this shortcut would be invalid and each permutation
# would cost a 3-hour deconvolution.)
#
# The null preserves every gene-gene and lineage-lineage correlation. That
# matters here more than usual: the five lineage estimates for one gene come
# from ONE bulk measurement on ONE subject, so they are strongly dependent by
# construction and no independence-based reference would be honest.
#
# COST. 200 permutations x 5 lineages x 2 contrasts. Run on a fixed random
# subsample of genes, drawn once with the config seed and reused for every
# permutation and for the observed statistic, so the comparison is like-for-
# like. Subsampling changes the precision of the estimate, not its validity.
# ==============================================================================

suppressPackageStartupMessages({library(limma)})
HERE <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/de_analysis"
setwd(HERE); set.seed(481)
RES <- "results"
p  <- readRDS("data/prep.rds")
bp <- readRDS(file.path(RES, "bmind_profiles.rds"))
A <- bp$A; SE <- bp$SE; LINS <- bp$lineages

# The 200 permutations are the pre-specified number and are NOT reduced.
# The gene subsample is what absorbs the compute budget instead, because
# subsampling changes the PRECISION of each statistic and not its validity,
# whereas cutting permutations would coarsen the empirical p-value itself.
# At 800 genes the standard error of pi0 is about 0.035 and of the discordance
# rate about 0.8 percentage points -- small against the differences being
# judged here.
NPERM  <- 200
NSUB   <- 800
NOMP   <- 0.05     # nominal p for "this lineage moved"
THRESH <- 10       # pre-specified discordance threshold, percent

pi0_storey <- function(pv, l = 0.5) min(1, mean(pv > l, na.rm = TRUE) / (1 - l))

make_weights <- function(S) {
  W <- 1 / (S^2); fin <- is.finite(W) & W > 0
  if (any(!fin)) W[!fin] <- min(W[fin])
  W / mean(W) }

design_for <- function(contrast) {
  s <- p$sel[[contrast]]
  mm <- p$meta[s$keep, , drop = FALSE]; mm$grp <- s$grp
  mm$plate <- droplevels(factor(mm$plate)); mm$sex <- droplevels(factor(mm$sex))
  if ("group" %in% s$cov) mm$group <- droplevels(factor(mm$group))
  cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]])) > 1, logical(1))]
  X <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse = " + "))), data = mm)
  ne <- nonEstimable(X); if (!is.null(ne)) X <- X[, setdiff(colnames(X), ne), drop = FALSE]
  list(X = X, samples = mm$title) }

# One pass over the five lineages: per-lineage pi0 and DEG count, plus the
# cross-lineage sign-discordance rate.
stats_for <- function(X, genes, samples, Wl) {
  P <- S <- matrix(NA_real_, length(genes), length(LINS), dimnames = list(genes, LINS))
  for (ct in LINS) {
    fit <- eBayes(lmFit(A[genes, ct, samples], X, weights = Wl[[ct]][genes, ]), trend = TRUE)
    tt <- topTable(fit, coef = "grpcase", number = Inf, sort.by = "none")
    P[rownames(tt), ct] <- tt$P.Value; S[rownames(tt), ct] <- sign(tt$logFC)
  }
  hit <- P < NOMP
  up <- rowSums(hit & S > 0, na.rm = TRUE); dn <- rowSums(hit & S < 0, na.rm = TRUE)
  list(pi0 = apply(P, 2, pi0_storey),
       deg = apply(P, 2, function(v) sum(p.adjust(v, "BH") < 0.05, na.rm = TRUE)),
       disc = 100 * mean(up > 0 & dn > 0)) }

ALL_PI0 <- list(); ALL_SUM <- list(); MECH <- NULL
for (cn in c("BPD", "LI")) {
  d <- design_for(cn); SAMP <- d$samples
  # Weights depend only on SE, never on labels -> identical in every permutation.
  Wl <- lapply(LINS, function(ct) make_weights(SE[, ct, SAMP])); names(Wl) <- LINS

  wbtt <- read.csv(file.path(RES, sprintf("WB_%s_raw_toptable.csv", cn)), stringsAsFactors = FALSE)
  bulk_null <- intersect(wbtt$gene[wbtt$adj.P.Val >= 0.05], dimnames(A)[[1]])
  sub <- sort(sample(bulk_null, min(NSUB, length(bulk_null))))

  cat(sprintf("\n================ %s  (n=%d) ================\n", cn, length(SAMP)))
  cat(sprintf("bulk-null genes %d of %d | subsample %d | %d permutations\n",
              length(bulk_null), nrow(wbtt), length(sub), NPERM))

  t0 <- Sys.time(); o <- stats_for(d$X, sub, SAMP, Wl)
  per1 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("observed: pi0 %s | disc %.2f%%  [%.0f s per fit set, ETA %.0f min]\n",
              paste(sprintf("%s=%.3f", LINS, o$pi0), collapse = " "), o$disc,
              per1, per1 * NPERM / 60)); flush.console()

  NP <- matrix(NA_real_, NPERM, length(LINS), dimnames = list(NULL, LINS))
  ND <- matrix(NA_real_, NPERM, length(LINS), dimnames = list(NULL, LINS))
  NDI <- numeric(NPERM)
  for (b in seq_len(NPERM)) {
    Xp <- d$X; Xp[, "grpcase"] <- Xp[sample(nrow(Xp)), "grpcase"]
    r <- stats_for(Xp, sub, SAMP, Wl)
    NP[b, ] <- r$pi0; ND[b, ] <- r$deg; NDI[b] <- r$disc
    if (b %% 25 == 0) {
      el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      cat(sprintf("  %3d/%d  elapsed %5.1f min  ETA %5.1f min\n", b, NPERM, el, el / b * (NPERM - b)))
      flush.console() }
  }

  cat(sprintf("\n--- %s: is pi0 calibrated? (null should sit near 1.0) ---\n", cn))
  S1 <- do.call(rbind, lapply(LINS, function(ct) data.frame(
    contrast = cn, lineage = ct, mean_frac = round(mean(p$lineage[, ct]), 4),
    obs_pi0 = round(o$pi0[ct], 4),
    null_pi0_median = round(median(NP[, ct]), 4),
    null_pi0_q05 = round(quantile(NP[, ct], 0.05, names = FALSE), 4),
    null_pi0_q95 = round(quantile(NP[, ct], 0.95, names = FALSE), 4),
    p_pi0 = round((sum(NP[, ct] <= o$pi0[ct]) + 1) / (NPERM + 1), 4),
    obs_deg = o$deg[ct], null_deg_median = median(ND[, ct]),
    null_deg_q95 = quantile(ND[, ct], 0.95, names = FALSE),
    p_deg = round((sum(ND[, ct] >= o$deg[ct]) + 1) / (NPERM + 1), 4),
    stringsAsFactors = FALSE)))
  print(S1[, c("lineage", "mean_frac", "obs_pi0", "null_pi0_median", "null_pi0_q05",
               "p_pi0", "obs_deg", "null_deg_q95", "p_deg")], row.names = FALSE)
  ALL_SUM[[cn]] <- S1; ALL_PI0[[cn]] <- NP

  if (cn == "BPD") {
    p_emp <- (sum(NDI >= o$disc) + 1) / (NPERM + 1)
    qn <- quantile(NDI, c(0.025, 0.5, 0.975), names = FALSE)
    cat("\n--- BPD: pre-specified MECHANISM test (sign discordance) ---\n")
    cat(sprintf("  observed      : %.3f%% of bulk-null genes\n", o$disc))
    cat(sprintf("  null          : median %.3f%%, 95%% range [%.3f%%, %.3f%%]\n", qn[2], qn[1], qn[3]))
    cat(sprintf("  empirical p   : %.4f | pre-specified threshold >=%d%%\n", p_emp, THRESH))
    v <- {
      if (o$disc >= THRESH && p_emp < 0.05) "SUPPORTED"
      else if (o$disc >= THRESH) "threshold met but not above the null -- NOT SUPPORTED"
      else if (p_emp < 0.05) "above the null but below the threshold -- NOT SUPPORTED"
      else "NOT SUPPORTED" }
    cat(sprintf("  VERDICT       : %s\n", v))
    MECH <- data.frame(
      statistic = "pct bulk-null genes with opposite-sign lineage effects, nominal p<0.05",
      n_bulk_null = length(bulk_null), n_subsample = length(sub), n_perm = NPERM,
      obs_pct = round(o$disc, 3), null_median = round(qn[2], 3),
      null_lo = round(qn[1], 3), null_hi = round(qn[3], 3),
      p_emp = round(p_emp, 4), threshold_pct = THRESH, verdict = v, stringsAsFactors = FALSE)
    write.csv(data.frame(perm = seq_len(NPERM), pct_discordant = round(NDI, 4)),
              file.path(RES, "celltype_mechanism_null.csv"), row.names = FALSE)
    write.csv(MECH, file.path(RES, "celltype_mechanism_test.csv"), row.names = FALSE)
  }
}

SUM <- do.call(rbind, ALL_SUM)
write.csv(SUM, file.path(RES, "celltype_permutation_calibration.csv"), row.names = FALSE)
saveRDS(ALL_PI0, file.path(RES, "celltype_permnull_pi0.rds"))

cat("\n\n================== CALIBRATION VERDICT ==================\n")
SUM$pi0_trustworthy <- SUM$null_pi0_median > 0.95
for (i in seq_len(nrow(SUM))) {
  r <- SUM[i, ]
  cat(sprintf("  %-4s %-5s obs pi0 %.3f | null median %.3f -> %s\n", r$contrast, r$lineage,
              r$obs_pi0, r$null_pi0_median,
              if (r$pi0_trustworthy) "calibrated, pi0 interpretable"
              else "NOT CALIBRATED -- pi0 uninterpretable at this level"))
}
write.csv(SUM, file.path(RES, "celltype_permutation_calibration.csv"), row.names = FALSE)
cat("\nwrote results/celltype_permutation_calibration.csv and celltype_mechanism_*.csv\n")
