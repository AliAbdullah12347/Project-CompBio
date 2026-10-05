#!/usr/bin/env Rscript
# ==============================================================================
# 01_refit_and_mtc.R -- re-fit all four comparisons from scratch and apply the
# full multiple-testing panel. No bands, no verdicts.
#
# WHAT CHANGED FROM de_analysis/ (v1), AND WHY
#
# v1 converted each result into SUPPORTED / INCONCLUSIVE / REFUTED using bands
# at 0.5% and 2% of genes. Those bands are gone. They were defensible as a
# pre-registration device -- they stopped a threshold being chosen after seeing
# the answer -- but they have two costs that outweigh that here:
#
#   1. They discard information. "0.881% of genes" became "INCONCLUSIVE", which
#      is strictly less than the number it replaced.
#   2. A band is a second, hidden multiple-testing decision layered on top of
#      the first, and it has no error guarantee attached to it at all.
#
# So this version reports the quantities and their uncertainty, and stops there.
# Thresholds appear only as columns in a table (how many genes at FDR 0.01, 0.05,
# 0.10, FWER 0.05), never as a verdict.
#
# WHAT IS RE-RUN AND WHAT IS REUSED
#
# Re-run from scratch: every model fit, all four comparisons, both adjustments.
# Reused as INPUT (not as results): de_analysis/data/prep.rds (the filtered,
# normalised matrices and the ILR balances) and de_analysis/results/
# bmind_profiles.rds (3.1 hours of deconvolution; it does not depend on any
# threshold or on the group labels, so re-running it would change nothing).
# No result file from v1 is read.
#
# THE FOUR COMPARISONS
#   WB_LI   whole blood, lithium users vs non-users within bipolar I   n=226
#   WB_BPD  whole blood, bipolar I off lithium vs healthy controls     n=308
#   CT_LI   per cell lineage, same contrast as WB_LI                   n=226
#   CT_BPD  per cell lineage, same contrast as WB_BPD                  n=308
# Each fitted twice: unadjusted, and adjusted for the four ILR cell-composition
# balances. Both are reported because for the cell-type level the adjusted fit
# is the circularity control (bMIND's posterior for a lineage depends on that
# sample's fraction for that lineage, and both contrasts move those fractions).
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
HERE <- file.path(ROOT, "de_v2")
setwd(HERE); set.seed(481)
source("scripts/mtc.R")
RES <- "results"; dir.create(RES, showWarnings = FALSE, recursive = TRUE)

p  <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
bp <- readRDS(file.path(ROOT, "de_analysis/results/bmind_profiles.rds"))
A <- bp$A; SE <- bp$SE; LINS <- bp$lineages
NG <- nrow(p$counts)
cat(sprintf("== 01_refit_and_mtc ==\n%d genes | %d samples | lineages %s\n\n",
            NG, ncol(p$counts), paste(LINS, collapse = ",")))

## ---- design -----------------------------------------------------------------
# Identical covariate set to v1, and identical across levels, so that whole
# blood and cell type differ only in what is being measured.
design_for <- function(contrast, adjust) {
  s <- p$sel[[contrast]]
  mm <- p$meta[s$keep, , drop = FALSE]; mm$grp <- s$grp
  mm$plate <- droplevels(factor(mm$plate)); mm$sex <- droplevels(factor(mm$sex))
  if ("group" %in% s$cov) mm$group <- droplevels(factor(mm$group))
  cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]])) > 1, logical(1))]
  X <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse = " + "))), data = mm)
  if (adjust == "ilr") X <- cbind(X, p$ILR[mm$title, , drop = FALSE])
  ne <- nonEstimable(X); if (!is.null(ne)) X <- X[, setdiff(colnames(X), ne), drop = FALSE]
  stopifnot(qr(X)$rank == ncol(X))        # refuse to fit a rank-deficient design
  list(X = X, samples = mm$title, dropped = setdiff(s$cov, cv))
}

# 1/SE^2 precision weights from bMIND's posterior, with non-finite values
# repaired rather than silently dropped (there are none in this dataset, but the
# guard is kept so the script cannot fail quietly on a re-run).
make_weights <- function(S) {
  W <- 1 / (S^2); fin <- is.finite(W) & W > 0
  if (any(!fin)) { hi <- max(W[fin]); lo <- min(W[fin])
    W[!is.finite(W) & !is.na(S) & S == 0] <- hi; W[!is.finite(W) | W <= 0] <- lo }
  W / mean(W) }

## ---- fits -------------------------------------------------------------------
fit_wb <- function(contrast, adjust) {
  d <- design_for(contrast, adjust)
  dge <- calcNormFactors(DGEList(p$counts[, d$samples, drop = FALSE]), method = "TMM")
  f <- eBayes(lmFit(voom(dge, d$X), d$X))
  tt <- topTable(f, coef = "grpcase", number = Inf, sort.by = "none")
  data.frame(gene = rownames(tt), lineage = NA_character_, logFC = tt$logFC,
             AveExpr = tt$AveExpr, t = tt$t, P.Value = tt$P.Value,
             n = length(d$samples), stringsAsFactors = FALSE)
}

fit_ct <- function(contrast, adjust) {
  d <- design_for(contrast, adjust)
  do.call(rbind, lapply(LINS, function(ct) {
    W <- make_weights(SE[, ct, d$samples])
    f <- eBayes(lmFit(A[, ct, d$samples], d$X, weights = W), trend = TRUE)
    tt <- topTable(f, coef = "grpcase", number = Inf, sort.by = "none")
    data.frame(gene = rownames(tt), lineage = ct, logFC = tt$logFC,
               AveExpr = tt$AveExpr, t = tt$t, P.Value = tt$P.Value,
               n = length(d$samples), stringsAsFactors = FALSE) }))
}

JOBS <- expand.grid(contrast = c("LI", "BPD"), level = c("WB", "CT"),
                    adjust = c("raw", "ilr"), stringsAsFactors = FALSE)
JOBS$key <- sprintf("%s_%s_%s", JOBS$level, JOBS$contrast, JOBS$adjust)

FITS <- list()
for (i in seq_len(nrow(JOBS))) {
  j <- JOBS[i, ]
  t0 <- Sys.time()
  FITS[[j$key]] <- if (j$level == "WB") fit_wb(j$contrast, j$adjust) else fit_ct(j$contrast, j$adjust)
  cat(sprintf("  %-14s %6d tests  n=%3d  [%.1f min]\n", j$key, nrow(FITS[[j$key]]),
              FITS[[j$key]]$n[1], as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  flush.console()
}
saveRDS(FITS, file.path(RES, "fits.rds"))

## ---- multiple testing -------------------------------------------------------
# Every method is applied to EVERY analysis, so no result depends on which
# correction a reader prefers. Hommel is skipped above 20,000 tests: its
# algorithm is quadratic and would take hours on the 61,840-test cell-type
# families for a negligible gain over Holm.
cat("\napplying the correction panel\n")
HOMMEL_MAX <- 20000
ALL <- list(); SUMM <- list()
for (k in names(FITS)) {
  f <- FITS[[k]]; pv <- f$P.Value; m <- length(pv)
  a <- data.frame(
    bonferroni = p.adjust(pv, "bonferroni"),
    holm       = p.adjust(pv, "holm"),
    BH         = p.adjust(pv, "BH"),
    BY         = p.adjust(pv, "BY"),
    stringsAsFactors = FALSE)
  a$hommel <- if (m <= HOMMEL_MAX) p.adjust(pv, "hommel") else NA_real_
  s_sp <- pi0_spline(pv); s_bs <- pi0_bootstrap(pv)
  a$q_storey_spline <- qvalue_storey(pv, pi0 = s_sp)
  a$q_storey_boot   <- qvalue_storey(pv, pi0 = s_bs)
  out <- cbind(f, a)
  ALL[[k]] <- out
  cnt <- function(col, thr) sum(out[[col]] < thr, na.rm = TRUE)
  SUMM[[k]] <- data.frame(
    analysis = k, level = substr(k, 1, 2),
    contrast = f$gene[0][1],   # placeholder, overwritten below
    n_tests = m, n_subjects = f$n[1],
    pi0_spline = round(s_sp, 4), pi0_boot = round(s_bs, 4),
    min_p = signif(min(pv), 3), max_abs_logFC = round(max(abs(f$logFC)), 3),
    BH_01 = cnt("BH", 0.01), BH_05 = cnt("BH", 0.05), BH_10 = cnt("BH", 0.10),
    BY_05 = cnt("BY", 0.05),
    storey_spline_05 = cnt("q_storey_spline", 0.05),
    storey_boot_05 = cnt("q_storey_boot", 0.05),
    holm_05 = cnt("holm", 0.05), bonferroni_05 = cnt("bonferroni", 0.05),
    hommel_05 = if (m <= HOMMEL_MAX) cnt("hommel", 0.05) else NA_integer_,
    stringsAsFactors = FALSE)
  cat(sprintf("  %-14s pi0 %.3f/%.3f | BH05 %5d | BY05 %5d | Storey05 %5d | Holm05 %4d\n",
              k, s_sp, s_bs, cnt("BH", 0.05), cnt("BY", 0.05),
              cnt("q_storey_spline", 0.05), cnt("holm", 0.05)))
  flush.console()
}
S <- do.call(rbind, SUMM)
S$contrast <- ifelse(grepl("_LI_", S$analysis), "LI", "BPD")
S$adjust <- ifelse(grepl("_ilr$", S$analysis), "ilr", "raw")
S <- S[, c("analysis", "level", "contrast", "adjust", setdiff(names(S), c("analysis", "level", "contrast", "adjust")))]

saveRDS(ALL, file.path(RES, "adjusted.rds"))
write.csv(S, file.path(RES, "summary_by_analysis.csv"), row.names = FALSE)
for (k in names(ALL))
  write.csv(ALL[[k]][order(ALL[[k]]$P.Value), ],
            gzfile(file.path(RES, sprintf("tests_%s.csv.gz", k))), row.names = FALSE)

## ---- independent filtering ---------------------------------------------------
# Bourgon et al. (2010): drop the genes least likely to be detectable using a
# statistic independent of the p-value under the null. Mean expression qualifies
# because it is computed across all samples ignoring group. This is reported as
# a diagnostic, NOT applied to the headline numbers -- our gene filter already
# removed low-expression genes, so the question is whether more filtering would
# still have bought power. If the best quantile is 0, it would not.
cat("\nindependent filtering (is further expression filtering worth it?)\n")
IF <- do.call(rbind, lapply(names(ALL), function(k) {
  f <- ALL[[k]]
  r <- independent_filter(f$P.Value, f$AveExpr, alpha = 0.05)
  r$analysis <- k; r }))
best <- IF[IF$is_best, ]
best <- best[!duplicated(best$analysis), ]
for (i in seq_len(nrow(best)))
  cat(sprintf("  %-14s best filter quantile %.2f -> %d rejections (unfiltered %d)\n",
              best$analysis[i], best$filter_quantile[i], best$n_rejected[i],
              IF$n_rejected[IF$analysis == best$analysis[i] & IF$filter_quantile == 0]))
write.csv(IF, file.path(RES, "independent_filtering.csv"), row.names = FALSE)

cat("\nwrote results/summary_by_analysis.csv, tests_*.csv.gz, fits.rds, adjusted.rds\n")
