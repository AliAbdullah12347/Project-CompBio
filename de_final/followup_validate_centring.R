#!/usr/bin/env Rscript
# ==============================================================================
# followup_validate_centring.R -- is the recentred result real, or did the
# correction manufacture it?
#
# WHAT HAPPENED. Per-sample median centring moved the cell-type counts from
#   CT_LI   161 -> 978 gene x lineage hits
#   CT_BPD    0 -> 348, of which 345 are in T cells alone
# and the T-cell jump is isolated to the one lineage with the largest offset
# (median logFC -0.066). A correction that turns 0 into 345 in a single lineage
# is either a real unmasking or a false-positive machine, and the two look
# identical from the count alone.
#
# WHY IT COULD BE AN ARTIFACT. Centring subtracts each sample's median across
# genes. That removes a between-sample axis of variation, which shrinks the
# residual variance, which shrinks every standard error, which inflates every
# t-statistic -- for null genes as much as for real ones. It also induces a
# weak negative dependence between genes within a sample, because they now
# share a subtracted quantity.
#
# THE TEST. Permute the group labels and run the identical pipeline. Centring
# is group-blind -- each sample's own median across genes -- so it is legitimate
# to apply it inside a permutation. If the centred analysis returns hundreds of
# genes on scrambled labels, the result is an artifact. If it returns ~0, the
# correction is sound and the genes are real.
#
# The uncentred pipeline is run through the same permutations as a control, so
# the two can be compared on equal terms.
#
# COST CONTROL. Centring uses all 12,368 genes to compute each sample's median,
# because that is what the correction is. Only the TESTING is restricted to a
# fixed random subsample, which changes the precision of the count, not its
# validity. The observed counts are reported on the same subsample so the
# comparison is like-for-like.
# ==============================================================================

suppressPackageStartupMessages({library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "de_final")); set.seed(481); options(width = 175)
RES <- "results"
FDR <- 0.05; NPERM <- 100; NSUB <- 2000

p  <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
bp <- readRDS(file.path(ROOT, "de_analysis/results/bmind_profiles.rds"))
A <- bp$A; SE <- bp$SE; LINS <- bp$lineages
GENES <- dimnames(A)[[1]]

design_for <- function(contrast) {
  s <- p$sel[[contrast]]
  mm <- p$meta[s$keep, , drop = FALSE]; mm$grp <- s$grp
  mm$plate <- droplevels(factor(mm$plate)); mm$sex <- droplevels(factor(mm$sex))
  if ("group" %in% s$cov) mm$group <- droplevels(factor(mm$group))
  cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]])) > 1, logical(1))]
  X <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse = " + "))), data = mm)
  list(X = X, samples = mm$title) }

mkw <- function(S) { W <- 1/(S^2); f <- is.finite(W) & W > 0
  if (any(!f)) W[!f] <- min(W[f]); W/mean(W) }

cat(sprintf("== followup_validate_centring ==\n%d permutations | testing on a fixed %d-gene subsample\n",
            NPERM, NSUB))

OUT <- list()
for (CN in c("LI", "BPD")) {
  d <- design_for(CN)
  sub <- sort(sample(GENES, NSUB))
  # Pre-compute per lineage: the centred and uncentred matrices, and the weights.
  # None of these depend on the group labels, so they are computed once and
  # reused across every permutation -- which is also why permuting is valid.
  PRE <- lapply(LINS, function(ct) {
    Y  <- A[, ct, d$samples]
    Yc <- sweep(Y, 2, apply(Y, 2, median), "-")      # centring uses ALL genes
    list(raw = Y[sub, ], cen = Yc[sub, ], W = mkw(SE[, ct, d$samples])[sub, ]) })
  names(PRE) <- LINS

  count_hits <- function(X, which) vapply(LINS, function(ct) {
    f <- eBayes(lmFit(PRE[[ct]][[which]], X, weights = PRE[[ct]]$W), trend = TRUE)
    tt <- topTable(f, coef = "grpcase", number = Inf, sort.by = "none")
    sum(p.adjust(tt$P.Value, "BH") < FDR) }, integer(1))

  obs_raw <- count_hits(d$X, "raw")
  obs_cen <- count_hits(d$X, "cen")

  t0 <- Sys.time()
  NR <- matrix(NA_integer_, NPERM, length(LINS), dimnames = list(NULL, LINS))
  NC <- NR
  for (b in seq_len(NPERM)) {
    Xp <- d$X; Xp[, "grpcase"] <- Xp[sample(nrow(Xp)), "grpcase"]
    NR[b, ] <- count_hits(Xp, "raw")
    NC[b, ] <- count_hits(Xp, "cen")
    if (b %% 25 == 0) { el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      cat(sprintf("  %s %3d/%d  %.1f min, ETA %.1f\n", CN, b, NPERM, el, el/b*(NPERM-b)))
      flush.console() } }

  cat(sprintf("\n=== %s: genes called on a %d-gene subsample, BH %.2f ===\n", CN, NSUB, FDR))
  R <- do.call(rbind, lapply(LINS, function(ct) data.frame(
    contrast = CN, lineage = ct,
    obs_raw = obs_raw[ct],
    null_raw_median = median(NR[, ct]), null_raw_q95 = quantile(NR[, ct], .95, names = FALSE),
    null_raw_max = max(NR[, ct]),
    p_raw = (sum(NR[, ct] >= obs_raw[ct]) + 1) / (NPERM + 1),
    obs_cen = obs_cen[ct],
    null_cen_median = median(NC[, ct]), null_cen_q95 = quantile(NC[, ct], .95, names = FALSE),
    null_cen_max = max(NC[, ct]),
    p_cen = (sum(NC[, ct] >= obs_cen[ct]) + 1) / (NPERM + 1),
    stringsAsFactors = FALSE)))
  print(R[, c("lineage", "obs_raw", "null_raw_median", "null_raw_q95", "p_raw",
              "obs_cen", "null_cen_median", "null_cen_q95", "null_cen_max", "p_cen")],
        row.names = FALSE)

  cat("\n  verdict per lineage:\n")
  for (i in seq_len(nrow(R))) {
    r <- R[i, ]
    v <- if (r$null_cen_median > 5)
           "CENTRING INFLATES FALSE POSITIVES -- the null itself returns genes"
         else if (r$obs_cen > 0 && r$p_cen < 0.05) "centred result exceeds its own null -- real"
         else if (r$obs_cen > 0) "centred result NOT above its null -- artifact"
         else "nothing either way"
    cat(sprintf("    %-5s obs %4d vs null median %5.1f  ->  %s\n", r$lineage, r$obs_cen,
                r$null_cen_median, v)) }
  OUT[[CN]] <- R
}
V <- do.call(rbind, OUT)
write.csv(V, file.path(RES, "followup_A_validation.csv"), row.names = FALSE)
cat("\nwrote results/followup_A_validation.csv\n")
