#!/usr/bin/env Rscript
# ==============================================================================
# method_sweep_null.R -- calibration for the method sweep.
#
# method_sweep.R answers "do the methods give different counts?". A different
# count is only interesting if we can say which count to believe, and a DEG
# total on its own cannot say that: a method can lead the table because it is
# more powerful or because it is anti-conservative, and those look identical.
# Two checks separate them.
#
#   P  PERMUTATION NULL. The lithium column is permuted and every covariate
#      stays attached to its own sample, so the null is "lithium carries no
#      signal", not "the samples are exchangeable". A calibrated method should
#      return ~0 DEG at FDR 0.05 and a flat p-value distribution. All methods
#      see the SAME permutations, so the comparison is paired.
#
#   S  SPLIT-HALF REPRODUCIBILITY. The 226 BP1 samples are split in two,
#      stratified on lithium, and each method is fitted to each half. A method
#      that calls more genes because it is more powerful will also agree with
#      itself better across halves; a method that calls more genes because it
#      is anti-conservative will not. The gene filter re-runs inside each half
#      (it is a cross-sample operation), so the halves are compared on their
#      intersected gene lists.
#
# Usage:  Rscript scripts/method_sweep_null.R <stage: P|S|PS> <B> [methods,csv]
# Writes only inside experimentation/.
# ==============================================================================

source("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation/scripts/common.R")

OUT <- file.path(EXP_HOME, "runs", "method_sweep")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
log <- function(...) { cat(sprintf(...), file = stderr()); flush(stderr()) }

args   <- commandArgs(trailingOnly = TRUE)
stage  <- if (length(args) >= 1) args[1] else "PS"
B      <- if (length(args) >= 2) as.integer(args[2]) else 10L
MTHS   <- if (length(args) >= 3) strsplit(args[3], ",")[[1]] else
          c("voom", "voomWQW", "trend", "QLF", "LRT")

SEED <- 20261218                    # the project deadline; one seed, threaded
p <- load_prep()

lambda_gc <- function(pv) median(qchisq(1 - pv, 1), na.rm = TRUE) / qchisq(0.5, 1)
signed_stat <- function(tt, method) {
  if (method %in% c("voom", "voomWQW", "trend")) tt$t
  else if (method == "QLF") sign(tt$logFC) * sqrt(pmax(tt$F,  0))
  else                      sign(tt$logFC) * sqrt(pmax(tt$LR, 0))
}

## ===========================================================================
## P. permutation null
## ===========================================================================
if (grepl("P", stage)) {
  d0 <- build_design(p, "lithium")
  n  <- d0$n
  set.seed(SEED)
  PERMS <- lapply(seq_len(B), function(b) sample.int(n))
  # A permutation that lands close to the observed assignment is not a null
  # draw in any useful sense; record the overlap so a degenerate draw is
  # visible rather than silently diluting the null.
  obs <- d0$X[, d0$coef]
  perm_overlap <- vapply(PERMS, function(ix) mean(obs[ix] == obs), numeric(1))
  log("[P] %d permutations, exposure agreement with observed: %.3f-%.3f (chance %.3f)\n",
      B, min(perm_overlap), max(perm_overlap),
      mean(obs)^2 + (1 - mean(obs))^2)

  rows <- list(); k <- 0
  for (mth in MTHS) for (b in seq_len(B)) {
    t0 <- proc.time()
    f  <- de_fit(p, "lithium", method = mth, permute_exposure = PERMS[[b]])
    el <- proc.time() - t0
    tt <- f$tt
    k <- k + 1
    rows[[k]] <- data.frame(
      method = mth, perm = b, n_genes = f$n_genes,
      deg_fdr05 = n_deg(tt, 0.05), deg_fdr10 = n_deg(tt, 0.10),
      deg_fdr20 = n_deg(tt, 0.20),
      frac_p_lt_05 = round(mean(tt$P.Value < 0.05, na.rm = TRUE), 5),
      frac_p_lt_01 = round(mean(tt$P.Value < 0.01, na.rm = TRUE), 5),
      frac_p_lt_001 = round(mean(tt$P.Value < 0.001, na.rm = TRUE), 5),
      pi0 = round(pi0_storey(tt$P.Value), 4),
      lambda_gc = round(lambda_gc(tt$P.Value), 4),
      min_p = signif(min(tt$P.Value, na.rm = TRUE), 3),
      ks_D = round(suppressWarnings(ks.test(tt$P.Value, "punif"))$statistic, 5),
      cpu = round(el[["user.self"]] + el[["sys.self"]], 1),
      stringsAsFactors = FALSE)
    log("[P] %-8s perm %2d  %4d DEG@05  p<.05 %.4f  lambda %.3f  pi0 %.3f\n",
        mth, b, rows[[k]]$deg_fdr05, rows[[k]]$frac_p_lt_05,
        rows[[k]]$lambda_gc, rows[[k]]$pi0)
  }
  # One file per method, so the stage can be run in pieces (voomWQW costs an
  # order of magnitude more than the rest and is run separately) without a
  # later invocation destroying an earlier one's results. The summary is then
  # rebuilt from every per-method file present.
  perm <- do.call(rbind, rows)
  for (mth in unique(perm$method))
    write.csv(perm[perm$method == mth, ],
              file.path(OUT, sprintf("permutation_null_raw_%s.csv", mth)),
              row.names = FALSE)
  have <- list.files(OUT, "^permutation_null_raw_.*\\.csv$", full.names = TRUE)
  perm <- do.call(rbind, lapply(have, read.csv, stringsAsFactors = FALSE))
  write.csv(perm, file.path(OUT, "permutation_null_raw.csv"), row.names = FALSE)

  agg <- do.call(rbind, lapply(split(perm, perm$method), function(x) data.frame(
    method = x$method[1], B = nrow(x),
    deg05_mean = round(mean(x$deg_fdr05), 2), deg05_max = max(x$deg_fdr05),
    deg05_n_nonzero = sum(x$deg_fdr05 > 0),
    deg10_mean = round(mean(x$deg_fdr10), 2),
    # Under a correct null this is 0.05 by construction. Anything above it is
    # the method's excess type-I error at the nominal level.
    frac_p05_mean = round(mean(x$frac_p_lt_05), 5),
    frac_p01_mean = round(mean(x$frac_p_lt_01), 5),
    frac_p001_mean = round(mean(x$frac_p_lt_001), 5),
    lambda_mean = round(mean(x$lambda_gc), 4),
    pi0_mean = round(mean(x$pi0), 4),
    ksD_mean = round(mean(x$ks_D), 5),
    cpu_mean = round(mean(x$cpu), 1),
    stringsAsFactors = FALSE)))
  agg$excess_type1_at_05 <- round(agg$frac_p05_mean / 0.05, 3)
  agg$excess_type1_at_01 <- round(agg$frac_p01_mean / 0.01, 3)
  agg <- agg[order(-agg$frac_p05_mean), ]
  log("\n== P: permutation null, averaged over %d permutations ==\n", B)
  print(agg, row.names = FALSE)
  write.csv(agg, file.path(OUT, "permutation_null_summary.csv"), row.names = FALSE)
}

## ===========================================================================
## S. split-half reproducibility
## ===========================================================================
if (grepl("S", stage)) {
  bp1 <- which(p$meta$dx == "BP1")
  lit <- p$meta$lithium[bp1]
  set.seed(SEED + 1)
  rows <- list(); k <- 0
  SPLITS <- lapply(seq_len(B), function(r) {
    # Stratified on lithium so both halves carry the same exposure ratio;
    # an unstratified split can hand one half far fewer non-users and the
    # reproducibility drop would then be a sample-size artefact.
    h1 <- unlist(lapply(split(seq_along(bp1), lit), function(ix)
      sample(ix, floor(length(ix) / 2))))
    a <- rep(FALSE, nrow(p$meta)); a[bp1[h1]] <- TRUE
    b <- rep(FALSE, nrow(p$meta)); b[bp1[setdiff(seq_along(bp1), h1)]] <- TRUE
    list(a = a, b = b)
  })

  for (mth in MTHS) for (r in seq_len(B)) {
    s  <- SPLITS[[r]]
    fa <- de_fit(p, "lithium", method = mth, samples = s$a)
    fb <- de_fit(p, "lithium", method = mth, samples = s$b)
    g  <- intersect(fa$tt$gene, fb$tt$gene)
    ia <- match(g, fa$tt$gene); ib <- match(g, fb$tt$gene)
    sa <- signed_stat(fa$tt, mth)[ia]; sb <- signed_stat(fb$tt, mth)[ib]
    pa <- fa$tt$P.Value[ia];           pb <- fb$tt$P.Value[ib]
    da <- fa$tt$adj.P.Val[ia] < 0.05;  db <- fb$tt$adj.P.Val[ib] < 0.05
    topA <- g[order(pa)][seq_len(500)]; topB <- g[order(pb)][seq_len(500)]
    k <- k + 1
    rows[[k]] <- data.frame(
      method = mth, rep = r, n_a = fa$n_samples, n_b = fb$n_samples,
      genes_a = fa$n_genes, genes_b = fb$n_genes, genes_shared = length(g),
      deg_a = sum(da, na.rm = TRUE), deg_b = sum(db, na.rm = TRUE),
      deg_both = sum(da & db, na.rm = TRUE),
      # Of the genes one half calls, how many does the other half also call?
      # This is the quantity that separates power from anti-conservatism.
      repro = round(if (sum(da, na.rm = TRUE) + sum(db, na.rm = TRUE) == 0) NA_real_
                    else 2 * sum(da & db, na.rm = TRUE) /
                         (sum(da, na.rm = TRUE) + sum(db, na.rm = TRUE)), 4),
      spearman = round(cor(sa, sb, method = "spearman", use = "complete.obs"), 4),
      top500_overlap = round(length(intersect(topA, topB)) / 500, 4),
      # Direction agreement among the genes the first half calls significant.
      sign_agree = round(if (sum(da, na.rm = TRUE) == 0) NA_real_
                         else mean(sign(sa[da]) == sign(sb[da]), na.rm = TRUE), 4),
      stringsAsFactors = FALSE)
    log("[S] %-8s rep %2d  deg %4d/%4d both %4d  repro %.3f  rho %.3f  top500 %.3f\n",
        mth, r, rows[[k]]$deg_a, rows[[k]]$deg_b, rows[[k]]$deg_both,
        rows[[k]]$repro, rows[[k]]$spearman, rows[[k]]$top500_overlap)
  }
  sh <- do.call(rbind, rows)
  for (mth in unique(sh$method))                     # same per-method policy
    write.csv(sh[sh$method == mth, ],
              file.path(OUT, sprintf("splithalf_raw_%s.csv", mth)), row.names = FALSE)
  have <- list.files(OUT, "^splithalf_raw_.*\\.csv$", full.names = TRUE)
  sh <- do.call(rbind, lapply(have, read.csv, stringsAsFactors = FALSE))
  write.csv(sh, file.path(OUT, "splithalf_raw.csv"), row.names = FALSE)

  shagg <- do.call(rbind, lapply(split(sh, sh$method), function(x) data.frame(
    method = x$method[1], reps = nrow(x),
    deg_half_mean = round(mean(c(x$deg_a, x$deg_b)), 1),
    deg_both_mean = round(mean(x$deg_both), 1),
    repro_mean = round(mean(x$repro, na.rm = TRUE), 4),
    repro_sd = round(sd(x$repro, na.rm = TRUE), 4),
    spearman_mean = round(mean(x$spearman), 4),
    top500_mean = round(mean(x$top500_overlap), 4),
    sign_agree_mean = round(mean(x$sign_agree, na.rm = TRUE), 4),
    stringsAsFactors = FALSE)))
  shagg <- shagg[order(-shagg$repro_mean), ]
  log("\n== S: split-half reproducibility over %d splits ==\n", B)
  print(shagg, row.names = FALSE)
  write.csv(shagg, file.path(OUT, "splithalf_summary.csv"), row.names = FALSE)
}
log("\n[null] done\n")
