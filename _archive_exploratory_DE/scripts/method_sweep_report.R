#!/usr/bin/env Rscript
# ==============================================================================
# method_sweep_report.R -- assemble the one table that answers the question.
#
# Everything the sweep produced is per-stage. This merges the pieces into a
# single row per configuration: how many genes it called, how far it sits from
# the other estimators, whether it is calibrated under a permutation null, and
# whether it agrees with itself across a split half. Those four columns
# together are what separates "more powerful" from "anti-conservative", which a
# DEG count on its own cannot do.
#
# Reads only files this experiment wrote. Writes runs/method_sweep/summary.csv.
# ==============================================================================

source("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation/scripts/common.R")
OUT <- file.path(EXP_HOME, "runs", "method_sweep")
log <- function(...) { cat(sprintf(...), file = stderr()); flush(stderr()) }

rd <- function(f) {
  p <- file.path(OUT, f)
  if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else NULL
}
pm <- rd("per_method.csv"); stopifnot(!is.null(pm))
od <- rd("oddness.csv"); pn <- rd("permutation_null_summary.csv")
sh <- rd("splithalf_summary.csv")

## --------------------------------------------------------------- paired null
# The five estimators were given different numbers of permutations, because
# voomWQW costs ~50x voom per fit and QLF/LRT ~15x. The draws are nested: the
# script seeds once and takes permutations in order, so method A's permutation
# k IS method B's permutation k for every k both of them ran. The comparison
# across methods must therefore be made on the permutations they SHARE, paired,
# which also removes the large between-permutation variance (voom's own null
# fraction ranges 0.011-0.109 across draws) that would otherwise swamp any
# difference between methods.
pr <- rd("permutation_null_raw.csv")
if (!is.null(pr)) {
  common <- Reduce(intersect, split(pr$perm, pr$method))
  sub <- pr[pr$perm %in% common, ]
  W <- reshape(sub[, c("perm", "method", "frac_p_lt_05")], idvar = "perm",
               timevar = "method", direction = "wide")
  names(W) <- sub("frac_p_lt_05\\.", "", names(W))
  mths <- setdiff(names(W), "perm")
  base <- if ("voom" %in% mths) "voom" else mths[1]
  paired <- data.frame(
    method = mths,
    n_shared_perms = length(common),
    null_frac_p05_shared = round(sapply(mths, function(m) mean(W[[m]])), 5),
    diff_vs_voom = round(sapply(mths, function(m) mean(W[[m]] - W[[base]])), 5),
    # Paired across the shared permutations. With 3-5 pairs this is a direction
    # and a magnitude, not a significance test, and is labelled as such.
    sd_of_paired_diff = round(sapply(mths, function(m)
      if (length(common) > 1) sd(W[[m]] - W[[base]]) else NA_real_), 5),
    n_perms_where_higher_than_voom = sapply(mths, function(m) sum(W[[m]] > W[[base]])),
    stringsAsFactors = FALSE)
  paired$ratio_vs_voom <- round(paired$null_frac_p05_shared /
                                paired$null_frac_p05_shared[paired$method == base], 3)
  paired <- paired[order(-paired$null_frac_p05_shared), ]
  log("\n== paired permutation null, on the %d permutation(s) all methods share ==\n",
      length(common))
  print(paired, row.names = FALSE)
  write.csv(paired, file.path(OUT, "permutation_null_paired.csv"), row.names = FALSE)
  write.csv(W, file.path(OUT, "permutation_null_paired_wide.csv"), row.names = FALSE)
}

## --------------------------------------------- scaling with an intercept too
# scaling_extrapolated.csv forces the fit through the origin, which is right
# for voom/trend (cost is per gene) but wrong for QLF/LRT, where estimateDisp
# has a large fixed cost. Refit with an intercept so the two readings sit side
# by side and neither is quietly presented as the answer.
sr <- rd("scaling_raw.csv")
if (!is.null(sr)) {
  sc <- do.call(rbind, lapply(split(sr, sr$method), function(x) {
    f0 <- lm(cpu ~ 0 + n_genes, data = x); f1 <- lm(cpu ~ n_genes, data = x)
    data.frame(method = x$method[1],
               fixed_cost_s = round(coef(f1)[[1]], 2),
               per_1000_genes_s = round(coef(f1)[[2]] * 1000, 3),
               pred_12173_with_intercept = round(predict(f1, data.frame(n_genes = 12173)), 1),
               pred_12173_through_origin = round(coef(f0)[[1]] * 12173, 1),
               r2 = round(summary(f1)$r.squared, 4), stringsAsFactors = FALSE)
  }))
  sc <- sc[order(-sc$pred_12173_with_intercept), ]
  log("\n== estimator cost, refitted with an intercept ==\n"); print(sc, row.names = FALSE)
  write.csv(sc, file.path(OUT, "scaling_with_intercept.csv"), row.names = FALSE)
}

m <- pm[, c("id", "method", "robust", "n_genes", "n_samples", "cpu_seconds",
            "deg_fdr05", "deg_fdr01", "deg_fdr05_lfc0.2", "up_fdr05",
            "down_fdr05", "pi0_l50", "lambda_gc")]
if (!is.null(od)) m <- merge(m, od, by = "id", all.x = TRUE)
# The null and split-half stages were run on the plain estimators only, so
# their columns attach by method and are NA for the robust rows rather than
# being silently copied across.
if (!is.null(pn)) {
  pn2 <- pn[, c("method", "deg05_mean", "deg05_max", "frac_p05_mean",
                "excess_type1_at_05", "excess_type1_at_01", "lambda_mean")]
  names(pn2) <- c("method", "null_deg05_mean", "null_deg05_max",
                  "null_frac_p05", "null_excess_t1_at05", "null_excess_t1_at01",
                  "null_lambda")
  m$.k <- ifelse(m$robust, NA, m$method)
  m <- merge(m, pn2, by.x = ".k", by.y = "method", all.x = TRUE); m$.k <- NULL
}
if (!is.null(sh)) {
  sh2 <- sh[, c("method", "deg_half_mean", "repro_mean", "spearman_mean",
                "top500_mean", "sign_agree_mean")]
  names(sh2) <- c("method", "sh_deg_half_mean", "sh_repro", "sh_spearman",
                  "sh_top500", "sh_sign_agree")
  m$.k <- ifelse(m$robust, NA, m$method)
  m <- merge(m, sh2, by.x = ".k", by.y = "method", all.x = TRUE); m$.k <- NULL
}
ord <- c("voom", "voomWQW", "trend", "QLF", "LRT")
m <- m[order(match(m$method, ord), m$robust), ]

notes <- readLines(file.path(OUT, "NOTES.md"), warn = FALSE)
notes <- if (length(notes)) notes else NULL
write_run("method_sweep", m, notes = notes)
log("\nmaster summary: %d rows x %d cols\n", nrow(m), ncol(m))
print(m[, intersect(c("id", "deg_fdr05", "mean_jaccard", "null_frac_p05",
                      "sh_repro"), names(m))], row.names = FALSE)
