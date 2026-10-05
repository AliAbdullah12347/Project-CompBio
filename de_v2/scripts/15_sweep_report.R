#!/usr/bin/env Rscript
# ==============================================================================
# 15_sweep_report.R -- assemble the filter sweeps into the tables worth reading.
#
# A DEG COUNT IS NOT COMPARABLE ACROSS FILTERS, so the count tables here are
# context, not the result. Changing the filter changes both the denominator and
# the number of tests the correction divides by, so two filters can disagree on
# a count while agreeing completely on which genes matter.
#
# The result is the STABILITY tables: how often each gene is called, how far its
# rank moves, and how much each filter's DEG set overlaps the baseline.
# ==============================================================================

ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "de_v2")); options(width = 185)
SW <- file.path("results", "sweep")
has <- function(f) file.exists(file.path(SW, f))
rd  <- function(f) read.csv(file.path(SW, f), stringsAsFactors = FALSE)
unpref <- function(x, p) sub(paste0(p, "\\."), "", x)

## ============================ WHOLE BLOOD ====================================
if (has("wb_sweep_summary.csv")) {
  s <- rd("wb_sweep_summary.csv")
  cat("================= WHOLE BLOOD =================\n")
  cat("\n--- genes called DE at BH 0.05, by filter ---\n")
  w <- reshape(s[, c("fid", "rule", "n_genes", "analysis", "BH_05")],
               idvar = c("fid", "rule", "n_genes"), timevar = "analysis", direction = "wide")
  names(w) <- unpref(names(w), "BH_05")
  w <- w[order(w$n_genes), ]
  w$pct_LI_raw <- round(100 * w$LI_raw / w$n_genes, 2)
  w$base <- ifelse(w$fid == "c10_p90", "<<", "")
  print(w[, c("fid", "rule", "n_genes", "LI_raw", "pct_LI_raw", "LI_ilr",
              "BPD_raw", "BPD_ilr", "base")], row.names = FALSE)

  cat("\n--- how the correction method interacts with the filter (LI_raw) ---\n")
  z <- s[s$analysis == "LI_raw", c("fid", "n_genes", "pi0", "BH_05", "storey_05", "BY_05", "holm_05")]
  print(z[order(z$n_genes), ], row.names = FALSE)

  cat("\n--- range across all 31 filters ---\n")
  for (a in unique(s$analysis)) {
    z <- s[s$analysis == a, ]
    cat(sprintf("  %-9s DEG %d to %d  |  %%genes %.2f to %.2f  |  pi0 %.3f to %.3f\n",
                a, min(z$BH_05), max(z$BH_05),
                min(100 * z$BH_05 / z$n_genes), max(100 * z$BH_05 / z$n_genes),
                min(z$pi0), max(z$pi0)))
  }
}

if (has("wb_sweep_overlap.csv")) {
  j <- rd("wb_sweep_overlap.csv")
  cat("\n--- overlap of each filter's DEG set with the baseline (LI_raw) ---\n")
  z <- j[j$analysis == "LI_raw", ]; z <- z[order(z$n_genes), ]
  print(z[, c("fid", "n_genes", "n_deg", "n_shared_with_baseline", "jaccard")], row.names = FALSE)
  cat(sprintf("\n  median Jaccard against baseline, LI_raw: %.3f\n",
              median(z$jaccard[z$fid != "c10_p90"], na.rm = TRUE)))
}

if (has("wb_gene_robustness.csv")) {
  r <- rd("wb_gene_robustness.csv")
  cat("\n--- per-gene stability ---\n")
  for (a in unique(r$analysis)) {
    z <- r[r$analysis == a, ]
    alw <- sum(z$n_settings_sig <- z$n_filters_sig == z$n_filters_tested)
    cat(sprintf("  %-9s %5d genes significant under >=1 filter | %4d under ALL (%.1f%%) | %4d under only one\n",
                a, nrow(z), sum(z$n_filters_sig == z$n_filters_tested),
                100 * sum(z$n_filters_sig == z$n_filters_tested) / nrow(z),
                sum(z$n_filters_sig == 1)))
  }
  cat("\n--- the 25 most robust lithium genes (whole blood, unadjusted) ---\n")
  z <- r[r$analysis == "LI_raw", ]
  z <- z[order(-z$frac_sig, z$median_rank), ]
  print(head(z[, c("symbol", "n_filters_sig", "n_filters_tested", "frac_sig",
                   "best_rank", "median_rank", "median_logFC", "min_BH")], 25), row.names = FALSE)

  cat("\n--- genes that appear at only one filter (least trustworthy), LI_raw ---\n")
  o <- z[z$n_filters_sig == 1, ]
  cat(sprintf("  %d such genes; median rank %.0f, median |logFC| %.3f\n",
              nrow(o), median(o$median_rank), median(abs(o$median_logFC))))
  if (nrow(o)) print(head(o[order(o$median_rank), c("symbol", "best_rank", "median_logFC", "min_BH")], 10),
                     row.names = FALSE)

  if (any(r$analysis == "BPD_raw")) {
    cat("\n--- bipolar, whole blood: every gene ever called, with its stability ---\n")
    z <- r[r$analysis == "BPD_raw", ]
    print(z[order(-z$frac_sig, z$median_rank),
            c("symbol", "n_filters_sig", "n_filters_tested", "frac_sig", "best_rank",
              "median_rank", "median_logFC", "min_BH")], row.names = FALSE)
  }
}

## ============================ CELL TYPE ======================================
if (has("ct_sweep_per_lineage.csv")) {
  per <- rd("ct_sweep_per_lineage.csv")
  cat("\n\n================= CELL TYPE =================\n")
  for (a in c("LI_raw", "BPD_raw")) {
    z <- per[per$analysis == a, ]
    if (!nrow(z)) next
    d <- reshape(z[, c("sid", "lineage", "deg_BH05")], idvar = "sid",
                 timevar = "lineage", direction = "wide")
    names(d) <- unpref(names(d), "deg_BH05")
    n <- reshape(z[, c("sid", "lineage", "n_genes")], idvar = "sid",
                 timevar = "lineage", direction = "wide")
    names(n) <- paste0("n_", unpref(names(n), "n_genes")); names(n)[1] <- "sid"
    lab <- unique(z[, c("sid", "kind", "label")])
    m <- merge(merge(lab, d, by = "sid"), n, by = "sid")
    m <- m[order(m$kind, m$sid), ]
    cat(sprintf("\n--- %s: genes DE per lineage (BH within lineage) ---\n", a))
    print(m[, c("sid", "kind", "label", "gran", "mono", "T", "NK", "B",
                "n_gran", "n_NK", "n_B")], row.names = FALSE)
  }
}

if (has("ct_sweep_pooled.csv")) {
  pool <- rd("ct_sweep_pooled.csv")
  cat("\n--- pooled across lineages, unique genes ---\n")
  pw <- reshape(pool[, c("sid", "analysis", "deg_unique_genes")], idvar = "sid",
                timevar = "analysis", direction = "wide")
  names(pw) <- unpref(names(pw), "deg_unique_genes")
  k <- unique(pool[, c("sid", "kind", "label")])
  print(merge(k, pw, by = "sid"), row.names = FALSE)
}

if (has("ct_gene_robustness.csv")) {
  r <- rd("ct_gene_robustness.csv")
  cat("\n--- per gene x lineage stability ---\n")
  for (a in unique(r$analysis)) {
    z <- r[r$analysis == a, ]
    cat(sprintf("  %-9s %4d gene x lineage pairs significant somewhere | %3d under every setting tested\n",
                a, nrow(z), sum(z$n_settings_sig == z$n_settings_tested)))
  }
  cat("\n--- the 20 most robust cell-type lithium findings ---\n")
  z <- r[r$analysis == "LI_raw", ]
  z <- z[order(-z$frac_sig, z$median_rank), ]
  print(head(z[, c("symbol", "lineage", "n_settings_sig", "n_settings_tested", "frac_sig",
                   "best_rank", "median_rank", "median_logFC", "min_BH")], 20), row.names = FALSE)
  cat("\n--- by lineage: how many robust findings each carries ---\n")
  z$always <- z$n_settings_sig == z$n_settings_tested
  print(as.data.frame(table(lineage = z$lineage, always_significant = z$always)), row.names = FALSE)
}
cat("\n")
