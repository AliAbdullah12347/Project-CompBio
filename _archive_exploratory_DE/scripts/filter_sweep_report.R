#!/usr/bin/env Rscript
# ==============================================================================
# filter_sweep_report.R -- final wrap-up.
#
# Routes the run through write_run() so runs/filter_sweep/ has the standard
# shape (summary.csv + NOTES.md) that every other experiment in this folder
# has, and prints the handful of numbers the write-up actually quotes, all
# recomputed here from the saved tables rather than transcribed.
#
# NOTES.md is authored as scripts/filter_sweep_NOTES.md and passed through, so
# the prose lives in version-controllable text and write_run only places it.
# ==============================================================================

source("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation/scripts/common.R")
RUN <- file.path(EXP_HOME, "runs", "filter_sweep")

S <- read.csv(file.path(RUN, "summary.csv"), stringsAsFactors = FALSE)
G <- S[S$ok & !is.na(S$min_count), ]

cat("== headline numbers (recomputed from summary.csv) ==\n")
cat(sprintf("grid cells fitted            %d of 36\n", nrow(G)))
cat(sprintf("genes tested                 %d .. %d  (%.2fx)\n",
            min(G$n_genes), max(G$n_genes), max(G$n_genes) / min(G$n_genes)))
cat(sprintf("DEG at FDR 0.05              %d .. %d  (%.2fx)\n",
            min(G$n_deg_05), max(G$n_deg_05), max(G$n_deg_05) / min(G$n_deg_05)))
cat(sprintf("DEG as %% of genes tested     %.2f%% .. %.2f%%  (%.2fx)\n",
            min(G$pct_deg_05), max(G$pct_deg_05), max(G$pct_deg_05) / min(G$pct_deg_05)))
cat(sprintf("pi0                          %.3f .. %.3f\n", min(G$pi0), max(G$pi0)))
cat(sprintf("implied non-null genes       %.0f .. %.0f\n",
            min((1 - G$pi0) * G$n_genes), max((1 - G$pi0) * G$n_genes)))
cat(sprintf("Jaccard vs canonical DEG set %.3f .. %.3f\n",
            min(G$jaccard_vs_canon), max(G$jaccard_vs_canon)))
cat(sprintf("fixed-universe DEG           %d .. %d\n",
            min(G$n_deg_fixedU), max(G$n_deg_fixedU)))
cat(sprintf("Spearman(t) vs canonical     %.4f .. %.4f\n",
            min(G$spearman_t_common, na.rm = TRUE), max(G$spearman_t_common, na.rm = TRUE)))
cat(sprintf("TMM factor corr vs canonical %.5f .. %.5f  (max abs diff %.4f)\n",
            min(G$tmm_cor_vs_canon, na.rm = TRUE), max(G$tmm_cor_vs_canon, na.rm = TRUE),
            max(G$tmm_maxabsdiff_vs_canon, na.rm = TRUE)))

sl <- coef(lm(log(n_deg_05) ~ log(n_genes), data = G[G$n_deg_05 > 0, ]))[2]
cat(sprintf("\nslope log(DEG) ~ log(genes)  %+.3f   [+1 = more tests give proportionally more hits; 0 = count independent of universe]\n", sl))
cat(sprintf("Spearman(n_genes, n_deg)     %+.3f\n", cor(G$n_genes, G$n_deg_05, method = "spearman")))
cat(sprintf("Spearman(n_genes, pct_deg)   %+.3f\n", cor(G$n_genes, G$pct_deg_05, method = "spearman")))
cat(sprintf("Spearman(n_genes, pi0)       %+.3f\n", cor(G$n_genes, G$pi0, method = "spearman")))

F <- S[is.na(S$min_count), ]
if (nrow(F)) {
  cat("\n== filterByExpr and no-filter ==\n")
  print(F[, c("tag", "n_genes", "n_deg_05", "pct_deg_05", "pi0",
              "jaccard_vs_canon", "frac_canon_deg_recovered", "n_deg_fixedU")],
        row.names = FALSE, digits = 4)
}


# ==============================================================================
# Is the 2-D grid actually 1-D?
#
# Reading the grid, cells reached by completely different (min_count, min_prop)
# pairs but retaining a similar number of genes give near-identical pi0. If that
# holds, min_count and min_prop are not two axes: they are two ways of turning
# one knob, and the knob is "how many genes survive". Test it by asking whether
# min_count and min_prop explain anything ABOUT pi0 and DEG% that n_genes has not
# already explained.
#
# Then the mechanism: if every rule removes genes in essentially the same order,
# a cell keeping n genes should be close to "the n most-expressed genes". Compare
# each cell's gene set to exactly that.
# ==============================================================================
cat("\n== is the grid 1-D? ==\n")
for (y in c("pi0", "pct_deg_05", "n_deg_05")) {
  f1 <- lm(reformulate("log(n_genes)", y), data = G)
  f2 <- lm(reformulate(c("log(n_genes)", "min_count", "min_prop"), y), data = G)
  f3 <- lm(reformulate(c("poly(log(n_genes), 2)"), y), data = G)
  a  <- anova(f1, f2)
  cat(sprintf("%-11s  R2(log n)=%.4f  R2(log n, quadratic)=%.4f  R2(+min_count+min_prop)=%.4f  added-terms p=%.3g\n",
              y, summary(f1)$r.squared, summary(f3)$r.squared,
              summary(f2)$r.squared, a$`Pr(>F)`[2]))
}

p <- load_prep()
d0 <- build_design(p, "lithium")
cnt <- p$counts[, d0$samples, drop = FALSE]
# Rank genes once by median CPM. "The n most-expressed genes" is then a single
# ordering, and every cell is compared against a prefix of it.
cpm_med <- apply(edgeR::cpm(cnt), 1, median)
ord <- names(sort(cpm_med, decreasing = TRUE))
J <- do.call(rbind, lapply(seq_len(nrow(G)), function(i) {
  f <- readRDS(file.path(RUN, "cache", paste0(G$tag[i], ".rds")))
  gs <- f$tt$gene; top <- ord[seq_len(length(gs))]
  data.frame(tag = G$tag[i], min_count = G$min_count[i], min_prop = G$min_prop[i],
             n_genes = length(gs),
             jaccard_vs_top_expressed = length(intersect(gs, top)) / length(union(gs, top)),
             stringsAsFactors = FALSE)
}))
J <- J[order(J$n_genes), ]
cat("\n== each cell's gene set vs simply taking the n most-expressed genes ==\n")
print(J, row.names = FALSE, digits = 4)
cat(sprintf("\nmedian Jaccard against a pure expression-rank cutoff: %.3f  (range %.3f-%.3f)\n",
            median(J$jaccard_vs_top_expressed), min(J$jaccard_vs_top_expressed),
            max(J$jaccard_vs_top_expressed)))
write.csv(J, file.path(RUN, "one_dimensional.csv"), row.names = FALSE)


# ==============================================================================
# Why filterByExpr's three parameterisations disagree so much.
#
# filterByExpr does not threshold raw counts. It converts min.count (10) to a
# CPM cutoff using the MEDIAN library size, then asks for that CPM in at least
# MinSampleSize samples -- and MinSampleSize is derived from whatever it is told
# about the groups:
#     group given   -> the smallest group's size
#     design given  -> 1 / max(hat(design))
#     neither       -> all samples
# then shrunk: values above large.n = 10 become 10 + (n - 10) * min.prop (0.7).
#
# A design carrying many nuisance columns (here: sequencing plate) has high
# leverage on its smallest cell, so 1/max(hat) is tiny and the filter becomes
# very lenient. Recomputing the constants makes that visible instead of leaving
# a three-fold difference in gene count unexplained.
# ==============================================================================
cat("\n== filterByExpr internals ==\n")
libs <- colSums(cnt)
cpm_cut <- 10 / median(libs) * 1e6
shrink <- function(n) if (n > 10) 10 + (n - 10) * 0.7 else n
h <- stats::hat(d0$X)
mss <- c(design = 1 / max(h),
         group = min(table(d0$meta$lithium)),
         nogroup = ncol(cnt))
CPM <- edgeR::cpm(cnt, lib.size = libs)
cat(sprintf("median library size %.0f -> CPM cutoff for min.count=10 is %.4f\n",
            median(libs), cpm_cut))
cat(sprintf("max leverage of the lithium design = %.4f (over %d design columns)\n",
            max(h), ncol(d0$X)))
# 1/max(hat) is set by the single most extreme sample in the covariate space,
# so one outlier decides the gene filter for all 226. Name it: the project
# briefing flags four QC-passing samples with RIN < 5, minimum 2.3, and that
# sample turns up here deciding how many genes get tested.
ho <- order(-h)[1:3]
cat("highest-leverage samples (these are what set MinSampleSize):\n")
print(data.frame(leverage = round(h[ho], 4), age = d0$meta$age[ho],
                 rin = d0$meta$rin[ho], seqpc2 = round(d0$meta$seqpc2[ho], 3),
                 lithium = d0$meta$lithium[ho]), row.names = FALSE)
cat(sprintf("mean leverage %.4f (= %d columns / %d rows); the top sample is %.1fx that\n",
            mean(h), ncol(d0$X), nrow(d0$X), max(h) / mean(h)))
FB <- do.call(rbind, lapply(names(mss), function(k) {
  m <- shrink(mss[[k]])
  keep <- rowSums(CPM >= cpm_cut) >= (m - 1e-14) & rowSums(cnt) >= 15 - 1e-14
  tg <- paste0("fbe_", k)
  f <- readRDS(file.path(RUN, "cache", paste0(tg, ".rds")))
  data.frame(tag = tg, raw_MinSampleSize = mss[[k]], shrunk = m,
             as_pct_of_samples = 100 * m / ncol(cnt),
             n_genes_recomputed = sum(keep), n_genes_in_fit = f$n_genes,
             n_deg = n_deg(f$tt), pi0 = pi0_storey(f$tt$P.Value),
             stringsAsFactors = FALSE)
}))
print(FB, row.names = FALSE, digits = 5)
cat(sprintf("\ncanonical filter for comparison: >10 counts in >=90%% of samples -> %d genes\n",
            S$n_genes[S$is_canon]))
write.csv(FB, file.path(RUN, "filterbyexpr_internals.csv"), row.names = FALSE)

# --- place NOTES.md and the standard summary.csv through the shared helper ---
notes_file <- file.path(EXP_HOME, "scripts", "filter_sweep_NOTES.md")
notes <- if (file.exists(notes_file)) readLines(notes_file, warn = FALSE) else NULL
write_run("filter_sweep", S, notes)
cat("runs/filter_sweep/ now holds:\n")
cat(paste0("  ", list.files(RUN), collapse = "\n"), "\n")
