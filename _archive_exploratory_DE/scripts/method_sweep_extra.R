#!/usr/bin/env Rscript
# ==============================================================================
# method_sweep_extra.R -- three follow-ups to method_sweep.R.
#
#   G  DOWNWEIGHTING vs EXCLUDING. voomWQW is the only estimator in the sweep
#      that reacts to the 6 sub-3M-read and 3 sub-RIN-5 samples the project
#      chose to keep. If its weights are doing the job an exclusion would have
#      done, then plain voom on the 217 clean samples should land near
#      voomWQW on all 226. If it is not, the two are doing different things and
#      the "we kept them, but downweighted" story does not hold.
#
#   H  EFFECTIVE SAMPLE SIZE. limma scales the array weights to mean 1, so
#      their spread is exactly the information lost. Kish's n_eff turns that
#      spread into a number of samples, which is the only interpretable unit.
#
#   I  HEADLINE DEPENDENCE. The one-line answer to "how much does the
#      estimator move the reported DEG count".
#
# Reads the cached objects method_sweep.R wrote. Writes only inside
# experimentation/.
# ==============================================================================

source("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation/scripts/common.R")

OUT <- file.path(EXP_HOME, "runs", "method_sweep")
log <- function(...) { cat(sprintf(...), file = stderr()); flush(stderr()) }

TB   <- readRDS(file.path(EXP_HOME, "data", "method_sweep_tables.rds"))
FITS <- readRDS(file.path(EXP_HOME, "data", "method_sweep_fits.rds"))
p    <- load_prep()
md   <- TB$weights            # 226 BP1 rows, fit order, with $weight
w    <- md$weight
stopifnot(length(w) == 226)

## ===========================================================================
## I. how method-dependent is the headline?
## ===========================================================================
s <- TB$summ
nonrob <- s[!s$robust, ]
hl <- data.frame(
  set = c("all 10 configurations", "5 estimators, robust off", "5 estimators, robust on"),
  n = c(nrow(s), sum(!s$robust), sum(s$robust)),
  min_deg = c(min(s$deg_fdr05), min(s$deg_fdr05[!s$robust]), min(s$deg_fdr05[s$robust])),
  max_deg = c(max(s$deg_fdr05), max(s$deg_fdr05[!s$robust]), max(s$deg_fdr05[s$robust])),
  stringsAsFactors = FALSE)
hl$fold_range <- round(hl$max_deg / hl$min_deg, 3)
hl$cv <- round(c(sd(s$deg_fdr05) / mean(s$deg_fdr05),
                 sd(s$deg_fdr05[!s$robust]) / mean(s$deg_fdr05[!s$robust]),
                 sd(s$deg_fdr05[s$robust]) / mean(s$deg_fdr05[s$robust])), 4)
hl$median_deg <- c(median(s$deg_fdr05), median(s$deg_fdr05[!s$robust]),
                   median(s$deg_fdr05[s$robust]))
log("\n== I: headline dependence (voom baseline = %d) ==\n", s$deg_fdr05[s$id == "voom"])
print(hl, row.names = FALSE)
write.csv(hl, file.path(OUT, "headline_dependence.csv"), row.names = FALSE)

# The robust switch, held within estimator, so its effect is separated from
# the estimator's.
rob <- do.call(rbind, lapply(unique(s$method), function(m) {
  a <- s$deg_fdr05[s$method == m & !s$robust]; b <- s$deg_fdr05[s$method == m & s$robust]
  data.frame(method = m, deg_plain = a, deg_robust = b, delta = b - a,
             pct_change = round(100 * (b - a) / a, 2),
             jaccard = round(TB$jac[m, paste0(m, "_rob")], 4),
             spearman = round(TB$sp_stat[m, paste0(m, "_rob")], 5),
             stringsAsFactors = FALSE)
}))
log("\n== I: effect of robust=TRUE, within estimator ==\n"); print(rob, row.names = FALSE)
write.csv(rob, file.path(OUT, "robust_effect.csv"), row.names = FALSE)

## ===========================================================================
## H. effective sample size implied by the voomWQW weights
## ===========================================================================
kish <- function(x) sum(x)^2 / sum(x^2)
eff <- data.frame(
  group = c("all BP1", "lithium non-users", "lithium users",
            "the 9 flagged (low depth or RIN)", "the other 217"),
  n = c(226, sum(md$lithium == 0), sum(md$lithium == 1),
        sum(md$low_depth | md$low_rin), sum(!(md$low_depth | md$low_rin))),
  mean_weight = round(c(mean(w), mean(w[md$lithium == 0]), mean(w[md$lithium == 1]),
                        mean(w[md$low_depth | md$low_rin]),
                        mean(w[!(md$low_depth | md$low_rin)])), 4),
  stringsAsFactors = FALSE)
eff$n_eff <- round(c(kish(w), kish(w[md$lithium == 0]), kish(w[md$lithium == 1]),
                     kish(w[md$low_depth | md$low_rin]),
                     kish(w[!(md$low_depth | md$low_rin)])), 2)
eff$samples_lost <- round(eff$n - eff$n_eff, 2)
eff$pct_lost <- round(100 * eff$samples_lost / eff$n, 2)
log("\n== H: Kish effective sample size from the voomWQW weights ==\n")
log("   mean weight over all 226 = %.6f (limma scales array weights to mean 1)\n", mean(w))
print(eff, row.names = FALSE)
write.csv(eff, file.path(OUT, "effective_sample_size.csv"), row.names = FALSE)

## ===========================================================================
## G. downweighting vs excluding
## ===========================================================================
flagged <- md$low_depth | md$low_rin
drop_titles <- md$title[flagged]
keep_clean  <- !(p$meta$title %in% drop_titles)        # logical over all 474
# The 9 lowest-weighted samples, whatever they are -- so that "exclude the
# samples voomWQW distrusts" can be compared with "exclude the samples our
# pre-declared thresholds distrust". They need not be the same nine.
low9        <- md$title[order(w)][1:9]
keep_low9   <- !(p$meta$title %in% low9)

variants <- list(
  list(id = "voom_excl9_QC",       method = "voom",    samples = keep_clean),
  list(id = "voomWQW_excl9_QC",    method = "voomWQW", samples = keep_clean),
  list(id = "voom_excl9_lowest_w", method = "voom",    samples = keep_low9))
# The two all-226 reference fits already exist in the stage-A cache. Refitting
# them would cost another voomWQW run (~10 CPU minutes) to reproduce a number
# that is already on disk, and would risk the two differing.
res <- list(voom_all226 = FITS$voom, voomWQW_all226 = FITS$voomWQW)
for (v in variants) {
  f <- de_fit(p, "lithium", method = v$method, samples = v$samples)
  res[[v$id]] <- f
  log("[G] %-20s n=%3d genes=%5d DEG=%4d\n", v$id, f$n_samples, f$n_genes, n_deg(f$tt))
}
res <- res[c("voom_all226", "voomWQW_all226", "voom_excl9_QC",
             "voomWQW_excl9_QC", "voom_excl9_lowest_w")]
gtab <- do.call(rbind, lapply(names(res), function(id) {
  f <- res[[id]]
  data.frame(id = id, n_samples = f$n_samples, n_genes = f$n_genes,
             deg_fdr05 = n_deg(f$tt), deg_fdr01 = n_deg(f$tt, 0.01),
             pi0 = round(pi0_storey(f$tt$P.Value), 4), stringsAsFactors = FALSE)
}))
jg <- function(a, b) {
  A <- deg_ids(res[[a]]$tt); B <- deg_ids(res[[b]]$tt)
  round(length(intersect(A, B)) / length(union(A, B)), 4)
}
gtab$jaccard_vs_voom_all <- vapply(names(res), function(id) jg("voom_all226", id), numeric(1))
gtab$jaccard_vs_voomWQW_all <- vapply(names(res), function(id) jg("voomWQW_all226", id), numeric(1))
# Rank correlation on the genes every variant tested, so the comparison is not
# contaminated by the filter having kept different genes in each subset.
common <- Reduce(intersect, lapply(res, function(f) f$tt$gene))
S <- sapply(res, function(f) f$tt$t[match(common, f$tt$gene)])
sg <- round(cor(S, method = "spearman"), 4)
log("\n== G: downweighting vs excluding (%d genes common to all variants) ==\n", length(common))
print(gtab, row.names = FALSE)
log("\n== G: Spearman of t-statistics ==\n"); print(sg)
log("\n  QC-flagged nine : %s\n", paste(sort(drop_titles), collapse = " "))
log("  nine lowest-weight: %s\n", paste(sort(low9), collapse = " "))
log("  overlap between the two sets: %d of 9\n", length(intersect(drop_titles, low9)))
write.csv(gtab, file.path(OUT, "downweight_vs_exclude.csv"), row.names = FALSE)
write.csv(sg, file.path(OUT, "downweight_vs_exclude_spearman.csv"))
writeLines(c(paste("QC_flagged_nine:", paste(sort(drop_titles), collapse = " ")),
             paste("lowest_weight_nine:", paste(sort(low9), collapse = " ")),
             paste("overlap:", length(intersect(drop_titles, low9)))),
           file.path(OUT, "excluded_sample_sets.txt"))

## ===========================================================================
## K/L. what KIND of gene does each method call?
## ===========================================================================
# data/gene_annot.csv is produced by a sibling experiment in this folder. It is
# read defensively: if it is absent or has changed shape, these two stages are
# skipped rather than guessed at, because nothing else here depends on them.
ANN <- file.path(EXP_HOME, "data", "gene_annot.csv")
if (file.exists(ANN)) {
  ann <- read.csv(ANN, stringsAsFactors = FALSE)
  if (all(c("gene", "symbol", "biotype") %in% names(ann))) {
    rownames(ann) <- ann$gene
    G <- rownames(TB$PVAL)
    bt <- ann[G, "biotype"]
    bt[is.na(bt)] <- "unannotated"
    # Collapse to the distinction that matters for this question: a method
    # that is picking up low-count noise shows it as an excess of pseudogene
    # and non-coding calls relative to the tested background.
    cls <- ifelse(bt == "protein_coding", "protein_coding",
           ifelse(grepl("pseudogene", bt), "pseudogene", "other_noncoding"))
    bg  <- prop.table(table(cls))
    rows <- lapply(colnames(TB$ADJ), function(id) {
      sel <- TB$ADJ[, id] < 0.05
      tb  <- table(factor(cls[sel], levels = names(bg)))
      data.frame(id = id, n_deg = sum(sel),
                 pct_protein_coding = round(100 * tb[["protein_coding"]] / sum(tb), 2),
                 pct_pseudogene = round(100 * tb[["pseudogene"]] / sum(tb), 2),
                 pct_other_noncoding = round(100 * tb[["other_noncoding"]] / sum(tb), 2),
                 enrich_pseudogene = round((tb[["pseudogene"]] / sum(tb)) / bg[["pseudogene"]], 3),
                 stringsAsFactors = FALSE)
    })
    biot <- do.call(rbind, rows)
    log("\n== K: biotype of the called genes (background: %.1f%% coding, %.1f%% pseudogene, %.1f%% other) ==\n",
        100 * bg[["protein_coding"]], 100 * bg[["pseudogene"]], 100 * bg[["other_noncoding"]])
    print(biot, row.names = FALSE)
    write.csv(biot, file.path(OUT, "deg_biotype_composition.csv"), row.names = FALSE)

    # The genes the estimators disagree about most, named. Disagreement is
    # measured as the spread of -log10(p) across the five non-robust methods.
    nr  <- c("voom", "voomWQW", "trend", "QLF", "LRT")
    L10 <- -log10(TB$PVAL[, nr, drop = FALSE])
    spread <- apply(L10, 1, function(x) max(x) - min(x))
    anyd   <- rowSums(TB$ADJ[, nr, drop = FALSE] < 0.05, na.rm = TRUE)
    # Only genes at least one method calls: a huge spread among genes nobody
    # calls is arithmetically real and scientifically uninteresting.
    cand <- which(anyd >= 1 & anyd <= 4)
    ordc <- cand[order(-spread[cand])][1:min(30, length(cand))]
    disag <- data.frame(
      gene = G[ordc], symbol = ann[G[ordc], "symbol"],
      biotype = ann[G[ordc], "biotype"],
      aveLogCPM = round(TB$ave[ordc], 3),
      logFC_voom = round(TB$LFC[ordc, "voom"], 4),
      n_methods_calling = anyd[ordc],
      spread_neglog10p = round(spread[ordc], 3),
      stringsAsFactors = FALSE)
    for (m in nr) disag[[paste0("adjP_", m)]] <- signif(TB$ADJ[ordc, m], 3)
    log("\n== L: the 12 genes the five estimators disagree about most ==\n")
    print(head(disag[, 1:7], 12), row.names = FALSE)
    write.csv(disag, file.path(OUT, "most_disputed_genes.csv"), row.names = FALSE)
    log("  genes called by at least one but not all five: %d of %d (%.1f%%)\n",
        sum(anyd >= 1 & anyd <= 4), nrow(L10),
        100 * sum(anyd >= 1 & anyd <= 4) / nrow(L10))
    log("  genes called by all five: %d ; by none: %d\n",
        sum(anyd == 5), sum(anyd == 0))
  } else log("\n[K/L] gene_annot.csv present but unexpected columns; skipped\n")
} else log("\n[K/L] gene_annot.csv absent; skipped\n")

## ===========================================================================
## J. three diagnostic figures
## ===========================================================================
# Base graphics only -- no new dependency. These exist to make the two claims
# in the notes checkable by eye rather than only by table.
png(file.path(OUT, "fig_weight_vs_quality.png"), 1100, 520, res = 120)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.4, 1))
cols <- ifelse(md$low_depth, "#c0392b", ifelse(md$low_rin, "#e67e22", "#7f8c8d"))
pch  <- ifelse(md$low_depth | md$low_rin, 19, 1)
plot(md$lib_M, w, col = cols, pch = pch, log = "x",
     xlab = "library size (millions of assigned reads, log scale)",
     ylab = "voomWQW sample weight", main = "weight vs depth")
abline(v = 3, lty = 2, col = "#c0392b"); abline(h = 1, lty = 3)
legend("bottomright", c("< 3M reads", "RIN < 5", "passes both"),
       col = c("#c0392b", "#e67e22", "#7f8c8d"), pch = c(19, 19, 1), bty = "n", cex = .7)
plot(md$rin, w, col = cols, pch = pch, xlab = "RIN",
     ylab = "voomWQW sample weight", main = "weight vs RIN")
abline(v = 5, lty = 2, col = "#e67e22"); abline(h = 1, lty = 3)
dev.off()

png(file.path(OUT, "fig_pvalue_histograms.png"), 1400, 620, res = 120)
par(mfrow = c(2, 5), mar = c(3.6, 3.6, 2.4, 1), mgp = c(2.2, .7, 0))
for (id in colnames(TB$PVAL))
  hist(TB$PVAL[, id], breaks = 50, col = "#5b8bb5", border = NA,
       xlab = "p", main = sprintf("%s\n%d DEG", id, TB$summ$deg_fdr05[TB$summ$id == id]))
dev.off()

png(file.path(OUT, "fig_concordance.png"), 900, 820, res = 120)
par(mar = c(7, 7, 3, 2))
M <- TB$jac[nrow(TB$jac):1, ]
image(seq_len(ncol(M)), seq_len(nrow(M)), t(M), col = hcl.colors(32, "Blues", rev = TRUE),
      axes = FALSE, xlab = "", ylab = "", zlim = c(0, 1),
      main = "Jaccard of DEG sets, FDR 0.05")
axis(1, seq_len(ncol(M)), colnames(M), las = 2, cex.axis = .7)
axis(2, seq_len(nrow(M)), rownames(M), las = 2, cex.axis = .7)
for (i in seq_len(nrow(M))) for (j in seq_len(ncol(M)))
  text(j, i, sprintf("%.2f", t(M)[j, i]), cex = .55,
       col = if (t(M)[j, i] > .6) "white" else "black")
dev.off()
log("\n[extra] wrote 3 figures\n")
log("\n[extra] done\n")
