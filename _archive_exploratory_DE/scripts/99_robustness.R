#!/usr/bin/env Rscript
# ==============================================================================
# 99_robustness.R -- the remaining experiments lost to the session limit.
#
#   A. covariate leave-one-out, and each covariate as its own exposure
#   B. sample-set sensitivity (does the depth/RIN decision matter?)
#   C. tobacco multiple-imputation sensitivity across all 20 completed columns
#   D. p-value histogram shape -- the single most diagnostic check in DE
# ==============================================================================

source("scripts/common.R")
p <- load_prep()
set.seed(481)
OUT <- file.path(EXP_HOME, "runs", "robustness")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ===================== A. covariates =======================================
cat("================ A. covariate leave-one-out (lithium contrast) ============\n")
base <- de_fit(p, "lithium")
base_ids <- deg_ids(base$tt)
COV <- c("age", "sex", "tob", "rin", "plate", "seqpc1", "seqpc2", "seqpc3")
A <- do.call(rbind, lapply(c("<none dropped>", COV), function(cv) {
  r <- if (cv == "<none dropped>") base else de_fit(p, "lithium", drop_covars = cv)
  u <- intersect(base$tt$gene, r$tt$gene)
  a <- intersect(base_ids, u); b <- intersect(deg_ids(r$tt), u)
  data.frame(dropped = cv, n_deg = n_deg(r$tt),
             jaccard_vs_base = round(length(intersect(a, b)) / length(union(a, b)), 3),
             t_spearman = round(cor(base$tt$t[match(u, base$tt$gene)],
                                    r$tt$t[match(u, r$tt$gene)], method = "spearman"), 3),
             stringsAsFactors = FALSE)
}))
# Fully unadjusted, for scale.
A <- rbind(A, {
  r <- de_fit(p, "lithium", drop_covars = COV)
  u <- intersect(base$tt$gene, r$tt$gene)
  a <- intersect(base_ids, u); b <- intersect(deg_ids(r$tt), u)
  data.frame(dropped = "<ALL covariates>", n_deg = n_deg(r$tt),
             jaccard_vs_base = round(length(intersect(a, b)) / length(union(a, b)), 3),
             t_spearman = round(cor(base$tt$t[match(u, base$tt$gene)],
                                    r$tt$t[match(u, r$tt$gene)], method = "spearman"), 3))
})
print(A, row.names = FALSE)

cat("\n  how much transcriptome signal does each variable carry ON ITS OWN?\n")
cat("  (same 226 subjects, each variable as the exposure, same covariate set)\n")
m0 <- p$meta[p$meta$dx == "BP1", ]; m0$tob <- p$TOB[m0$title, "tobacco_imp_01"]
m0 <- m0[complete.cases(m0[, c(COV, "lithium")]), ]
m0$plate <- droplevels(m0$plate); m0$sex <- droplevels(m0$sex)
cnt0 <- p$counts[, m0$title, drop = FALSE]
cnt0 <- cnt0[filter_genes(cnt0, 10, 0.90), , drop = FALSE]
dge0 <- calcNormFactors(DGEList(cnt0), method = "TMM")
m0$b1 <- p$ILR[m0$title, "b1_myeloid_vs_lymphoid"]
solo <- do.call(rbind, lapply(c("lithium", "b1", COV), function(v) {
  others <- setdiff(c(COV, "lithium"), v)
  X <- model.matrix(as.formula(paste("~", v, "+", paste(others, collapse = " + "))), data = m0)
  ne <- nonEstimable(X); if (!is.null(ne)) X <- X[, setdiff(colnames(X), ne), drop = FALSE]
  cf <- grep(paste0("^", v), colnames(X), value = TRUE)[1]
  tt <- topTable(eBayes(lmFit(voom(dge0, X), X)), coef = cf, number = Inf, sort.by = "none")
  data.frame(exposure = v, n_deg = sum(tt$adj.P.Val < 0.05),
             pi0 = round(pi0_storey(tt$P.Value), 3), stringsAsFactors = FALSE)
}))
print(solo[order(-solo$n_deg), ], row.names = FALSE)
write.csv(A, file.path(OUT, "covariate_loo.csv"), row.names = FALSE)
write.csv(solo, file.path(OUT, "covariate_solo.csv"), row.names = FALSE)

## ===================== B. sample set =======================================
cat("\n================ B. sample-set sensitivity ================================\n")
scr <- p$meta$depth_ok & p$meta$rin_ok
B <- do.call(rbind, list(
  {r <- base; data.frame(set = "all BP1", n = r$n_samples, n_deg = n_deg(r$tt),
                         jaccard = 1, stringsAsFactors = FALSE)},
  {r <- de_fit(p, "lithium", samples = scr)
   u <- intersect(base$tt$gene, r$tt$gene)
   a <- intersect(base_ids, u); b <- intersect(deg_ids(r$tt), u)
   data.frame(set = "depth>=3M & RIN>=5", n = r$n_samples, n_deg = n_deg(r$tt),
              jaccard = round(length(intersect(a, b)) / length(union(a, b)), 3))}))
print(B, row.names = FALSE)
cat("\n  100 random subsamples at n=200, to show how much DEG count moves from sampling alone:\n")
sub <- vapply(1:100, function(i) {
  keep <- rep(FALSE, nrow(p$meta))
  idx <- which(p$meta$dx == "BP1")
  keep[sample(idx, 200)] <- TRUE
  n_deg(de_fit(p, "lithium", samples = keep)$tt)
}, numeric(1))
cat(sprintf("  n=200 subsamples: median %.0f  IQR %.0f-%.0f  range %.0f-%.0f\n",
            median(sub), quantile(sub, .25), quantile(sub, .75), min(sub), max(sub)))
write.csv(B, file.path(OUT, "sample_set.csv"), row.names = FALSE)
write.csv(data.frame(draw = 1:100, n_deg = sub), file.path(OUT, "subsample_n200.csv"),
          row.names = FALSE)

## ===================== C. tobacco imputations ==============================
cat("\n================ C. tobacco multiple imputation (casecon) =================\n")
cat("  the 30 imputed subjects are 29 controls + 1 BP2, so casecon is where they act\n")
mi <- lapply(colnames(p$TOB), function(tc) de_fit(p, "casecon", tobcol = tc)$tt)
names(mi) <- colnames(p$TOB)
cnts <- vapply(mi, n_deg, numeric(1))
cat(sprintf("  DEG per imputation: median %.0f  range %.0f-%.0f  SD %.1f\n",
            median(cnts), min(cnts), max(cnts), sd(cnts)))
allg <- Reduce(intersect, lapply(mi, function(t) t$gene))
degs <- lapply(mi, function(t) intersect(deg_ids(t), allg))
tab <- table(unlist(degs))
cat(sprintf("  DEG in ALL 20 imputations: %d | in >=1: %d | in exactly 1: %d\n",
            sum(tab == 20), length(tab), sum(tab == 1)))
# Rubin pooling at the gene level.
lfc <- sapply(mi, function(t) t$logFC[match(allg, t$gene)])
se  <- sapply(mi, function(t) (t$logFC / t$t)[match(allg, t$gene)])
qbar <- rowMeans(lfc); wbar <- rowMeans(se^2); bvar <- apply(lfc, 1, var)
tot <- wbar + (1 + 1/20) * bvar
pooled_p <- 2 * pnorm(-abs(qbar / sqrt(tot)))
naive_p  <- 2 * pnorm(-abs(qbar / sqrt(wbar)))
cat(sprintf("  Rubin-pooled DEG: %d | ignoring imputation uncertainty: %d\n",
            sum(p.adjust(pooled_p, "BH") < 0.05), sum(p.adjust(naive_p, "BH") < 0.05)))
cat(sprintf("  SEs are understated by a median of %.2f%% if imputation uncertainty is ignored\n",
            100 * (median(sqrt(tot / wbar)) - 1)))
obs <- de_fit(p, "casecon", tobcol = "obs")
cat(sprintf("  observed-tobacco-only (drops the 30): n=%d, DEG=%d\n",
            obs$n_samples, n_deg(obs$tt)))
write.csv(data.frame(imputation = names(cnts), n_deg = cnts),
          file.path(OUT, "tobacco_mi.csv"), row.names = FALSE)

## ===================== D. p-value histograms ===============================
cat("\n================ D. p-value histogram shape ===============================\n")
cat("  A well-specified model gives a flat null with a spike near 0. A rising\n")
cat("  right tail means conservative; a U shape means misspecified.\n\n")
shape <- do.call(rbind, lapply(c("lithium", "casecon", "lithium_krebs"), function(cn) {
  tt <- de_fit(p, cn)$tt
  h <- hist(tt$P.Value, breaks = seq(0, 1, 0.05), plot = FALSE)$counts
  h <- h / sum(h)
  data.frame(contrast = cn, bin_0_05 = round(h[1], 4),
             mid_mean = round(mean(h[8:13]), 4), bin_95_1 = round(h[20], 4),
             tail_over_mid = round(h[20] / mean(h[8:13]), 3),
             spike_over_mid = round(h[1] / mean(h[8:13]), 2),
             pi0 = round(pi0_storey(tt$P.Value), 3), stringsAsFactors = FALSE)
}))
print(shape, row.names = FALSE)
cat("\n  tail_over_mid near 1 = well behaved; >1.2 = conservative; <0.8 = anticonservative\n")
for (i in seq_len(nrow(shape))) {
  v <- shape$tail_over_mid[i]
  cat(sprintf("  %-14s %s\n", shape$contrast[i],
              if (v > 1.2) "CONSERVATIVE -- investigate"
              else if (v < 0.8) "ANTICONSERVATIVE -- investigate"
              else "well specified"))
}
write.csv(shape, file.path(OUT, "pvalue_histogram.csv"), row.names = FALSE)
cat("\nwrote runs/robustness/\n")
