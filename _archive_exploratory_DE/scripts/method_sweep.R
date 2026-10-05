#!/usr/bin/env Rscript
# ==============================================================================
# method_sweep.R -- vary the ESTIMATOR only, on the lithium contrast.
#
# Everything else is held at the project baseline: BP1 only (n=226), filter
# >10 counts in >=90% of samples, TMM, the full covariate set from
# CONTRASTS$lithium, tobacco completion column tobacco_imp_01.
#
# Question: how method-dependent is the headline DEG count, and is any method
# an outlier? Particular attention to voomWithQualityWeights, because this
# cohort deliberately retains 9 samples under 3M reads and 4 with RIN < 5 --
# voomWQW is the only estimator here that can notice that.
#
# Stages
#   A  fit the 10 configurations (5 methods x robust off/on), time each
#   B  concordance: Jaccard, Spearman, top-K overlap, clustering
#   C  voomWQW sample weights and what predicts them
#   D  mechanistic decomposition: weights vs re-estimated mean-variance trend
#   E  where along the abundance axis the methods disagree
#
# Writes only inside experimentation/.
# ==============================================================================

source("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation/scripts/common.R")

OUT   <- file.path(EXP_HOME, "runs", "method_sweep")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
CACHE <- file.path(EXP_HOME, "data", "method_sweep_fits.rds")

# stderr is unbuffered on Windows; stdout redirected to a file is not.
# Progress therefore goes to stderr so a long background run can be watched.
log <- function(...) { cat(sprintf(...), file = stderr()); flush(stderr()) }

p <- load_prep()

METHODS <- c("voom", "voomWQW", "trend", "QLF", "LRT")
GRID <- expand.grid(method = METHODS, robust = c(FALSE, TRUE),
                    stringsAsFactors = FALSE)
GRID$id <- sprintf("%s%s", GRID$method, ifelse(GRID$robust, "_rob", ""))

## ===========================================================================
## A. fit every configuration
## ===========================================================================
if (file.exists(CACHE)) {
  log("loading cached fits\n"); FITS <- readRDS(CACHE)
} else {
  FITS <- list()
  for (i in seq_len(nrow(GRID))) {
    id <- GRID$id[i]
    log("[A] %-12s ", id)
    t0 <- proc.time()
    f  <- de_fit(p, "lithium", method = GRID$method[i], robust = GRID$robust[i])
    el <- proc.time() - t0
    # Five other experiments are running on this machine, so wall clock is
    # contended and not a property of the estimator. CPU seconds is the
    # honest number; elapsed is kept only so the gap is visible.
    f$cpu <- el[["user.self"]] + el[["sys.self"]]
    f$elapsed <- el[["elapsed"]]
    f$id <- id; f$robust <- GRID$robust[i]
    FITS[[id]] <- f
    log("%5d genes  %4d DEG  cpu %6.1fs wall %6.1fs\n",
        f$n_genes, n_deg(f$tt), f$cpu, f$elapsed)
  }
  saveRDS(FITS, CACHE)
}
FITS <- FITS[GRID$id]

# The filter is a cross-sample operation, but the sample set is identical in
# every configuration, so the gene list must be too. Assert rather than assume:
# if it ever differed, every Jaccard below would be comparing different
# universes and would be silently meaningless.
GENES <- FITS[[1]]$tt$gene
for (f in FITS) stopifnot(identical(f$tt$gene, GENES))
NG <- length(GENES)
log("\ncommon gene universe: %d genes, %d samples\n", NG, FITS[[1]]$n_samples)

# A signed statistic, so that "these two rank genes the same way" and "these
# two agree on direction" collapse into one question. F and LR are unsigned by
# construction; sqrt() puts them on a t-like scale and logFC supplies the sign.
signed_stat <- function(tt, method) {
  if (method %in% c("voom", "voomWQW", "trend")) tt$t
  else if (method == "QLF") sign(tt$logFC) * sqrt(pmax(tt$F,  0))
  else                      sign(tt$logFC) * sqrt(pmax(tt$LR, 0))
}
STAT <- sapply(FITS, function(f) signed_stat(f$tt, f$method))
PVAL <- sapply(FITS, function(f) f$tt$P.Value)
LFC  <- sapply(FITS, function(f) f$tt$logFC)
ADJ  <- sapply(FITS, function(f) f$tt$adj.P.Val)
rownames(STAT) <- rownames(PVAL) <- rownames(LFC) <- rownames(ADJ) <- GENES

## ---- per-configuration summary --------------------------------------------
# lambda_GC: the genomic-inflation factor, median chi-square of the observed
# p-values over the median of a 1-df null. Under genuine signal it is
# legitimately > 1, so it means nothing in isolation -- it is only interpretable
# ACROSS methods on the same data, which is exactly this comparison.
lambda_gc <- function(pv) median(qchisq(1 - pv, 1), na.rm = TRUE) / qchisq(0.5, 1)

summ <- do.call(rbind, lapply(names(FITS), function(id) {
  f <- FITS[[id]]; tt <- f$tt
  data.frame(
    id = id, method = f$method, robust = f$robust,
    n_genes = f$n_genes, n_samples = f$n_samples,
    cpu_seconds = round(f$cpu, 2), wall_seconds = round(f$elapsed, 2),
    deg_fdr10 = n_deg(tt, 0.10), deg_fdr05 = n_deg(tt, 0.05),
    deg_fdr01 = n_deg(tt, 0.01), deg_fdr001 = n_deg(tt, 0.001),
    deg_fdr05_lfc0.1 = n_deg(tt, 0.05, 0.1),
    deg_fdr05_lfc0.2 = n_deg(tt, 0.05, 0.2),
    deg_fdr05_lfc0.5 = n_deg(tt, 0.05, 0.5),
    up_fdr05   = sum(tt$adj.P.Val < 0.05 & tt$logFC > 0, na.rm = TRUE),
    down_fdr05 = sum(tt$adj.P.Val < 0.05 & tt$logFC < 0, na.rm = TRUE),
    p_lt_05  = sum(tt$P.Value < 0.05,  na.rm = TRUE),
    p_lt_001 = sum(tt$P.Value < 0.001, na.rm = TRUE),
    pi0_l50 = round(pi0_storey(tt$P.Value, 0.5), 4),
    pi0_l80 = round(pi0_storey(tt$P.Value, 0.8), 4),
    lambda_gc = round(lambda_gc(tt$P.Value), 4),
    median_p = round(median(tt$P.Value, na.rm = TRUE), 4),
    max_abs_lfc = round(max(abs(tt$logFC)), 4),
    stringsAsFactors = FALSE)
}))
log("\n== A: per-configuration ==\n")
print(summ[, c("id", "cpu_seconds", "deg_fdr05", "deg_fdr01", "pi0_l50",
               "lambda_gc", "up_fdr05", "down_fdr05")], row.names = FALSE)
write.csv(summ, file.path(OUT, "per_method.csv"), row.names = FALSE)

## ===========================================================================
## B. concordance
## ===========================================================================
ids  <- names(FITS)
DEGS <- lapply(FITS, function(f) deg_ids(f$tt, 0.05)); names(DEGS) <- ids

pairmat <- function(fun) {
  M <- matrix(NA_real_, length(ids), length(ids), dimnames = list(ids, ids))
  for (a in ids) for (b in ids) M[a, b] <- fun(a, b)
  M
}
jac <- pairmat(function(a, b) {
  u <- length(union(DEGS[[a]], DEGS[[b]]))
  if (u == 0) NA else length(intersect(DEGS[[a]], DEGS[[b]])) / u })
# Containment: the fraction of the ROW method's DEGs that the COLUMN method
# also calls. Asymmetric on purpose -- it separates "disagrees with" from
# "is simply a smaller set than".
cont <- pairmat(function(a, b) {
  A <- DEGS[[a]]
  if (!length(A)) NA else length(intersect(A, DEGS[[b]])) / length(A) })

sp_stat <- cor(STAT, method = "spearman")
sp_p    <- cor(-log10(PVAL), method = "spearman")
pe_lfc  <- cor(LFC, method = "pearson")
sp_lfc  <- cor(LFC, method = "spearman")

topk <- function(K) pairmat(function(a, b)
  length(intersect(GENES[order(PVAL[, a])][seq_len(K)],
                   GENES[order(PVAL[, b])][seq_len(K)])) / K)
TOPK <- lapply(c(100, 500, 1000, 2000), topk)
names(TOPK) <- paste0("top", c(100, 500, 1000, 2000))

wr <- function(M, nm) write.csv(round(M, 5), file.path(OUT, nm))
wr(jac,     "concordance_jaccard_fdr05.csv")
wr(cont,    "concordance_containment_fdr05.csv")
wr(sp_stat, "concordance_spearman_statistic.csv")
wr(sp_p,    "concordance_spearman_neglog10p.csv")
wr(pe_lfc,  "concordance_pearson_logFC.csv")
wr(sp_lfc,  "concordance_spearman_logFC.csv")
for (nm in names(TOPK)) wr(TOPK[[nm]], sprintf("concordance_%s.csv", nm))

log("\n== B: Jaccard of DEG sets (FDR 0.05) ==\n");        print(round(jac, 3))
log("\n== B: Spearman of signed test statistic ==\n");     print(round(sp_stat, 4))
log("\n== B: Spearman of logFC ==\n");                     print(round(sp_lfc, 4))
log("\n== B: containment, row's DEGs also called by column ==\n"); print(round(cont, 3))

off <- function(M) { diag(M) <- NA; M }
oddness <- data.frame(
  id = ids,
  mean_jaccard  = round(rowMeans(off(jac),        na.rm = TRUE), 4),
  mean_spearman = round(rowMeans(off(sp_stat),    na.rm = TRUE), 4),
  mean_top500   = round(rowMeans(off(TOPK$top500), na.rm = TRUE), 4),
  stringsAsFactors = FALSE)
oddness <- oddness[order(oddness$mean_spearman), ]
log("\n== B: how far each configuration sits from the others ==\n")
print(oddness, row.names = FALSE)
write.csv(oddness, file.path(OUT, "oddness.csv"), row.names = FALSE)

hc <- hclust(as.dist(1 - sp_stat), method = "average")
log("\n== B: clustering on 1 - Spearman ==\n")
log("  order: %s\n", paste(hc$labels[hc$order], collapse = " "))
log("  heights: %s\n", paste(sprintf("%.4f", hc$height), collapse = " "))

## ===========================================================================
## C. voomWQW sample weights
## ===========================================================================
# Rebuilt with build_design()/filter_genes() from common.R, i.e. the exact
# objects de_fit() uses, so the weights belong to the fit reported above and
# not to a lookalike assembled a second way.
d   <- build_design(p, "lithium")
cnt <- p$counts[, d$samples, drop = FALSE]
g   <- filter_genes(cnt, 10, 0.90)
dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = cnt[g, , drop = FALSE]),
                              method = "TMM")
stopifnot(identical(rownames(dge), GENES))

vq <- limma::voomWithQualityWeights(dge, d$X, plot = FALSE)
vp <- limma::voom(dge, d$X, plot = FALSE)
stopifnot(!is.null(vq$targets$sample.weights))
w  <- as.numeric(vq$targets$sample.weights)
names(w) <- colnames(dge)

md <- d$meta                             # the 226 BP1 rows, in fit order
stopifnot(identical(md$title, names(w)))
md$weight   <- w
md$lib_raw  <- p$meta$lib[d$samples]     # total assigned reads, pre-filter
md$lib_M    <- md$lib_raw / 1e6
md$normfac  <- dge$samples$norm.factors
md$cs_cor   <- p$cs_qc[md$title, "Correlation"]
md$cs_rmse  <- p$cs_qc[md$title, "RMSE"]
md$cs_p     <- p$cs_qc[md$title, "P-value"]
md$low_depth <- !md$depth_ok
md$low_rin   <- !md$rin_ok

log("\n== C: voomWQW sample weights, n = %d ==\n", length(w))
log("  min %.4f  q05 %.4f  median %.4f  q95 %.4f  max %.4f  (ratio max/min %.1fx)\n",
    min(w), quantile(w, .05), median(w), quantile(w, .95), max(w), max(w) / min(w))
log("  n below 0.5: %d   below 0.25: %d   above 2: %d\n",
    sum(w < 0.5), sum(w < 0.25), sum(w > 2))
log("  low-depth (<3M) in this subset: %d   low-RIN (<5): %d\n",
    sum(md$low_depth), sum(md$low_rin))

cor_row <- function(v, nm) {
  ok <- is.finite(v) & is.finite(w)
  ct <- suppressWarnings(cor.test(v[ok], w[ok], method = "spearman"))
  data.frame(variable = nm, n = sum(ok), spearman_rho = round(unname(ct$estimate), 4),
             p = signif(ct$p.value, 3), stringsAsFactors = FALSE)
}
wcor <- rbind(
  cor_row(log10(md$lib_raw), "log10_library_size"),
  cor_row(md$lib_M,   "library_size_millions"),
  cor_row(md$rin,     "RIN"),
  cor_row(md$age,     "age"),
  cor_row(md$seqpc1,  "seq_metric_pc1"),
  cor_row(md$seqpc2,  "seq_metric_pc2"),
  cor_row(md$seqpc3,  "seq_metric_pc3"),
  cor_row(md$normfac, "TMM_norm_factor"),
  cor_row(md$cs_cor,  "CIBERSORTx_correlation"),
  cor_row(md$cs_rmse, "CIBERSORTx_RMSE"),
  cor_row(as.numeric(md$lithium), "lithium_exposure"),
  cor_row(as.numeric(md$sex == "M"), "sex_male"),
  cor_row(md$tob,     "tobacco"))
wcor <- wcor[order(-abs(wcor$spearman_rho)), ]
log("\n== C: what predicts the weight (Spearman) ==\n"); print(wcor, row.names = FALSE)
write.csv(wcor, file.path(OUT, "weight_correlates.csv"), row.names = FALSE)

grp_test <- function(flag, nm) {
  if (sum(flag) == 0 || sum(!flag) == 0)
    return(data.frame(group = nm, n_flagged = sum(flag), med_flagged = NA,
                      med_rest = NA, wilcox_p = NA, stringsAsFactors = FALSE))
  wt <- suppressWarnings(wilcox.test(w[flag], w[!flag]))
  data.frame(group = nm, n_flagged = sum(flag),
             med_flagged = round(median(w[flag]), 4),
             med_rest = round(median(w[!flag]), 4),
             wilcox_p = signif(wt$p.value, 3), stringsAsFactors = FALSE)
}
gt <- rbind(
  grp_test(md$low_depth,              "library < 3M reads"),
  grp_test(md$low_rin,                "RIN < 5"),
  grp_test(md$low_depth | md$low_rin, "either flag"),
  grp_test(md$lib_M < quantile(md$lib_M, 0.10), "bottom decile of depth"),
  grp_test(md$rin < quantile(md$rin, 0.10),     "bottom decile of RIN"),
  grp_test(md$lithium == 1,           "lithium user"),
  grp_test(md$sex == "M",             "male"))
log("\n== C: is the weight lower for the samples we deliberately kept? ==\n")
print(gt, row.names = FALSE)
write.csv(gt, file.path(OUT, "weight_group_tests.csv"), row.names = FALSE)

# How much of the weight is explained jointly, and by what. Log-weight because
# the weights are a multiplicative precision scale.
lmw <- lm(log(weight) ~ log10(lib_raw) + rin + seqpc1 + seqpc2 + seqpc3 +
            age + sex + tob + lithium, data = md)
log("\n== C: joint model for log(weight), adj R2 = %.4f ==\n",
    summary(lmw)$adj.r.squared)
print(round(summary(lmw)$coefficients, 4))
capture.output(summary(lmw), file = file.path(OUT, "weight_joint_model.txt"))

wt_tab <- md[order(md$weight),
             c("title", "weight", "lib_M", "rin", "lithium", "sex", "age",
               "plate", "low_depth", "low_rin", "cs_cor", "cs_rmse")]
wt_tab$weight <- round(wt_tab$weight, 4); wt_tab$lib_M <- round(wt_tab$lib_M, 3)
write.csv(wt_tab, file.path(OUT, "voomWQW_sample_weights.csv"), row.names = FALSE)
log("\n== C: 15 lowest-weighted samples ==\n"); print(head(wt_tab, 15), row.names = FALSE)
log("\n== C: 8 highest-weighted samples ==\n"); print(tail(wt_tab, 8), row.names = FALSE)

# Rank agreement between "voomWQW thinks this sample is bad" and the two
# quality axes we screened on but chose not to exclude.
log("\n== C: rank of the flagged samples among all 226 by weight (1 = lowest) ==\n")
rk <- rank(md$weight)
for (nm in c("low_depth", "low_rin")) {
  f <- md[[nm]]
  if (any(f)) log("  %-10s ranks: %s   (of 226)\n", nm,
                  paste(sort(round(rk[f])), collapse = " "))
}

## ===========================================================================
## D. decomposition: is voomWQW's difference the WEIGHTS or the re-fitted trend?
## ===========================================================================
# voomWithQualityWeights does two things at once: it estimates array weights,
# and it re-runs the voom mean-variance fit with those weights in place. Fitting
# plain voom's observation weights scaled by voomWQW's sample weights isolates
# the first from the second.
Wmix <- vp$weights * rep(w, each = nrow(vp$weights))
fit_mix  <- limma::eBayes(limma::lmFit(vp$E, d$X, weights = Wmix))
tt_mix   <- limma::topTable(fit_mix, coef = d$coef, number = Inf, sort.by = "none")
tt_mix$gene <- rownames(tt_mix)
stopifnot(identical(tt_mix$gene, GENES))

decomp <- data.frame(
  variant = c("voom (no sample weights)",
              "voom E + voomWQW weights (weights only)",
              "voomWQW (weights + re-fitted trend)"),
  deg_fdr05 = c(n_deg(FITS$voom$tt), n_deg(tt_mix), n_deg(FITS$voomWQW$tt)),
  deg_fdr01 = c(n_deg(FITS$voom$tt, .01), n_deg(tt_mix, .01), n_deg(FITS$voomWQW$tt, .01)),
  pi0 = round(c(pi0_storey(FITS$voom$tt$P.Value), pi0_storey(tt_mix$P.Value),
                pi0_storey(FITS$voomWQW$tt$P.Value)), 4),
  stringsAsFactors = FALSE)
decomp$jaccard_vs_voom <- round(c(
  1,
  length(intersect(deg_ids(tt_mix), DEGS$voom)) / length(union(deg_ids(tt_mix), DEGS$voom)),
  jac["voom", "voomWQW"]), 4)
log("\n== D: decomposing voomWQW ==\n"); print(decomp, row.names = FALSE)
write.csv(decomp, file.path(OUT, "voomWQW_decomposition.csv"), row.names = FALSE)

## ===========================================================================
## E. where on the abundance axis do the methods disagree?
## ===========================================================================
# Common binning variable for every method, taken from the voom baseline, so
# the bins are the same set of genes in every column of the table below.
ave <- FITS$voom$tt$AveExpr
bin <- cut(ave, breaks = quantile(ave, seq(0, 1, 0.1)), include.lowest = TRUE,
           labels = paste0("D", 1:10))
abund <- data.frame(decile = levels(bin),
                    median_logCPM = round(tapply(ave, bin, median), 3),
                    n = as.integer(table(bin)))
for (id in ids)
  abund[[id]] <- as.integer(tapply(ADJ[, id] < 0.05, bin, sum, na.rm = TRUE))
log("\n== E: DEG count by expression decile (D1 = lowest) ==\n")
print(abund, row.names = FALSE)
write.csv(abund, file.path(OUT, "deg_by_abundance_decile.csv"), row.names = FALSE)

# Same question for the genes the methods disagree about, relative to voom.
disc <- do.call(rbind, lapply(setdiff(ids, "voom"), function(id) {
  only_here <- setdiff(DEGS[[id]], DEGS$voom)
  only_voom <- setdiff(DEGS$voom, DEGS[[id]])
  data.frame(id = id,
             n_only_this = length(only_here), n_only_voom = length(only_voom),
             med_logCPM_only_this = round(median(ave[GENES %in% only_here]), 3),
             med_logCPM_only_voom = round(median(ave[GENES %in% only_voom]), 3),
             med_logCPM_shared = round(median(ave[GENES %in% intersect(DEGS[[id]], DEGS$voom)]), 3),
             med_absLFC_only_this = round(median(abs(LFC[GENES %in% only_here, id])), 4),
             med_absLFC_shared = round(median(abs(LFC[GENES %in% intersect(DEGS[[id]], DEGS$voom), id])), 4),
             stringsAsFactors = FALSE)
}))
log("\n== E: genes gained/lost relative to voom ==\n"); print(disc, row.names = FALSE)
write.csv(disc, file.path(OUT, "discordant_gene_profile.csv"), row.names = FALSE)
log("  (median logCPM over all %d tested genes: %.3f)\n", NG, median(ave))

## ===========================================================================
## F. the shrinkage parameters, so the robust effect is explained not just seen
## ===========================================================================
# `robust` does not change the model; it changes how hard the variance estimate
# is shrunk towards the prior and how much a genuinely-overdispersed gene is
# allowed to keep its own variance. Recording df.prior and s2.prior turns
# "robust moved the count" into a statement about which knob moved.
lcpm    <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
lim_fits <- list(
  voom    = limma::lmFit(vp, d$X),
  voomWQW = limma::lmFit(vq, d$X),
  trend   = limma::lmFit(lcpm, d$X))
prior_rows <- list()
for (nm in names(lim_fits)) for (rb in c(FALSE, TRUE)) {
  eb <- limma::eBayes(lim_fits[[nm]], trend = (nm == "trend"), robust = rb)
  dp <- eb$df.prior
  prior_rows[[length(prior_rows) + 1]] <- data.frame(
    id = paste0(nm, if (rb) "_rob" else ""), method = nm, robust = rb,
    df_prior_min = round(min(dp), 3), df_prior_med = round(median(dp), 3),
    df_prior_max = round(max(dp), 3),
    n_df_prior_downweighted = if (length(dp) > 1) sum(dp < median(dp) * 0.5) else NA_integer_,
    df_residual = eb$df.residual[1],
    s2_prior_med = signif(median(eb$s2.prior), 4),
    stringsAsFactors = FALSE)
}
# edgeR's equivalent knob. estimateDisp is the slow step, so it is run once per
# robust setting and both tests read off the same dispersion object.
for (rb in c(FALSE, TRUE)) {
  dd <- edgeR::estimateDisp(dge, d$X, robust = rb)
  f  <- edgeR::glmQLFit(dge, d$X, dispersion = dd$trended.dispersion, robust = rb)
  prior_rows[[length(prior_rows) + 1]] <- data.frame(
    id = paste0("edgeR", if (rb) "_rob" else ""), method = "edgeR", robust = rb,
    df_prior_min = round(min(f$df.prior), 3),
    df_prior_med = round(median(f$df.prior), 3),
    df_prior_max = round(max(f$df.prior), 3),
    n_df_prior_downweighted = sum(f$df.prior < median(f$df.prior) * 0.5),
    df_residual = round(median(dd$df.residual %||% NA_real_), 3),
    s2_prior_med = signif(dd$common.dispersion, 4),
    stringsAsFactors = FALSE)
  log("  edgeR robust=%-5s common.dispersion %.5f  prior.df %.3f  BCV %.4f\n",
      rb, dd$common.dispersion, dd$prior.df, sqrt(dd$common.dispersion))
}
priors <- do.call(rbind, prior_rows)
log("\n== F: shrinkage parameters (s2_prior_med is common.dispersion for edgeR) ==\n")
print(priors, row.names = FALSE)
write.csv(priors, file.path(OUT, "shrinkage_parameters.csv"), row.names = FALSE)

saveRDS(list(STAT = STAT, PVAL = PVAL, LFC = LFC, ADJ = ADJ, summ = summ,
             DEGS = DEGS, grid = GRID, jac = jac, cont = cont,
             sp_stat = sp_stat, sp_p = sp_p, pe_lfc = pe_lfc, sp_lfc = sp_lfc,
             TOPK = TOPK, hc = hc, weights = md, ave = ave,
             decomp = decomp, abund = abund, disc = disc, priors = priors),
        file.path(EXP_HOME, "data", "method_sweep_tables.rds"))
log("\n[A-F] done\n")
