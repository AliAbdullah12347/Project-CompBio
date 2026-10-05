#!/usr/bin/env Rscript
# ==============================================================================
# norm_sweep_followup.R -- chase down what the main sweep turned up.
#
# The sweep produced three things that need explaining rather than reporting:
#
#   1. norm="none" gives 4013 DEG against TMM's 1382, but the logFC difference
#      between the two runs has sd 0.0034 around a median of -0.0713. That is
#      not a redistribution of signal, it is a near-perfect CONSTANT. Section F1
#      tests whether the constant equals the lithium difference in log2 scaling
#      factor, which would mean the 2631 extra genes are one number wearing
#      2631 hats.
#
#   2. Putting log2(TMM factor) into the design as a covariate drops DEG from
#      1382 to 158. That is a big claim -- a single scalar absorbing 89% of the
#      result -- and it is only interesting if it is SPECIFIC. Section F2 runs
#      the same test with negative controls (a random covariate, the same
#      factor permuted) and with rival global covariates, so the claim can be
#      graded rather than asserted.
#
#   3. The factor/lithium association might be an artefact of this contrast.
#      Section F6 checks the other two contrasts, and F5/F7 check whether the
#      association survives changes to how the factor is computed at all.
#
# Also extends the sweep past scaling factors into DISTRIBUTIONAL
# normalisation (quantile, cyclic loess), which the main sweep could not reach
# because de_fit() only exposes edgeR's scaling methods.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})
source("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation/scripts/common.R")

OUTDIR <- file.path(EXP_HOME, "runs", "norm_sweep")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
set.seed(20261218)
p   <- load_prep()
OBJ <- readRDS(file.path(OUTDIR, "norm_sweep_objects.rds"))
NF  <- OBJ$NF

gmap <- read.delim(file.path(EXP_HOME, "data", "gene_map.tsv"), header = FALSE,
                   col.names = c("gene", "symbol", "biotype"))
sym  <- setNames(gmap$symbol, gmap$gene)

d    <- build_design(p, "lithium")
cnt  <- p$counts[, d$samples, drop = FALSE]
gid  <- filter_genes(cnt, 10, 0.90)
cntf <- cnt[gid, , drop = FALSE]
meta <- d$meta
lith <- meta$lithium
Xcov <- d$X[, setdiff(colnames(d$X), "lithium"), drop = FALSE]
jacc <- function(a, b) length(intersect(a, b)) / length(union(a, b))
ref_tt  <- OBJ$fits_tt$voom_TMM
ref_deg <- ref_tt$gene[ref_tt$adj.P.Val < 0.05]

cat("=========================================================\n")
cat("norm_sweep_followup\n")
cat("=========================================================\n")

## ===========================================================================
## F1 -- is the norm="none" excess one constant offset?
## ===========================================================================
cat("\n### F1  the constant-offset decomposition ###\n")
a <- ref_tt[order(ref_tt$gene), ]
b <- OBJ$fits_tt$voom_none_CPM[order(OBJ$fits_tt$voom_none_CPM$gene), ]
stopifnot(identical(a$gene, b$gene))
shift <- b$logFC - a$logFC

# The quantity the offset SHOULD equal: the lithium coefficient when log2 of
# the TMM factor is itself the outcome, adjusted for the same covariates. If
# normalisation is doing nothing but removing a global lithium-aligned scale
# difference, dropping it must add back exactly that number to every gene.
nf_beta <- coef(lm(log2(NF[, "TMM"]) ~ lith + Xcov - 1))["lith"]
f1 <- data.frame(
  quantity = c("median(logFC_none - logFC_TMM)", "mean(same)", "sd(same)",
               "IQR(same)", "lithium coef on log2(TMM factor), adjusted",
               "ratio median_shift / nf_coef",
               "R^2 of shift on a constant (1 - var/var of logFC_TMM)"),
  value = c(median(shift), mean(shift), sd(shift), IQR(shift),
            -unname(nf_beta), median(shift) / -unname(nf_beta),
            1 - var(shift) / var(a$logFC)))
print(format(f1, digits = 5), row.names = FALSE)

# Direction and abundance of the genes that only norm="none" calls. A constant
# downward offset makes the most precisely measured (highest expressed, lowest
# residual variance) genes significant first, regardless of biology.
only_none <- setdiff(b$gene[b$adj.P.Val < 0.05], ref_deg)
shared    <- intersect(b$gene[b$adj.P.Val < 0.05], ref_deg)
ave <- rowMeans(cpm(calcNormFactors(DGEList(cntf)), log = TRUE))
f1b <- data.frame(
  set = c("DEG in both", "DEG only under norm=none", "never DEG"),
  n = c(length(shared), length(only_none),
        nrow(a) - length(union(shared, only_none))),
  mean_logCPM = c(mean(ave[shared]), mean(ave[only_none]),
                  mean(ave[setdiff(a$gene, union(shared, only_none))])),
  pct_down_in_users = c(
    100 * mean(b$logFC[match(shared, b$gene)] < 0),
    100 * mean(b$logFC[match(only_none, b$gene)] < 0),
    NA),
  median_abs_logFC_TMM = c(
    median(abs(a$logFC[match(shared, a$gene)])),
    median(abs(a$logFC[match(only_none, a$gene)])), NA))
print(format(f1b, digits = 4), row.names = FALSE)
write.csv(f1,  file.path(OUTDIR, "F1_offset_decomposition.csv"), row.names = FALSE)
write.csv(f1b, file.path(OUTDIR, "F1_extra_deg_profile.csv"),    row.names = FALSE)

## ===========================================================================
## F2 -- is the "factor as covariate" collapse specific?
## ===========================================================================
cat("\n### F2  specificity of the norm-factor covariate ###\n")
# PCs of the log-CPM matrix. Flagged as circular: these are computed from the
# same matrix that supplies the outcome, so a collapse here is partly guaranteed
# and is a comparator for the magnitude, not a competing explanation.
lcpm <- cpm(calcNormFactors(DGEList(cntf)), log = TRUE, prior.count = 3)
pcs  <- prcomp(t(lcpm), scale. = FALSE)$x[, 1:3]
colnames(pcs) <- paste0("exprPC", 1:3)
stopifnot(identical(rownames(pcs), meta$title))

cpm_raw <- sweep(cntf, 2, colSums(cntf), "/")
top50   <- apply(cpm_raw, 2, function(v) sum(sort(v, TRUE)[1:50]))

mk <- function(v, nm) { m <- matrix(v, ncol = 1,
  dimnames = list(meta$title, nm)); m }
CAND <- list(
  "none (baseline)"              = NULL,
  "log2 TMM factor"              = mk(log2(NF[, "TMM"]), "nf"),
  "log2 TMM factor, PERMUTED"    = mk(log2(NF[sample.int(nrow(NF)), "TMM"]), "nfperm"),
  "random N(0,1)"                = mk(rnorm(nrow(NF)), "rnd"),
  "log2 library size"            = mk(log2(colSums(cntf)), "lib"),
  "top-50 count share"           = mk(top50, "top50"),
  "granulocyte fraction"         = mk(p$lineage[meta$title, "gran"], "gran"),
  "ILR b1 (myeloid vs lymphoid)" = mk(p$ILR[meta$title, 1], "b1"),
  "all 4 ILR balances"           = p$ILR[meta$title, , drop = FALSE],
  "expression PC1 (circular)"    = pcs[, 1, drop = FALSE],
  "expression PC1-3 (circular)"  = pcs[, 1:3, drop = FALSE],
  "CIBERSORTx RMSE"              = mk(p$cs_qc[meta$title, "RMSE"], "rmse")
)
f2 <- do.call(rbind, lapply(names(CAND), function(nm) {
  A <- CAND[[nm]]
  f <- de_fit(p, "lithium", adjust = A, norm = "TMM", method = "voom")
  dg <- deg_ids(f$tt)
  # how much of the exposure's own variance the added covariate eats: a model
  # that kills DEG only by making lithium unidentifiable is not evidence of
  # mediation, it is evidence of collinearity, and these two columns separate
  # the cases
  r2l <- if (is.null(A)) 0 else summary(lm(lith ~ as.matrix(A) + Xcov - 1))$r.squared -
                                summary(lm(lith ~ Xcov - 1))$r.squared
  data.frame(covariate = nm, k = if (is.null(A)) 0L else ncol(A),
             deg_fdr05 = length(dg), pi0 = pi0_storey(f$tt$P.Value),
             pct_of_baseline = 100 * length(dg) / length(ref_deg),
             jaccard_vs_baseline = jacc(dg, ref_deg),
             extra_R2_on_lithium = r2l,
             se_inflation = {
               s <- sqrt(diag(solve(crossprod(f$design$X))))[f$coef]
               s0 <- sqrt(diag(solve(crossprod(d$X))))[["lithium"]]
               unname(s / s0) },
             stringsAsFactors = FALSE)
}))
print(format(f2, digits = 4), row.names = FALSE)
write.csv(f2, file.path(OUTDIR, "F2_covariate_specificity.csv"), row.names = FALSE)

## ===========================================================================
## F3 -- distributional normalisation (beyond scaling factors)
## ===========================================================================
cat("\n### F3  quantile and cyclic-loess normalisation ###\n")
X <- d$X
run_custom <- function(label, E, wts = NULL) {
  fit <- eBayes(lmFit(E, X, weights = wts))
  tt  <- topTable(fit, coef = "lithium", number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  tt
}
dge_tmm  <- calcNormFactors(DGEList(cntf))
dge_none <- DGEList(cntf)
CUSTOM <- list()
v <- voom(dge_tmm,  X, normalize.method = "quantile")
CUSTOM[["voom_TMM_quantile"]]  <- run_custom("", v$E, v$weights)
v <- voom(dge_none, X, normalize.method = "quantile")
CUSTOM[["voom_none_quantile"]] <- run_custom("", v$E, v$weights)
lc <- cpm(dge_none, log = TRUE, prior.count = 3)
CUSTOM[["trend_none_quantile"]]   <- { fit <- eBayes(lmFit(normalizeBetweenArrays(lc, "quantile"), X), trend = TRUE)
  tt <- topTable(fit, coef = "lithium", number = Inf, sort.by = "none"); tt$gene <- rownames(tt); tt }
CUSTOM[["trend_none_cyclicloess"]] <- { fit <- eBayes(lmFit(normalizeBetweenArrays(lc, "cyclicloess"), X), trend = TRUE)
  tt <- topTable(fit, coef = "lithium", number = Inf, sort.by = "none"); tt$gene <- rownames(tt); tt }
f3 <- do.call(rbind, lapply(names(CUSTOM), function(nm) {
  tt <- CUSTOM[[nm]]; dg <- deg_ids(tt)
  o  <- order(tt$gene)
  data.frame(config = nm, deg_fdr05 = length(dg), pi0 = pi0_storey(tt$P.Value),
             jaccard_vs_voomTMM = jacc(dg, ref_deg),
             spearman_t_vs_voomTMM = cor(tt$t[o], a$t, method = "spearman"),
             median_lfc_shift = median(tt$logFC[o] - a$logFC),
             stringsAsFactors = FALSE)
}))
print(format(f3, digits = 4), row.names = FALSE)
write.csv(f3, file.path(OUTDIR, "F3_distributional_norm.csv"), row.names = FALSE)

## ===========================================================================
## F4 -- does the factor/lithium link survive changes to how it is computed?
## ===========================================================================
cat("\n### F4  robustness of the factor/lithium association ###\n")
assoc_p <- function(y) {
  y <- as.numeric(y)
  m0 <- lm(y ~ Xcov - 1); m1 <- lm(y ~ lith + Xcov - 1)
  c(beta = unname(coef(m1)["lith"]), p = anova(m0, m1)$`Pr(>F)`[2],
    r = cor(residuals(m0), residuals(lm(lith ~ Xcov - 1))))
}
# (a) filtered vs unfiltered gene set. edgeR's own advice is to filter first,
# but plenty of pipelines do not, and the factor is a different number either
# way.
nf_unf <- calcNormFactors(DGEList(cnt))$samples$norm.factors
# (b) does the choice of TMM reference column matter? TMM picks one sample as
# the reference by default; if the association is an artefact of that choice it
# should wander as the reference moves.
refsweep <- t(vapply(seq_len(ncol(cntf)), function(i)
  assoc_p(log2(calcNormFactors(DGEList(cntf), method = "TMM", refColumn = i)$samples$norm.factors)),
  numeric(3)))
f4 <- rbind(
  data.frame(variant = "TMM on filtered genes (as used)",  t(assoc_p(log2(NF[, "TMM"])))),
  data.frame(variant = "TMM on all 57,773 genes",          t(assoc_p(log2(nf_unf)))),
  data.frame(variant = "TMM, reference column = median over all 226 choices",
             beta = median(refsweep[, 1]), p = median(refsweep[, 2]), r = median(refsweep[, 3])),
  data.frame(variant = "TMM, reference column = WORST of 226 choices",
             beta = refsweep[which.max(refsweep[, 2]), 1],
             p = max(refsweep[, 2]), r = refsweep[which.max(refsweep[, 2]), 3]))
print(format(f4, digits = 4), row.names = FALSE)
cat(sprintf("reference-column sweep: p < 0.05 in %d of %d choices; max p = %.4g\n",
            sum(refsweep[, 2] < 0.05), nrow(refsweep), max(refsweep[, 2])))
write.csv(f4, file.path(OUTDIR, "F4_factor_association_robustness.csv"), row.names = FALSE)

## ===========================================================================
## F5 -- is this lithium-specific, or does every exposure do it?
## ===========================================================================
cat("\n### F5  the same test in the other two contrasts ###\n")
f5 <- do.call(rbind, lapply(c("lithium", "casecon", "lithium_krebs"), function(ct) {
  dd <- build_design(p, ct)
  cc <- p$counts[, dd$samples, drop = FALSE]
  cc <- cc[filter_genes(cc, 10, 0.90), , drop = FALSE]
  nf <- log2(calcNormFactors(DGEList(cc))$samples$norm.factors)
  ex <- dd$X[, dd$coef]
  XC <- dd$X[, setdiff(colnames(dd$X), dd$coef), drop = FALSE]
  m0 <- lm(nf ~ XC - 1); m1 <- lm(nf ~ ex + XC - 1)
  gr <- p$lineage[dd$meta$title, "gran"]
  data.frame(contrast = ct, n = dd$n, exposure = dd$coef,
             beta_log2nf = unname(coef(m1)["ex"]),
             p_adj = anova(m0, m1)$`Pr(>F)`[2],
             partial_r = cor(residuals(m0), residuals(lm(ex ~ XC - 1))),
             r_nf_gran = cor(nf, gr),
             gran_beta = unname(coef(lm(gr ~ ex + XC - 1))["ex"]),
             gran_p = anova(lm(gr ~ XC - 1), lm(gr ~ ex + XC - 1))$`Pr(>F)`[2],
             stringsAsFactors = FALSE)
}))
print(format(f5, digits = 4), row.names = FALSE)
write.csv(f5, file.path(OUTDIR, "F5_other_contrasts.csv"), row.names = FALSE)

## ===========================================================================
## F6 -- the mechanism: composition -> dominance -> factor
## ===========================================================================
cat("\n### F6  mediation chain, granulocytes -> dominance -> factor ###\n")
gran <- p$lineage[meta$title, "gran"]
y    <- log2(NF[, "TMM"])
# Which transcripts actually set the factor? Correlate each gene's share of the
# library with the factor, then look at what sits at the top. If the mechanism
# story is right these should be neutrophil transcripts.
share <- t(apply(cpm_raw, 1, function(v) v))
r_gene <- apply(cpm_raw, 1, function(v) cor(v, y, method = "spearman"))
top_drivers <- head(sort(r_gene), 25)
cat("genes whose library share most strongly ANTI-correlates with the TMM factor\n")
cat("(these are the transcripts whose dominance pushes the factor down):\n")
print(data.frame(gene = names(top_drivers), symbol = unname(sym[names(top_drivers)]),
                 spearman_vs_log2nf = unname(top_drivers),
                 mean_pct_of_library = 100 * rowMeans(cpm_raw[names(top_drivers), ]),
                 r_with_gran = apply(cpm_raw[names(top_drivers), ], 1, function(v) cor(v, gran)),
                 row.names = NULL), digits = 3)
f6 <- data.frame(
  link = c("gran fraction ~ lithium (adj)", "top50 dominance ~ gran",
           "log2 TMM factor ~ top50 dominance", "log2 TMM factor ~ gran",
           "log2 TMM factor ~ lithium (adj)"),
  statistic = c("partial r", "pearson r", "pearson r", "pearson r", "partial r"),
  value = c(assoc_p(gran)["r"], cor(top50, gran), cor(y, top50), cor(y, gran),
            assoc_p(y)["r"]))
print(format(f6, digits = 4), row.names = FALSE)
write.csv(data.frame(gene = names(top_drivers), symbol = unname(sym[names(top_drivers)]),
                     spearman_vs_log2nf = unname(top_drivers)),
          file.path(OUTDIR, "F6_factor_driver_genes.csv"), row.names = FALSE)
write.csv(f6, file.path(OUTDIR, "F6_mechanism_chain.csv"), row.names = FALSE)

## ---- top DE genes with symbols, for the record ----------------------------
topg <- a[order(a$P.Value), ][1:30, ]
topg$symbol <- sym[topg$gene]
write.csv(topg[, c("gene", "symbol", "logFC", "t", "P.Value", "adj.P.Val")],
          file.path(OUTDIR, "top30_voom_TMM.csv"), row.names = FALSE)
cat("\ntop 10 lithium DEGs (voom+TMM):\n")
print(head(topg[, c("gene", "symbol", "logFC", "P.Value")], 10), row.names = FALSE)
cat("\ndone.\n")
