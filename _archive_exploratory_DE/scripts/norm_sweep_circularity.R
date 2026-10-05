#!/usr/bin/env Rscript
# ==============================================================================
# norm_sweep_circularity.R -- the control the F2 result needs.
#
# F2 found that adding log2(TMM factor) to the design drops the lithium DEG
# count from 1382 to 158, with negative controls (a random covariate, the same
# factor permuted, library size) leaving it untouched. The obvious objection is
# CIRCULARITY: the factor is computed FROM the expression matrix that supplies
# the outcome, so of course conditioning on it removes expression variance.
#
# A gene split settles it. Compute the factor from half the genes, then test the
# OTHER half. No gene in the test set contributed to its own covariate, so a
# collapse there cannot be per-gene circularity -- it has to be a real
# sample-level property that the factor indexes. Ten random splits, because one
# split is an anecdote.
#
# The same split is applied to expression PC1, which is the natural rival
# explanation ("it is just the dominant axis of the matrix").
#
# Second question, separate from circularity: TMM assumes most genes are not
# differentially expressed. If lithium really does shift the cell mixture, that
# assumption is FALSE here by construction, and TMM is not correcting an
# artefact -- it is deleting the effect of interest. Nothing internal to the
# data can adjudicate that without an absolute anchor, but a HOUSEKEEPING-gene
# normalisation is a third opinion with a different failure mode, so it is
# worth knowing which side it lands on.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})
source("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation/scripts/common.R")

OUTDIR <- file.path(EXP_HOME, "runs", "norm_sweep")
set.seed(20261218)
p <- load_prep()

d    <- build_design(p, "lithium")
cnt  <- p$counts[, d$samples, drop = FALSE]
gid  <- filter_genes(cnt, 10, 0.90)
cntf <- cnt[gid, , drop = FALSE]
meta <- d$meta; lith <- meta$lithium
Xcov <- d$X[, setdiff(colnames(d$X), "lithium"), drop = FALSE]
G    <- nrow(cntf)
gmap <- read.delim(file.path(EXP_HOME, "data", "gene_map.tsv"), header = FALSE,
                   col.names = c("gene", "symbol", "biotype"))
sym  <- setNames(gmap$symbol, gmap$gene)

cat("=========================================================\n")
cat("norm_sweep_circularity\n")
cat("=========================================================\n")

## ---------------------------------------------------------------------------
## helper: DEG count restricted to a gene subset, with BH recomputed inside it
## ---------------------------------------------------------------------------
# Recomputing BH within the subset matters. Carrying over q-values that were
# adjusted across all 12,173 genes would make the two halves depend on each
# other through the multiplicity correction, which is the very coupling this
# script is trying to break.
deg_in <- function(tt, genes, fdr = 0.05) {
  s <- tt[tt$gene %in% genes, ]
  sum(p.adjust(s$P.Value, "BH") < fdr)
}
run <- function(adjust = NULL) {
  f <- de_fit(p, "lithium", adjust = adjust, norm = "TMM", method = "voom")
  f$tt
}
base_tt <- run(NULL)
cat(sprintf("baseline: %d DEG over all %d genes\n", n_deg(base_tt), G))

## ===========================================================================
## S1 -- gene-split test of the norm-factor covariate
## ===========================================================================
cat("\n### S1  split-gene control (10 random half-splits) ###\n")
NSPLIT <- 10
rows <- list()
for (k in seq_len(NSPLIT)) {
  A <- sample(gid, floor(G / 2)); B <- setdiff(gid, A)
  # factor computed from half A ONLY
  nfA <- calcNormFactors(DGEList(cntf[A, , drop = FALSE]))$samples$norm.factors
  covA <- matrix(log2(nfA), ncol = 1, dimnames = list(meta$title, "nfA"))
  # expression PC1 from half A only, the rival explanation
  lA <- cpm(calcNormFactors(DGEList(cntf[A, , drop = FALSE])), log = TRUE, prior.count = 3)
  pcA <- matrix(prcomp(t(lA))$x[, 1], ncol = 1, dimnames = list(meta$title, "pcA"))

  tt_nf <- run(covA); tt_pc <- run(pcA)
  # does the half-A factor still track lithium? if not, nothing below means much
  m0 <- lm(log2(nfA) ~ Xcov - 1); m1 <- lm(log2(nfA) ~ lith + Xcov - 1)
  rows[[k]] <- data.frame(
    split = k,
    r_nfA_vs_nfFull = cor(log2(nfA), log2(calcNormFactors(DGEList(cntf))$samples$norm.factors)),
    p_nfA_lithium   = anova(m0, m1)$`Pr(>F)`[2],
    base_A = deg_in(base_tt, A), base_B = deg_in(base_tt, B),
    nfA_A  = deg_in(tt_nf, A),   nfA_B  = deg_in(tt_nf, B),
    pcA_A  = deg_in(tt_pc, A),   pcA_B  = deg_in(tt_pc, B),
    stringsAsFactors = FALSE)
  cat(sprintf("  split %2d  nfA~lithium p=%.2e | DEG in held-out half B: base %4d -> +nfA %3d (%.0f%%), +pcA %3d (%.0f%%)\n",
              k, rows[[k]]$p_nfA_lithium, rows[[k]]$base_B, rows[[k]]$nfA_B,
              100 * rows[[k]]$nfA_B / rows[[k]]$base_B, rows[[k]]$pcA_B,
              100 * rows[[k]]$pcA_B / rows[[k]]$base_B))
  flush(stdout())
}
S1 <- do.call(rbind, rows)
S1$pct_retained_B_nfA <- 100 * S1$nfA_B / S1$base_B
S1$pct_retained_A_nfA <- 100 * S1$nfA_A / S1$base_A
S1$pct_retained_B_pcA <- 100 * S1$pcA_B / S1$base_B
cat(sprintf("\nheld-out half B, mean DEG retained: +nfA %.1f%% (sd %.1f), +pcA %.1f%%\n",
            mean(S1$pct_retained_B_nfA), sd(S1$pct_retained_B_nfA), mean(S1$pct_retained_B_pcA)))
cat(sprintf("training half A, mean DEG retained: +nfA %.1f%%\n", mean(S1$pct_retained_A_nfA)))
cat("if held-out and training retention are similar, the collapse is NOT circularity\n")
write.csv(S1, file.path(OUTDIR, "S1_gene_split_control.csv"), row.names = FALSE)

## ===========================================================================
## S2 -- housekeeping normalisation: a third opinion
## ===========================================================================
cat("\n### S2  housekeeping-gene normalisation ###\n")
# A conventional panel, chosen by name before looking at any result rather than
# picked for stability in this dataset (which would be circular). The caveat is
# real and stated: "housekeeping" genes in whole blood are not cell-type
# neutral -- B2M and ACTB genuinely differ between lymphocytes and
# granulocytes -- so this is a third opinion, not an arbiter.
HK <- c("ACTB","GAPDH","B2M","PPIA","TBP","RPL13A","HPRT1","UBC","YWHAZ",
        "SDHA","RPLP0","PGK1","TFRC","GUSB","HMBS","POLR2A","PSMB4","VPS29")
hk_ids <- names(sym)[sym %in% HK]
hk_ids <- intersect(hk_ids, gid)
cat(sprintf("housekeeping genes present in filtered set: %d of %d panel members (%s)\n",
            length(hk_ids), length(HK), paste(sort(unique(sym[hk_ids])), collapse = ",")))

lib  <- colSums(cntf)
hk   <- colSums(cntf[hk_ids, , drop = FALSE])
# norm factor that makes effective library size proportional to housekeeping
# signal, centred so the geometric mean is 1 (edgeR's convention)
nf_hk <- (hk / lib); nf_hk <- nf_hk / exp(mean(log(nf_hk)))

X <- d$X
fit_with_nf <- function(nf, label) {
  dge <- DGEList(cntf); dge$samples$norm.factors <- nf
  v   <- voom(dge, X)
  fit <- eBayes(lmFit(v, X))
  tt  <- topTable(fit, coef = "lithium", number = Inf, sort.by = "none")
  tt$gene <- rownames(tt); tt
}
tt_hk  <- fit_with_nf(nf_hk, "housekeeping")
nf_tmm <- calcNormFactors(DGEList(cntf))$samples$norm.factors
ref_deg <- deg_ids(base_tt)
none_deg <- { O <- readRDS(file.path(OUTDIR, "norm_sweep_objects.rds"))
              t <- O$fits_tt$voom_none_CPM; t$gene[t$adj.P.Val < 0.05] }
jacc <- function(a, b) length(intersect(a, b)) / length(union(a, b))
hkd <- deg_ids(tt_hk)
m0 <- lm(log2(nf_hk) ~ Xcov - 1); m1 <- lm(log2(nf_hk) ~ lith + Xcov - 1)
S2 <- data.frame(
  quantity = c("DEG under housekeeping normalisation", "DEG under TMM", "DEG under none",
               "Jaccard housekeeping vs TMM", "Jaccard housekeeping vs none",
               "corr log2(hk factor) vs log2(TMM factor)",
               "lithium coef on log2(hk factor)", "p for that coef (adjusted)",
               "corr log2(hk factor) vs granulocyte fraction"),
  value = c(length(hkd), length(ref_deg), length(none_deg),
            jacc(hkd, ref_deg), jacc(hkd, none_deg),
            cor(log2(nf_hk), log2(nf_tmm)),
            unname(coef(m1)["lith"]), anova(m0, m1)$`Pr(>F)`[2],
            cor(log2(nf_hk), p$lineage[meta$title, "gran"])))
print(format(S2, digits = 4), row.names = FALSE)
write.csv(S2, file.path(OUTDIR, "S2_housekeeping_norm.csv"), row.names = FALSE)

## ===========================================================================
## S3 -- the spread of defensible answers
## ===========================================================================
# Every row here is a pipeline a competent analyst could write down and defend
# in a methods section. The range is the honest uncertainty attributable to
# normalisation and composition handling, before a single biological claim is
# made.
cat("\n### S3  the defensible range ###\n")
gran <- matrix(p$lineage[meta$title, "gran"], ncol = 1,
               dimnames = list(meta$title, "gran"))
S3 <- data.frame(
  pipeline = c("voom + TMM (project baseline)",
               "voom + RLE", "voom + upperquartile", "voom + TMM + quantile",
               "voom, no normalisation (log-CPM only)",
               "voom + housekeeping-gene normalisation",
               "voom + TMM + granulocyte fraction",
               "voom + TMM + 4 ILR balances",
               "voom + TMM + log2 TMM factor"),
  deg_fdr05 = c(length(ref_deg),
                { O <- readRDS(file.path(OUTDIR, "norm_sweep_objects.rds"))
                  c(sum(O$fits_tt$voom_RLE$adj.P.Val < 0.05),
                    sum(O$fits_tt$voom_upperquartile$adj.P.Val < 0.05)) },
                read.csv(file.path(OUTDIR, "F3_distributional_norm.csv"))$deg_fdr05[1],
                length(none_deg), length(hkd),
                read.csv(file.path(OUTDIR, "F2_covariate_specificity.csv"))$deg_fdr05[c(7, 9, 2)]),
  stringsAsFactors = FALSE)
S3$fold_vs_baseline <- round(S3$deg_fdr05 / length(ref_deg), 2)
print(S3, row.names = FALSE)
cat(sprintf("\nrange across defensible pipelines: %d to %d DEG (%.1f-fold)\n",
            min(S3$deg_fdr05), max(S3$deg_fdr05), max(S3$deg_fdr05) / min(S3$deg_fdr05)))
write.csv(S3, file.path(OUTDIR, "S3_defensible_range.csv"), row.names = FALSE)
cat("\ndone.\n")
