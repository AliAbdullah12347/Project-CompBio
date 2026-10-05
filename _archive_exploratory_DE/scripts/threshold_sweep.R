#!/usr/bin/env Rscript
# ==============================================================================
# threshold_sweep.R
#
# ONE canonical lithium fit. Nothing about the model, the filter, the
# normalisation or the covariates changes anywhere in this script. The only
# thing that varies is the rule used to call a gene significant.
#
# The question behind it: the baseline reports 1382 DEGs for lithium at
# FDR 0.05. A count is not a result. Is that 1382 genes with large expression
# changes, or 1382 genes whose changes are tiny but measured precisely enough
# (n = 226) to be distinguishable from zero? Those two worlds look identical in
# a DEG count and completely different in a volcano plot.
#
# Sections:
#   A  canonical fit + reproduction check against the stated baseline
#   B  FDR x |logFC| grid
#   C  Bonferroni, minimum adjusted p, significance floor
#   D  effect-size distribution, DEG vs non-DEG, in log2 and in fold change
#   E  within-fit yardstick: the sex coefficient from the SAME model
#   F  precision vs magnitude: standard errors, partial R^2, MDE
#   G  TREAT vs naive post-hoc |logFC| filtering
#   H  permutation null run through the identical grid
#
# Writes only into experimentation/runs/threshold_sweep/.
# ==============================================================================

source("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation/scripts/common.R")
set.seed(481)

OUT <- file.path(EXP_HOME, "runs", "threshold_sweep")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
wr <- function(df, f) write.csv(df, file.path(OUT, f), row.names = FALSE)
hr <- function(s) cat(sprintf("\n===== %s =====\n", s))

p <- load_prep()

## ============================================================== A  canonical
hr("A  canonical lithium fit")

fit0 <- de_fit(p, "lithium")            # every argument at its default
tt   <- fit0$tt
cat(sprintf("genes %d   samples %d   coef %s   method %s/%s\n",
            fit0$n_genes, fit0$n_samples, fit0$coef, fit0$method, fit0$norm))
cat(sprintf("DEG @ FDR 0.05            %d\n", n_deg(tt, 0.05)))
cat(sprintf("pi0 (Storey)              %.4f\n", pi0_storey(tt$P.Value)))

# The baseline the experiment is anchored to. If this does not reproduce,
# nothing below is interpretable, so stop rather than report a sweep of the
# wrong fit.
stopifnot(fit0$n_genes == 12173, fit0$n_samples == 226, n_deg(tt, 0.05) == 1382)
cat("reproduction check        PASS (12173 genes / 226 samples / 1382 DEG)\n")

## ================================================================= B  grid
hr("B  FDR x |logFC| grid")

FDRS <- c(0.001, 0.01, 0.05, 0.10, 0.20)
LFCS <- c(0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0)

grid <- do.call(rbind, lapply(FDRS, function(f) do.call(rbind, lapply(LFCS, function(l) {
  k <- n_deg(tt, f, l)
  data.frame(fdr = f, lfc = l, n_deg = k,
             pct_of_tested = 100 * k / nrow(tt),
             # what fraction of the genes that pass THIS fdr also clear THIS lfc
             pct_of_fdr_set = 100 * k / max(1, n_deg(tt, f, 0)))
}))))
wr(grid, "grid_fdr_x_lfc.csv")

# printed as the matrix a reader actually wants to look at
mat <- matrix(grid$n_deg, nrow = length(FDRS), byrow = TRUE,
              dimnames = list(paste0("FDR<", FDRS), paste0("|lfc|>=", LFCS)))
print(mat)

## ========================================== C  Bonferroni and the floor
hr("C  Bonferroni / extremes / significance floor")

bonf <- p.adjust(tt$P.Value, "bonferroni")
sig  <- tt$adj.P.Val < 0.05

# The largest raw p-value that still clears FDR 0.05 -- the operative alpha.
p_cut <- if (any(sig)) max(tt$P.Value[sig]) else NA_real_

extremes <- data.frame(
  metric = c("n_tested", "min_raw_p", "min_adj_p_BH", "min_adj_p_bonferroni",
             "n_bonferroni_0.05", "n_bonferroni_0.01", "n_BH_0.05",
             "operative_raw_p_cutoff_at_BH0.05",
             "min_absLFC_among_BH0.05", "max_absLFC_among_BH0.05",
             "max_absLFC_all_genes", "pi0_storey_lambda0.5"),
  value = c(nrow(tt), min(tt$P.Value), min(tt$adj.P.Val), min(bonf),
            sum(bonf < 0.05), sum(bonf < 0.01), sum(sig),
            p_cut,
            min(abs(tt$logFC[sig])), max(abs(tt$logFC[sig])),
            max(abs(tt$logFC)), pi0_storey(tt$P.Value)))
print(extremes, digits = 6)
wr(extremes, "extremes.csv")

# Storey's pi0 implies a number of genuinely non-null genes. Compare that to
# how many the FDR rule actually returns: the gap is the power deficit.
m_tot <- nrow(tt)
cat(sprintf("\nimplied non-null genes (1-pi0)*m = %.0f ; BH 0.05 returns %d (%.0f%% recovered)\n",
            (1 - pi0_storey(tt$P.Value)) * m_tot, sum(sig),
            100 * sum(sig) / ((1 - pi0_storey(tt$P.Value)) * m_tot)))

## ====================================================== D  effect sizes
hr("D  effect-size distribution, DEG vs non-DEG")

qs <- c(0, .05, .25, .50, .75, .95, 1)
dist_tbl <- rbind(
  data.frame(set = "DEG_BH0.05",  n = sum(sig),
             t(quantile(abs(tt$logFC[sig]), qs)),  mean_abs = mean(abs(tt$logFC[sig])),
             check.names = FALSE),
  data.frame(set = "nonDEG",      n = sum(!sig),
             t(quantile(abs(tt$logFC[!sig]), qs)), mean_abs = mean(abs(tt$logFC[!sig])),
             check.names = FALSE),
  data.frame(set = "all",         n = nrow(tt),
             t(quantile(abs(tt$logFC), qs)),       mean_abs = mean(abs(tt$logFC)),
             check.names = FALSE))
names(dist_tbl)[3:9] <- paste0("q", c("00","05","25","50","75","95","100"))
print(dist_tbl, digits = 4)
wr(dist_tbl, "logfc_distribution.csv")

# Same numbers as multiplicative fold change, which is the scale a biologist
# reads. 2^0.1 = 1.072 -- a 7% change.
fc_tbl <- data.frame(
  quantity = c("median |logFC| of DEGs", "median fold change of DEGs",
               "median % change of DEGs", "90th pct |logFC| of DEGs",
               "90th pct fold change of DEGs",
               "median |logFC| of non-DEGs", "max |logFC| of DEGs",
               "max fold change of DEGs"),
  value = c(median(abs(tt$logFC[sig])), 2^median(abs(tt$logFC[sig])),
            100 * (2^median(abs(tt$logFC[sig])) - 1),
            quantile(abs(tt$logFC[sig]), .90), 2^quantile(abs(tt$logFC[sig]), .90),
            median(abs(tt$logFC[!sig])), max(abs(tt$logFC[sig])),
            2^max(abs(tt$logFC[sig]))))
print(fc_tbl, digits = 5)
wr(fc_tbl, "fold_change_summary.csv")

# Conventional biological filters, applied post hoc to the FDR 0.05 set.
surv <- data.frame(
  rule = c("BH0.05 only", "BH0.05 & |lfc|>0.10", "BH0.05 & |lfc|>0.20",
           "BH0.05 & |lfc|>0.32 (1.25x)", "BH0.05 & |lfc|>0.50",
           "BH0.05 & |lfc|>0.58 (1.5x)", "BH0.05 & |lfc|>1.00 (2x)"),
  lfc  = c(0, 0.10, 0.20, log2(1.25), 0.50, log2(1.5), 1.00))
surv$n <- sapply(surv$lfc, function(l) sum(sig & abs(tt$logFC) > l))
surv$pct_of_1382 <- 100 * surv$n / sum(sig)
print(surv, digits = 4)
wr(surv, "conventional_lfc_survivors.csv")

# Are the strongest p-values also the strongest effects? If the two rankings
# disagree, "top DEG" means different genes depending on which column is sorted.
cat(sprintf("\nSpearman rho( -log10 p , |logFC| ) all genes = %.3f | DEGs only = %.3f\n",
            cor(-log10(tt$P.Value), abs(tt$logFC), method = "spearman"),
            cor(-log10(tt$P.Value[sig]), abs(tt$logFC[sig]), method = "spearman")))
top_p   <- tt$gene[order(tt$P.Value)][1:100]
top_lfc <- tt$gene[order(-abs(tt$logFC))][1:100]
cat(sprintf("overlap of top-100 by p and top-100 by |logFC| = %d genes\n",
            length(intersect(top_p, top_lfc))))

# Does |logFC| track expression level? Large apparent effects concentrated in
# weakly expressed genes are the classic low-count artefact.
cat(sprintf("Spearman rho( AveExpr , |logFC| ) = %.3f  (all genes)\n",
            cor(tt$AveExpr, abs(tt$logFC), method = "spearman")))
ab <- cut(tt$AveExpr, quantile(tt$AveExpr, 0:5/5), include.lowest = TRUE,
          labels = paste0("expr_Q", 1:5))
expr_tbl <- data.frame(
  bin = levels(ab),
  n = as.integer(table(ab)),
  median_AveExpr = tapply(tt$AveExpr, ab, median),
  median_absLFC  = tapply(abs(tt$logFC), ab, median),
  n_DEG          = tapply(sig, ab, sum),
  pct_DEG        = 100 * tapply(sig, ab, mean),
  n_DEG_lfc0.5   = tapply(sig & abs(tt$logFC) > 0.5, ab, sum))
print(expr_tbl, digits = 4, row.names = FALSE)
wr(expr_tbl, "effect_by_expression_bin.csv")

## ================================================ E  within-fit yardstick
hr("E  sex coefficient from the SAME fit -- an effect-size yardstick")

# Sex is already a covariate in the canonical design. Pulling its coefficient
# out of the identical model gives a reference effect measured on the same
# genes, same samples, same normalisation, same residual variance. Anything
# about sequencing depth or the filter that inflates or deflates lithium's
# logFC does the same to sex's. Sex-chromosome genes are a known large,
# unambiguous biological contrast, so this calibrates the scale.
sexcoef <- grep("^sex", colnames(fit0$design$X), value = TRUE)[1]
fit_sex <- de_fit(p, "lithium", coef_override = sexcoef)
tts <- fit_sex$tt
stopifnot(identical(tts$gene, tt$gene))   # same genes, same fit, other column

yard <- data.frame(
  coefficient = c("lithium", sexcoef),
  n_BH0.05      = c(sum(sig), sum(tts$adj.P.Val < 0.05)),
  n_BH0.05_lfc0.5 = c(sum(sig & abs(tt$logFC) > 0.5),
                      sum(tts$adj.P.Val < 0.05 & abs(tts$logFC) > 0.5)),
  median_absLFC_of_sig = c(median(abs(tt$logFC[sig])),
                           median(abs(tts$logFC[tts$adj.P.Val < 0.05]))),
  max_absLFC   = c(max(abs(tt$logFC)), max(abs(tts$logFC))),
  min_adj_p    = c(min(tt$adj.P.Val), min(tts$adj.P.Val)))
print(yard, digits = 4)
wr(yard, "yardstick_sex_vs_lithium.csv")

## =========================================== F  precision vs magnitude
hr("F  precision vs magnitude")

# Under limma, adj.P is driven by logFC / SE. Recovering SE from t lets the
# question "is 1382 a magnitude result or a precision result?" be answered
# numerically instead of asserted.
se  <- abs(tt$logFC) / abs(tt$t)
dfr <- fit0$n_samples - ncol(fit0$design$X)
# squared partial correlation of the exposure with each gene, the share of
# residual variance lithium accounts for once covariates are removed
pr2 <- tt$t^2 / (tt$t^2 + dfr)

prec <- data.frame(
  metric = c("residual df", "median SE of logFC", "median SE among DEGs",
             "median SE among non-DEGs",
             "median partial R^2 all genes", "median partial R^2 DEGs",
             "max partial R^2", "90th pct partial R^2 DEGs",
             "smallest |logFC| called significant",
             "MDE at 80% power, BH-operative alpha, median SE",
             "MDE as fold change"),
  value = c(dfr, median(se), median(se[sig]), median(se[!sig]),
            median(pr2), median(pr2[sig]), max(pr2), quantile(pr2[sig], .90),
            min(abs(tt$logFC[sig])),
            NA, NA))
# minimum detectable effect: the |logFC| a gene of median precision would need
# for 80% power at the alpha the BH rule is actually operating at
tcrit <- qt(1 - p_cut / 2, dfr)
mde   <- (tcrit + qt(0.80, dfr)) * median(se)
prec$value[prec$metric == "MDE at 80% power, BH-operative alpha, median SE"] <- mde
prec$value[prec$metric == "MDE as fold change"] <- 2^mde
print(prec, digits = 5)
wr(prec, "precision_vs_magnitude.csv")

## ================================================== G  TREAT vs naive
hr("G  TREAT vs naive post-hoc |logFC| filtering")

# TREAT tests H0: |true logFC| <= tau. Naive filtering tests H0: logFC = 0 and
# then discards genes whose POINT ESTIMATE is small. The second has no error
# control at the threshold: a gene with logFC 0.51 and a wide interval passes
# the naive filter while TREAT correctly refuses it. They are different
# hypotheses, so they are allowed to disagree; the size of the disagreement is
# the result.
#
# The design and the count matrix come from the engine (build_design +
# filter_genes + calcNormFactors + voom), so this is the canonical fit with
# eBayes swapped for treat and nothing else. Equivalence to de_fit is asserted.
d  <- build_design(p, "lithium")
cnt <- p$counts[, d$samples, drop = FALSE]
cnt <- cnt[filter_genes(cnt, 10, 0.90), , drop = FALSE]
dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = cnt), method = "TMM")
v   <- limma::voom(dge, d$X)
lf  <- limma::lmFit(v, d$X)

chk <- limma::topTable(limma::eBayes(lf), coef = d$coef, number = Inf, sort.by = "none")
stopifnot(nrow(chk) == nrow(tt), identical(rownames(chk), tt$gene),
          max(abs(chk$logFC - tt$logFC)) < 1e-10,
          max(abs(chk$adj.P.Val - tt$adj.P.Val)) < 1e-10)
cat("local refit is byte-identical to de_fit canonical  PASS\n")

TAUS <- c(0.1, 0.2, 0.5)
treat_rows <- do.call(rbind, lapply(TAUS, function(tau) {
  tr <- limma::treat(lf, lfc = tau)
  tp <- limma::topTreat(tr, coef = d$coef, number = Inf, sort.by = "none")
  tg <- rownames(tp)[tp$adj.P.Val < 0.05]
  ng <- tt$gene[sig & abs(tt$logFC) > tau]      # naive: BH0.05 then discard
  data.frame(tau = tau,
             n_treat_BH0.05 = length(tg),
             n_naive_BH0.05_then_lfc = length(ng),
             n_both = length(intersect(tg, ng)),
             n_treat_only = length(setdiff(tg, ng)),
             n_naive_only = length(setdiff(ng, tg)),
             jaccard = length(intersect(tg, ng)) / max(1, length(union(tg, ng))),
             treat_min_absLFC = if (length(tg)) min(abs(tp$logFC[match(tg, rownames(tp))])) else NA,
             naive_min_absLFC = if (length(ng)) min(abs(tt$logFC[match(ng, tt$gene)])) else NA)
}))
print(treat_rows, digits = 4)
wr(treat_rows, "treat_vs_naive.csv")

# TREAT at several FDRs as well, so the comparison is not anchored to one alpha
treat_grid <- do.call(rbind, lapply(TAUS, function(tau) {
  tp <- limma::topTreat(limma::treat(lf, lfc = tau), coef = d$coef,
                        number = Inf, sort.by = "none")
  do.call(rbind, lapply(FDRS, function(f)
    data.frame(tau = tau, fdr = f,
               n_treat = sum(tp$adj.P.Val < f),
               n_naive = sum(tt$adj.P.Val < f & abs(tt$logFC) > tau))))
}))
wr(treat_grid, "treat_grid.csv")
print(treat_grid)

## ================================================== H  permutation null
hr("H  permutation null through the identical grid")

# The grid above says how many genes a rule returns. It does not say how many a
# rule returns when there is nothing to find. Permuting the lithium column
# only -- every covariate stays attached to its own sample -- gives the null
# for "lithium carries no signal" rather than the null of exchangeable samples.
NPERM <- 20
n_i <- fit0$n_samples
perm <- do.call(rbind, lapply(seq_len(NPERM), function(i) {
  pi_ <- sample.int(n_i)
  tp  <- de_fit(p, "lithium", permute_exposure = pi_)$tt
  cbind(perm = i,
        do.call(rbind, lapply(FDRS, function(f) do.call(rbind, lapply(LFCS, function(l)
          data.frame(fdr = f, lfc = l, n_deg = n_deg(tp, f, l)))))),
        pi0 = pi0_storey(tp$P.Value), min_adj_p = min(tp$adj.P.Val),
        n_bonf = sum(p.adjust(tp$P.Value, "bonferroni") < 0.05))
}))
wr(perm, "permutation_null_raw.csv")

null_summ <- aggregate(n_deg ~ fdr + lfc, perm, function(x)
  c(mean = mean(x), max = max(x), median = median(x)))
null_summ <- data.frame(null_summ[, 1:2], as.data.frame(null_summ$n_deg))
names(null_summ)[3:5] <- c("null_mean", "null_max", "null_median")
obs <- grid[, c("fdr", "lfc", "n_deg")]
cmp <- merge(obs, null_summ, by = c("fdr", "lfc"))
cmp$enrichment <- cmp$n_deg / pmax(cmp$null_mean, 0.5)   # 0.5 floor: no /0
cmp <- cmp[order(cmp$fdr, cmp$lfc), ]
print(cmp, digits = 4, row.names = FALSE)
wr(cmp, "observed_vs_null.csv")

perm1 <- perm[!duplicated(perm$perm), ]
cat(sprintf("\nnull pi0: median %.3f (range %.3f-%.3f)   observed pi0 %.3f\n",
            median(perm1$pi0), min(perm1$pi0), max(perm1$pi0), pi0_storey(tt$P.Value)))
cat(sprintf("null Bonferroni hits: median %g, max %g   observed %d\n",
            median(perm1$n_bonf), max(perm1$n_bonf), sum(bonf < 0.05)))

## =================================================== main summary table
hr("summary")

summ <- grid
summ$null_mean <- cmp$null_mean[match(paste(summ$fdr, summ$lfc),
                                      paste(cmp$fdr, cmp$lfc))]
summ$null_max  <- cmp$null_max[match(paste(summ$fdr, summ$lfc),
                                     paste(cmp$fdr, cmp$lfc))]
summ$treat_n <- NA_integer_
for (i in seq_len(nrow(summ))) {
  hit <- treat_grid$tau == summ$lfc[i] & treat_grid$fdr == summ$fdr[i]
  if (any(hit)) summ$treat_n[i] <- treat_grid$n_treat[hit][1]
}
write_run("threshold_sweep", summ, NULL)   # NOTES.md written separately

saveRDS(list(tt = tt, tt_sex = tts, grid = grid, cmp = cmp,
             treat_rows = treat_rows, treat_grid = treat_grid,
             extremes = extremes, se = se, pr2 = pr2),
        file.path(OUT, "objects.rds"))

cat("\n---- HEADLINE ----\n")
cat(sprintf("1382 DEGs at FDR 0.05; median |logFC| = %.4f (fold change %.3f, %.1f%%)\n",
            median(abs(tt$logFC[sig])), 2^median(abs(tt$logFC[sig])),
            100 * (2^median(abs(tt$logFC[sig])) - 1)))
cat(sprintf("survivors at |logFC|>0.5: %d (%.2f%% of the 1382)\n",
            sum(sig & abs(tt$logFC) > 0.5), 100 * sum(sig & abs(tt$logFC) > 0.5) / sum(sig)))
cat(sprintf("TREAT tau=0.1 at FDR 0.05: %d genes; naive BH0.05 & |lfc|>0.1: %d\n",
            treat_rows$n_treat_BH0.05[treat_rows$tau == 0.1],
            treat_rows$n_naive_BH0.05_then_lfc[treat_rows$tau == 0.1]))
cat("done.\n")
