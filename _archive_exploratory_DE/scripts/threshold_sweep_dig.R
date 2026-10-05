#!/usr/bin/env Rscript
# ==============================================================================
# threshold_sweep_dig.R -- follow-ups to threshold_sweep.R
#
# Three results in the main sweep were surprising enough to be worth chasing,
# and one claim in it needs a stronger form of evidence:
#
#   I    Not one of the 12,173 genes reaches |logFC| > 1. Is that biology, or
#        did the abundance filter remove every gene capable of showing a large
#        effect before the model was ever fitted?
#   II   The permutation null is zero in 19 of 20 replicates at FDR 0.20 and
#        1,755 in the twentieth. A null that bimodal is not a null you can
#        summarise with a mean.
#   III  DEG rate climbs monotonically with expression (7.4% -> 18.4%) while
#        median |logFC| stays flat. That is a precision gradient, not a
#        magnitude gradient -- shown here with standard errors per bin.
#   IV   TREAT is the formal answer to "how big is this effect". The intuitive
#        answer is a confidence interval. Both are computed so they can be
#        checked against each other.
#
# Everything reuses the canonical fit from threshold_sweep.R (objects.rds).
# Section I deliberately steps outside the canonical filter -- it is labelled a
# diagnostic, its numbers never replace a canonical number, and it exists only
# to say whether the canonical ceiling is real or imposed.
#
# Writes only into experimentation/runs/threshold_sweep/.
# ==============================================================================

source("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation/scripts/common.R")
set.seed(481)

OUT <- file.path(EXP_HOME, "runs", "threshold_sweep")
wr <- function(df, f) write.csv(df, file.path(OUT, f), row.names = FALSE)
hr <- function(s) cat(sprintf("\n===== %s =====\n", s))

p   <- load_prep()
OBJ <- readRDS(file.path(OUT, "objects.rds"))
tt  <- OBJ$tt
sig <- tt$adj.P.Val < 0.05

## ================================ I  is the |logFC| ceiling real or imposed?
hr("I  DIAGNOSTIC -- does the abundance filter cap the observable effect size?")

# The canonical filter is "more than 10 counts in at least 90% of samples".
# A gene that is genuinely OFF in one group and ON in the other fails that test
# by construction: it has <=10 counts in the OFF group, and if the OFF group is
# more than 10% of the cohort the gene is deleted. So the filter systematically
# removes the class of gene that produces the largest log fold changes. Whether
# that matters here is an empirical question, and the sex coefficient answers
# it: sex-chromosome genes are the textbook on/off case and sex is already a
# covariate in the canonical design, so its coefficient is available for free.

SEXGENES <- c(XIST = "ENSG00000229807", RPS4Y1 = "ENSG00000129824",
              DDX3Y = "ENSG00000067048", KDM5D = "ENSG00000012817",
              UTY   = "ENSG00000183878", USP9Y = "ENSG00000114374",
              EIF1AY= "ENSG00000198692", NLGN4Y= "ENSG00000165246",
              TXLNGY= "ENSG00000131002")
present <- SEXGENES %in% rownames(p$counts)
cat(sprintf("sex-marker genes present in the count matrix: %d/%d\n",
            sum(present), length(SEXGENES)))
cat(sprintf("of those, surviving the CANONICAL filter:     %d\n",
            sum(SEXGENES %in% tt$gene)))
print(data.frame(symbol = names(SEXGENES), ensembl = unname(SEXGENES),
                 in_matrix = present,
                 passes_canonical_filter = SEXGENES %in% tt$gene), row.names = FALSE)

# A ladder of filters, from the canonical one down to almost nothing. For each,
# the largest effect the model is capable of reporting, for lithium and for sex.
LADDER <- list(
  list(lab = "canonical  >10 in >=90%",  mc = 10, mp = 0.90),
  list(lab = ">10 in >=50%",             mc = 10, mp = 0.50),
  list(lab = ">0  in >=50%",             mc =  0, mp = 0.50),
  list(lab = ">0  in >=20%",             mc =  0, mp = 0.20))

sexcoef <- grep("^sex", colnames(de_fit(p, "lithium")$design$X), value = TRUE)[1]

lad <- do.call(rbind, lapply(LADDER, function(s) {
  fl <- de_fit(p, "lithium", min_count = s$mc, min_prop = s$mp)
  fs <- de_fit(p, "lithium", min_count = s$mc, min_prop = s$mp,
               coef_override = sexcoef)
  a <- fl$tt; b <- fs$tt
  sx <- SEXGENES[SEXGENES %in% b$gene]
  data.frame(filter = s$lab, n_genes = fl$n_genes,
             li_nDEG          = n_deg(a, 0.05),
             li_max_absLFC    = max(abs(a$logFC)),
             li_n_lfc_gt_0.5  = n_deg(a, 0.05, 0.5),
             li_n_lfc_gt_1    = n_deg(a, 0.05, 1.0),
             sex_nDEG         = n_deg(b, 0.05),
             sex_max_absLFC   = max(abs(b$logFC)),
             sex_n_lfc_gt_1   = n_deg(b, 0.05, 1.0),
             n_sexmarkers_kept = length(sx),
             sexmarker_max_absLFC = if (length(sx))
               max(abs(b$logFC[match(sx, b$gene)])) else NA_real_)
}))
print(lad, digits = 4, row.names = FALSE)
wr(lad, "diag_filter_effect_ceiling.csv")

## =========================== II  the bimodal permutation null
hr("II  why one permutation in twenty behaves nothing like the other nineteen")

# Hypothesis: the outlier permutation is simply one that landed close to the
# true assignment. Testing that with MORE random permutations is the wrong
# instrument -- a random permutation preserves the 74/152 split, so its
# agreement with the truth concentrates near 126/226 and rarely strays far.
# Random draws sample the middle of the axis densely and the interesting ends
# almost never.
#
# So instead of resampling, the axis is walked deliberately. Starting from the
# identity permutation (agreement 226/226, i.e. the truth) and applying r
# random 0<->1 transpositions gives a permutation with agreement exactly
# 226-2r. That constructs a dose-response: DEG count as a function of how much
# of the real exposure survives. A random permutation sits at r ~ 49.8 on this
# ladder (expected agreement 74*74/226 + 152*152/226 = 126.4), so the ladder
# brackets the random null on both sides.
#
# The 20 random permutations from the main script are reused for the null
# distribution itself rather than recomputed -- same seed, same engine.
d0  <- build_design(p, "lithium")
lit <- d0$X[, d0$coef]
n_i <- d0$n
i1  <- which(lit == 1); i0 <- which(lit == 0)
stopifnot(length(i1) == 152, length(i0) == 74)

swap_perm <- function(r) {          # r transpositions -> agreement n - 2r
  pm <- seq_len(n_i)
  a <- sample(i1, r); b <- sample(i0, r)
  pm[a] <- b; pm[b] <- a
  pm
}
RS <- c(74, 60, 50, 45, 40, 30, 20, 10, 5)
ladder2 <- do.call(rbind, lapply(RS, function(r) {
  pm <- swap_perm(r)
  tp <- de_fit(p, "lithium", permute_exposure = pm)$tt
  data.frame(n_swaps = r,
             agreement = sum(lit[pm] == lit),
             pct_agree = 100 * sum(lit[pm] == lit) / n_i,
             cor_with_truth = cor(lit[pm], lit),
             pi0 = pi0_storey(tp$P.Value),
             deg_005 = n_deg(tp, 0.05), deg_010 = n_deg(tp, 0.10),
             deg_020 = n_deg(tp, 0.20),
             min_adj_p = min(tp$adj.P.Val),
             n_bonf = sum(p.adjust(tp$P.Value, "bonferroni") < 0.05))
}))
cat("dose-response: DEG count vs how much of the true exposure survives\n")
cat(sprintf("(a random permutation sits at agreement ~%.1f, i.e. r ~ %.1f)\n",
            length(i0)^2/n_i + length(i1)^2/n_i,
            (n_i - (length(i0)^2/n_i + length(i1)^2/n_i))/2))
print(ladder2, digits = 4, row.names = FALSE)
wr(ladder2, "diag_null_doseresponse.csv")
cat(sprintf("\nSpearman rho( agreement , deg_020 ) across the ladder = %.3f\n",
            cor(ladder2$agreement, ladder2$deg_020, method = "spearman")))

# null distribution: the 20 random permutations already computed and saved
raw <- read.csv(file.path(OUT, "permutation_null_raw.csv"))
w <- function(f) raw$n_deg[raw$fdr == f & raw$lfc == 0]
nullq <- data.frame(
  fdr    = c(0.05, 0.10, 0.20),
  mean   = sapply(c(0.05,0.10,0.20), function(f) mean(w(f))),
  median = sapply(c(0.05,0.10,0.20), function(f) median(w(f))),
  q90    = sapply(c(0.05,0.10,0.20), function(f) unname(quantile(w(f), .90))),
  q95    = sapply(c(0.05,0.10,0.20), function(f) unname(quantile(w(f), .95))),
  max    = sapply(c(0.05,0.10,0.20), function(f) max(w(f))),
  n_perm_zero = sapply(c(0.05,0.10,0.20), function(f) sum(w(f) == 0)),
  n_perm      = sapply(c(0.05,0.10,0.20), function(f) length(w(f))),
  observed = c(n_deg(tt, 0.05), n_deg(tt, 0.10), n_deg(tt, 0.20)))
print(nullq, digits = 4, row.names = FALSE)
wr(nullq, "diag_null_quantiles.csv")
cat(sprintf("\npermutation p (deg_005 >= observed): (%d+1)/(%d+1) = %.4f\n",
            sum(w(0.05) >= n_deg(tt, 0.05)), length(w(0.05)),
            (sum(w(0.05) >= n_deg(tt, 0.05)) + 1) / (length(w(0.05)) + 1)))

## ================= III  precision gradient, not magnitude gradient
hr("III  is the expression-level DEG gradient precision or magnitude?")

se  <- abs(tt$logFC) / abs(tt$t)
ab  <- cut(tt$AveExpr, quantile(tt$AveExpr, 0:5/5), include.lowest = TRUE,
           labels = paste0("expr_Q", 1:5))
g <- data.frame(
  bin            = levels(ab),
  n              = as.integer(table(ab)),
  median_AveExpr = tapply(tt$AveExpr, ab, median),
  median_absLFC  = tapply(abs(tt$logFC), ab, median),
  median_SE      = tapply(se, ab, median),
  median_absT    = tapply(abs(tt$t), ab, median),
  pct_DEG        = 100 * tapply(sig, ab, mean),
  # what |logFC| a gene in this bin needs for 80% power at the operative alpha
  MDE_80         = NA_real_)
p_cut <- max(tt$P.Value[sig])
dfr   <- d0$n - ncol(d0$X)          # computed, not copied from the other script
g$MDE_80 <- (qt(1 - p_cut/2, dfr) + qt(0.80, dfr)) * g$median_SE
print(g, digits = 4, row.names = FALSE)
wr(g, "diag_precision_gradient.csv")
cat(sprintf("\nSE falls %.2fx from Q1 to Q5 while median |logFC| changes %.2fx\n",
            g$median_SE[1] / g$median_SE[5], g$median_absLFC[1] / g$median_absLFC[5]))

## ============================ IV  confidence intervals on the effect
hr("IV  95% CI on logFC -- the interval form of the TREAT question")

d  <- build_design(p, "lithium")
cn <- p$counts[, d$samples, drop = FALSE]
cn <- cn[filter_genes(cn, 10, 0.90), , drop = FALSE]
dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = cn), method = "TMM")
fit <- limma::eBayes(limma::lmFit(limma::voom(dge, d$X), d$X))
ci  <- limma::topTable(fit, coef = d$coef, number = Inf, sort.by = "none",
                       confint = 0.95)
stopifnot(identical(rownames(ci), tt$gene), max(abs(ci$logFC - tt$logFC)) < 1e-10)

# lower bound of |effect|: 0 if the interval straddles zero, else the nearer end
lo <- ifelse(ci$CI.L > 0, ci$CI.L, ifelse(ci$CI.R < 0, -ci$CI.R, 0))
TAUS <- c(0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0)
cit <- data.frame(
  tau = TAUS,
  n_CIlower_above_tau      = sapply(TAUS, function(x) sum(lo > x)),
  n_CIlower_above_tau_sig  = sapply(TAUS, function(x) sum(lo > x & sig)),
  n_pointest_above_tau_sig = sapply(TAUS, function(x) sum(abs(tt$logFC) > x & sig)))
cit$ratio_point_to_CI <- with(cit, n_pointest_above_tau_sig /
                                pmax(1, n_CIlower_above_tau_sig))
print(cit, digits = 4, row.names = FALSE)
wr(cit, "diag_confint.csv")

cat(sprintf("\nmedian CI width among DEGs: %.3f log2 units (%.2fx multiplicative)\n",
            median((ci$CI.R - ci$CI.L)[sig]), 2^median((ci$CI.R - ci$CI.L)[sig])))
cat(sprintf("widest DEG point estimate %.3f, its CI [%.3f, %.3f]\n",
            tt$logFC[which.max(abs(tt$logFC))],
            ci$CI.L[which.max(abs(tt$logFC))], ci$CI.R[which.max(abs(tt$logFC))]))

# top 15 by |logFC| with intervals -- the genes anyone would quote
top <- order(-abs(tt$logFC))[1:15]
tops <- data.frame(gene = tt$gene[top], logFC = tt$logFC[top],
                   CI.L = ci$CI.L[top], CI.R = ci$CI.R[top],
                   AveExpr = tt$AveExpr[top], adj.P = tt$adj.P.Val[top],
                   foldchange = 2^abs(tt$logFC[top]))
print(tops, digits = 4, row.names = FALSE)
wr(tops, "diag_top15_by_effect.csv")

## ================================================ V  direction balance
hr("V  direction of the significant effects")

# If lithium's transcriptome signal were mostly a shift in cell mixture -- more
# neutrophils, fewer lymphocytes -- the DEG set would be strongly polarised,
# because a mixture shift moves whole marker programmes in one direction. A
# roughly 50/50 up/down split argues against a single dominant mixture axis and
# is worth recording even though this script cannot test it directly.
dirtab <- do.call(rbind, lapply(c(0, 0.1, 0.2, 0.5), function(l) {
  s <- sig & abs(tt$logFC) > l
  data.frame(lfc_cut = l, n = sum(s),
             n_up = sum(s & tt$logFC > 0), n_down = sum(s & tt$logFC < 0),
             pct_up = 100 * sum(s & tt$logFC > 0) / max(1, sum(s)),
             binom_p = if (sum(s) > 0)
               binom.test(sum(s & tt$logFC > 0), sum(s), 0.5)$p.value else NA)
}))
print(dirtab, digits = 4, row.names = FALSE)
wr(dirtab, "diag_direction.csv")

cat("\ndig complete.\n")
