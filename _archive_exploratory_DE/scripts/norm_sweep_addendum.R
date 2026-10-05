#!/usr/bin/env Rscript
# ==============================================================================
# norm_sweep_addendum.R -- three loose ends from the sweep, made explicit.
#
# A1. F2 compared the TMM factor against "expression PC1" as a rival
#     explanation and the factor won by a wide margin (158 vs 470 DEG). That
#     comparison was UNFAIR and it is worth saying so in code rather than in a
#     footnote: PC1 there was computed on TMM-NORMALISED log-CPM, so the scale
#     axis had already been removed from it before it was asked to compete for
#     the scale axis. The fair rival is PC1 of UN-normalised log-CPM.
#
# A2. Which covariates change the factor/lithium association, and whether they
#     do it by removing confounding (beta moves) or by removing noise (SE
#     moves). These are different claims and the distinction is invisible in a
#     p-value.
#
# A3. The relationship between each method's global logFC offset and the
#     directional composition of its gene list, saved rather than only printed.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})
source("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation/scripts/common.R")

OUTDIR <- file.path(EXP_HOME, "runs", "norm_sweep")
set.seed(20261218)
p   <- load_prep()
OBJ <- readRDS(file.path(OUTDIR, "norm_sweep_objects.rds"))
NF  <- OBJ$NF

d    <- build_design(p, "lithium")
cnt  <- p$counts[, d$samples, drop = FALSE]
cntf <- cnt[filter_genes(cnt, 10, 0.90), , drop = FALSE]
meta <- d$meta; lith <- meta$lithium
Xcov <- d$X[, setdiff(colnames(d$X), "lithium"), drop = FALSE]
ref_deg <- { t <- OBJ$fits_tt$voom_TMM; t$gene[t$adj.P.Val < 0.05] }
jacc <- function(a, b) length(intersect(a, b)) / length(union(a, b))

## ---------------------------------------------------------------- A1 -------
cat("### A1  the fair rival: PC1 of UN-normalised log-CPM ###\n")
lc_none <- cpm(DGEList(cntf), log = TRUE, prior.count = 3)                 # no TMM
lc_tmm  <- cpm(calcNormFactors(DGEList(cntf)), log = TRUE, prior.count = 3) # with TMM
pc_none <- prcomp(t(lc_none))$x[, 1:2]
pc_tmm  <- prcomp(t(lc_tmm))$x[, 1:2]
y <- log2(NF[, "TMM"])
cat(sprintf("|r| between log2 TMM factor and PC1 of un-normalised log-CPM: %.3f\n",
            abs(cor(y, pc_none[, 1]))))
cat(sprintf("|r| between log2 TMM factor and PC1 of TMM-normalised log-CPM: %.3f\n",
            abs(cor(y, pc_tmm[, 1]))))

mk <- function(v, nm) matrix(v, ncol = 1, dimnames = list(meta$title, nm))
CAND <- list("none (baseline)"               = NULL,
             "log2 TMM factor"               = mk(y, "nf"),
             "PC1 of UN-normalised log-CPM"  = mk(pc_none[, 1], "pc1u"),
             "PC1 of TMM-normalised log-CPM" = mk(pc_tmm[, 1], "pc1t"))
A1 <- do.call(rbind, lapply(names(CAND), function(nm) {
  f  <- de_fit(p, "lithium", adjust = CAND[[nm]], norm = "TMM", method = "voom")
  dg <- deg_ids(f$tt)
  data.frame(covariate = nm, deg_fdr05 = length(dg), pi0 = pi0_storey(f$tt$P.Value),
             pct_of_baseline = 100 * length(dg) / length(ref_deg),
             jaccard_vs_baseline = jacc(dg, ref_deg),
             abs_r_with_log2nf = if (is.null(CAND[[nm]])) NA else abs(cor(y, CAND[[nm]][, 1])),
             stringsAsFactors = FALSE)
}))
print(format(A1, digits = 4), row.names = FALSE)
write.csv(A1, file.path(OUTDIR, "A1_fair_pc1_rival.csv"), row.names = FALSE)

## ---------------------------------------------------------------- A2 -------
cat("\n### A2  confounding vs noise in the factor/lithium association ###\n")
sets <- list(none = NULL, age = "age", sex = "sex", tobacco = "tob", rin = "rin",
             plate = "plate", seq_pcs = "seqpc1+seqpc2+seqpc3",
             full_DE_covariates = "age+sex+tob+rin+plate+seqpc1+seqpc2+seqpc3")
dat <- cbind(meta, y = y, lith = lith)
A2 <- do.call(rbind, lapply(names(sets), function(nm) {
  f  <- if (is.null(sets[[nm]])) "y ~ lith" else paste("y ~ lith +", sets[[nm]])
  s  <- summary(lm(as.formula(f), data = dat))$coefficients["lith", ]
  data.frame(covariates = nm, beta = s[1], se = s[2], t = s[3], p = s[4],
             stringsAsFactors = FALSE)
}))
A2$beta_shift_vs_crude <- 100 * (A2$beta / A2$beta[1] - 1)
A2$se_shift_vs_crude   <- 100 * (A2$se   / A2$se[1]   - 1)
print(format(A2, digits = 4), row.names = FALSE)
cat("beta moving = confounding removed; se moving = noise removed. Here it is mostly noise.\n")
write.csv(A2, file.path(OUTDIR, "A2_covariate_decomposition.csv"), row.names = FALSE)

## ---------------------------------------------------------------- A3 -------
cat("\n### A3  the offset determines the gene list ###\n")
s <- read.csv(file.path(OUTDIR, "summary.csv"))
A3 <- s[, c("config", "norm", "method", "median_lfc_shift", "sd_lfc_shift",
            "deg_up", "deg_dn", "deg_fdr05", "jaccard_vs_TMM")]
A3$log2_up_over_dn <- log2(A3$deg_up / A3$deg_dn)
A3 <- A3[order(-A3$median_lfc_shift), ]
print(format(A3, digits = 4), row.names = FALSE)
cat(sprintf("\nPearson r(offset, log2 up/down) = %.4f over %d configurations\n",
            cor(A3$median_lfc_shift, A3$log2_up_over_dn), nrow(A3)))
cat(sprintf("linear fit: log2(up/dn) = %.3f + %.2f * offset\n",
            coef(lm(log2_up_over_dn ~ median_lfc_shift, A3))[1],
            coef(lm(log2_up_over_dn ~ median_lfc_shift, A3))[2]))
write.csv(A3, file.path(OUTDIR, "A3_offset_vs_direction.csv"), row.names = FALSE)

## ---------------------------------------------------------------- A4 -------
# How much of lithium's effect on the scaling factor runs through the cell
# mixture? Plain product-of-coefficients on the same covariate block. This is a
# DESCRIPTIVE decomposition, not a causal mediation estimate: no sensitivity
# analysis for unmeasured confounding of the mediator-outcome path, and the
# mediator is CIBERSORTx output estimated from the same matrix, so it inherits
# the circularity caveat in CLAUDE.md 4.2. Reported as a share, not a claim.
cat("\n### A4  how much of lithium -> scaling factor runs through cell mixture ###\n")
gran  <- p$lineage[meta$title, "gran"]
cpmr  <- sweep(cntf, 2, colSums(cntf), "/")
top50 <- apply(cpmr, 2, function(v) sum(sort(v, TRUE)[1:50]))
path <- function(M, nm) {
  a <- coef(lm(M    ~ lith + Xcov - 1))["lith"]
  f <- lm(y ~ M + lith + Xcov - 1)
  b <- coef(f)["M"]; dir <- coef(f)["lith"]
  tot <- coef(lm(y ~ lith + Xcov - 1))["lith"]
  data.frame(mediator = nm, a_lithium_to_mediator = unname(a),
             b_mediator_to_factor = unname(b), indirect = unname(a * b),
             direct = unname(dir), total = unname(tot),
             pct_of_total = unname(100 * a * b / tot), stringsAsFactors = FALSE)
}
A4 <- rbind(path(gran, "granulocyte lineage fraction"),
            path(top50, "top-50-gene share of library"),
            path(p$ILR[meta$title, 1], "ILR b1 (myeloid vs lymphoid)"))
print(format(A4, digits = 4), row.names = FALSE)
cat(sprintf("\ngranulocyte fraction: non-users %.3f -> users %.3f (adjusted shift %+.1f percentage points)\n",
            mean(gran[lith == 0]), mean(gran[lith == 1]),
            100 * coef(lm(gran ~ lith + Xcov - 1))["lith"]))
cat(sprintf("gran vs top-50 share correlate only r = %.3f, so they are largely SEPARATE routes into the factor\n",
            cor(gran, top50)))
write.csv(A4, file.path(OUTDIR, "A4_path_decomposition.csv"), row.names = FALSE)
cat("\ndone.\n")
