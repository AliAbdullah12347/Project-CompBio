#!/usr/bin/env Rscript
# ==============================================================================
# cellcomp.R -- Is the lithium expression signature a cell-composition shift?
#
# THE central question of the project. Lithium raises neutrophil proportion.
# Bulk RNA from whole blood is a weighted average over cell types, so a shift in
# the weights moves thousands of genes without any cell changing its behaviour.
# If the lithium signature is composition, adjusting for composition should
# delete it. If it is transcriptional, adjustment should leave it standing.
#
# Five models on the lithium contrast (BP1 only, n=226, 74 non-users / 152 users):
#   a  canonical                       -- no composition term
#   b  + 4 ILR balances                -- the pre-specified compositional model
#   c  + 4 of 5 raw lineage proportions-- "minus one" because they sum to 1
#   d  + first 5 PCs of the 22 fractions
#   e  + the single myeloid/lymphoid balance
#   e2 + logit(granulocyte)            -- the "routine blood count" version
#
# Plus a NEGATIVE CONTROL that the brief did not ask for but without which none
# of the above means anything: adding 4 covariates costs 4 residual degrees of
# freedom and shrinks the DEG count even when the covariates are pure noise.
# Model f permutes the rows of the ILR matrix, preserving its marginal
# distribution and its 4 df while destroying its link to the samples. Any DEG
# loss beyond model f's is the part attributable to real confounding.
#
# Writes only to experimentation/runs/cellcomp/.
# ==============================================================================


source("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation/scripts/common.R")

set.seed(20261218)                      # project seed, fixed before any result
OUT <- file.path(EXP_HOME, "runs", "cellcomp")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
FDR <- 0.05

p <- load_prep()
cat("== cellcomp ==\n")

## ---------------------------------------------------------------------------
## Adjustment matrices. All are built on all 474 samples and indexed by title
## inside build_design, so no experiment can accidentally align them by
## position. All are per-sample transforms of a per-sample quantity, so none
## of them leaks information across samples.
## ---------------------------------------------------------------------------

ILR4 <- p$ILR

# One lineage must be omitted: the five proportions sum to 1 exactly, so the
# full set is rank-deficient. B cells are the reference (smallest lineage).
# Uses the zero-replaced closure so this is the same composition the ILR sees.
LIN4 <- p$lineage_z[, c("gran", "mono", "T", "NK"), drop = FALSE]
colnames(LIN4) <- paste0("lin_", colnames(LIN4))

# PCA on the 22 fractions. center only, NOT scale: four LM22 types are exactly
# zero in all 474 samples of this CIBERSORTx b-mode run, so scaling divides by
# zero. Those four columns contribute nothing and are dropped first.
frac_nz <- p$frac[, apply(p$frac, 2, var) > 0, drop = FALSE]
pca <- prcomp(frac_nz, center = TRUE, scale. = FALSE)
ve <- pca$sdev^2 / sum(pca$sdev^2)
PC5 <- pca$x[, 1:5, drop = FALSE]
colnames(PC5) <- paste0("fracPC", 1:5)
cat(sprintf("fraction PCA: %d non-zero types, PC1-5 explain %.1f%% of variance\n",
            ncol(frac_nz), 100 * sum(ve[1:5])))

B1 <- p$ILR[, 1, drop = FALSE]

# logit of the granulocyte share. This is the quantity a haematology analyser
# actually reports, so it is the composition adjustment a clinician could make
# without any sequencing at all.
g <- p$lineage_z[, "gran"]
GRAN <- matrix(log(g / (1 - g)), ncol = 1,
               dimnames = list(names(g), "logit_gran"))

## ---------------------------------------------------------------------------
## Fit the models. Every argument except `adjust` is held at the canonical
## value, so any difference between rows is the composition term and nothing
## else.
## ---------------------------------------------------------------------------

MODELS <- list(
  a_canonical   = NULL,
  b_ILR4        = ILR4,
  c_lineage4    = LIN4,
  d_fracPC5     = PC5,
  e_ILR1        = B1,
  e2_logit_gran = GRAN
)

fits <- lapply(names(MODELS), function(nm) {
  cat(sprintf("  fitting %s ...\n", nm))
  de_fit(p, "lithium", adjust = MODELS[[nm]])
})
names(fits) <- names(MODELS)

# The gene set and sample set must be identical across models, or "how many
# canonical DEGs survive" is comparing different populations of genes. The
# adjustment columns have no missing values, so complete.cases is unchanged.
gset <- rownames(fits$a_canonical$tt)
for (nm in names(fits)) {
  stopifnot(identical(rownames(fits[[nm]]$tt), gset),
            fits[[nm]]$n_samples == fits$a_canonical$n_samples)
}
cat(sprintf("all models: %d genes x %d samples (identical sets)\n",
            length(gset), fits$a_canonical$n_samples))

## ---------------------------------------------------------------------------
## Comparison table
## ---------------------------------------------------------------------------

can <- fits$a_canonical$tt
can_deg <- deg_ids(can, FDR)

compare <- function(tt, nm) {
  deg <- deg_ids(tt, FDR)
  surv <- intersect(can_deg, deg)
  new  <- setdiff(deg, can_deg)
  lost <- setdiff(can_deg, deg)

  i_can <- match(can_deg, can$gene)
  i_adj <- match(can_deg, tt$gene)
  fc_c <- can$logFC[i_can]; fc_a <- tt$logFC[i_adj]

  is  <- match(surv, can$gene); ia <- match(surv, tt$gene)
  il  <- match(lost, can$gene); ial <- match(lost, tt$gene)

  data.frame(
    model          = nm,
    n_extra_covar  = tt_ncov[[nm]],
    n_genes        = length(tt$gene),
    n_deg          = length(deg),
    pi0            = round(pi0_storey(tt$P.Value), 4),
    canon_survive  = length(surv),
    canon_lost     = length(lost),
    prop_lost      = round(length(lost) / length(can_deg), 4),
    n_new          = length(new),
    # correlation over every gene, not just the DEGs: a high value over DEGs
    # alone would be guaranteed by the selection, not informative.
    r_logFC_all    = round(cor(can$logFC, tt$logFC), 4),
    r_logFC_degs   = round(cor(fc_c, fc_a), 4),
    # shrinkage: 1 means unchanged, <1 means attenuated. Ratio of medians of
    # |logFC| rather than median of ratios, which blows up near zero.
    shrink_degs    = round(median(abs(fc_a)) / median(abs(fc_c)), 4),
    shrink_surv    = round(median(abs(tt$logFC[ia])) / median(abs(can$logFC[is])), 4),
    shrink_lost    = if (length(lost)) round(median(abs(tt$logFC[ial])) /
                                             median(abs(can$logFC[il])), 4) else NA,
    # slope of adjusted on canonical through the origin, over canonical DEGs.
    # A slope near 1 with genes dropping out means power was lost, not effect.
    slope_degs     = round(sum(fc_c * fc_a) / sum(fc_c^2), 4),
    med_absFC_can  = round(median(abs(fc_c)), 4),
    med_absFC_adj  = round(median(abs(fc_a)), 4),
    n_signflip     = sum(sign(fc_c) != sign(fc_a)),
    stringsAsFactors = FALSE)
}

tt_ncov <- list(a_canonical = 0, b_ILR4 = 4, c_lineage4 = 4,
                d_fracPC5 = 5, e_ILR1 = 1, e2_logit_gran = 1)

tab <- do.call(rbind, lapply(names(fits), function(nm) compare(fits[[nm]]$tt, nm)))
cat("\n--- model comparison ---\n"); print(tab, row.names = FALSE)

## ---------------------------------------------------------------------------
## NEGATIVE CONTROL. Permuted ILR: same 4 columns, same marginals, same 4 df,
## no true association with any sample. Whatever DEG loss this produces is the
## price of the degrees of freedom alone.
## ---------------------------------------------------------------------------

NPERM <- 50
cat(sprintf("\n--- negative control: %d permuted-ILR fits ---\n", NPERM))
perm <- t(vapply(seq_len(NPERM), function(i) {
  A <- ILR4[sample(nrow(ILR4)), , drop = FALSE]
  rownames(A) <- rownames(ILR4)                   # relabel: breaks the link only
  f <- de_fit(p, "lithium", adjust = A)
  d <- deg_ids(f$tt, FDR)
  c(n_deg = length(d),
    surv  = length(intersect(can_deg, d)),
    r     = cor(can$logFC, f$tt$logFC),
    shrink = median(abs(f$tt$logFC[match(can_deg, f$tt$gene)])) /
             median(abs(can$logFC[match(can_deg, can$gene)])))
}, numeric(4)))

pnull <- data.frame(
  model = "f_ILR4_permuted_null",
  n_deg_mean = round(mean(perm[, "n_deg"]), 1),
  n_deg_sd   = round(sd(perm[, "n_deg"]), 1),
  n_deg_min  = min(perm[, "n_deg"]), n_deg_max = max(perm[, "n_deg"]),
  surv_mean  = round(mean(perm[, "surv"]), 1),
  prop_lost_mean = round(1 - mean(perm[, "surv"]) / length(can_deg), 4),
  r_mean     = round(mean(perm[, "r"]), 4),
  shrink_mean = round(mean(perm[, "shrink"]), 4))
print(pnull, row.names = FALSE)

obs_lost <- tab$canon_lost[tab$model == "b_ILR4"]
null_lost <- length(can_deg) - perm[, "surv"]
cat(sprintf("\nreal ILR loses %d canonical DEGs; permuted ILR loses %.1f (sd %.1f).\n",
            obs_lost, mean(null_lost), sd(null_lost)))
cat(sprintf("excess loss attributable to genuine confounding: %d genes (%.1f%% of canonical)\n",
            obs_lost - round(mean(null_lost)),
            100 * (obs_lost - mean(null_lost)) / length(can_deg)))
cat(sprintf("permutation p for observed loss: %.4f  (z = %.1f)\n",
            (sum(null_lost >= obs_lost) + 1) / (NPERM + 1),
            (obs_lost - mean(null_lost)) / sd(null_lost)))

## ---------------------------------------------------------------------------
## MEDIATION PATH a. If the mechanism is real, lithium must move the mediator.
## Regress each compositional coordinate on lithium + the same covariates the
## DE model uses. Linear model, not limma -- there is one outcome, not 12,000.
## ---------------------------------------------------------------------------

d_lit <- build_design(p, "lithium")
Xa <- d_lit$X
mt <- d_lit$meta$title

path_a <- function(y, label) {
  y <- y[mt]
  # partial t for the lithium column from the full QR fit
  f <- lm(y ~ Xa[, setdiff(colnames(Xa), "(Intercept)"), drop = FALSE])
  cf <- summary(f)$coefficients
  row <- grep("lithium", rownames(cf), value = TRUE)[1]
  # partial R2 of lithium: comparing to the model without it
  f0 <- lm(y ~ Xa[, setdiff(colnames(Xa), c("(Intercept)", "lithium")), drop = FALSE])
  r2p <- (sum(resid(f0)^2) - sum(resid(f)^2)) / sum(resid(f0)^2)
  data.frame(mediator = label, beta_a = cf[row, 1], se = cf[row, 2],
             t = cf[row, 3], p = cf[row, 4], partial_R2 = r2p,
             sd_mediator = sd(y),
             beta_in_SD = cf[row, 1] / sd(y), stringsAsFactors = FALSE)
}

meds <- c(
  setNames(lapply(colnames(p$ILR), function(k) p$ILR[, k]), colnames(p$ILR)),
  setNames(lapply(colnames(p$lineage), function(k) p$lineage[, k]),
           paste0("prop_", colnames(p$lineage))),
  list(logit_gran = GRAN[, 1],
       neutrophil_raw = p$frac[, "Neutrophils"],
       lymphocyte_raw = rowSums(p$frac[, c(p$LIN$T, p$LIN$NK, p$LIN$B)]),
       monocyte_raw   = p$frac[, "Monocytes"],
       NLR_log = log((p$lineage_z[, "gran"]) /
                     rowSums(p$lineage_z[, c("T", "NK", "B")]))))

pa <- do.call(rbind, lapply(names(meds), function(k) path_a(meds[[k]], k)))
pa$fdr <- p.adjust(pa$p, "BH")
pa[, c("beta_a", "se", "t", "partial_R2", "sd_mediator", "beta_in_SD")] <-
  round(pa[, c("beta_a", "se", "t", "partial_R2", "sd_mediator", "beta_in_SD")], 4)
pa$p <- signif(pa$p, 3); pa$fdr <- signif(pa$fdr, 3)
cat("\n--- mediation path a: lithium -> composition (n=226, same covariates) ---\n")
print(pa, row.names = FALSE)

## ---------------------------------------------------------------------------
## INTERNAL CONSISTENCY. For a linear model the difference-of-coefficients and
## the product-of-coefficients definitions of mediation are algebraically the
## same thing. So for every gene,
##        logFC_canonical - logFC_adjusted  ==  sum_k a_k * b_gk
## where a_k is the lithium->balance_k coefficient above and b_gk is the
## balance_k coefficient in the adjusted gene model. They will not match to
## machine precision because voom reweights, but if they do not match closely
## something is wrong with one of the two fits. This is a check, not a result.
## ---------------------------------------------------------------------------

a_k <- setNames(pa$beta_a[match(colnames(p$ILR), pa$mediator)], colnames(p$ILR))
b_gk <- sapply(colnames(p$ILR), function(k)
  de_fit(p, "lithium", adjust = ILR4, coef_override = make.names(k))$tt$logFC)
rownames(b_gk) <- gset
pred_med <- as.vector(b_gk %*% a_k)
obs_diff <- can$logFC - fits$b_ILR4$tt$logFC
cat(sprintf("\nconsistency check  cor(indirect via a*b, logFC drop) = %.4f  slope %.4f\n",
            cor(pred_med, obs_diff), sum(pred_med * obs_diff) / sum(pred_med^2)))

## ---------------------------------------------------------------------------
## WHICH genes are lost? If composition is the mechanism, the casualties should
## be exactly the genes whose expression tracks the myeloid/lymphoid balance.
## Bin all genes by how strongly they load on b1 in the adjusted model and
## report the survival rate of canonical DEGs in each bin.
## ---------------------------------------------------------------------------

adj <- fits$b_ILR4$tt
lost_flag <- can$gene %in% setdiff(can_deg, deg_ids(adj, FDR))
is_deg <- can$gene %in% can_deg
b1_load <- abs(b_gk[, 1])
qb <- cut(b1_load, quantile(b1_load, 0:5 / 5), include.lowest = TRUE,
          labels = paste0("Q", 1:5))
bin <- do.call(rbind, lapply(levels(qb), function(q) {
  s <- qb == q
  data.frame(b1_loading_quintile = q,
             median_abs_b1 = round(median(b1_load[s]), 4),
             n_canon_deg = sum(is_deg & s),
             n_lost = sum(lost_flag & s),
             pct_lost = round(100 * sum(lost_flag & s) / max(1, sum(is_deg & s)), 1))
}))
cat("\n--- canonical DEG loss by strength of association with the myeloid/lymphoid balance ---\n")
print(bin, row.names = FALSE)

## ---------------------------------------------------------------------------
## Is composition a bigger nuisance than lithium? Test the balances themselves
## as exposures in the adjusted model and count genes.
## ---------------------------------------------------------------------------

comp_deg <- sapply(colnames(p$ILR), function(k) {
  t2 <- de_fit(p, "lithium", adjust = ILR4, coef_override = make.names(k))$tt
  c(n_deg = n_deg(t2, FDR), pi0 = round(pi0_storey(t2$P.Value), 3))
})
cat("\n--- genes associated with each balance, lithium held constant ---\n")
print(t(comp_deg))
cat(sprintf("lithium in the same model: %d DEG, pi0 %.3f\n",
            tab$n_deg[tab$model == "b_ILR4"], tab$pi0[tab$model == "b_ILR4"]))

## ---------------------------------------------------------------------------
## Effect-size floor. If adjustment only removes small effects, the survivors
## are the biologically meaningful ones. Recount at |logFC| thresholds.
## ---------------------------------------------------------------------------

lfc_tab <- do.call(rbind, lapply(c(0, 0.1, 0.2, 0.3, 0.5), function(L)
  data.frame(lfc = L,
             canonical = n_deg(can, FDR, L),
             ILR4      = n_deg(fits$b_ILR4$tt, FDR, L),
             ILR1      = n_deg(fits$e_ILR1$tt, FDR, L),
             logit_gran= n_deg(fits$e2_logit_gran$tt, FDR, L),
             retained  = round(n_deg(fits$b_ILR4$tt, FDR, L) /
                               max(1, n_deg(can, FDR, L)), 3))))
cat("\n--- DEG count vs effect-size floor ---\n"); print(lfc_tab, row.names = FALSE)

## ---------------------------------------------------------------------------
## SPECIFICITY CONTROL. Does composition adjustment hit every signature, or
## just lithium's? Run the same b-vs-a comparison on the case/control contrast,
## where lithium is not the exposure.
## ---------------------------------------------------------------------------

cc_a <- de_fit(p, "casecon")
cc_b <- de_fit(p, "casecon", adjust = ILR4)
cc_can <- deg_ids(cc_a$tt, FDR)
cc <- data.frame(
  contrast = c("lithium (BP1)", "case vs control (474)"),
  n_canon_deg = c(length(can_deg), length(cc_can)),
  n_deg_ILR4  = c(tab$n_deg[tab$model == "b_ILR4"], n_deg(cc_b$tt, FDR)),
  prop_lost   = c(tab$prop_lost[tab$model == "b_ILR4"],
                  round(length(setdiff(cc_can, deg_ids(cc_b$tt, FDR))) /
                        length(cc_can), 4)),
  r_logFC     = c(tab$r_logFC_all[tab$model == "b_ILR4"],
                  round(cor(cc_a$tt$logFC, cc_b$tt$logFC), 4)))
cat("\n--- specificity: same adjustment applied to a different exposure ---\n")
print(cc, row.names = FALSE)

## ---------------------------------------------------------------------------
## Write everything out.
## ---------------------------------------------------------------------------

write.csv(pa,      file.path(OUT, "path_a_mediator_models.csv"), row.names = FALSE)
write.csv(pnull,   file.path(OUT, "permuted_ILR_null.csv"),      row.names = FALSE)
write.csv(bin,     file.path(OUT, "loss_by_b1_loading.csv"),     row.names = FALSE)
write.csv(lfc_tab, file.path(OUT, "deg_by_lfc_floor.csv"),       row.names = FALSE)
write.csv(cc,      file.path(OUT, "specificity_casecon.csv"),    row.names = FALSE)
write.csv(data.frame(perm), file.path(OUT, "permutation_draws.csv"), row.names = FALSE)
write.csv(t(comp_deg), file.path(OUT, "balance_as_exposure.csv"))

gene_tab <- data.frame(
  gene = can$gene,
  logFC_canonical = round(can$logFC, 5),
  adjP_canonical  = signif(can$adj.P.Val, 4),
  logFC_ILR4      = round(adj$logFC, 5),
  adjP_ILR4       = signif(adj$adj.P.Val, 4),
  b1_coef         = round(b_gk[, 1], 5),
  indirect_pred   = round(pred_med, 5),
  status = ifelse(is_deg & lost_flag, "lost",
           ifelse(is_deg, "survives",
           ifelse(can$gene %in% deg_ids(adj, FDR), "new", "never"))))
write.csv(gene_tab, file.path(OUT, "gene_level.csv"), row.names = FALSE)

write_run("cellcomp", tab, notes = NULL)   # NOTES.md written by hand, not clobbered
cat("\ndone.\n")
