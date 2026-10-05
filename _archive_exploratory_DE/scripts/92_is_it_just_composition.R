#!/usr/bin/env Rscript
# ==============================================================================
# 92_is_it_just_composition.R -- the aggregate mediation test, done in closed form.
#
# If lithium changes expression ONLY by shifting the cell mixture, then for every
# gene the algebra is forced. Write b for the myeloid/lymphoid balance:
#
#     expression_g  =  alpha_g + gamma_g * b + noise        (gene responds to mixture)
#     b             =  delta * lithium + noise              (lithium shifts mixture)
#  => marginal lithium effect on gene g  =  gamma_g * delta
#
# So under pure mediation, the vector of lithium logFCs across genes must be
# PROPORTIONAL to the vector of gene-vs-balance slopes, with the constant of
# proportionality equal to delta -- a quantity measured independently.
#
# That gives three checkable predictions rather than one:
#   1. beta_lithium and gamma should be strongly correlated across genes.
#   2. The fitted slope of beta on gamma should equal the separately measured delta.
#   3. Residual scatter identifies the genes lithium acts on DIRECTLY -- the
#      ones a compositional explanation cannot account for.
#
# Prediction 2 is the one that makes this a real test: a correlation alone could
# arise from any shared structure, but the slope matching an independently
# measured effect size is hard to get by accident.
# ==============================================================================

source("scripts/common.R")
p <- load_prep()
set.seed(481)

OUT <- file.path(EXP_HOME, "runs", "just_composition")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

m <- p$meta[p$meta$dx == "BP1", ]
m$tob <- p$TOB[m$title, "tobacco_imp_01"]
m$b1  <- p$ILR[m$title, "b1_myeloid_vs_lymphoid"]
covars <- c("age", "sex", "tob", "rin", "plate", "seqpc1", "seqpc2", "seqpc3")
m <- m[complete.cases(m[, c("lithium", "b1", covars)]), ]
m$plate <- droplevels(m$plate); m$sex <- droplevels(m$sex)
cat(sprintf("bipolar I: n=%d (%d non-users, %d users)\n",
            nrow(m), sum(m$lithium == 0), sum(m$lithium == 1)))

cnt <- p$counts[, m$title, drop = FALSE]
cnt <- cnt[filter_genes(cnt, 10, 0.90), , drop = FALSE]
dge <- calcNormFactors(DGEList(cnt), method = "TMM")
cat(sprintf("genes tested: %d\n", nrow(cnt)))

## ---- path a: lithium -> balance -------------------------------------------
fa <- lm(as.formula(paste("b1 ~ lithium +", paste(covars, collapse = " + "))), data = m)
delta <- coef(fa)["lithium"]; se_d <- summary(fa)$coefficients["lithium", 2]
pa <- summary(fa)$coefficients["lithium", 4]
cat(sprintf("\nPATH a  lithium -> myeloid/lymphoid balance\n"))
cat(sprintf("  delta = %+.4f  (SE %.4f, p = %.3g)  partial R2 = %.4f\n",
            delta, se_d, pa,
            summary(fa)$r.squared - summary(update(fa, . ~ . - lithium))$r.squared))

## ---- beta: lithium -> expression (marginal, NOT adjusted for balance) -----
Xb <- model.matrix(as.formula(paste("~ lithium +", paste(covars, collapse = " + "))), data = m)
vb <- voom(dge, Xb); fb <- eBayes(lmFit(vb, Xb))
beta <- fb$coefficients[, "lithium"]
tt_b <- topTable(fb, coef = "lithium", number = Inf, sort.by = "none")

## ---- gamma: balance -> expression (adjusted for lithium and covariates) ---
# Adjusted for lithium so gamma is the mixture's effect holding exposure fixed,
# which is the b-path of the product of coefficients.
Xg <- model.matrix(as.formula(paste("~ b1 + lithium +", paste(covars, collapse = " + "))), data = m)
vg <- voom(dge, Xg); fg <- eBayes(lmFit(vg, Xg))
gamma <- fg$coefficients[, "b1"]
cprime <- fg$coefficients[, "lithium"]          # direct effect
tt_g <- topTable(fg, coef = "b1", number = Inf, sort.by = "none")

stopifnot(identical(names(beta), names(gamma)))

## ---- the three predictions -------------------------------------------------
r_all <- cor(beta, gamma)
fitline <- lm(beta ~ gamma)
slope <- coef(fitline)["gamma"]; ci <- confint(fitline)["gamma", ]
cat(sprintf("\nPREDICTION 1  beta vs gamma across %d genes\n", length(beta)))
cat(sprintf("  Pearson r = %+.3f   Spearman = %+.3f   R2 = %.3f\n",
            r_all, cor(beta, gamma, method = "spearman"), summary(fitline)$r.squared))
cat(sprintf("\nPREDICTION 2  fitted slope should equal delta = %+.4f\n", delta))
cat(sprintf("  slope = %+.4f  95%% CI [%+.4f, %+.4f]\n", slope, ci[1], ci[2]))
cat(sprintf("  delta inside the slope CI? %s\n",
            if (delta >= ci[1] && delta <= ci[2]) "YES -- consistent with pure mediation"
            else "NO -- a purely compositional account is rejected"))

## ---- prediction 3: who is left over ---------------------------------------
pred <- gamma * delta
resid <- beta - pred
mediated_share <- 1 - var(resid) / var(beta)
cat(sprintf("\nPREDICTION 3  variance of the lithium effect explained by composition\n"))
cat(sprintf("  share of var(beta) accounted for by gamma*delta : %.1f%%\n", 100 * mediated_share))

deg <- rownames(tt_b)[tt_b$adj.P.Val < 0.05]
cat(sprintf("\n  among the %d lithium DEGs: r(beta,gamma) = %+.3f, share = %.1f%%\n",
            length(deg), cor(beta[deg], gamma[deg]),
            100 * (1 - var(beta[deg] - pred[deg]) / var(beta[deg]))))

# Genes whose lithium effect composition cannot explain: large |residual| and
# a direct effect c' that stays significant with the balance in the model.
tt_c <- topTable(fg, coef = "lithium", number = Inf, sort.by = "none")
direct <- rownames(tt_c)[tt_c$adj.P.Val < 0.05]
cat(sprintf("  genes with a DIRECT effect surviving balance adjustment: %d\n", length(direct)))
cat(sprintf("  of which also marginally DE: %d\n", length(intersect(direct, deg))))

out <- data.frame(gene = names(beta), beta_lithium = beta, gamma_balance = gamma,
                  c_prime_direct = cprime, predicted_from_composition = pred,
                  residual = resid,
                  p_marginal = tt_b$P.Value, fdr_marginal = tt_b$adj.P.Val,
                  p_direct = tt_c$P.Value, fdr_direct = tt_c$adj.P.Val,
                  stringsAsFactors = FALSE)
out <- out[order(-abs(out$residual)), ]
write.csv(out, file.path(OUT, "per_gene.csv"), row.names = FALSE)
cat("\n  10 genes least explained by composition (largest |residual|):\n")
print(head(out[, c("gene", "beta_lithium", "predicted_from_composition", "residual", "fdr_direct")], 10),
      row.names = FALSE, digits = 3)

summ <- data.frame(
  delta_path_a = delta, delta_p = pa,
  r_beta_gamma = r_all, slope = slope, slope_lo = ci[1], slope_hi = ci[2],
  delta_in_slope_ci = delta >= ci[1] && delta <= ci[2],
  var_share_mediated = mediated_share,
  n_genes = length(beta), n_deg_marginal = length(deg), n_direct = length(direct),
  stringsAsFactors = FALSE)
write.csv(summ, file.path(OUT, "summary.csv"), row.names = FALSE)
saveRDS(list(beta = beta, gamma = gamma, delta = delta, summ = summ, out = out),
        file.path(OUT, "fits.rds"))
cat("\nwrote runs/just_composition/\n")
