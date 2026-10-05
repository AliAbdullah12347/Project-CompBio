#!/usr/bin/env Rscript
# ==============================================================================
# cellcomp_power.R -- calibrating the headline number from cellcomp.R.
#
# cellcomp.R found that adjusting the lithium contrast for cell composition
# destroys ~78% of the differentially expressed genes. That number on its own
# is ambiguous, because adjustment removes DEGs for two quite different reasons:
#
#   1. CONFOUNDING -- the gene never had a transcriptional lithium effect; it
#      moved because the cell mixture moved. This is the finding we care about.
#   2. POWER       -- lithium is correlated with composition, so conditioning on
#      composition removes variance from the exposure itself. A gene with a
#      genuine, purely transcriptional lithium effect can drop below FDR purely
#      because its standard error grew. This is collateral damage.
#
# The permuted-ILR control in cellcomp.R prices the degrees of freedom but NOT
# reason 2, because permuting destroys the very lithium-composition correlation
# that causes the variance inflation. So here we calibrate directly, by spiking
# known effects into the real expression matrix:
#
#   TRANS  genes: a purely transcriptional lithium effect, orthogonal to
#                 composition by construction. Survival rate = the ceiling.
#                 Anything below this is what adjustment costs an honest signal.
#   COMPOS genes: an effect applied through the myeloid/lymphoid balance, so
#                 the gene only looks lithium-associated because lithium shifts
#                 the balance. Survival rate = the floor. Adjustment SHOULD
#                 annihilate these, and if it does not, the adjustment is weak.
#
# The real 78% is then read against those two anchors.
#
# Writes only to experimentation/runs/cellcomp/.
# ==============================================================================

source("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation/scripts/common.R")
set.seed(20261218)
OUT <- file.path(EXP_HOME, "runs", "cellcomp"); dir.create(OUT, FALSE, TRUE)
FDR <- 0.05
p <- load_prep()
cat("== cellcomp_power ==\n")

## ---------------------------------------------------------------------------
## Rebuild the canonical fit by hand. This is the same sequence de_fit() runs;
## it is repeated here only because the spike-in needs the voom object, which
## de_fit does not return. The gene count is asserted against de_fit's so the
## two cannot silently drift apart.
## ---------------------------------------------------------------------------
d  <- build_design(p, "lithium")
Xa <- d$X
mt <- d$meta$title
cnt <- p$counts[, d$samples, drop = FALSE]
cnt <- cnt[filter_genes(cnt, 10, 0.90), , drop = FALSE]
dge <- calcNormFactors(DGEList(counts = cnt), method = "TMM")
v   <- voom(dge, Xa)
ref <- de_fit(p, "lithium")
stopifnot(identical(rownames(v$E), rownames(ref$tt)), nrow(v$E) == ref$n_genes)
cat(sprintf("canonical rebuilt: %d genes x %d samples (matches de_fit)\n",
            nrow(v$E), ncol(v$E)))

ILR4 <- p$ILR[mt, , drop = FALSE]
Xb   <- cbind(Xa, ILR4)
colnames(Xb) <- make.names(colnames(Xb))
lith <- Xa[, d$coef]
b1   <- ILR4[, 1]

## ---------------------------------------------------------------------------
## COLLINEARITY. How much of lithium's usable variance does conditioning on
## composition actually consume? This is the mechanical source of power loss,
## and it is a property of the design matrix alone -- no genes involved.
## ---------------------------------------------------------------------------
oth  <- setdiff(colnames(Xa), c("(Intercept)", d$coef))
r_a  <- resid(lm(lith ~ Xa[, oth]))                      # lithium | covariates
r_b  <- resid(lm(lith ~ cbind(Xa[, oth], ILR4)))         # lithium | covariates + ILR
vif  <- var(r_a) / var(r_b)
r2   <- 1 - var(r_b) / var(r_a)
collin <- data.frame(
  resid_var_canonical = round(var(r_a), 5),
  resid_var_adjusted  = round(var(r_b), 5),
  R2_lithium_on_ILR_given_covars = round(r2, 4),
  variance_inflation  = round(vif, 4),
  se_inflation        = round(sqrt(vif), 4),
  # residual df also falls by 4, a second, smaller penalty
  df_canonical = ncol(v$E) - ncol(Xa), df_adjusted = ncol(v$E) - ncol(Xb))
cat("\n--- collinearity between lithium and composition ---\n")
print(collin, row.names = FALSE)
cat(sprintf("=> every standard error grows by x%.3f before any gene is considered.\n",
            sqrt(vif)))

## ---------------------------------------------------------------------------
## SPIKE-IN. Effects are drawn from the observed canonical DEG |logFC|
## distribution so the synthetic signal is the same size as the real one --
## calibrating against an unrealistically large effect would flatter the
## adjustment.
## ---------------------------------------------------------------------------
can <- ref$tt
can_deg <- deg_ids(can, FDR)
eff_pool <- abs(can$logFC[match(can_deg, can$gene)])
cat(sprintf("\neffect pool: %d canonical DEGs, |logFC| median %.3f IQR %.3f-%.3f\n",
            length(eff_pool), median(eff_pool),
            quantile(eff_pool, .25), quantile(eff_pool, .75)))

# Spike only into genes with no detectable canonical lithium signal, so the
# planted effect is the whole effect and is not sitting on top of a real one.
null_pool <- which(can$P.Value > 0.5)
NSPIKE <- 800
idx <- sample(null_pool, 2 * NSPIKE)
i_trans <- idx[1:NSPIKE]; i_comp <- idx[(NSPIKE + 1):(2 * NSPIKE)]

delta <- sample(eff_pool, 2 * NSPIKE, replace = TRUE) * sample(c(-1, 1), 2 * NSPIKE, TRUE)
d_tr <- delta[1:NSPIKE]; d_cp <- delta[(NSPIKE + 1):(2 * NSPIKE)]

# Path a for b1, needed to size the compositional spike so that the effect it
# INDUCES on the lithium coefficient matches d_cp. Recomputed here rather than
# read from the other script's CSV, so this file stands alone.
a1 <- coef(summary(lm(b1 ~ Xa[, oth] + lith)))
a1 <- a1[grep("lith", rownames(a1)), 1]
cat(sprintf("path a (lithium -> b1) = %.4f\n", a1))

Es <- v$E
# TRANS: effect enters directly on lithium, carrying no composition component.
Es[i_trans, ] <- Es[i_trans, ] + outer(d_tr, lith - mean(lith))
# COMPOS: effect enters only through the balance. The gene has no lithium term
# at all; its apparent lithium association is entirely inherited from a1.
Es[i_comp, ]  <- Es[i_comp, ]  + outer(d_cp / a1, b1 - mean(b1))

spike_fit <- function(X, E) {
  tt <- topTable(eBayes(lmFit(E, X, weights = v$weights)),
                 coef = make.names(d$coef), number = Inf, sort.by = "none")
  tt$gene <- rownames(tt); tt
}
colnames(Xa) <- make.names(colnames(Xa))
s_a <- spike_fit(Xa, Es); s_b <- spike_fit(Xb, Es)

hit <- function(tt, i) sum(tt$adj.P.Val[i] < FDR)
unspiked <- setdiff(seq_len(nrow(Es)), c(i_trans, i_comp))

spike <- data.frame(
  gene_class = c("TRANS (true transcriptional)", "COMPOS (composition-driven)",
                 "unspiked background"),
  n = c(NSPIKE, NSPIKE, length(unspiked)),
  detected_canonical = c(hit(s_a, i_trans), hit(s_a, i_comp), hit(s_a, unspiked)),
  detected_adjusted  = c(hit(s_b, i_trans), hit(s_b, i_comp), hit(s_b, unspiked)))
spike$survival_rate <- round(spike$detected_adjusted /
                             pmax(1, spike$detected_canonical), 4)
spike$median_absFC_canonical <- round(c(
  median(abs(s_a$logFC[i_trans])), median(abs(s_a$logFC[i_comp])),
  median(abs(s_a$logFC[unspiked]))), 4)
spike$median_absFC_adjusted <- round(c(
  median(abs(s_b$logFC[i_trans])), median(abs(s_b$logFC[i_comp])),
  median(abs(s_b$logFC[unspiked]))), 4)
spike$shrinkage <- round(spike$median_absFC_adjusted /
                         spike$median_absFC_canonical, 4)
cat("\n--- spike-in calibration ---\n"); print(spike, row.names = FALSE)

rs <- length(intersect(can_deg, deg_ids(de_fit(p, "lithium", adjust = p$ILR)$tt, FDR))) /
      length(can_deg)
cat(sprintf("\nsurvival:  TRANS %.3f (ceiling)   REAL %.3f   COMPOS %.3f (floor)\n",
            spike$survival_rate[1], rs, spike$survival_rate[2]))
# Where does the real signature sit between the two anchors? 0 = fully
# compositional, 1 = fully transcriptional.
pos <- (rs - spike$survival_rate[2]) /
       (spike$survival_rate[1] - spike$survival_rate[2])
cat(sprintf("real signature sits %.1f%% of the way from 'pure composition' to 'pure transcription'\n",
            100 * pos))

## ---------------------------------------------------------------------------
## Does the conclusion survive a change of estimator? If the 78% is an artifact
## of voom's weighting it will move a lot here.
## ---------------------------------------------------------------------------
est <- do.call(rbind, lapply(c("voom", "voomWQW", "trend", "QLF"), function(mm) {
  fa <- de_fit(p, "lithium", method = mm)
  fb <- de_fit(p, "lithium", method = mm, adjust = p$ILR)
  da <- deg_ids(fa$tt, FDR); db <- deg_ids(fb$tt, FDR)
  data.frame(method = mm, n_deg_canonical = length(da), n_deg_adjusted = length(db),
             prop_lost = round(length(setdiff(da, db)) / max(1, length(da)), 4),
             r_logFC = round(cor(fa$tt$logFC, fb$tt$logFC), 4))
}))
cat("\n--- estimator sensitivity ---\n"); print(est, row.names = FALSE)

## ---------------------------------------------------------------------------
## The survivors. These are the candidate genuinely transcriptional lithium
## targets -- significant with composition held constant, and with the least
## attenuation. Reported for inspection, not as a validated list.
## ---------------------------------------------------------------------------
fb <- de_fit(p, "lithium", adjust = p$ILR)$tt
surv <- intersect(can_deg, deg_ids(fb, FDR))
si <- match(surv, fb$gene); ci <- match(surv, can$gene)
top <- data.frame(gene = surv,
                  logFC_canonical = round(can$logFC[ci], 4),
                  logFC_adjusted  = round(fb$logFC[si], 4),
                  retained = round(fb$logFC[si] / can$logFC[ci], 3),
                  adjP_adjusted = signif(fb$adj.P.Val[si], 3))
top <- top[order(fb$adj.P.Val[si]), ]
cat("\n--- 25 strongest composition-independent lithium genes ---\n")
print(head(top, 25), row.names = FALSE)
cat(sprintf("\nof %d survivors, %d retain >=90%% of their canonical effect, %d retain >=100%%\n",
            nrow(top), sum(top$retained >= 0.9), sum(top$retained >= 1)))

write.csv(collin, file.path(OUT, "collinearity.csv"), row.names = FALSE)
write.csv(spike,  file.path(OUT, "spikein_calibration.csv"), row.names = FALSE)
write.csv(est,    file.path(OUT, "estimator_sensitivity.csv"), row.names = FALSE)
write.csv(top,    file.path(OUT, "survivor_genes.csv"), row.names = FALSE)
cat("\ndone.\n")
