#!/usr/bin/env Rscript
# ==============================================================================
# 98_interaction_surrogate.R -- two experiments the workflow lost to the session
# limit, rewritten directly.
#
# PART A -- exposure x mediator interaction.
#   The product-of-coefficients decomposition the project pre-specifies is only
#   valid if lithium's effect on a gene does NOT depend on the cell mixture. The
#   pre-proposal fixes the decision rule in advance: test lithium x ILR first,
#   and if any interaction survives BH at FDR 0.05, abandon the product of
#   coefficients for the counterfactual formulation of Valeri & VanderWeele.
#   This runs that test. It is a gate, not an exploration.
#
# PART B -- hidden structure.
#   sva is not installed, so surrogate variables are built by hand: fit the
#   canonical model, take residuals, PCA them, and feed leading residual PCs
#   back in as covariates. If the lithium signal collapses once a few
#   unmodelled axes are absorbed, it was never robust. If a residual PC turns
#   out to correlate strongly with the cell balances, that names composition as
#   the dominant unmodelled axis directly.
# ==============================================================================

source("scripts/common.R")
p <- load_prep()
set.seed(481)
OUT <- file.path(EXP_HOME, "runs", "interaction_surrogate")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

m <- p$meta[p$meta$dx == "BP1", ]
m$tob <- p$TOB[m$title, "tobacco_imp_01"]
covars <- c("age", "sex", "tob", "rin", "plate", "seqpc1", "seqpc2", "seqpc3")
m <- m[complete.cases(m[, c("lithium", covars)]), ]
m$plate <- droplevels(m$plate); m$sex <- droplevels(m$sex)
B <- p$ILR[m$title, , drop = FALSE]
cnt <- p$counts[, m$title, drop = FALSE]
cnt <- cnt[filter_genes(cnt, 10, 0.90), , drop = FALSE]
dge <- calcNormFactors(DGEList(cnt), method = "TMM")
cat(sprintf("n = %d bipolar I, %d genes\n", nrow(m), nrow(cnt)))

base_rhs <- paste(covars, collapse = " + ")

## ===================== PART A -- interaction gate ==========================
cat("\n================ PART A: exposure x covariate interactions ================\n")
inter_test <- function(mod_var, label) {
  d <- cbind(m, v = mod_var)
  X <- model.matrix(as.formula(paste("~ lithium * v +", base_rhs)), data = d)
  ne <- nonEstimable(X); if (!is.null(ne)) X <- X[, setdiff(colnames(X), ne), drop = FALSE]
  cf <- grep("^lithium:", colnames(X), value = TRUE)
  if (!length(cf)) return(NULL)
  v <- voom(dge, X); f <- eBayes(lmFit(v, X))
  tt <- topTable(f, coef = cf[1], number = Inf, sort.by = "none")
  data.frame(term = label, n_fdr05 = sum(tt$adj.P.Val < 0.05),
             n_fdr10 = sum(tt$adj.P.Val < 0.10),
             min_fdr = min(tt$adj.P.Val), pi0 = pi0_storey(tt$P.Value),
             frac_p05 = mean(tt$P.Value < 0.05), stringsAsFactors = FALSE)
}
A <- do.call(rbind, list(
  inter_test(B[, "b1_myeloid_vs_lymphoid"], "lithium x b1 (myeloid/lymphoid)"),
  inter_test(B[, "b2_gran_vs_mono"],        "lithium x b2 (gran/mono)"),
  inter_test(m$age,                          "lithium x age"),
  inter_test(as.numeric(m$sex == levels(m$sex)[2]), "lithium x sex")))
print(A, row.names = FALSE, digits = 3)
cat(sprintf("\n  Under the null, frac_p05 should sit near 0.05 and pi0 near 1.\n"))
gate <- sum(A$n_fdr05) == 0
cat(sprintf("  GATE: any interaction surviving BH at FDR 0.05? %s\n",
            if (gate) "NO -- product of coefficients remains valid"
            else "YES -- switch to the counterfactual formulation"))
write.csv(A, file.path(OUT, "interaction.csv"), row.names = FALSE)

## ===================== PART B -- surrogate variables =======================
cat("\n================ PART B: hidden structure =================================\n")
X0 <- model.matrix(as.formula(paste("~ lithium +", base_rhs)), data = m)
v0 <- voom(dge, X0); f0 <- lmFit(v0, X0)
R  <- v0$E - f0$coefficients %*% t(X0)           # residual matrix, genes x samples
pca <- prcomp(t(R), center = TRUE, scale. = FALSE)
ve <- pca$sdev^2 / sum(pca$sdev^2)
cat(sprintf("residual PCs: variance explained %s\n",
            paste(sprintf("PC%d=%.1f%%", 1:6, 100 * ve[1:6]), collapse = " ")))

cat("\n  what each residual PC correlates with (|r| shown, max over plate levels for factors):\n")
targets <- list(lithium = m$lithium, age = m$age, rin = m$rin,
                depth = m$lib, seqpc1 = m$seqpc1,
                b1_myeloid_lymphoid = B[, "b1_myeloid_vs_lymphoid"],
                b2_gran_mono = B[, "b2_gran_vs_mono"],
                gran_prop = p$lineage[m$title, "gran"])
cormat <- sapply(targets, function(t) abs(cor(pca$x[, 1:5], t)))
rownames(cormat) <- paste0("PC", 1:5)
print(round(cormat, 3))

sv_res <- do.call(rbind, lapply(c(0, 1, 2, 3, 5, 10), function(k) {
  Xk <- if (k == 0) X0 else cbind(X0, pca$x[, seq_len(k), drop = FALSE])
  vk <- voom(dge, Xk); fk <- eBayes(lmFit(vk, Xk))
  tt <- topTable(fk, coef = "lithium", number = Inf, sort.by = "none")
  data.frame(n_sv = k, n_deg = sum(tt$adj.P.Val < 0.05),
             pi0 = pi0_storey(tt$P.Value),
             mean_abs_logFC = mean(abs(tt$logFC)), stringsAsFactors = FALSE)
}))
cat("\n  lithium DEG count as residual PCs are absorbed:\n")
print(sv_res, row.names = FALSE, digits = 3)
write.csv(sv_res, file.path(OUT, "surrogate.csv"), row.names = FALSE)
write.csv(as.data.frame(cormat), file.path(OUT, "residual_pc_correlations.csv"))

## ---- PCA of the full matrix, not just residuals ---------------------------
cat("\n  PCA of the full log-CPM matrix (all 226), |r| with each variable:\n")
lc <- cpm(dge, log = TRUE, prior.count = 3)
pf <- prcomp(t(lc), center = TRUE, scale. = FALSE)
vf <- pf$sdev^2 / sum(pf$sdev^2)
cm2 <- sapply(targets, function(t) abs(cor(pf$x[, 1:5], t)))
rownames(cm2) <- sprintf("PC%d(%.0f%%)", 1:5, 100 * vf[1:5])
print(round(cm2, 3))
write.csv(as.data.frame(cm2), file.path(OUT, "full_pc_correlations.csv"))
cat("\nwrote runs/interaction_surrogate/\n")
