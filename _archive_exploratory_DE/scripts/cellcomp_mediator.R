#!/usr/bin/env Rscript
# ==============================================================================
# cellcomp_mediator.R -- interrogating the mediator itself.
#
# Everything else in this experiment treats the cell fractions as a given. They
# are not: they are ESTIMATED, from the same expression matrix that supplies the
# outcome. Three things therefore need checking before the mediation story can
# be believed, and none of them involves a differential expression fit.
#
#   1. WHICH cell types move? The lineage-level path a says "myeloid up", but
#      lithium's documented haematological effect is specifically neutrophilia.
#      If the shift is spread evenly over unrelated types that is a red flag
#      for a deconvolution artifact rather than a drug effect.
#   2. Is the DECONVOLUTION ITSELF biased by lithium? CIBERSORTx returns a fit
#      p-value, correlation and RMSE per sample. If lithium users are fitted
#      systematically worse, their fractions are systematically wrong, and an
#      apparent composition shift could be a fitting artifact.
#   3. Is the shift SPECIFIC to lithium, or does bipolar illness itself move
#      the blood composition? If illness moves it too, composition is a general
#      nuisance in this cohort rather than a lithium mechanism.
#
# Fast: no DE fits, only per-sample linear models.
# Writes only to experimentation/runs/cellcomp/.
# ==============================================================================

source("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation/scripts/common.R")
OUT <- file.path(EXP_HOME, "runs", "cellcomp"); dir.create(OUT, FALSE, TRUE)
p <- load_prep()
cat("== cellcomp_mediator ==\n")

# Same design the DE models use, so the covariate adjustment is identical.
d  <- build_design(p, "lithium"); Xa <- d$X; mt <- d$meta$title
oth <- setdiff(colnames(Xa), c("(Intercept)", d$coef))
lith <- Xa[, d$coef]

assoc <- function(y, exposure, X_other, label) {
  f  <- lm(y ~ X_other + exposure)
  f0 <- lm(y ~ X_other)
  cf <- coef(summary(f)); r <- grep("^exposure$", rownames(cf))
  data.frame(what = label, beta = cf[r, 1], se = cf[r, 2], t = cf[r, 3],
             p = cf[r, 4],
             partial_R2 = (sum(resid(f0)^2) - sum(resid(f)^2)) / sum(resid(f0)^2),
             stringsAsFactors = FALSE)
}

## ---- 1. which of the 22 types move with lithium ---------------------------
fr <- p$frac[mt, , drop = FALSE]
keep <- apply(fr, 2, var) > 0
cells <- do.call(rbind, lapply(colnames(fr)[keep], function(k)
  assoc(fr[, k], lith, Xa[, oth], k)))
cells$fdr <- p.adjust(cells$p, "BH")
cells$mean_nonuser <- sapply(colnames(fr)[keep], function(k) mean(fr[lith == 0, k]))
cells$mean_user    <- sapply(colnames(fr)[keep], function(k) mean(fr[lith == 1, k]))
cells$diff_pct_pts <- 100 * (cells$mean_user - cells$mean_nonuser)
cells <- cells[order(cells$p), ]
num <- c("beta", "se", "t", "partial_R2", "mean_nonuser", "mean_user", "diff_pct_pts")
cells[num] <- lapply(cells[num], function(x) round(x, 5))
cells$p <- signif(cells$p, 3); cells$fdr <- signif(cells$fdr, 3)
cat("\n--- lithium vs each LM22 cell type (n=226, covariate-adjusted) ---\n")
print(cells[, c("what", "mean_nonuser", "mean_user", "diff_pct_pts", "t", "p", "fdr",
                "partial_R2")], row.names = FALSE)
cat(sprintf("\n%d of %d testable types significant at FDR 0.05\n",
            sum(cells$fdr < 0.05), nrow(cells)))

## ---- 2. is the deconvolution itself biased by lithium? --------------------
qc <- p$cs_qc[mt, , drop = FALSE]
qcs <- do.call(rbind, lapply(names(qc), function(k)
  assoc(as.numeric(qc[[k]]), lith, Xa[, oth], paste0("CIBERSORTx_", k))))
qcs[c("beta", "se", "t", "partial_R2")] <-
  lapply(qcs[c("beta", "se", "t", "partial_R2")], function(x) signif(x, 4))
qcs$p <- signif(qcs$p, 3)
qcs$mean_nonuser <- sapply(names(qc), function(k) round(mean(as.numeric(qc[[k]])[lith == 0]), 4))
qcs$mean_user    <- sapply(names(qc), function(k) round(mean(as.numeric(qc[[k]])[lith == 1]), 4))
cat("\n--- deconvolution fit quality vs lithium (a biased mediator would show here) ---\n")
print(qcs, row.names = FALSE)

# Library depth is the other way a fraction estimate can go wrong.
dep <- assoc(log10(p$meta$lib[d$samples]), lith, Xa[, oth], "log10 library size")
cat(sprintf("log10 library size vs lithium: beta %.4f  p %.3f\n", dep$beta, dep$p))

## ---- 3. does illness move composition too? --------------------------------
dc  <- build_design(p, "casecon"); Xc <- dc$X; mtc <- dc$meta$title
othc <- setdiff(colnames(Xc), c("(Intercept)", dc$coef))
caseX <- Xc[, dc$coef]
comp_names <- c(colnames(p$ILR), "gran_prop")
comp_mat <- cbind(p$ILR, gran_prop = p$lineage[, "gran"])
ill <- do.call(rbind, lapply(comp_names, function(k)
  assoc(comp_mat[mtc, k], caseX, Xc[, othc], k)))
ill$contrast <- "case vs control (n=474)"
lit <- do.call(rbind, lapply(comp_names, function(k)
  assoc(comp_mat[mt, k], lith, Xa[, oth], k)))
lit$contrast <- "lithium within BP1 (n=226)"
both <- rbind(lit, ill)
both[c("beta", "se", "t", "partial_R2")] <-
  lapply(both[c("beta", "se", "t", "partial_R2")], function(x) round(x, 4))
both$p <- signif(both$p, 3)
cat("\n--- is the composition shift specific to lithium, or is it illness? ---\n")
print(both[, c("contrast", "what", "beta", "t", "p", "partial_R2")], row.names = FALSE)

## ---- extra: does composition differ between unmedicated BD and controls? --
# The cleanest version of question 3. Controls vs bipolar I who are NOT on
# lithium: any composition difference here cannot be a lithium effect.
mm <- p$meta
sel <- (mm$dx == "Control") | (mm$dx == "BP1" & mm$lithium == 0)
d3 <- build_design(p, "casecon", samples = sel)
X3 <- d3$X; oth3 <- setdiff(colnames(X3), c("(Intercept)", d3$coef))
off <- do.call(rbind, lapply(comp_names, function(k)
  assoc(comp_mat[d3$meta$title, k], X3[, d3$coef], X3[, oth3], k)))
off[c("beta", "se", "t", "partial_R2")] <-
  lapply(off[c("beta", "se", "t", "partial_R2")], function(x) round(x, 4))
off$p <- signif(off$p, 3)
cat(sprintf("\n--- unmedicated BP1 (n=%d) vs controls: composition ---\n",
            sum(sel & mm$dx == "BP1")))
print(off, row.names = FALSE)

write.csv(cells, file.path(OUT, "lithium_vs_celltypes.csv"), row.names = FALSE)
write.csv(qcs,   file.path(OUT, "deconv_qc_vs_lithium.csv"), row.names = FALSE)
write.csv(both,  file.path(OUT, "composition_lithium_vs_illness.csv"), row.names = FALSE)
write.csv(off,   file.path(OUT, "composition_unmedicated_vs_control.csv"), row.names = FALSE)
cat("\ndone.\n")
