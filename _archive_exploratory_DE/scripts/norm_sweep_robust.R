#!/usr/bin/env Rscript
# ==============================================================================
# norm_sweep_robust.R -- is the factor/lithium association fragile?
#
# The whole interpretation of this experiment rests on one association: log2 of
# the TMM scaling factor is lower in lithium users (adjusted beta -0.0718,
# p = 3.1e-4, partial r = -0.244). Before that is written down anywhere, it
# should survive the obvious ways a correlation of that size fails: a handful
# of influential samples, one bad sequencing plate, or an unlucky draw.
#
# No model fitting here beyond lm, so this is cheap and can run alongside the
# heavier scripts.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})
source("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation/scripts/common.R")

OUTDIR <- file.path(EXP_HOME, "runs", "norm_sweep")
set.seed(20261218)
p   <- load_prep()
OBJ <- readRDS(file.path(OUTDIR, "norm_sweep_objects.rds"))

d    <- build_design(p, "lithium")
meta <- d$meta; lith <- meta$lithium
y    <- log2(OBJ$NF[, "TMM"])
Xcov <- d$X[, setdiff(colnames(d$X), "lithium"), drop = FALSE]

est <- function(yy, ll, XX) {
  m0 <- lm(yy ~ XX - 1); m1 <- lm(yy ~ ll + XX - 1)
  c(beta = unname(coef(m1)["ll"]), p = anova(m0, m1)$`Pr(>F)`[2],
    r = cor(residuals(m0), residuals(lm(ll ~ XX - 1))))
}
base <- est(y, lith, Xcov)
cat(sprintf("baseline: beta %.5f  p %.3e  partial r %.4f\n", base[1], base[2], base[3]))

## ---- R1  leave-one-out influence ------------------------------------------
loo <- t(vapply(seq_along(y), function(i)
  est(y[-i], lith[-i], Xcov[-i, , drop = FALSE]), numeric(3)))
cat(sprintf("\nR1 leave-one-out over %d samples:\n", length(y)))
cat(sprintf("   beta range      %.5f to %.5f (baseline %.5f)\n",
            min(loo[, 1]), max(loo[, 1]), base[1]))
cat(sprintf("   p range         %.3e to %.3e\n", min(loo[, 2]), max(loo[, 2])))
cat(sprintf("   p < 0.05 in %d of %d single-sample deletions\n",
            sum(loo[, 2] < 0.05), nrow(loo)))
infl <- order(abs(loo[, 1] - base[1]), decreasing = TRUE)[1:5]
cat("   five most influential samples (delta beta):\n")
for (i in infl) cat(sprintf("     %s  lithium=%d  log2nf=%+.3f  delta=%+.5f\n",
                            meta$title[i], lith[i], y[i], loo[i, 1] - base[1]))

## ---- R2  trimmed: drop the most extreme normalisation factors -------------
cat("\nR2 trimming the most extreme scaling factors:\n")
R2 <- do.call(rbind, lapply(c(0, 0.01, 0.025, 0.05, 0.10), function(tr) {
  lo <- quantile(y, tr); hi <- quantile(y, 1 - tr)
  k  <- y >= lo & y <= hi
  e  <- est(y[k], lith[k], Xcov[k, , drop = FALSE])
  data.frame(trim_each_tail = tr, n = sum(k), beta = e[1], p = e[2], partial_r = e[3])
}))
print(format(R2, digits = 4), row.names = FALSE)

## ---- R3  within-plate ------------------------------------------------------
# Plate is already in the model, but a single aberrant plate can still drive an
# adjusted coefficient. Refitting with each plate removed in turn shows whether
# any one of them is carrying the result.
cat("\nR3 dropping one sequencing plate at a time:\n")
pl <- droplevels(meta$plate)
R3 <- do.call(rbind, lapply(levels(pl), function(L) {
  k <- pl != L
  if (length(unique(lith[k])) < 2 || sum(k) < 40) return(NULL)
  XX <- Xcov[k, , drop = FALSE]
  # dropping a plate makes its own dummy column all-zero, hence constant; strip
  # constants first, then let limma catch anything still aliased
  keepcol <- apply(XX, 2, function(v) length(unique(v)) > 1) | colnames(XX) == "(Intercept)"
  XX <- XX[, keepcol, drop = FALSE]
  ne <- limma::nonEstimable(XX)
  if (!is.null(ne)) XX <- XX[, setdiff(colnames(XX), ne), drop = FALSE]
  e <- est(y[k], lith[k], XX)
  data.frame(plate_dropped = L, n_removed = sum(!k), n = sum(k),
             beta = e[1], p = e[2], partial_r = e[3])
}))
print(format(R3, digits = 4), row.names = FALSE)
cat(sprintf("   p < 0.05 with every plate dropped in turn: %s\n",
            if (all(R3$p < 0.05)) "YES" else
              paste("NO -- fails when dropping", paste(R3$plate_dropped[R3$p >= 0.05], collapse = ", "))))

## ---- R4  bootstrap over subjects ------------------------------------------
cat("\nR4 bootstrap over subjects (2000 resamples):\n")
B  <- 2000
bs <- replicate(B, {
  i <- sample.int(length(y), replace = TRUE)
  XX <- Xcov[i, , drop = FALSE]
  ne <- limma::nonEstimable(XX)
  if (!is.null(ne)) XX <- XX[, setdiff(colnames(XX), ne), drop = FALSE]
  ll <- lith[i]
  if (length(unique(ll)) < 2) return(NA_real_)
  tryCatch(unname(coef(lm(y[i] ~ ll + XX - 1))["ll"]), error = function(e) NA_real_)
})
bs <- bs[!is.na(bs)]
ci <- quantile(bs, c(.025, .975))
cat(sprintf("   beta %.5f, bootstrap 95%% CI [%.5f, %.5f], %d/%d resamples negative\n",
            base[1], ci[1], ci[2], sum(bs < 0), length(bs)))

R  <- rbind(
  data.frame(check = "baseline", detail = "full sample, full covariate block",
             beta = base[1], p = base[2], partial_r = base[3]),
  data.frame(check = "leave-one-out", detail = sprintf("worst of %d deletions", nrow(loo)),
             beta = loo[which.max(loo[, 2]), 1], p = max(loo[, 2]),
             partial_r = loo[which.max(loo[, 2]), 3]),
  data.frame(check = "trimmed 10% each tail", detail = sprintf("n = %d", R2$n[R2$trim_each_tail == .10]),
             beta = R2$beta[R2$trim_each_tail == .10], p = R2$p[R2$trim_each_tail == .10],
             partial_r = R2$partial_r[R2$trim_each_tail == .10]),
  data.frame(check = "drop one plate", detail = sprintf("worst of %d plates", nrow(R3)),
             beta = R3$beta[which.max(R3$p)], p = max(R3$p),
             partial_r = R3$partial_r[which.max(R3$p)]),
  data.frame(check = "bootstrap", detail = sprintf("95%% CI [%.4f, %.4f]", ci[1], ci[2]),
             beta = base[1], p = NA, partial_r = base[3]))
write.csv(R,  file.path(OUTDIR, "R_robustness_summary.csv"), row.names = FALSE)
write.csv(R2, file.path(OUTDIR, "R2_trimmed.csv"), row.names = FALSE)
write.csv(R3, file.path(OUTDIR, "R3_drop_plate.csv"), row.names = FALSE)
write.csv(data.frame(title = meta$title, lithium = lith, log2nf = y,
                     loo_beta = loo[, 1], loo_p = loo[, 2]),
          file.path(OUTDIR, "R1_leave_one_out.csv"), row.names = FALSE)
cat("\n-- robustness summary --\n"); print(format(R, digits = 4), row.names = FALSE)
cat("\ndone.\n")
