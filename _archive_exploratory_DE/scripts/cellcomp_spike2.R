#!/usr/bin/env Rscript
# ==============================================================================
# cellcomp_spike2.R -- a proper calibration curve for the spike-in.
#
# WHY THIS EXISTS. The first pass (cellcomp_power.R) gave a clean TRANS ceiling
# from 583 detected genes but a COMPOS floor resting on only 14, which is not
# enough to quote. The reason is a real and interesting asymmetry rather than a
# coding error: a composition-driven gene carries the mediator's variance in
# its residual under the canonical model, because the canonical model does not
# condition on the mediator. Injecting an effect of size d through the balance
# therefore requires a coefficient of d/a1 = d/0.156 = 6.4d on the balance, and
# that inflates the gene's residual SD far more than an effect injected
# directly on the exposure does. Same nominal logFC, much worse t statistic.
#
# So the two arms are not comparable at one effect size. This script sweeps the
# effect multiplier for both arms and reports survival as a function of how
# strongly each class was detected in the first place, which is the comparison
# that actually controls for power.
#
# Writes only to experimentation/runs/cellcomp/.
# ==============================================================================

source("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation/scripts/common.R")
set.seed(20261218)
OUT <- file.path(EXP_HOME, "runs", "cellcomp"); dir.create(OUT, FALSE, TRUE)
FDR <- 0.05
p <- load_prep()
cat("== cellcomp_spike2 ==\n")

d  <- build_design(p, "lithium"); Xa <- d$X; mt <- d$meta$title
cnt <- p$counts[, d$samples, drop = FALSE]
cnt <- cnt[filter_genes(cnt, 10, 0.90), , drop = FALSE]
v   <- voom(calcNormFactors(DGEList(counts = cnt), method = "TMM"), Xa)
ILR4 <- p$ILR[mt, , drop = FALSE]
Xb <- cbind(Xa, ILR4); colnames(Xb) <- make.names(colnames(Xb))
lith <- Xa[, d$coef]; b1 <- ILR4[, 1]
oth <- setdiff(colnames(Xa), c("(Intercept)", d$coef))
a1 <- coef(summary(lm(b1 ~ Xa[, oth] + lith)))
a1 <- a1[grep("lith", rownames(a1)), 1]
colnames(Xa) <- make.names(colnames(Xa))
CO <- make.names(d$coef)

ref <- de_fit(p, "lithium")
can <- ref$tt
eff_pool <- abs(can$logFC[match(deg_ids(can, FDR), can$gene)])
null_pool <- which(can$P.Value > 0.5)
cat(sprintf("%d genes, %d null-pool candidates, a1 = %.4f\n",
            nrow(v$E), length(null_pool), a1))

fitboth <- function(E) {
  ta <- topTable(eBayes(lmFit(E, Xa, weights = v$weights)), coef = CO,
                 number = Inf, sort.by = "none")
  tb <- topTable(eBayes(lmFit(E, Xb, weights = v$weights)), coef = CO,
                 number = Inf, sort.by = "none")
  list(a = ta, b = tb)
}

N <- 1200
MULT <- c(1, 1.5, 2, 3, 4, 6)
rows <- list()
for (arm in c("TRANS", "COMPOS")) {
  for (mu in MULT) {
    i <- sample(null_pool, N)
    dd <- sample(eff_pool, N, replace = TRUE) * sample(c(-1, 1), N, TRUE) * mu
    E <- v$E
    if (arm == "TRANS") E[i, ] <- E[i, ] + outer(dd, lith - mean(lith))
    else                E[i, ] <- E[i, ] + outer(dd / a1, b1 - mean(b1))
    f <- fitboth(E)
    det_a <- which(f$a$adj.P.Val[i] < FDR)
    det_b <- sum(f$b$adj.P.Val[i] < FDR)
    surv  <- sum(f$b$adj.P.Val[i][det_a] < FDR)
    rows[[length(rows) + 1]] <- data.frame(
      arm = arm, multiplier = mu,
      target_absFC = round(mean(abs(dd)), 4),
      n_spiked = N,
      detected_canonical = length(det_a),
      power_canonical = round(length(det_a) / N, 4),
      detected_adjusted = det_b,
      survivors = surv,
      survival_rate = round(surv / max(1, length(det_a)), 4),
      # attenuation among the genes that WERE detected canonically
      shrinkage = round(median(abs(f$b$logFC[i][det_a])) /
                        median(abs(f$a$logFC[i][det_a])), 4),
      median_t_canonical = round(median(abs(f$a$t[i][det_a])), 3))
    cat(sprintf("  %-7s x%-4.1f  canon %4d (power %.2f)  survive %4d (%.3f)  shrink %.3f\n",
                arm, mu, length(det_a), length(det_a) / N, surv,
                surv / max(1, length(det_a)),
                rows[[length(rows)]]$shrinkage))
  }
}
cal <- do.call(rbind, rows)

## ---------------------------------------------------------------------------
## Read the real signature against the curve, matched on canonical power.
## ---------------------------------------------------------------------------
real_n <- length(deg_ids(can, FDR))
real_surv <- length(intersect(deg_ids(can, FDR),
                              deg_ids(de_fit(p, "lithium", adjust = p$ILR)$tt, FDR))) / real_n
cat(sprintf("\nREAL: %d canonical DEGs, survival %.4f\n", real_n, real_surv))

cat("\n--- survival by arm, pooled over multipliers weighted by detections ---\n")
pooled <- do.call(rbind, lapply(c("TRANS", "COMPOS"), function(a) {
  s <- cal[cal$arm == a, ]
  data.frame(arm = a, total_detected = sum(s$detected_canonical),
             total_survivors = sum(s$survivors),
             survival_rate = round(sum(s$survivors) / sum(s$detected_canonical), 4),
             mean_shrinkage = round(mean(s$shrinkage), 4))
}))
print(pooled, row.names = FALSE)

tr <- pooled$survival_rate[pooled$arm == "TRANS"]
cp <- pooled$survival_rate[pooled$arm == "COMPOS"]
cat(sprintf("\nanchors: COMPOS %.3f  ..  REAL %.3f  ..  TRANS %.3f\n", cp, real_surv, tr))
cat(sprintf("the real lithium signature behaves like a signal that is %.0f%% compositional.\n",
            100 * (1 - (real_surv - cp) / (tr - cp))))

write.csv(cal, file.path(OUT, "spikein_calibration_curve.csv"), row.names = FALSE)
write.csv(pooled, file.path(OUT, "spikein_anchors.csv"), row.names = FALSE)
cat("\ndone.\n")
