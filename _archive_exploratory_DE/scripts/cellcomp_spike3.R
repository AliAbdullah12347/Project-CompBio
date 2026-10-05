#!/usr/bin/env Rscript
# ==============================================================================
# cellcomp_spike3.R -- the corrected compositional spike-in.
#
# WHY A THIRD VERSION. spike2 sweeps the effect size and still detects almost
# no COMPOS genes (2 of 1200 at 1x, 4 at 1.5x). That is not a power problem to
# be solved with a bigger multiplier -- it is a mis-specification, and worth
# recording because it is easy to get wrong.
#
# spike2 ADDS gamma*b1 on top of the gene's existing residual. The gene's total
# noise therefore grows by gamma*sd(b1), and since gamma = d/a1 = 6.4d the
# added noise dwarfs the induced effect. The simulated gene is a noisy gene
# that happens to track composition.
#
# A real composition-driven gene is not like that. Its expression IS the
# mixture-weighted average of the cell types present, so composition explains
# part of the variance it already had -- it does not bolt extra variance on.
# The correct construction takes the compositional signal OUT OF the gene's
# existing noise budget:
#
#     E_new = (covariate fit, lithium removed) + gamma*(b1 - b1bar) + eps
#     eps   = original residual rescaled so total residual variance is unchanged
#
# The gene then has no direct lithium term at all -- every bit of its apparent
# lithium association is inherited through a1 -- while its detectability stays
# comparable to a TRANS gene of the same induced effect size. That is the
# apples-to-apples floor.
#
# Both arms get the same variance-preserving treatment so the comparison is
# symmetric.
#
# Writes only to experimentation/runs/cellcomp/.
# ==============================================================================

source("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation/scripts/common.R")
set.seed(20261218)
OUT <- file.path(EXP_HOME, "runs", "cellcomp"); dir.create(OUT, FALSE, TRUE)
FDR <- 0.05
p <- load_prep()
cat("== cellcomp_spike3 ==\n")

d  <- build_design(p, "lithium"); Xa <- d$X; mt <- d$meta$title
cnt <- p$counts[, d$samples, drop = FALSE]
cnt <- cnt[filter_genes(cnt, 10, 0.90), , drop = FALSE]
v   <- voom(calcNormFactors(DGEList(counts = cnt), method = "TMM"), Xa)
ILR4 <- p$ILR[mt, , drop = FALSE]
Xb <- cbind(Xa, ILR4); colnames(Xb) <- make.names(colnames(Xb))
lith <- Xa[, d$coef]; b1 <- ILR4[, 1]
oth  <- setdiff(colnames(Xa), c("(Intercept)", d$coef))
X0   <- Xa[, c("(Intercept)", oth), drop = FALSE]     # design WITHOUT lithium
a1 <- coef(summary(lm(b1 ~ Xa[, oth] + lith)))
a1 <- a1[grep("lith", rownames(a1)), 1]
colnames(Xa) <- make.names(colnames(Xa)); CO <- make.names(d$coef)

ref <- de_fit(p, "lithium"); can <- ref$tt
eff_pool  <- abs(can$logFC[match(deg_ids(can, FDR), can$gene)])
null_pool <- which(can$P.Value > 0.5)

# Baseline decomposition of every gene against the lithium-free design. Doing
# it once keeps the loop cheap.
H0   <- X0 %*% solve(crossprod(X0), t(X0))
FIT0 <- v$E %*% t(H0)                     # covariate-explained part
RES0 <- v$E - FIT0                        # everything else, incl. real lithium
sdR  <- sqrt(rowMeans(RES0^2))
cat(sprintf("%d genes; sd(b1) = %.4f; a1 = %.4f; median residual sd = %.4f\n",
            nrow(v$E), sd(b1), a1, median(sdR)))

b1c <- b1 - mean(b1)
fitboth <- function(E) list(
  a = topTable(eBayes(lmFit(E, Xa, weights = v$weights)), coef = CO,
               number = Inf, sort.by = "none"),
  b = topTable(eBayes(lmFit(E, Xb, weights = v$weights)), coef = CO,
               number = Inf, sort.by = "none"))

N <- 1200
MULT <- c(0.5, 1, 1.5, 2, 3)
rows <- list()
for (arm in c("TRANS", "COMPOS")) {
  for (mu in MULT) {
    i  <- sample(null_pool, N)
    dd <- sample(eff_pool, N, replace = TRUE) * sample(c(-1, 1), N, TRUE) * mu
    # the driver the effect is injected on, and the coefficient needed to make
    # the induced lithium logFC equal dd
    drv  <- if (arm == "TRANS") lith - mean(lith) else b1c
    gam  <- if (arm == "TRANS") dd else dd / a1
    sig  <- outer(gam, drv)
    # variance preservation: the injected signal comes out of the gene's own
    # residual budget rather than being added to it
    vs   <- rowMeans(sig^2)
    frac <- 1 - vs / (sdR[i]^2)
    capped <- sum(frac < 0)
    eps  <- RES0[i, ] * sqrt(pmax(frac, 0))
    E <- v$E; E[i, ] <- FIT0[i, ] + sig + eps

    f <- fitboth(E)
    det_a <- which(f$a$adj.P.Val[i] < FDR)
    surv  <- sum(f$b$adj.P.Val[i][det_a] < FDR)
    rows[[length(rows) + 1]] <- data.frame(
      arm = arm, multiplier = mu, n_spiked = N, n_capped = capped,
      mean_induced_absFC = round(mean(abs(f$a$logFC[i])), 4),
      detected_canonical = length(det_a),
      power_canonical = round(length(det_a) / N, 4),
      detected_adjusted = sum(f$b$adj.P.Val[i] < FDR),
      survivors = surv,
      survival_rate = round(surv / max(1, length(det_a)), 4),
      shrinkage = round(median(abs(f$b$logFC[i][det_a])) /
                        median(abs(f$a$logFC[i][det_a])), 4))
    r <- rows[[length(rows)]]
    cat(sprintf("  %-7s x%-4.1f  capped %4d  inducedFC %.3f  canon %4d (%.2f)  surv %4d (%.3f)  shrink %.3f\n",
                arm, mu, capped, r$mean_induced_absFC, r$detected_canonical,
                r$power_canonical, surv, r$survival_rate, r$shrinkage))
  }
}
cal <- do.call(rbind, rows)

real_deg  <- deg_ids(can, FDR)
real_surv <- length(intersect(real_deg,
               deg_ids(de_fit(p, "lithium", adjust = p$ILR)$tt, FDR))) / length(real_deg)

pooled <- do.call(rbind, lapply(c("TRANS", "COMPOS"), function(a) {
  s <- cal[cal$arm == a, ]
  data.frame(arm = a, detected = sum(s$detected_canonical),
             survivors = sum(s$survivors),
             survival_rate = round(sum(s$survivors) / sum(s$detected_canonical), 4),
             mean_shrinkage = round(mean(s$shrinkage), 4))
}))
cat("\n--- pooled anchors (variance-preserving spike) ---\n"); print(pooled, row.names = FALSE)

tr <- pooled$survival_rate[pooled$arm == "TRANS"]
cp <- pooled$survival_rate[pooled$arm == "COMPOS"]
cat(sprintf("\nCOMPOS floor %.3f  |  REAL %.3f  |  TRANS ceiling %.3f\n",
            cp, real_surv, tr))
cat(sprintf("=> the real lithium signature behaves as %.0f%% composition-driven.\n",
            100 * (1 - (real_surv - cp) / (tr - cp))))

write.csv(cal,    file.path(OUT, "spikein_varpreserving_curve.csv"), row.names = FALSE)
write.csv(pooled, file.path(OUT, "spikein_varpreserving_anchors.csv"), row.names = FALSE)
cat("\ndone.\n")
