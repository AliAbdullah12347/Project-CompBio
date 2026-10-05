#!/usr/bin/env Rscript
# ==============================================================================
# cellcomp_nonlinear.R -- two loose ends in the headline number.
#
# 1. LINEARITY. Every adjustment so far enters the design linearly. If the true
#    relation between a gene and the cell mixture is curved, a linear term
#    under-corrects, and 78% would be a LOWER bound rather than an estimate.
#    Tested with a natural spline on logit(granulocyte) and, separately, with
#    quintile dummies -- which assume nothing about shape at all.
#
# 2. HOW MUCH IS ENOUGH. Adding composition dimensions must eventually stop
#    removing genes. Sweeping the number of fraction PCs from 1 to 17 gives the
#    asymptote: the most the signature can be reduced by composition alone,
#    and the point past which extra dimensions only cost degrees of freedom.
#
# Raw fraction PCs are used for the sweep rather than ILR because ILR needs a
# sequential binary partition and there is no principled 18-part tree; PCs need
# no tree and span the same space.
#
# Writes only to experimentation/runs/cellcomp/.
# ==============================================================================

suppressPackageStartupMessages(library(splines))
source("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation/scripts/common.R")
set.seed(20261218)
OUT <- file.path(EXP_HOME, "runs", "cellcomp"); dir.create(OUT, FALSE, TRUE)
FDR <- 0.05
p <- load_prep()
cat("== cellcomp_nonlinear ==\n")

g <- p$lineage_z[, "gran"]
lg <- log(g / (1 - g))
frac_nz <- p$frac[, apply(p$frac, 2, var) > 0, drop = FALSE]
pcs <- prcomp(frac_nz, center = TRUE, scale. = FALSE)$x
ve  <- (prcomp(frac_nz, center = TRUE, scale. = FALSE)$sdev)^2
ve  <- cumsum(ve) / sum(ve)
cat(sprintf("%d non-zero cell types -> %d PCs\n", ncol(frac_nz), ncol(pcs)))

mk <- function(M, nm) { M <- as.matrix(M); rownames(M) <- rownames(p$frac)
                        colnames(M) <- paste0(nm, seq_len(ncol(M))); M }

# Quintile dummies: 4 indicator columns, no functional form assumed.
qt <- cut(lg, quantile(lg, 0:5 / 5), include.lowest = TRUE)
QD <- mk(model.matrix(~ qt)[, -1, drop = FALSE], "granQ")
SP <- mk(ns(lg, df = 4), "granSpline")
LN <- mk(matrix(lg, ncol = 1), "granLinear")

ADJ <- c(
  list(`00_canonical` = NULL,
       `01_gran_linear` = LN,
       `02_gran_spline_df4` = SP,
       `03_gran_quintiles` = QD,
       `04_ILR4` = p$ILR),
  setNames(lapply(c(1, 2, 3, 5, 8, 12, 17),
                  function(k) mk(pcs[, 1:k, drop = FALSE], "fPC")),
           sprintf("05_fracPC%02d", c(1, 2, 3, 5, 8, 12, 17))))

can <- de_fit(p, "lithium")$tt
can_deg <- deg_ids(can, FDR)
cat(sprintf("canonical: %d DEG\n\n", length(can_deg)))

res <- do.call(rbind, lapply(names(ADJ), function(nm) {
  f <- de_fit(p, "lithium", adjust = ADJ[[nm]])
  dg <- deg_ids(f$tt, FDR)
  i <- match(can_deg, f$tt$gene); j <- match(can_deg, can$gene)
  r <- data.frame(
    adjustment = nm,
    df_spent = if (is.null(ADJ[[nm]])) 0L else ncol(ADJ[[nm]]),
    var_explained = if (grepl("fracPC", nm))
      round(ve[ncol(ADJ[[nm]])], 4) else NA_real_,
    n_deg = length(dg),
    survive = length(intersect(can_deg, dg)),
    prop_lost = round(1 - length(intersect(can_deg, dg)) / length(can_deg), 4),
    pi0 = round(pi0_storey(f$tt$P.Value), 4),
    shrink = round(median(abs(f$tt$logFC[i])) / median(abs(can$logFC[j])), 4),
    r_logFC = round(cor(can$logFC, f$tt$logFC), 4))
  cat(sprintf("  %-20s df %2d  DEG %4d  lost %.3f  pi0 %.3f  shrink %.3f\n",
              nm, r$df_spent, r$n_deg, r$prop_lost, r$pi0, r$shrink))
  r
}))

lin <- res$prop_lost[res$adjustment == "01_gran_linear"]
spl <- res$prop_lost[res$adjustment == "02_gran_spline_df4"]
qdm <- res$prop_lost[res$adjustment == "03_gran_quintiles"]
cat(sprintf("\nlinearity: linear %.3f | spline(df4) %.3f | quintiles %.3f\n",
            lin, spl, qdm))
cat(if (max(spl, qdm) - lin > 0.02)
  "=> a curved adjustment removes MORE. The linear figure under-corrects, so it is a lower bound.\n"
  else
  "=> curved and linear agree. The linear adjustment is adequate; no shape effect to report.\n")

mx <- res[which.max(res$prop_lost), ]
cat(sprintf("\nmost aggressive adjustment: %s, %.1f%% of canonical DEGs removed, %d survive\n",
            mx$adjustment, 100 * mx$prop_lost, mx$survive))

write.csv(res, file.path(OUT, "adjustment_dose_response.csv"), row.names = FALSE)
cat("\ndone.\n")
