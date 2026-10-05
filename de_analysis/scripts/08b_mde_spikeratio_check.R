#!/usr/bin/env Rscript
# ==============================================================================
# 08b_mde_spikeratio_check.R -- is the granulocyte MDE an artifact of how many
# genes 08 spiked?
#
# THE PROBLEM 08 REVEALED
#
# 08 returned a false-positive count among its UNSPIKED control genes of:
#
#     gran 26/200 (13.0%)   mono 2/200 (1.0%)   T 0/200 (0.0%)
#
# 13% of nulls called significant at BH FDR 0.05 is too many, and it appears in
# exactly the lineage whose result matters most -- granulocytes are where the
# lithium effect lives.
#
# A LIKELY CAUSE, AND IT IS OUR DESIGN, NOT bMIND
#
# 08 spikes three blocks of 200 genes and keeps one block of 200 as controls,
# so 600 of 800 genes (75%) carry a real effect. limma's eBayes borrows
# information ACROSS GENES to form its variance prior, and `trend = TRUE` fits
# that prior against average expression. When three quarters of the genes
# carry a spike, the prior is estimated mostly from spiked genes, and the few
# genuine nulls are judged against a reference that does not describe them.
#
# The whole-blood MDE in 03 did not have this problem: it spiked 1,200 of
# 12,368 genes, under 10%.
#
# THE CHECK
#
# Re-measure granulocytes with the ratio INVERTED -- 200 spiked against 600
# controls (25% spiked) -- one run per effect size so the blocks stay disjoint.
# If the false-positive rate falls and the power curve holds, the MDE from 08
# stands and only the FP figure was an artifact. If the power curve ALSO moves,
# then 08's granulocyte MDE is biased and the corrected one here replaces it.
#
# Either way the honest number is the one measured with nulls in the majority,
# because that is the situation a real analysis is in.
# ==============================================================================

suppressPackageStartupMessages({library(MIND); library(limma)})
HERE <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/de_analysis"
setwd(HERE); set.seed(481)
RES <- "results"
p <- readRDS("data/prep.rds")

CT      <- "gran"
NSPIKE  <- 200
NCTRL   <- 600
NCORE   <- 4
LEVELS  <- c(0.25, 0.50, 0.80, 1.00)
PARTS   <- file.path(RES, "mde_spikeratio_parts"); dir.create(PARTS, showWarnings = FALSE, recursive = TRUE)

s  <- p$sel$LI
mm <- p$meta[s$keep, , drop = FALSE]; mm$grp <- s$grp
mm$plate <- droplevels(factor(mm$plate)); mm$sex <- droplevels(factor(mm$sex))
cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]])) > 1, logical(1))]
X <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse = " + "))), data = mm)
X[, "grpcase"] <- X[sample(nrow(X)), "grpcase"]          # permute: measure the spike alone
y <- X[, "grpcase"]

bulk0 <- p$logtpm[, mm$title, drop = FALSE]
frac  <- p$lineage[mm$title, , drop = FALSE]
allg  <- sample(rownames(bulk0))
spk   <- allg[1:NSPIKE]
ctrl  <- allg[(NSPIKE + 1):(NSPIKE + NCTRL)]
keepg <- c(spk, ctrl)

cat(sprintf("== 08b_mde_spikeratio_check ==\nlineage %s | %d spiked vs %d controls (%.0f%% spiked)\n",
            CT, NSPIKE, NCTRL, 100 * NSPIKE / (NSPIKE + NCTRL)))
cat(sprintf("compare against 08, which spiked 600 of 800 (75%%)\n\n"))

run <- function(L) {
  B <- bulk0
  add <- L * frac[, CT] * y
  B[spk, ] <- B[spk, ] + matrix(rep(add, each = NSPIKE), nrow = NSPIKE)
  r <- bMIND(B[keepg, , drop = FALSE], frac = frac, ncore = NCORE)
  Y <- r$A[, CT, ]; W <- 1 / (r$SE[, CT, ]^2)
  fin <- is.finite(W) & W > 0; if (any(!fin)) W[!fin] <- min(W[fin]); W <- W / mean(W)
  tt <- topTable(eBayes(lmFit(Y[, mm$title], X, weights = W[, mm$title]), trend = TRUE),
                 coef = "grpcase", number = Inf, sort.by = "none")
  q <- setNames(p.adjust(tt$P.Value, "BH"), rownames(tt))
  data.frame(lineage = CT, log2FC = L, n_spiked = NSPIKE, n_control = NCTRL,
             detected = sum(q[spk] < 0.05, na.rm = TRUE),
             power = round(sum(q[spk] < 0.05, na.rm = TRUE) / NSPIKE, 3),
             false_pos = sum(q[ctrl] < 0.05, na.rm = TRUE),
             fp_rate = round(sum(q[ctrl] < 0.05, na.rm = TRUE) / NCTRL, 4),
             stringsAsFactors = FALSE) }

OUT <- list()
for (L in LEVELS) {
  f <- file.path(PARTS, sprintf("L_%03.0f.rds", L * 100))
  if (file.exists(f)) { OUT[[length(OUT) + 1]] <- readRDS(f); cat(sprintf("  L=%.2f from checkpoint\n", L)); next }
  t0 <- Sys.time(); r <- run(L); saveRDS(r, f); OUT[[length(OUT) + 1]] <- r
  cat(sprintf("  L=%.2f  power %.3f  FP %d/%d (%.2f%%)  [%.1f min]\n", L, r$power, r$false_pos,
              NCTRL, 100 * r$fp_rate, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  flush.console()
}
D <- do.call(rbind, OUT)
cat("\n=== granulocytes, 25% spiked ===\n"); print(D, row.names = FALSE)

mde80 <- if (all(D$power < 0.8)) NA_real_ else
         if (all(D$power >= 0.8)) min(D$log2FC) else
         approx(D$power, D$log2FC, xout = 0.8, ties = "ordered")$y
old <- read.csv(file.path(RES, "celltype_mde.csv"))
og <- old[old$lineage == CT, ]
mde80_old <- if (all(og$power < 0.8)) NA_real_ else approx(og$power, og$log2FC, xout = 0.8, ties = "ordered")$y

cat(sprintf("\n  MDE80 at 25%% spiked : %s log2FC\n", if (is.na(mde80)) "> largest tested" else sprintf("%.3f", mde80)))
cat(sprintf("  MDE80 at 75%% spiked : %s log2FC   (08)\n", if (is.na(mde80_old)) "> largest tested" else sprintf("%.3f", mde80_old)))
cat(sprintf("  mean FP rate at 25%% spiked: %.2f%%   (08 at 75%%: %.2f%%)\n",
            100 * mean(D$fp_rate), 100 * mean(og$false_pos / og$n_control)))
cat(sprintf("\n  abundance-only prediction: %.2f log2FC\n", 0.203 / mean(frac[, CT])))
cat(sprintf("  config meaningfulness floor: 0.5 log2FC -> granulocyte null %s\n",
            if (!is.na(mde80) && mde80 <= 0.5) "INTERPRETABLE" else "NOT INTERPRETABLE (underpowered)"))

write.csv(D, file.path(RES, "mde_spikeratio_check.csv"), row.names = FALSE)
cat("\nwrote results/mde_spikeratio_check.csv\n")
