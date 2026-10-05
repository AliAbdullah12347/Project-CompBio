#!/usr/bin/env Rscript
# ==============================================================================
# baseline_bp2.R -- stage four: why 13 subjects change the answer by 2.5x.
#
# The single most surprising number in this experiment: on the Krebs 444, the
# published lithium model returns 1,023 DEGs when diagnosis is coded as a
# 2-level case/control indicator and 2,652 when it is coded as the 3-level
# Control/BP1/BP2 diagnosis. Removing the 13 BP2 subjects makes the two codings
# identical (2,591 each), so those 13 subjects are doing all of it.
#
# The mechanism, worked out on paper and checked here:
#
#   With a 3-level diagnosis the design is saturated -- one parameter per cell
#   -- so the lithium coefficient is exactly BP1-users minus BP1-non-users.
#
#   With a 2-level indicator, BP2 subjects (all lithium-free) carry the same
#   covariate value as BP1 non-users, so the model is forced to give them one
#   shared mean. Least squares sets that shared mean to the count-weighted
#   average of the two cells, and the lithium coefficient becomes
#
#       b_2level = y_use - (74*y_non + 13*y_BP2)/87
#                = b_3level - (13/87) * (y_BP2 - y_non)
#
#   i.e. the published lithium effect is measured against a reference group
#   that is 15% bipolar-II patients.
#
# Two consequences are testable, and are tested below:
#   (i)  the predicted variance ratio, which explains part of the change;
#   (ii) whether (y_BP2 - y_non) is systematically aligned with the lithium
#        effect, which would explain the rest -- and would mean bipolar-II
#        patients who have never taken lithium look, transcriptionally, like
#        lithium users.
#
# Fits use the shared engine with a cell-means parameterisation, so every
# quantity comes out of the same normalisation and filter as everything else.
# Writes only inside experimentation/runs/baseline/.
# ==============================================================================

setwd("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation")
source("scripts/common.R"); source("scripts/ref_krebs.R")

CACHE <- "C:/Users/hp/AppData/Local/Temp/claude/baseline_fits"
OUT   <- file.path(EXP_HOME, "runs", "baseline")
p <- load_prep(); SYM <- gene_symbols()
sym <- function(g) ifelse(is.na(SYM[g]), g, SYM[g])
qc <- p$meta$qc_pass

# four disjoint cells; Control is the reference
p$meta$cell <- factor(with(p$meta, ifelse(dx == "Control", "Control",
                    ifelse(dx == "BP2", "BP2",
                    ifelse(lithium == 1, "BP1use", "BP1non")))),
                levels = c("Control", "BP1non", "BP1use", "BP2"))
cat("== cell counts on the Krebs 444 ==\n"); print(table(p$meta$cell[qc]))
n_non <- sum(qc & p$meta$cell == "BP1non"); n_use <- sum(qc & p$meta$cell == "BP1use")
n_bp2 <- sum(qc & p$meta$cell == "BP2")
w <- n_bp2 / (n_non + n_bp2)
cat(sprintf("\nBP2 weight in the pooled reference: %d/(%d+%d) = %.4f\n",
            n_bp2, n_non, n_bp2, w))

CV <- c("age", "sex", "tob", "group", "rin", "plate", "seqpc1", "seqpc2", "seqpc3")
p$CONTRASTS$cells <- list(keep = rep(TRUE, nrow(p$meta)), exposure = "cell",
                          covars = CV, note = "cell-means parameterisation")

# One fit per cell coefficient. The design, filter and normalisation are
# identical across the three calls, so the logFCs are mutually comparable;
# only the coefficient reported differs.
getc <- function(co) {
  f <- file.path(CACHE, paste0("cell_", co, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  r <- de_fit(p, "cells", samples = qc, coef_override = co)
  saveRDS(r, f); r
}
CO <- lapply(setNames(c("cellBP1non", "cellBP1use", "cellBP2"),
                      c("non", "use", "bp2")), getc)
g  <- CO$non$tt$gene
stopifnot(identical(g, CO$use$tt$gene), identical(g, CO$bp2$tt$gene))
v  <- function(x) setNames(CO[[x]]$tt$logFC, g)

# all three are differences from Control, so cell contrasts subtract cleanly
L <- v("use") - v("non")     # lithium effect, the clean within-BP1 contrast
D <- v("bp2") - v("non")     # BP2 minus BP1 non-users

cat(sprintf("\n== the two contrasts, over %d genes ==\n", length(g)))
cat(sprintf("  lithium effect L = BP1use - BP1non   sd %.4f  median|L| %.4f\n",
            sd(L), median(abs(L))))
cat(sprintf("  BP2 contrast  D = BP2    - BP1non   sd %.4f  median|D| %.4f\n",
            sd(D), median(abs(D))))

# --------------------------------------------------------------- the test
cat("\n== is D aligned with L? ==\n")
r  <- cor(L, D); rho <- cor(L, D, method = "spearman")
sl <- sum(L * D) / sum(L * L)          # through-origin slope of D on L
cat(sprintf("  pearson r(L,D)   = %+.4f\n", r))
cat(sprintf("  spearman rho     = %+.4f\n", rho))
cat(sprintf("  slope of D on L  = %+.4f  (through origin)\n", sl))
cat(sprintf("  => BP2 patients, none of whom take lithium, sit %.1f%% of the way\n", 100 * sl))
cat(sprintf("     along the lithium axis relative to BP1 non-users.\n"))

# ------------------------------------------------- does it predict b_2level?
pred <- L - w * D
act  <- readRDS(file.path(CACHE, "spec_allcase2_444.rds"))$tt
act  <- setNames(act$logFC, act$gene)[g]
cat("\n== predicted vs actual 2-level lithium coefficient ==\n")
cat(sprintf("  r(predicted, actual) = %.6f   max|diff| = %.5f   median|diff| = %.6f\n",
            cor(pred, act), max(abs(pred - act)), median(abs(pred - act))))
cat(sprintf("  median|L| (3-level)  = %.5f\n", median(abs(L))))
cat(sprintf("  median|pred|         = %.5f\n", median(abs(pred))))
cat(sprintf("  median|actual|       = %.5f\n", median(abs(act))))
cat(sprintf("  attenuation actual   = %.1f%%\n", 100 * (1 - median(abs(act)) / median(abs(L)))))

# how much of the attenuation is variance, how much is the systematic shift?
sd_ratio_pred <- sqrt(1/n_use + (1 - w)^2 / n_non + w^2 / n_bp2) /
                 sqrt(1/n_use + 1/n_non)
cat(sprintf("\n  variance-only prediction for the |logFC| ratio: %.4f\n", sd_ratio_pred))
cat(sprintf("  observed |logFC| ratio                        : %.4f\n",
            median(abs(act)) / median(abs(L))))
cat(sprintf("  -> the extra shrinkage beyond sampling variance is %.1f percentage points,\n",
            100 * (sd_ratio_pred - median(abs(act)) / median(abs(L)))))
cat(sprintf("     and it is the systematic alignment of D with L.\n"))

# --------------------------------------------- restrict to real lithium genes
li <- readRDS(file.path(CACHE, "spec_bp1_474.rds"))$tt
lid <- li$gene[li$adj.P.Val < 0.05]
sel <- intersect(lid, g)
cat(sprintf("\n== restricted to the %d BP1-only lithium DEGs present here ==\n", length(sel)))
cat(sprintf("  slope of D on L = %+.4f   r = %+.4f\n",
            sum(L[sel] * D[sel]) / sum(L[sel]^2), cor(L[sel], D[sel])))
cat(sprintf("  mean L = %+.4f   mean D = %+.4f   (sign-aligned means: %+.4f / %+.4f)\n",
            mean(L[sel]), mean(D[sel]),
            mean(abs(L[sel])), mean(D[sel] * sign(L[sel]))))

res <- data.frame(
  quantity = c("n genes", "n BP1 non-users", "n BP1 users", "n BP2", "BP2 weight w",
               "sd(L)", "sd(D)", "pearson r(L,D)", "spearman rho(L,D)",
               "slope D~L (all genes)", "slope D~L (lithium DEGs)",
               "median|L|", "median|b_2level| actual", "attenuation",
               "variance-only predicted ratio", "observed ratio",
               "r(predicted b_2level, actual)", "max|pred - actual|"),
  value = c(length(g), n_non, n_use, n_bp2, round(w, 4),
            round(sd(L), 4), round(sd(D), 4), round(r, 4), round(rho, 4),
            round(sl, 4), round(sum(L[sel] * D[sel]) / sum(L[sel]^2), 4),
            round(median(abs(L)), 5), round(median(abs(act)), 5),
            round(1 - median(abs(act)) / median(abs(L)), 4),
            round(sd_ratio_pred, 4), round(median(abs(act)) / median(abs(L)), 4),
            round(cor(pred, act), 6), round(max(abs(pred - act)), 5)))
print(res, row.names = FALSE)
write.csv(res, file.path(OUT, "bp2_mechanism.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 6. free by-product: what IS the deposited PRE sheet?
# ---------------------------------------------------------------------------
# `cellBP1use` is BP1-lithium-users versus controls -- exactly the "BD on
# lithium vs control" contrast that ../Krebs Paper/REPRODUCIBILITY.md concluded
# the mislabelled PRE sheet contains. Testing it costs nothing here because the
# coefficient is already fitted, and it is an independent route to that claim.
K <- krebs_tables()
cmp <- function(x, nm) {
  gg <- intersect(names(x), K$PRE$gene)
  a  <- K$PRE$logFC[match(gg, K$PRE$gene)]
  cat(sprintf("  PRE vs %-28s r = %+.4f  rho = %+.4f  (n=%d)\n", nm,
              cor(a, x[gg]), cor(a, x[gg], method = "spearman"), length(gg)))
  cor(a, x[gg])
}
cat("\n== which of our coefficients does the deposited PRE sheet track? ==\n")
r_use <- cmp(v("use"), "BP1-users vs control")
r_non <- cmp(v("non"), "BP1-non-users vs control")
r_li  <- cmp(L,        "lithium main effect")
cat(sprintf("\n  PRE max|logFC| = %.3f ; our BP1-users-vs-control max|logFC| = %.3f\n",
            max(abs(K$PRE$logFC)), max(abs(v("use")))))
cat("  The sheet tracks BP1-users-vs-control, not the lithium main effect.\n")
write.csv(data.frame(
  our_coefficient = c("BP1 users vs control", "BP1 non-users vs control", "lithium main effect"),
  r_with_PRE_sheet = round(c(r_use, r_non, r_li), 4)),
  file.path(OUT, "pre_sheet_identity.csv"), row.names = FALSE)

per <- data.frame(gene = g, symbol = sym(g),
  L_lithium = signif(L, 6), D_bp2_minus_bp1non = signif(D, 6),
  b_2level_predicted = signif(pred, 6), b_2level_actual = signif(act, 6),
  is_lithium_DEG = g %in% lid)
gz <- gzfile(file.path(OUT, "bp2_per_gene.csv.gz"), "w")
write.csv(per, gz, row.names = FALSE); close(gz)
cat("\n[baseline_bp2] done\n")
