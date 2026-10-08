#!/usr/bin/env Rscript
# ==============================================================================
# 05_sensitivity.R -- the pre-declared sensitivity analyses S1, S2, S5, S6, S9,
# S12 from PRESPEC.md section 2. Each re-runs the primary analysis with one
# thing changed, so the question is always "does the conclusion move?".
#
#   S1  the 350 LM22 signature genes, analysed separately   (circularity)
#   S2  a second, independent deconvolution as the mediator (Krebs' CIBERSORT)
#   S5  all 20 tobacco imputations, pooled by Rubin's rules
#   S6  zero-replacement delta swept over five values
#   S9  unadjusted, and age-only-adjusted
#   S12 unweighted outcome models, where total = direct + indirect is exact
#
# S3, S4 and S10 (calibration, E-value, power) are in 06_calibration.R.
# S7 and S8 (correction panel, product-of-coefficients) were produced in
# 02_effects.R and are summarised in the final report.
# ==============================================================================

ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "mediation")); set.seed(481)
suppressPackageStartupMessages({library(edgeR); library(limma)})
source("scripts/med.R")
options(width = 160)
say <- function(...) cat(sprintf(...))
rule <- function(t) cat(sprintf("\n%s\n%s\n", t, strrep("=", nchar(t))))

prep <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
dl <- readRDS("data/med_input.rds")        # lineage,  4 balances
dp <- readRDS("data/med_input_proposal.rds")  # proposal, 5 balances
OUT <- list()

## ---- one helper: run the whole primary on any (counts, M, C, x) -------------
run_variant <- function(CNT, M, C, x, weighted = TRUE, label = "") {
  cM <- centre(M); cC <- centre(C)
  Xd <- med_design_full(x, cM$X, cC$X)
  if (qr(Xd)$rank < ncol(Xd)) stop("rank deficient: ", label)
  dge <- normLibSizes(DGEList(CNT), method = "TMM")
  v <- voom(dge, Xd)
  W <- if (weighted) v$weights else matrix(1, nrow(v$E), ncol(v$E))
  mf <- mediator_fit(cM$X, x, cC$X)
  of <- outcome_fit(v$E, W, Xd, se = TRUE, moderate = TRUE)
  tb <- med_table(med_effects(of, mf), rownames(v$E))
  tb$NIE_q_BH <- p.adjust(tb$NIE_p, "BH"); tb$TE_q_BH <- p.adjust(tb$TE_p, "BH")
  sig <- tb$TE_q_BH < 0.05
  list(tb = tb, mf = mf,
       summ = data.frame(analysis = label, n_genes = nrow(tb),
         n_TE = sum(sig), n_NIE = sum(tb$NIE_q_BH < 0.05),
         median_PM = if (any(sig)) median(tb$PM[sig]) else NA_real_,
         median_NIE = median(tb$NIE), median_TE = median(tb$TE),
         b1_1 = mf$b1[1], b1_1_t = (mf$b1 / mf$b1_se)[1],
         stringsAsFactors = FALSE))
}

CNTl <- prep$counts[dl$gene_main, dl$meta$title, drop = FALSE]
CNTm <- prep$counts[dl$gene_mark, dl$meta$title, drop = FALSE]

## ---- S1: the circularity check ---------------------------------------------
rule("S1  LM22 signature genes vs the rest -- how much circularity is there?")
# The cell fractions were estimated FROM the expression matrix using these 547
# genes (350 survive our filter). For them the mediator is close to a
# deterministic function of the outcome, so a mediated proportion is partly
# built in. The gap between the two sets MEASURES that, rather than assuming it
# away.
s1_main <- run_variant(CNTl, dl$M, dl$C, dl$lithium, label = "primary (12,018 non-marker)")
s1_mark <- run_variant(CNTm, dl$M, dl$C, dl$lithium, label = "LM22 markers (350)")
S1 <- rbind(s1_main$summ, s1_mark$summ)
print(S1[, c("analysis", "n_genes", "n_TE", "n_NIE", "median_PM", "median_NIE")], row.names = FALSE)
say("\n  median proportion mediated: non-markers %.3f vs markers %.3f  (ratio %.2fx)\n",
    S1$median_PM[1], S1$median_PM[2], S1$median_PM[2] / S1$median_PM[1])
say("  -> %s\n", if (S1$median_PM[2] > S1$median_PM[1])
  "markers show MORE mediation, as circularity predicts; the 12,018 figure is the defensible one" else
  "markers do NOT show more mediation; circularity is not visibly inflating the estimate")
OUT$S1 <- S1

## ---- S2: a second, independent deconvolution -------------------------------
rule("S2  does the result survive a different deconvolution of the same samples?")
# Krebs et al. deposited CIBERSORT fractions for these same subjects, produced
# by different software, different settings and different operators. Running
# them through the identical lineage aggregation and ILR construction gives a
# genuinely separate measurement of the mediator.
rel <- readRDS("results/reliability.rds")
ILRk <- rel$ILR_kr; rownames(ILRk) <- prep$meta$title
Mk <- ILRk[dl$meta$title, , drop = FALSE]; colnames(Mk) <- colnames(dl$M)
s2 <- run_variant(CNTl, Mk, dl$C, dl$lithium, label = "Krebs CIBERSORT mediator")
S2 <- rbind(s1_main$summ, s2$summ)
print(S2[, c("analysis", "n_TE", "n_NIE", "median_PM", "b1_1", "b1_1_t")], row.names = FALSE)
say("\n  per-gene NIE correlation between the two deconvolutions: %.3f\n",
    cor(s1_main$tb$NIE, s2$tb$NIE))
say("  median PM: ours %.3f vs Krebs %.3f\n", S2$median_PM[1], S2$median_PM[2])
OUT$S2 <- S2

## ---- S6: zero replacement ---------------------------------------------------
rule("S6  does the zero-replacement constant change anything?")
# The feedback asked directly: "CIBERSORT returns exact zeros, and log-ratios
# are undefined at zero. How will you replace zeros, and does that choice change
# the ratios?" For the lineage mediator the answer is provable: there are no
# zeros, so the rule never fires. For the proposal's finer tree it fires on 15%
# of samples, so it is swept.
nz_lin <- sum(rel$frac_cx == 0)
say("  lineage mediator: %d zeros among %d cells -> replacement never fires, delta is irrelevant\n",
    nz_lin, length(rel$frac_cx))
zrepl <- function(X, delta) {
  lim <- apply(X, 2, function(v) { v <- v[v > 0]; if (!length(v)) 1e-8 else min(v) }) * delta
  z <- X <= 0; ins <- rowSums(sweep(z, 2, lim, "*"))
  X[z] <- matrix(rep(lim, each = nrow(X)), nrow(X))[z]
  X[!z] <- (X * (1 - ins))[!z]; X / rowSums(X)
}
S6 <- NULL
for (dlt in c(0.1, 0.3, 0.5, 0.65, 0.8)) {
  Mz <- log(zrepl(dp$parts, dlt)) %*% t(dp$psi)
  colnames(Mz) <- colnames(dp$M); rownames(Mz) <- rownames(dp$parts)
  r <- run_variant(prep$counts[dp$gene_main, dp$meta$title, drop = FALSE],
                   Mz[dp$meta$title, , drop = FALSE], dp$C, dp$lithium,
                   label = sprintf("proposal tree, delta = %.2f", dlt))
  S6 <- rbind(S6, r$summ)
}
print(S6[, c("analysis", "n_TE", "n_NIE", "median_PM", "median_NIE")], row.names = FALSE)
say("\n  median PM across delta: %.4f to %.4f (spread %.4f)\n",
    min(S6$median_PM), max(S6$median_PM), diff(range(S6$median_PM)))
say("  -> the zero-replacement constant is %s\n",
    if (diff(range(S6$median_PM)) < 0.02) "NOT load-bearing" else "load-bearing; report the sweep")
OUT$S6 <- S6

## ---- S9: covariate adjustment ----------------------------------------------
rule("S9  how much do the covariates matter? (age is the imbalanced one, p = 0.0007)")
S9 <- s1_main$summ; S9$analysis <- "full adjustment (primary)"
Cn <- dl$C[, "age", drop = FALSE]
S9 <- rbind(S9, run_variant(CNTl, dl$M, Cn, dl$lithium, label = "age only")$summ)
# unadjusted needs at least one column in C for the design builder; use a
# constant-free zero-column trick by passing a single centred dummy that is
# exactly collinear with nothing -- simplest is to pass age centred and then
# also report it, so instead we fit with NO covariates via a rank-safe route
Cz <- matrix(0, nrow(dl$C), 0)
S9 <- rbind(S9, tryCatch(run_variant(CNTl, dl$M, Cz, dl$lithium, label = "no covariates")$summ,
                         error = function(e) { say("  (unadjusted fit failed: %s)\n", conditionMessage(e)); NULL }))
print(S9[, c("analysis", "n_TE", "n_NIE", "median_PM", "b1_1", "b1_1_t")], row.names = FALSE)
OUT$S9 <- S9

## ---- S12: unweighted outcome models ----------------------------------------
rule("S12  unweighted outcome models, where total = direct + indirect is EXACT")
# The mediator models are unweighted (voom weights describe count noise in gene
# expression, and the lithium-to-cell-mix effect is one gene-independent fact).
# The outcome models are weighted, for comparability with the DE arm. The two
# use different projections, so the textbook identity total = direct + indirect
# holds only to about 5% here rather than exactly. Dropping the weights makes it
# exact; if the conclusions are the same either way, the mismatch is harmless.
s12 <- run_variant(CNTl, dl$M, dl$C, dl$lithium, weighted = FALSE,
                   label = "unweighted outcome models")
S12 <- rbind(s1_main$summ, s12$summ)
print(S12[, c("analysis", "n_TE", "n_NIE", "median_PM", "median_NIE")], row.names = FALSE)
say("\n  per-gene NIE correlation weighted vs unweighted: %.4f\n",
    cor(s1_main$tb$NIE, s12$tb$NIE))
OUT$S12 <- S12

## ---- S5: all 20 tobacco imputations, pooled by Rubin's rules ---------------
rule("S5  all 20 tobacco imputations, pooled by Rubin's rules")
# 30 subjects in the 474-sample cohort lack tobacco status and were multiply
# imputed (20 datasets). The primary analysis uses imputation 1. Rubin's rules
# combine all 20 so the imputation uncertainty enters the standard error instead
# of being ignored.
md <- read.csv(file.path(ROOT, "data/cohort_474/metadata_474_imputed.csv"),
               check.names = FALSE, stringsAsFactors = FALSE)
rownames(md) <- md$title
imp_cols <- sprintf("tobacco_imp_%02d", 1:20)
stopifnot(all(imp_cols %in% names(md)))
subj <- dl$meta$title
n_imputed_here <- sum(as.logical(md[subj, "tobacco_was_imputed"]))
say("  subjects in this analysis whose tobacco status was imputed: %d of %d\n",
    n_imputed_here, length(subj))

Qs <- Us <- NULL
for (mi in seq_along(imp_cols)) {
  Ci <- dl$C; Ci[, "tobacco"] <- as.numeric(md[subj, imp_cols[mi]])
  r <- run_variant(CNTl, dl$M, Ci, dl$lithium, label = sprintf("imp %02d", mi))
  Qs <- cbind(Qs, r$tb$NIE); Us <- cbind(Us, r$tb$NIE_se^2)
  if (mi %% 5 == 0) say("    ... imputation %d of 20\n", mi)
}
m <- ncol(Qs)
Qbar <- rowMeans(Qs); Ubar <- rowMeans(Us)
Bvar <- apply(Qs, 1, var)
Tvar <- Ubar + (1 + 1 / m) * Bvar
lambda <- (1 + 1 / m) * Bvar / Tvar              # fraction of missing information
nu_old <- (m - 1) / pmax(lambda, 1e-12)^2
dfc <- 226 - ncol(med_design_full(dl$lithium, centre(dl$M)$X, centre(dl$C)$X))
nu_obs <- (dfc + 1) / (dfc + 3) * dfc * (1 - lambda)
nu_BR <- 1 / (1 / nu_old + 1 / nu_obs)           # Barnard & Rubin (1999)
tstat <- Qbar / sqrt(Tvar)
pval <- 2 * pt(-abs(tstat), df = nu_BR)
S5 <- data.frame(gene = dl$gene_main, NIE_pooled = Qbar, se_pooled = sqrt(Tvar),
                 fmi = lambda, df = nu_BR, p = pval, q_BH = p.adjust(pval, "BH"),
                 stringsAsFactors = FALSE)
say("\n  pooled vs imputation-1: NIE correlation %.5f | median |difference| %.5f\n",
    cor(S5$NIE_pooled, s1_main$tb$NIE), median(abs(S5$NIE_pooled - s1_main$tb$NIE)))
say("  fraction of missing information: median %.4f, max %.4f\n",
    median(S5$fmi), max(S5$fmi))
say("  SE inflation from imputation uncertainty: median %.3f%%\n",
    100 * (median(S5$se_pooled / s1_main$tb$NIE_se) - 1))
say("  genes with pooled NIE at BH 5%%: %d (imputation 1 alone gave %d)\n",
    sum(S5$q_BH < 0.05), sum(s1_main$tb$NIE_q_BH < 0.05))
write.csv(S5, "results/sens_S5_rubin_pooled.csv", row.names = FALSE)
OUT$S5 <- S5

## ---- save -------------------------------------------------------------------
dir.create("results/sens", showWarnings = FALSE, recursive = TRUE)
for (nm in c("S1", "S2", "S6", "S9", "S12"))
  write.csv(OUT[[nm]], sprintf("results/sens/%s.csv", nm), row.names = FALSE)
saveRDS(list(OUT = OUT, s1_main = s1_main$tb, s1_mark = s1_mark$tb,
             s2 = s2$tb, s12 = s12$tb), "results/sens/sensitivity.rds")
cat("\n05_sensitivity complete.\n")
