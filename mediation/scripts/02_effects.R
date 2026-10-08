#!/usr/bin/env Rscript
# ==============================================================================
# 02_effects.R -- the primary mediation result.
#
# Runs the pre-specified model (PRESPEC.md section 1) on BOTH mediator
# definitions, because the written proposal and the implemented pipeline
# disagree about the tree and Ali asked for both:
#
#   lineage   5 lineages -> 4 balances   (what the pipeline built)
#   proposal  6 parts    -> 5 balances   (what the proposal text specifies)
#
# For each: the mediator models, the per-gene outcome models, the natural
# direct and indirect effects with delta-method standard errors, the four-way
# decomposition, the joint test of exposure-mediator interaction, and multiple
# testing correction across genes.
#
# Nothing here is chosen after seeing a result. Every threshold is in PRESPEC.md,
# which was written before this script ran.
# ==============================================================================

ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
HERE <- file.path(ROOT, "mediation")
setwd(HERE); set.seed(481)
suppressPackageStartupMessages({library(limma); library(edgeR)})
source("scripts/med.R")
source(file.path(ROOT, "de_v2/scripts/mtc.R"))
options(width = 150)
say <- function(...) cat(sprintf(...))
rule <- function(t) cat(sprintf("\n%s\n%s\n", t, strrep("=", nchar(t))))

VARIANTS <- list(
  lineage  = list(rds = "data/med_input.rds",
                  lab = "5 lineages -> 4 balances (implemented pipeline)"),
  proposal = list(rds = "data/med_input_proposal.rds",
                  lab = "6 parts -> 5 balances (proposal text)"))

dir.create("results/primary", showWarnings = FALSE, recursive = TRUE)
ALL <- list()

for (vn in names(VARIANTS)) {
  V <- VARIANTS[[vn]]
  rule(sprintf("VARIANT: %s  --  %s", vn, V$lab))
  d <- readRDS(V$rds)
  K <- ncol(d$M)
  cM <- centre(d$M); cC <- centre(d$C)
  x  <- d$lithium

  ## ---- stage 1: lithium -> cell mixture -------------------------------------
  # The proposal predicts lithium raises myeloid relative to lymphoid ONLY.
  # The pre-proposal feedback challenged the word "only", so all balances are
  # tested and all are reported.
  mf <- mediator_fit(cM$X, x, cC$X)
  zb <- mf$b1 / mf$b1_se
  pb <- 2 * pt(-abs(zb), df = mf$df)
  med_tab <- data.frame(
    balance  = colnames(d$M),
    contrast = if (!is.null(d$summary)) d$summary$contrast else
               c("gran+mono vs T+NK+B", "gran vs mono", "T vs NK+B", "NK vs B"),
    b1 = mf$b1, se = mf$b1_se, t = zb, p = pb,
    q_BH = p.adjust(pb, "BH"), stringsAsFactors = FALSE)
  cat("\nSTAGE 1  effect of lithium on each balance (the mediator models)\n")
  print(transform(med_tab, b1 = round(b1, 4), se = round(se, 4),
                  t = round(t, 2), p = signif(p, 3), q_BH = signif(q_BH, 3)),
        row.names = FALSE)
  say("\n  balances moved by lithium at BH 5%%: %s\n",
      paste(med_tab$balance[med_tab$q_BH < 0.05], collapse = ", "))

  ## ---- stage 2: the per-gene outcome models ---------------------------------
  Y <- d$v_main$E; W <- d$v_main$weights
  Xd <- med_design_full(x, cM$X, cC$X)
  stopifnot(qr(Xd)$rank == ncol(Xd))
  t0 <- Sys.time()
  of <- outcome_fit(Y, W, Xd, se = TRUE, moderate = TRUE)
  say("\nSTAGE 2  %d outcome models, design %d x %d, fitted in %.1f s\n",
      nrow(Y), nrow(Xd), ncol(Xd), as.numeric(difftime(Sys.time(), t0, units = "secs")))
  say("  residual df %d + empirical-Bayes prior df %.1f = %.1f total\n",
      of$df, of$df_prior, of$df_total)

  ## ---- the estimands --------------------------------------------------------
  eff <- med_effects(of, mf)
  tb  <- med_table(eff, rownames(Y))

  # multiple testing across genes, on the PRIMARY estimand
  adj <- adjust_panel(tb$NIE_p)
  tb$NIE_q_BH  <- adj$BH;  tb$NIE_q_BY <- adj$BY
  tb$NIE_q_bonf <- adj$bonferroni; tb$NIE_q_holm <- adj$holm
  tb$NIE_q_storey <- qvalue_storey(tb$NIE_p)
  tb$TE_q_BH  <- p.adjust(tb$TE_p,  "BH")
  tb$NDE_q_BH <- p.adjust(tb$NDE_p, "BH")

  ## ---- the interaction test, which the feedback demanded --------------------
  it <- interaction_test(Y, W, Xd, K)
  tb$INT_F <- it$F; tb$INT_p <- it$p; tb$INT_q_BH <- p.adjust(it$p, "BH")

  ## ---- proportion mediated, only where the total effect is real -------------
  sigTE <- tb$TE_q_BH < 0.05
  tb$PM_reported <- ifelse(sigTE, tb$PM, NA_real_)

  ## ---- report ---------------------------------------------------------------
  cat("\nSTAGE 3  gene-level counts\n")
  cnt <- data.frame(
    quantity = c("total effect TE", "natural indirect NIE", "natural direct NDE",
                 "interaction (joint 4/5 df)"),
    n_BH05 = c(sum(tb$TE_q_BH < 0.05), sum(tb$NIE_q_BH < 0.05),
               sum(tb$NDE_q_BH < 0.05), sum(tb$INT_q_BH < 0.05)),
    pct = NA_real_, stringsAsFactors = FALSE)
  cnt$pct <- round(100 * cnt$n_BH05 / nrow(tb), 1)
  print(cnt, row.names = FALSE)

  say("\n  genes with significant TE: %d\n", sum(sigTE))
  if (sum(sigTE) > 0) {
    pm <- tb$PM[sigTE]
    say("  proportion mediated among them: median %.3f  IQR %.3f-%.3f\n",
        median(pm), quantile(pm, .25), quantile(pm, .75))
    say("  share with PM > 0.50 (the pre-declared 'substantial' threshold): %.1f%%\n",
        100 * mean(pm > 0.5))
    say("  share with PM > 1 (over-mediation, direct effect opposing): %.1f%%\n",
        100 * mean(pm > 1))
    say("  share with PM < 0 (mediated and direct effects opposite in sign): %.1f%%\n",
        100 * mean(pm < 0))
  }

  cat("\n  four-way decomposition, mean across genes with significant TE:\n")
  fw <- data.frame(
    component = c("CDE  (neither mediation nor interaction)",
                  "INTref (interaction only)",
                  "INTmed (mediation and interaction)",
                  "PIE  (pure mediation)"),
    mean_abs = c(mean(abs(tb$CDE[sigTE])), mean(abs(tb$INTref[sigTE])),
                 mean(abs(tb$INTmed[sigTE])), mean(abs(tb$PIE[sigTE]))),
    stringsAsFactors = FALSE)
  fw$pct_of_total <- round(100 * fw$mean_abs / sum(fw$mean_abs), 1)
  fw$mean_abs <- round(fw$mean_abs, 5)
  print(fw, row.names = FALSE)

  ## S8 -- the proposal's original estimator, REFITTED without the interaction
  ## columns. Dropping t3 from the formula of the interaction model is NOT the
  ## same fit: it would merely reproduce PIE under a second name and make any
  ## "cost of assuming no interaction" tautological. See med.R prod_of_coef().
  pc <- prod_of_coef(Y, W, x, cM$X, cC$X, mf$b1)
  tb$NIE_prod <- pc$NIE_prod; tb$NDE_prod <- pc$NDE_prod
  cat("\n  S8  product-of-coefficients (the proposal's original estimator) vs ours:\n")
  say("    median NIE, interaction-allowing : %.5f\n", median(tb$NIE))
  say("    median NIE, product-of-coeffs    : %.5f\n", median(tb$NIE_prod))
  say("    median absolute difference       : %.5f (%.1f%% of median |NIE|)\n",
      median(abs(tb$NIE - tb$NIE_prod)),
      100 * median(abs(tb$NIE - tb$NIE_prod)) / median(abs(tb$NIE)))
  say("    genes where the two disagree in SIGN: %d (%.2f%%)\n",
      sum(sign(tb$NIE) != sign(tb$NIE_prod)),
      100 * mean(sign(tb$NIE) != sign(tb$NIE_prod)))
  say("    NOTE |NIE - PIE| = INTmed by construction (%.5f); that is a within-model\n",
      median(abs(tb$NIE - tb$PIE)))
  say("    decomposition term, NOT a comparison of two estimators.\n")

  ## ---- save -----------------------------------------------------------------
  ord <- order(tb$NIE_p)
  write.csv(tb[ord, ], sprintf("results/primary/%s_gene_effects.csv", vn), row.names = FALSE)
  write.csv(med_tab, sprintf("results/primary/%s_mediator_models.csv", vn), row.names = FALSE)
  saveRDS(list(tb = tb, mf = mf, of_s2 = of$s2_use, df_total = of$df_total,
               med_tab = med_tab, K = K, sigTE = sigTE, label = V$lab),
          sprintf("results/primary/%s_fit.rds", vn))
  ALL[[vn]] <- tb
  say("\n  wrote results/primary/%s_gene_effects.csv (%d genes x %d columns)\n",
      vn, nrow(tb), ncol(tb))
}

## ---- do the two mediator definitions agree? ---------------------------------
rule("AGREEMENT BETWEEN THE TWO MEDIATOR DEFINITIONS")
a <- ALL$lineage; b <- ALL$proposal
stopifnot(identical(a$gene, b$gene))
say("  correlation of per-gene NIE : Pearson %.3f | Spearman %.3f\n",
    cor(a$NIE, b$NIE), cor(a$NIE, b$NIE, method = "spearman"))
say("  correlation of per-gene TE  : Pearson %.3f\n", cor(a$TE, b$TE))
say("  NIE agrees in sign for %.1f%% of genes\n", 100 * mean(sign(a$NIE) == sign(b$NIE)))
sa <- a$NIE_q_BH < 0.05; sb <- b$NIE_q_BH < 0.05
say("  significant NIE: lineage %d | proposal %d | both %d | Jaccard %.3f\n",
    sum(sa), sum(sb), sum(sa & sb), sum(sa & sb) / sum(sa | sb))
pa <- a$PM[a$TE_q_BH < 0.05]; pb2 <- b$PM[b$TE_q_BH < 0.05]
say("  median proportion mediated: lineage %.3f | proposal %.3f\n", median(pa), median(pb2))

saveRDS(ALL, "results/primary/all_variants.rds")
cat("\n02_effects complete.\n")
