#!/usr/bin/env Rscript
# ==============================================================================
# 02c_identity.R -- permanent record of the exact algebraic check behind the one
# thing in this arm that initially looked like a bug.
#
# THE PUZZLE
#
# In a linear model the decomposition total = direct + indirect is an EXACT
# algebraic identity, so the mediation total effect should equal the lithium
# coefficient from the plain marginal model Y ~ lithium + C to machine
# precision. It did not: the gap reached 4.3e-02.
#
# THE RESOLUTION
#
# The identity requires the mediator model and the outcome model to use the SAME
# projection. Ours deliberately do not:
#
#   outcome models   WEIGHTED least squares, with voom precision weights, for
#                    comparability with the differential-expression arm.
#   mediator models  UNWEIGHTED. voom weights describe count-level sampling
#                    noise in GENE EXPRESSION. They say nothing about how
#                    precisely a subject's cell fractions were estimated, and
#                    the lithium-to-cell-mix effect is a single, gene-independent
#                    quantity. Weighting it per gene would claim that lithium
#                    shifts the myeloid/lymphoid balance by 0.1558 when you look
#                    at one gene and 0.1949 when you look at another, which is
#                    not a coherent statement about blood.
#
# This script proves that is the entire explanation: refit the mediator model
# with each gene's weights, so the projections match, and the identity returns
# to machine precision. The unweighted choice is kept, and S12 in
# 05_sensitivity.R confirms the conclusions do not depend on it (per-gene NIE
# correlation 0.9984, median proportion mediated 0.396 vs 0.400).
# ==============================================================================

ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "mediation")); set.seed(481)
source("scripts/med.R"); options(width = 150)
say <- function(...) cat(sprintf(...))

d <- readRDS("data/med_input.rds"); K <- ncol(d$M)
cM <- centre(d$M); cC <- centre(d$C); x <- d$lithium
Y <- d$v_main$E; W <- d$v_main$weights
G <- 400; gs <- sort(sample.int(nrow(Y), G))
Xni <- cbind(1, x, cM$X, cC$X)      # outcome model WITHOUT the interaction
Xmg <- cbind(1, x, cC$X)            # marginal model
Dm  <- cbind(1, x, cC$X)            # mediator design
Bu  <- qr.coef(qr(Dm), cM$X)        # unweighted mediator fit (what we use)

gapU <- gapW <- numeric(G)
for (j in seq_len(G)) {
  g <- gs[j]; sw <- sqrt(W[g, ])
  bn <- qr.coef(qr(Xni * sw), Y[g, ] * sw)
  bm <- qr.coef(qr(Xmg * sw), Y[g, ] * sw)
  te_marginal <- bm[2]
  gapU[j] <- (bn[2] + sum(bn[2 + (1:K)] * Bu[2, ])) - te_marginal
  Bw <- qr.coef(qr(Dm * sw), cM$X * sw)          # mediator fit with THIS gene's weights
  gapW[j] <- (bn[2] + sum(bn[2 + (1:K)] * Bw[2, ])) - te_marginal
}
say("== 02c_identity ==  %d sampled genes, no-interaction outcome model\n\n", G)
say("  mediator model UNWEIGHTED (our choice)          max |gap| = %.3e\n", max(abs(gapU)))
say("  mediator model weighted per gene (projections match) max |gap| = %.3e\n", max(abs(gapW)))
say("\n  verdict: %s\n", if (max(abs(gapW)) < 1e-10)
  "the identity holds EXACTLY once the projections match. The implementation is\n           correct; the gap is the deliberate weighting choice, documented above." else
  "**the identity does not recover -- there is a genuine bug**")
stopifnot(max(abs(gapW)) < 1e-10)
say("\n  unweighted mediator b1 : %s\n", paste(sprintf("%+.4f", Bu[2, ]), collapse = "  "))
Bw1 <- qr.coef(qr(Dm * sqrt(W[gs[1], ])), cM$X * sqrt(W[gs[1], ]))
say("  gene-weighted b1       : %s   <- differs per gene, which is why we do not do this\n",
    paste(sprintf("%+.4f", Bw1[2, ]), collapse = "  "))
writeLines(sprintf("unweighted gap %.3e | matched-projection gap %.3e | PASS",
                   max(abs(gapU)), max(abs(gapW))), "results/identity_check.txt")
