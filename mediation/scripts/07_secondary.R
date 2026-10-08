#!/usr/bin/env Rscript
# ==============================================================================
# 07_secondary.R -- T1 and T2 from PRESPEC.md section 3. EXPLORATORY, and
# labelled as such everywhere they appear.
#
# The proposal's rule: "We specify a single primary model here and designate
# every other configuration exploratory, to be reported as such and never as a
# confirmatory result." These are those configurations, and the list was closed
# in PRESPEC.md before any of them ran.
#
#   T1  age as the exposure, cell mix as the mediator, same outcome set
#   T2  moderated mediation: does age's mediated effect differ between lithium
#       users and non-users? (the proposal's hypothesis, after Salarda et al.
#       2021 on lithium and ageing)
#
# T2 is reported WITH its power, because the feedback asked directly: "With 74
# non-users, do you have power to detect a difference in mediated effects?"
# ==============================================================================

ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "mediation")); set.seed(481)
suppressPackageStartupMessages({library(limma); library(edgeR)})
source("scripts/med.R"); options(width = 160)
say <- function(...) cat(sprintf(...))
rule <- function(t) cat(sprintf("\n%s\n%s\n", t, strrep("=", nchar(t))))

prep <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
d <- readRDS("data/med_input.rds"); K <- ncol(d$M)
DELTA <- 10                                  # contrast: 10 years of age

## ---- natural effects for a CONTINUOUS exposure ------------------------------
# For exposure x compared with x* = 0 (age is centred, so 0 is the mean age):
#   NDE = D (t1 + sum_k t3_k b0_k)
#   NIE = D sum_k (t2_k + D t3_k) b1_k
# which reduces to the binary formulas at D = 1. Gradients follow directly and
# the same block-diagonal delta method applies.
cont_effects <- function(of, mf, D) {
  Kk <- length(mf$b1); ix <- med_idx(Kk)
  t1 <- of$theta[, ix$iX]
  T2 <- of$theta[, ix$iM, drop = FALSE]; T3 <- of$theta[, ix$iI, drop = FALSE]
  NDE <- D * (t1 + T3 %*% mf$b0)
  NIE <- D * ((T2 + D * T3) %*% mf$b1)
  mkA <- function(xx = 0, m = numeric(Kk), i = numeric(Kk)) {
    a <- numeric(of$nsub); a[ix$jX] <- xx; a[ix$jM] <- m; a[ix$jI] <- i; a }
  qt <- function(a) as.vector(of$Vsub %*% as.vector(outer(a, a))) * of$s2_use
  Vb <- mf$Vbeta
  V11 <- Vb[1:Kk, 1:Kk, drop = FALSE]; V22 <- Vb[Kk + 1:Kk, Kk + 1:Kk, drop = FALSE]
  qb <- function(B, V) rowSums((B %*% V) * B)
  se_NDE <- sqrt(qt(mkA(xx = D, i = D * mf$b0)) + qb(D * T3, V11))
  Sb <- D * (T2 + D * T3)
  se_NIE <- sqrt(qt(mkA(m = D * mf$b1, i = D^2 * mf$b1)) + qb(Sb, V22))
  data.frame(NDE = as.vector(NDE), NDE_se = se_NDE,
             NIE = as.vector(NIE), NIE_se = se_NIE,
             TE = as.vector(NDE + NIE), stringsAsFactors = FALSE)
}

run_age <- function(keep, label) {
  mm <- d$meta[keep, , drop = FALSE]
  age <- mm$age
  # exposure = age; lithium moves INTO the covariate set (except within strata,
  # where it is constant and must be dropped)
  Cdf <- d$C[keep, setdiff(colnames(d$C), "age"), drop = FALSE]
  lith <- d$lithium[keep]
  if (length(unique(lith)) > 1) Cdf <- cbind(Cdf, lithium = lith)
  keepcol <- apply(Cdf, 2, function(v) length(unique(v)) > 1)
  Cdf <- Cdf[, keepcol, drop = FALSE]
  M <- d$M[keep, , drop = FALSE]
  cA <- age - mean(age); cM <- centre(M); cC <- centre(Cdf)
  Xd <- med_design_full(cA, cM$X, cC$X)
  if (qr(Xd)$rank < ncol(Xd)) stop("rank deficient in ", label)
  CNT <- prep$counts[d$gene_main, mm$title, drop = FALSE]
  v <- voom(normLibSizes(DGEList(CNT), method = "TMM"), Xd)
  mf <- mediator_fit(cM$X, cA, cC$X)
  of <- outcome_fit(v$E, v$weights, Xd, se = TRUE, moderate = TRUE)
  e <- cont_effects(of, mf, DELTA)
  e$gene <- rownames(v$E)
  e$NIE_p <- 2 * pnorm(-abs(e$NIE / e$NIE_se))
  e$NIE_q <- p.adjust(e$NIE_p, "BH")
  e$TE_p <- 2 * pnorm(-abs(e$TE / sqrt(e$NDE_se^2 + e$NIE_se^2)))
  e$TE_q <- p.adjust(e$TE_p, "BH")
  e$PM <- e$NIE / e$TE
  list(e = e, mf = mf, n = sum(keep), label = label)
}

## =============================================================================
## T1  age as the exposure
## =============================================================================
rule("T1  EXPLORATORY -- age as exposure, cell mix as mediator (per 10 years)")
t1 <- run_age(rep(TRUE, nrow(d$M)), "all 226 bipolar I")
say("  n = %d, contrast = %d years of age, lithium adjusted as a covariate\n", t1$n, DELTA)
say("\n  effect of age on each balance (the mediator models):\n")
tt <- t1$mf$b1 / t1$mf$b1_se
for (k in seq_len(K))
  say("    %-4s b1 = %+.5f per year  t = %+.2f  p = %.4f\n", colnames(d$M)[k],
      t1$mf$b1[k], tt[k], 2 * pt(-abs(tt[k]), t1$mf$df))
sig <- t1$e$TE_q < 0.05
say("\n  genes with a significant total effect of age: %d of %d\n", sum(sig), nrow(t1$e))
say("  genes with NIE at BH 5%%: %d\n", sum(t1$e$NIE_q < 0.05))
if (any(sig)) say("  median proportion mediated among them: %.3f (IQR %.3f-%.3f)\n",
                  median(t1$e$PM[sig]), quantile(t1$e$PM[sig], .25), quantile(t1$e$PM[sig], .75))
say("  median |NIE| per 10 years: %.5f | median |NDE|: %.5f\n",
    median(abs(t1$e$NIE)), median(abs(t1$e$NDE)))
write.csv(t1$e, "results/sens/T1_age_exposure.csv", row.names = FALSE)

## =============================================================================
## T2  moderated mediation
## =============================================================================
rule("T2  EXPLORATORY -- does age's mediated effect differ by lithium status?")
# The proposal hypothesises that "the cell-mediated effect of age on gene
# expression increases in lithium-users". Fitting the age mediation separately
# within each stratum keeps the comparison clean: the two groups are disjoint,
# so the difference in NIE has variance equal to the sum of the two variances.
tu <- run_age(d$lithium == 1, "lithium users")
tn <- run_age(d$lithium == 0, "lithium non-users")
say("  lithium users n = %d | non-users n = %d\n", tu$n, tn$n)

say("\n  effect of age on balance 1 (myeloid vs lymphoid), by stratum:\n")
for (r in list(tu, tn)) {
  tt <- r$mf$b1[1] / r$mf$b1_se[1]
  say("    %-18s b1 = %+.5f per year  t = %+.2f\n", r$label, r$mf$b1[1], tt)
}

stopifnot(identical(tu$e$gene, tn$e$gene))
dif <- tu$e$NIE - tn$e$NIE
sed <- sqrt(tu$e$NIE_se^2 + tn$e$NIE_se^2)
z   <- dif / sed; pz <- 2 * pnorm(-abs(z)); qz <- p.adjust(pz, "BH")
say("\n  difference in mediated effect (users minus non-users), per 10 years:\n")
say("    median difference %+.5f | median SE %.5f | median |z| %.2f\n",
    median(dif), median(sed), median(abs(z)))
say("    genes with a significant difference at BH 5%%: %d\n", sum(qz < 0.05))
say("    share where the users' mediated effect is LARGER in absolute size: %.1f%%\n",
    100 * mean(abs(tu$e$NIE) > abs(tn$e$NIE)))

rule("T2b  the power question the feedback asked")
z95 <- qnorm(0.975); z80 <- qnorm(0.80)
mde <- (z95 + z80) * median(sed)
mde_bh <- (qnorm(1 - 0.05 / (2 * nrow(t1$e))) + z80) * median(sed)
say("  smallest detectable difference in mediated effect, 80%% power:\n")
say("    at an uncorrected 5%% threshold : %.5f\n", mde)
say("    at a Bonferroni threshold      : %.5f\n", mde_bh)
say("  observed median |difference|     : %.5f\n", median(abs(dif)))
say("  -> the design can detect a difference only if it is %.1fx larger than the\n",
    mde / median(abs(dif)))
say("     one observed. With 74 non-users this comparison is UNDERPOWERED, which\n")
say("     is the honest answer to the question, and it was knowable in advance.\n")
say("     A null here is uninformative and is reported as such, not as evidence\n")
say("     that the mediated effect is equal across strata.\n")

T2 <- data.frame(gene = tu$e$gene, NIE_users = tu$e$NIE, NIE_nonusers = tn$e$NIE,
                 difference = dif, se = sed, z = z, p = pz, q_BH = qz,
                 stringsAsFactors = FALSE)
write.csv(T2, "results/sens/T2_moderated_mediation.csv", row.names = FALSE)
saveRDS(list(T1 = t1, T2 = T2, users = tu, nonusers = tn,
             mde = mde, mde_bh = mde_bh), "results/sens/secondary.rds")
cat("\n07_secondary complete.\n")
