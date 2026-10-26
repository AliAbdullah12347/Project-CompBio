#!/usr/bin/env Rscript
# ==============================================================================
# 12_diversity_alpha.R -- alpha diversity: does lithium change it, and is the
# change mediated by the cell mix? (PRESPEC.md section 8, tests A1-A3, D-B, D-C)
#
# Alpha diversity is ONE number per person, so the arm's validated estimator
# (med.R) applies unchanged with a single outcome. Ordinary least squares, no
# voom weights: diversity is a per-person summary, not a count.
#
#   A1  total effect of the exposure on diversity          marginal t-test
#   A2  natural indirect effect through the cell mix       bootstrap CI (B=10,000)
#                                                          + permutation (B=10,000)
#   A3  exposure x cell-mix interaction                     F-test, K df
#
# Outcomes are modelled on the LOG scale of the Hill numbers, so an effect of d
# means the effective number of genes changes by a factor exp(d):
#   log Hill0 (richness), H = log Hill1 (Shannon, PRIMARY), log Hill2.
#
# Configurations: lithium within bipolar I on both trees (primary + sensitivity);
# bipolar-off-lithium vs control on the lineage tree (EXPLORATORY).
# ==============================================================================

ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "mediation")); set.seed(481)
source("scripts/med.R"); options(width = 170)
say  <- function(...) cat(sprintf(...))
rule <- function(t) cat(sprintf("\n%s\n%s\n", t, strrep("=", nchar(t))))
B <- 10000L

prep <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
dv   <- readRDS("data/diversity_input.rds")
A    <- dv$alpha; rownames(A) <- A$title

## ---- outcomes ----------------------------------------------------------------
OUT <- c(
  H_primary          = "PRIMARY: Shannon, protein-coding non-LM22, TPM",
  logHill0           = "richness (all protein-coding non-LM22 genes)",
  logHill2           = "inverse Simpson",
  H_with_markers     = "sens: LM22 markers included",
  H_no_mito          = "sens: mitochondrial genes excluded",
  H_counts           = "sens: read-count proportions",
  H_all_filtered     = "sens: all filtered genes incl. ncRNA (expected to be technical)",
  H_pw_hallmark      = "PATHWAY primary: Hallmark, fractional counting",
  H_pw_hallmark_noLM22 = "pathway sens: Hallmark without LM22 markers",
  H_pw_kegg          = "pathway sens: KEGG legacy")
A$logHill0 <- log(A$Hill0_richness)
A$logHill2 <- log(A$Hill2_primary)

## ---- configurations ----------------------------------------------------------
build_cfg <- function(sel, tree, label) {
  s <- prep$sel[[sel]]; mm <- prep$meta[s$keep, , drop = FALSE]
  x <- as.integer(s$grp == "case")
  if (sel == "LI") {
    d <- readRDS(if (tree == "lineage") "data/med_input.rds" else "data/med_input_proposal.rds")
    stopifnot(identical(rownames(d$M), mm$title))
    M <- d$M; C <- d$C
  } else {
    stopifnot(tree == "lineage")
    M <- prep$ILR[mm$title, , drop = FALSE]
    mm$sex <- droplevels(factor(mm$sex)); mm$plate <- droplevels(factor(mm$plate))
    mm$group <- droplevels(factor(mm$group))
    C <- model.matrix(~ age + sex + tobacco + rin + plate + seqpc1 + seqpc2 + seqpc3 + group,
                      data = mm)[, -1, drop = FALSE]
  }
  C <- cbind(C, log_libsize = A[mm$title, "log_libsize"])
  stopifnot(!any(is.na(C)), !any(is.na(M)), qr(cbind(1, x, C))$rank == ncol(C) + 2)
  list(label = label, sel = sel, tree = tree, titles = mm$title, x = x, M = M, C = C)
}
CFG <- list(
  build_cfg("LI",  "lineage",  "lithium | lineage tree (PRIMARY)"),
  build_cfg("LI",  "proposal", "lithium | proposal tree"),
  build_cfg("BPD", "lineage",  "bipolar vs control | lineage tree (EXPLORATORY)"))

## ---- fast fit of NIE / NDE / TE for many outcomes at once --------------------
fast_eff <- function(x, M, C, Y) {
  cM <- sweep(M, 2, colMeans(M)); cC <- sweep(C, 2, colMeans(C)); K <- ncol(M)
  D <- cbind(1, x, cC); Bm <- qr.coef(qr(D), cM)
  Xd <- cbind(1, x, cM, cM * x, cC); Th <- qr.coef(qr(Xd), Y)
  t1 <- Th[2, ]; T2 <- Th[2 + seq_len(K), , drop = FALSE]; T3 <- Th[2 + K + seq_len(K), , drop = FALSE]
  NIE <- colSums((T2 + T3) * Bm[2, ]); NDE <- t1 + colSums(T3 * Bm[1, ])
  rbind(NIE = NIE, NDE = NDE, TE = NIE + NDE)
}

RES <- list(); BOOT <- list()
for (cfg in CFG) {
  rule(cfg$label)
  Y <- as.matrix(A[cfg$titles, names(OUT)]); n <- nrow(Y); K <- ncol(cfg$M)
  say("  n = %d (%d exposed) | mediator %d balances | %d covariate columns\n",
      n, sum(cfg$x), K, ncol(cfg$C))

  ## point estimates and delta-method SEs from the validated estimator
  cM <- centre(cfg$M); cC <- centre(cfg$C)
  mf <- mediator_fit(cM$X, cfg$x, cC$X)
  Xd <- med_design_full(cfg$x, cM$X, cC$X)
  of <- outcome_fit(t(Y), matrix(1, ncol(Y), n), Xd, se = TRUE, moderate = FALSE)
  ef <- med_effects(of, mf)
  # cross-check: the fast path must reproduce med.R exactly
  fe <- fast_eff(cfg$x, cfg$M, cfg$C, Y)
  stopifnot(max(abs(fe["NIE", ] - ef$est$NIE)) < 1e-10, max(abs(fe["TE", ] - ef$est$TE)) < 1e-10)

  ## A1: total effect from the marginal model
  a1 <- t(apply(Y, 2, function(y) {
    f <- summary(lm(y ~ cfg$x + cfg$C)); c(est = f$coef[2, 1], se = f$coef[2, 2], p = f$coef[2, 4]) }))

  ## A3: interaction F-test
  a3 <- apply(Y, 2, function(y) {
    f0 <- lm(y ~ cfg$x + cM$X + cC$X); f1 <- lm(y ~ cfg$x + cM$X + I(cM$X * cfg$x) + cC$X)
    anova(f0, f1)$`Pr(>F)`[2] })

  ## paths for joint-significance reading: a = exposure -> mediator (multivariate),
  ## b = mediator -> outcome given exposure (K-df F-test)
  pa <- anova(lm(cM$X ~ cC$X), lm(cM$X ~ cfg$x + cC$X), test = "Pillai")$`Pr(>F)`[2]
  pb <- apply(Y, 2, function(y)
    anova(lm(y ~ cfg$x + cC$X), lm(y ~ cfg$x + cM$X + cC$X))$`Pr(>F)`[2])

  ## A2: bootstrap and permutation of NIE (and NDE, TE, PM from the same draws)
  set.seed(481)
  bt <- array(NA_real_, c(B, 3, ncol(Y)), dimnames = list(NULL, c("NIE", "NDE", "TE"), colnames(Y)))
  for (b in seq_len(B)) {
    i <- sample.int(n, n, replace = TRUE)
    bt[b, , ] <- fast_eff(cfg$x[i], cfg$M[i, , drop = FALSE], cfg$C[i, , drop = FALSE], Y[i, , drop = FALSE])
  }
  set.seed(4812)
  pn <- matrix(NA_real_, B, ncol(Y), dimnames = list(NULL, colnames(Y)))
  for (b in seq_len(B)) pn[b, ] <- fast_eff(sample(cfg$x), cfg$M, cfg$C, Y)["NIE", ]

  q <- function(v, p) quantile(v, p, names = FALSE)
  r <- data.frame(config = cfg$label, outcome = names(OUT), description = unname(OUT),
    n = n, sd_outcome = apply(Y, 2, sd),
    TE_marginal = a1[, "est"], TE_p = a1[, "p"],
    TE = ef$est$TE, TE_lo = apply(bt[, "TE", ], 2, q, .025), TE_hi = apply(bt[, "TE", ], 2, q, .975),
    NDE = ef$est$NDE, NDE_lo = apply(bt[, "NDE", ], 2, q, .025), NDE_hi = apply(bt[, "NDE", ], 2, q, .975),
    NIE = ef$est$NIE, NIE_se_delta = ef$se$NIE,
    NIE_lo = apply(bt[, "NIE", ], 2, q, .025), NIE_hi = apply(bt[, "NIE", ], 2, q, .975),
    NIE_p_perm = sapply(seq_len(ncol(Y)), function(j) (1 + sum(abs(pn[, j]) >= abs(ef$est$NIE[j]))) / (1 + B)),
    p_path_a = pa, p_path_b = pb, p_interaction = a3,
    stringsAsFactors = FALSE)
  r$PM <- ifelse(r$TE_p < 0.05, r$NIE / r$TE, NA_real_)
  pmb <- bt[, "NIE", ] / bt[, "TE", ]
  r$PM_lo <- ifelse(r$TE_p < 0.05, apply(pmb, 2, q, .025), NA_real_)
  r$PM_hi <- ifelse(r$TE_p < 0.05, apply(pmb, 2, q, .975), NA_real_)
  r$TE_pct_effective_genes <- 100 * (exp(r$TE) - 1)       # log-scale outcomes only
  r$NIE_ci_excludes_0 <- r$NIE_lo > 0 | r$NIE_hi < 0
  r$joint_significance <- r$p_path_a < 0.05 & r$p_path_b < 0.05

  show <- transform(r[, c("outcome", "TE", "TE_p", "NIE", "NIE_lo", "NIE_hi", "NIE_p_perm",
                          "p_path_a", "p_path_b", "PM", "p_interaction")],
                    TE = signif(TE, 3), TE_p = signif(TE_p, 3), NIE = signif(NIE, 3),
                    NIE_lo = signif(NIE_lo, 3), NIE_hi = signif(NIE_hi, 3),
                    NIE_p_perm = signif(NIE_p_perm, 3), p_path_a = signif(p_path_a, 3),
                    p_path_b = signif(p_path_b, 3), PM = round(PM, 3), p_interaction = signif(p_interaction, 3))
  print(show, row.names = FALSE)
  RES[[cfg$label]] <- r; BOOT[[cfg$label]] <- list(boot = bt, perm = pn)
}

## =============================================================================
## D-B: where does lithium's change in bulk diversity come from?
## =============================================================================
rule("D-B  decomposition of lithium's effect on bulk Shannon diversity (lineage, lithium)")
cfg <- CFG[[1]]; ti <- cfg$titles
lz <- prep$lineage_z[ti, , drop = FALSE]
Hk <- as.matrix(A[ti, paste0("H_lineage_", colnames(lz))])
colnames(Hk) <- colnames(lz)        # so Hk[, "gran"] etc. index by lineage name
fbar <- colMeans(lz); Hbar <- colMeans(Hk)
comp <- data.frame(
  H_bulk        = A[ti, "H_primary"],
  reweighting   = as.vector(lz %*% Hbar),                         # mix changes, cells fixed
  within_lineage = as.vector(Hk %*% fbar) ,                        # cells change, mix fixed
  cross         = rowSums((lz - rep(fbar, each = nrow(lz))) * (Hk - rep(Hbar, each = nrow(Hk)))),
  between_JSD   = A[ti, "JSD_between"],
  residual      = A[ti, "resid_decomp"])
# H_within = reweighting + within_lineage + cross - sum(fbar * Hbar), exactly
stopifnot(max(abs(A[ti, "H_within"] - (comp$reweighting + comp$within_lineage + comp$cross - sum(fbar * Hbar)))) < 1e-10)
coef_x <- function(y, x, C) unname(coef(lm(y ~ x + C))[2])
est <- sapply(comp, coef_x, x = cfg$x, C = cfg$C)
say("  identity check: effect on H_bulk = sum of the five component effects: %.3g vs %.3g\n",
    est["H_bulk"], sum(est[-1]))
stopifnot(abs(est["H_bulk"] - sum(est[-1])) < 1e-10)
set.seed(481)
bd <- t(replicate(B, { i <- sample.int(nrow(comp), replace = TRUE)
  sapply(comp[i, ], coef_x, x = cfg$x[i], C = cfg$C[i, , drop = FALSE]) }))
DB <- data.frame(component = names(est),
  meaning = c("bulk diversity (total)", "MIXTURE: lineage weights change, lineage diversity held at average",
              "CELLS: each lineage's own diversity changes, weights held at average",
              "weights and lineage diversity changing together",
              "MIXTURE: divergence between lineages, weighted by the mix",
              "what the bMIND mixture does not reproduce"),
  effect_of_lithium = est,
  lo = apply(bd, 2, quantile, .025), hi = apply(bd, 2, quantile, .975),
  share_of_total = est / est["H_bulk"], stringsAsFactors = FALSE)
print(transform(DB[, -2], effect_of_lithium = signif(effect_of_lithium, 3), lo = signif(lo, 3),
                hi = signif(hi, 3), share_of_total = round(share_of_total, 3)), row.names = FALSE)
lin_eff <- t(sapply(colnames(lz), function(k) {
  f <- summary(lm(Hk[, k] ~ cfg$x + cfg$C)); c(effect = f$coef[2, 1], p = f$coef[2, 4]) }))
say("\n  lithium's effect on each lineage's OWN diversity (bMIND profiles):\n")
print(signif(as.data.frame(lin_eff), 3))

## =============================================================================
## D-C: the cell mix's own diversity, and a 1-D mediation through it
## =============================================================================
rule("D-C  diversity of the cell mix itself")
DC <- do.call(rbind, lapply(CFG, function(cfg) {
  hcol <- if (cfg$tree == "proposal") "H_cellmix6" else "H_cellmix5"
  h <- A[cfg$titles, hcol]
  f <- summary(lm(h ~ cfg$x + cfg$C))
  # 1-D mediation: exposure -> Shannon of the cell mix -> transcriptome Shannon
  Mh <- matrix(h, ncol = 1, dimnames = list(NULL, "H_cellmix"))
  fe <- fast_eff(cfg$x, Mh, cfg$C, matrix(A[cfg$titles, "H_primary"], ncol = 1))
  set.seed(481)
  bb <- replicate(2000, { i <- sample.int(length(h), replace = TRUE)
    fast_eff(cfg$x[i], Mh[i, , drop = FALSE], cfg$C[i, , drop = FALSE],
             matrix(A[cfg$titles[i], "H_primary"], ncol = 1))["NIE", 1] })
  data.frame(config = cfg$label, cellmix_diversity = hcol, mean = mean(h),
             effect_of_exposure = f$coef[2, 1], p = f$coef[2, 4],
             NIE_via_cellmix_diversity = fe["NIE", 1],
             NIE_lo = quantile(bb, .025, names = FALSE), NIE_hi = quantile(bb, .975, names = FALSE),
             stringsAsFactors = FALSE)
}))
print(transform(DC, mean = round(mean, 4), effect_of_exposure = signif(effect_of_exposure, 3),
                p = signif(p, 3), NIE_via_cellmix_diversity = signif(NIE_via_cellmix_diversity, 3),
                NIE_lo = signif(NIE_lo, 3), NIE_hi = signif(NIE_hi, 3)), row.names = FALSE)

## ---- post-hoc sensitivity, motivated by 11_diversity_prep.R ------------------
# The primary H correlates -0.50 with the share of structural ncRNA, a library
# quality marker. NOT pre-declared: added because that diagnostic appeared.
rule("POST HOC  primary alpha mediation adjusting also for the structural-ncRNA share")
cfg <- CFG[[1]]
C2 <- cbind(cfg$C, ncRNA_share = A[cfg$titles, "share_struct_ncRNA"])
fe <- fast_eff(cfg$x, cfg$M, C2, matrix(A[cfg$titles, "H_primary"], ncol = 1))
set.seed(481)
bb <- t(replicate(B, { i <- sample.int(length(cfg$x), replace = TRUE)
  fast_eff(cfg$x[i], cfg$M[i, , drop = FALSE], C2[i, , drop = FALSE],
           matrix(A[cfg$titles[i], "H_primary"], ncol = 1))[, 1] }))
PH <- data.frame(analysis = "H_primary, + ncRNA share covariate (post hoc)",
                 TE = fe["TE", 1], NIE = fe["NIE", 1],
                 NIE_lo = quantile(bb[, "NIE"], .025, names = FALSE),
                 NIE_hi = quantile(bb[, "NIE"], .975, names = FALSE),
                 PM = fe["NIE", 1] / fe["TE", 1], stringsAsFactors = FALSE)
print(transform(PH, TE = signif(TE, 3), NIE = signif(NIE, 3), NIE_lo = signif(NIE_lo, 3),
                NIE_hi = signif(NIE_hi, 3), PM = round(PM, 3)), row.names = FALSE)

## ---- save --------------------------------------------------------------------
dir.create("results/diversity", showWarnings = FALSE)
AR <- do.call(rbind, RES)
write.csv(AR, "results/diversity/alpha_results.csv", row.names = FALSE)
write.csv(DB, "results/diversity/alpha_decomposition_DB.csv", row.names = FALSE)
write.csv(DC, "results/diversity/cellmix_diversity_DC.csv", row.names = FALSE)
write.csv(PH, "results/diversity/alpha_posthoc_ncRNA.csv", row.names = FALSE)
saveRDS(list(results = AR, decomposition = DB, lineage_effects = lin_eff, cellmix = DC,
             posthoc = PH, boot = BOOT, B = B), "results/diversity/alpha.rds")
cat("\n12_diversity_alpha complete.\n")
