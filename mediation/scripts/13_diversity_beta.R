#!/usr/bin/env Rscript
# ==============================================================================
# 13_diversity_beta.R -- beta diversity: does the exposure shift whole-
# transcriptome structure, change its heterogeneity, and is the shift mediated
# by the cell mix? (PRESPEC.md section 8, tests B1-B4, plus D-C beta.)
#
#   B1  location   PERMANOVA pseudo-F for the exposure given covariates
#   B2  dispersion PERMDISP (Anderson 2006), as declared: standard, unadjusted;
#                  a covariate-adjusted version is reported beside it
#   B3  mediation  (i)   partial dbRDA: exposure given covariates AND balances
#                  (ii)  share of the exposure's adjusted R^2 shared with the
#                        cell mix, with a bootstrap CI (Euclidean layers)
#                  (iii) joint significance of path a and path b
#   B4  only if B2 is significant on the primary layer
#   +   leading principal axes as single outcomes through med.R (descriptive)
#
# EUCLIDEAN LAYERS ARE FITTED AS RDA. Aitchison distance is Euclidean distance on
# CLR data, so PERMANOVA on it is exactly redundancy analysis on the CLR matrix,
# and that matrix can be replaced by its principal-component scores without
# changing any distance. This script verifies that equivalence numerically
# before relying on it. Jensen-Shannon and Bray-Curtis are not Euclidean and
# go through dbRDA on the distance matrix.
# ==============================================================================

ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "mediation")); set.seed(481)
suppressPackageStartupMessages({library(vegan); library(permute)})
source("scripts/med.R"); options(width = 170)
say  <- function(...) cat(sprintf(...))
rule <- function(t) cat(sprintf("\n%s\n%s\n", t, strrep("=", nchar(t))))
NPERM <- 9999L; BBOOT <- 1000L; PAR <- 3L

prep <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
dv   <- readRDS("data/diversity_input.rds")
A    <- dv$alpha; rownames(A) <- A$title
mats <- c(list(gene_aitchison = dv$CLR$primary, gene_aitchison_markers = dv$CLR$with_markers),
          setNames(dv$PCLR, paste0("pw_", names(dv$PCLR), "_aitchison")))

## ---- configurations (same as the alpha script) -------------------------------
build_cfg <- function(sel, tree, label) {
  s <- prep$sel[[sel]]; mm <- prep$meta[s$keep, , drop = FALSE]
  x <- as.integer(s$grp == "case")
  if (sel == "LI") {
    d <- readRDS(if (tree == "lineage") "data/med_input.rds" else "data/med_input_proposal.rds")
    M <- d$M; C <- d$C
  } else {
    M <- prep$ILR[mm$title, , drop = FALSE]
    mm$sex <- droplevels(factor(mm$sex)); mm$plate <- droplevels(factor(mm$plate))
    mm$group <- droplevels(factor(mm$group))
    C <- model.matrix(~ age + sex + tobacco + rin + plate + seqpc1 + seqpc2 + seqpc3 + group,
                      data = mm)[, -1, drop = FALSE]
  }
  C <- cbind(C, log_libsize = A[mm$title, "log_libsize"])
  colnames(C) <- make.names(colnames(C)); colnames(M) <- paste0("bal_", colnames(M))
  list(label = label, tree = tree, titles = mm$title, x = x, M = M, C = C)
}
CFG <- list(build_cfg("LI", "lineage", "lithium | lineage tree (PRIMARY)"),
            build_cfg("LI", "proposal", "lithium | proposal tree"),
            build_cfg("BPD", "lineage", "bipolar vs control | lineage tree (EXPLORATORY)"))

## ---- R^2 machinery for Euclidean layers --------------------------------------
pc_scores <- function(Y) {                      # exact Euclidean embedding
  Yc <- scale(Y, center = TRUE, scale = FALSE)
  s <- svd(Yc, nu = min(dim(Yc)) - 1, nv = 0)
  s$u %*% diag(s$d[seq_len(ncol(s$u))])
}
r2adj <- function(S, X) {                        # multivariate adjusted R^2
  n <- nrow(S); Sc <- scale(S, scale = FALSE); Xc <- scale(X, scale = FALSE)
  fit <- qr.fitted(qr(cbind(1, Xc)), Sc)
  r2 <- sum(fit^2) / sum(Sc^2); p <- qr(Xc)$rank
  1 - (1 - r2) * (n - 1) / (n - p - 1)
}
shared_fraction <- function(S, x, C, M) {
  ux_C  <- r2adj(S, cbind(x, C)) - r2adj(S, C)            # exposure | covariates
  ux_CM <- r2adj(S, cbind(x, C, M)) - r2adj(S, cbind(C, M))  # exposure | covariates + mix
  c(R2adj_exposure_given_C = ux_C, R2adj_exposure_given_C_M = ux_CM,
    shared_with_cellmix = 1 - ux_CM / ux_C)
}

ptest <- function(mod) { a <- anova(mod, permutations = how(nperm = NPERM), parallel = PAR)
  c(F = a$F[1], p = a$`Pr(>F)`[1]) }

RES <- list(); AXES <- list(); first <- TRUE
for (cfg in CFG) {
  rule(cfg$label)
  ti <- cfg$titles; x <- cfg$x; C <- cfg$C; M <- cfg$M; n <- length(x)
  df <- data.frame(x = x, C, M, check.names = TRUE)
  fC  <- paste(colnames(C), collapse = " + "); fM <- paste(colnames(M), collapse = " + ")

  ## path a, once per configuration: exposure -> cell mix (Euclidean on ILR =
  ## Aitchison distance between cell compositions)
  # PATH A AS DECLARED (PRESPEC 8.2 / DIVERSITY_RESEARCH 5.2): a MANOVA of the
  # balances on the exposure given covariates, Pillai's trace. Pillai uses the
  # balances' covariance, so a shift in a LOW-variance balance (gran vs mono, the
  # one lithium moves) counts in full. The RDA version below weights balances by
  # raw variance -- i.e. Aitchison distance on the cell mix -- and is reported
  # beside it, not instead of it.
  cMa <- scale(M, scale = FALSE)
  pa_pillai <- anova(lm(cMa ~ C), lm(cMa ~ x + C), test = "Pillai")$`Pr(>F)`[2]
  pa_rda <- ptest(rda(as.formula(paste("M ~ x + Condition(", fC, ")")), data = df))
  pa <- c(F = NA_real_, p = pa_pillai)
  disp_mix <- permutest(betadisper(dist(M), factor(x)), permutations = how(nperm = NPERM))
  say("  path a  exposure -> cell mix, MANOVA Pillai (declared): p = %.4f\n", pa_pillai)
  say("  path a  as RDA on ILR (Aitchison geometry, reported beside): F = %.3f, p = %.4f\n",
      pa_rda["F"], pa_rda["p"])
  say("  cell-mix dispersion (PERMDISP on ILR): F = %.3f, p = %.4f\n",
      disp_mix$tab$F[1], disp_mix$tab$`Pr(>F)`[1])

  layers <- c(names(mats), "gene_jsd", "gene_bray")
  if (cfg$tree == "proposal") layers <- c("gene_aitchison", "pw_hallmark_aitchison")   # trees differ only in M
  for (ly in layers) {
    t0 <- Sys.time()
    if (ly %in% names(mats)) {
      S <- pc_scores(mats[[ly]][ti, , drop = FALSE]); D <- NULL
      if (first) {    # verify RDA-on-scores == PERMANOVA on the Aitchison distance
        Dchk <- dist(mats[[ly]][ti, , drop = FALSE])
        # the F statistic does not depend on the permutations; 99 just satisfies adonis2
        a1 <- adonis2(as.formula(paste("Dchk ~ x +", fC)), data = df, by = "margin", permutations = 99)
        a2 <- anova(rda(as.formula(paste("S ~ x + Condition(", fC, ")")), data = df), permutations = 99)
        say("  EQUIVALENCE CHECK  PERMANOVA F = %.6f | RDA-on-scores F = %.6f -> %s\n",
            a1$F[1], a2$F[1], if (abs(a1$F[1] - a2$F[1]) < 1e-8) "IDENTICAL" else "**DIFFER**")
        stopifnot(abs(a1$F[1] - a2$F[1]) < 1e-8)
        first <- FALSE
      }
      mk <- function(f) rda(as.formula(f), data = df)
      B1  <- ptest(mk(paste("S ~ x + Condition(", fC, ")")))
      B3i <- ptest(mk(paste("S ~ x + Condition(", fC, "+", fM, ")")))
      Pb  <- ptest(mk(paste("S ~", fM, "+ Condition(x +", fC, ")")))
      Dd  <- dist(S)
      sf  <- shared_fraction(S, x, C, M)
      set.seed(481)
      bs <- replicate(BBOOT, { i <- sample.int(n, n, replace = TRUE)
        shared_fraction(S[i, , drop = FALSE], x[i], C[i, , drop = FALSE], M[i, , drop = FALSE])["shared_with_cellmix"] })
      sf_lo <- quantile(bs, .025, names = FALSE, na.rm = TRUE); sf_hi <- quantile(bs, .975, names = FALSE, na.rm = TRUE)
      # covariate-adjusted dispersion: residualise on covariates (not exposure)
      Sr <- qr.resid(qr(cbind(1, C)), S)
      B2adj <- permutest(betadisper(dist(Sr), factor(x)), permutations = how(nperm = NPERM))
      B2adj_F <- B2adj$tab$F[1]; B2adj_p <- B2adj$tab$`Pr(>F)`[1]
    } else {
      P <- dv$P_primary[, ti]
      Dd <- if (ly == "gene_bray") vegdist(t(P), method = "bray") else {
        Pt <- t(P); m <- nrow(Pt); J <- matrix(0, m, m)
        H  <- -rowSums(ifelse(Pt > 0, Pt * log(Pt), 0))
        for (a in seq_len(m - 1)) {
          b <- (a + 1):m
          Mx <- (matrix(Pt[a, ], length(b), ncol(Pt), byrow = TRUE) + Pt[b, , drop = FALSE]) / 2
          Hm <- -rowSums(ifelse(Mx > 0, Mx * log(Mx), 0))
          J[b, a] <- J[a, b] <- sqrt(pmax(Hm - (H[a] + H[b]) / 2, 0))
        }
        as.dist(J) }
      mk <- function(f) dbrda(as.formula(f), data = df)
      B1  <- ptest(mk(paste("Dd ~ x + Condition(", fC, ")")))
      B3i <- ptest(mk(paste("Dd ~ x + Condition(", fC, "+", fM, ")")))
      Pb  <- ptest(mk(paste("Dd ~", fM, "+ Condition(x +", fC, ")")))
      r2 <- function(f) RsquareAdj(mk(f))$adj.r.squared
      ux_C  <- r2(paste("Dd ~ x +", fC)) - r2(paste("Dd ~", fC))
      ux_CM <- r2(paste("Dd ~ x +", fC, "+", fM)) - r2(paste("Dd ~", fC, "+", fM))
      sf <- c(R2adj_exposure_given_C = ux_C, R2adj_exposure_given_C_M = ux_CM,
              shared_with_cellmix = 1 - ux_CM / ux_C)
      sf_lo <- sf_hi <- NA_real_; B2adj_F <- B2adj_p <- NA_real_
    }
    B2 <- permutest(betadisper(Dd, factor(x)), permutations = how(nperm = NPERM))
    bd <- betadisper(Dd, factor(x))
    r <- data.frame(config = cfg$label, layer = ly,
      B1_F = B1["F"], B1_p = B1["p"],
      B2_F = B2$tab$F[1], B2_p = B2$tab$`Pr(>F)`[1],
      B2_mean_dist_unexposed = mean(bd$distances[x == 0]), B2_mean_dist_exposed = mean(bd$distances[x == 1]),
      B2adj_F = B2adj_F, B2adj_p = B2adj_p,
      B3i_direct_F = B3i["F"], B3i_direct_p = B3i["p"],
      path_a_p = pa["p"], path_a_rda_p = pa_rda["p"], path_b_F = Pb["F"], path_b_p = Pb["p"],
      cellmix_dispersion_p = disp_mix$tab$`Pr(>F)`[1],
      R2adj_exposure_given_C = sf["R2adj_exposure_given_C"],
      R2adj_exposure_given_C_M = sf["R2adj_exposure_given_C_M"],
      shared_with_cellmix = sf["shared_with_cellmix"], shared_lo = sf_lo, shared_hi = sf_hi,
      stringsAsFactors = FALSE)
    r$joint_significance <- r$path_a_p < 0.05 & r$path_b_p < 0.05
    RES[[paste(cfg$label, ly)]] <- r
    say("  %-26s B1 p=%.4f | B2 p=%.4f | direct p=%.4f | shared %.3f [%s] | path b p=%.4f | %.0fs\n",
        ly, r$B1_p, r$B2_p, r$B3i_direct_p, r$shared_with_cellmix,
        if (is.na(sf_lo)) "no CI" else sprintf("%.3f, %.3f", sf_lo, sf_hi), r$path_b_p,
        as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }

  ## leading principal axes of the primary layer, each mediated through med.R
  S <- pc_scores(mats$gene_aitchison[ti, , drop = FALSE])
  ve <- colSums(S^2) / sum(S^2)
  cM <- centre(M); cC <- centre(C)
  mf <- mediator_fit(cM$X, x, cC$X); Xd <- med_design_full(x, cM$X, cC$X)
  Yax <- t(S[, 1:5, drop = FALSE])
  of <- outcome_fit(Yax, matrix(1, 5, n), Xd, se = TRUE, moderate = FALSE)
  ef <- med_effects(of, mf)
  ax <- data.frame(config = cfg$label, axis = paste0("PC", 1:5), var_explained = ve[1:5],
                   TE = ef$est$TE, TE_p = ef$se$TE, NIE = ef$est$NIE, NIE_se = ef$se$NIE,
                   stringsAsFactors = FALSE)
  ax$TE_p <- 2 * pnorm(-abs(ax$TE / ax$TE_p)); ax$NIE_p <- 2 * pnorm(-abs(ax$NIE / ax$NIE_se))
  ax$PM <- ifelse(ax$TE_p < 0.05, ax$NIE / ax$TE, NA_real_)
  ax$corr_with_granulocyte_fraction <- as.vector(cor(S[, 1:5], prep$lineage_z[ti, "gran"]))
  AXES[[cfg$label]] <- ax
  cat("\n"); print(transform(ax[, -1], var_explained = round(var_explained, 3), TE = signif(TE, 3),
                            TE_p = signif(TE_p, 3), NIE = signif(NIE, 3), NIE_p = signif(NIE_p, 3),
                            PM = round(PM, 3), corr_with_granulocyte_fraction = round(corr_with_granulocyte_fraction, 3)),
                   row.names = FALSE)
}

## ---- B4: only if dispersion differs on the PRIMARY layer ----------------------
prim <- RES[["lithium | lineage tree (PRIMARY) gene_aitchison"]]
rule("B4  dispersion mediation (pre-declared to run only if B2 is significant)")
if (prim$B2_p < 0.05) {
  say("  B2 was significant (p = %.4f): running B4\n", prim$B2_p)
  cfg <- CFG[[1]]; ti <- cfg$titles
  S <- pc_scores(mats$gene_aitchison[ti, , drop = FALSE])
  Sr <- qr.resid(qr(cbind(1, cfg$C)), S)
  ctr <- colMeans(Sr[cfg$x == 0, , drop = FALSE])
  ydist <- sqrt(rowSums(sweep(Sr, 2, ctr)^2))
  cM <- centre(cfg$M); cC <- centre(cfg$C)
  of <- outcome_fit(matrix(ydist, 1), matrix(1, 1, length(ydist)),
                    med_design_full(cfg$x, cM$X, cC$X), se = TRUE, moderate = FALSE)
  ef <- med_effects(of, mediator_fit(cM$X, cfg$x, cC$X))
  B4 <- data.frame(NIE = ef$est$NIE, NIE_se = ef$se$NIE, TE = ef$est$TE, TE_se = ef$se$TE)
  print(signif(B4, 3))
} else {
  say("  B2 was NOT significant (p = %.4f): B4 is not run, exactly as declared.\n", prim$B2_p)
  B4 <- NULL
}

BR <- do.call(rbind, RES); rownames(BR) <- NULL
AX <- do.call(rbind, AXES); rownames(AX) <- NULL
write.csv(BR, "results/diversity/beta_results.csv", row.names = FALSE)
write.csv(AX, "results/diversity/beta_principal_axes.csv", row.names = FALSE)
saveRDS(list(results = BR, axes = AX, B4 = B4, nperm = NPERM, nboot = BBOOT),
        "results/diversity/beta.rds")
cat("\n13_diversity_beta complete.\n")
