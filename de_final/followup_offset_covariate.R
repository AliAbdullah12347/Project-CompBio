#!/usr/bin/env Rscript
# ==============================================================================
# followup_offset_covariate.R -- the right way to remove the bMIND tilt.
#
# THE PROBLEM, RESTATED
#
# bMIND's estimates carry a per-sample, per-lineage additive offset. It shows up
# as a lineage-wide shift in the fold changes (granulocytes, lithium: median
# -0.118 against a gene-to-gene SD of 0.090 -- the artifact is larger than the
# signal it is burying), and it correlates with group, because bMIND's estimate
# for a lineage depends on that sample's fraction for that lineage and both
# contrasts move those fractions.
#
# WHY SUBTRACTING IT FAILED
#
# followup_validate_centring.R tested per-sample median centring against 100
# label permutations. The uncentred pipeline is well calibrated -- scrambled
# labels give 0-2 genes every time. The centred pipeline is not: for bipolar
# granulocytes the scrambled null reached a 95th percentile of 34.6 genes and a
# maximum of 298. Subtraction removes a between-sample axis of variation
# without telling the model it has done so, so the residual degrees of freedom
# are overstated, every standard error shrinks, and every t inflates -- null
# genes included.
#
# THE FIX
#
# Put the offset in the design matrix instead of subtracting it. The model then
# removes the same nuisance, but it SPENDS A DEGREE OF FREEDOM doing so and the
# standard errors stay honest. This is the ordinary treatment of a nuisance
# covariate; subtraction is the same operation with the accounting thrown away.
#
# THE COST, WHICH MUST BE MEASURED NOT ASSUMED
#
# If the offset is strongly correlated with group, adjusting for it also
# removes any genuine lineage-wide biological effect, and the test narrows to
# "does this gene deviate from its lineage's overall shift?". That is a
# defensible question but a different one, and the correlation determines how
# much it costs. It is reported per lineage below, before any count.
#
# All four treatments are compared on identical data, and the two defensible
# ones (raw, covariate) are validated against the same permutation null.
# ==============================================================================

suppressPackageStartupMessages({library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "de_final")); set.seed(481); options(width = 180)
RES <- "results"
FDR <- 0.05; NPERM <- 100; NSUB <- 2000

p  <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
bp <- readRDS(file.path(ROOT, "de_analysis/results/bmind_profiles.rds"))
A <- bp$A; SE <- bp$SE; LINS <- bp$lineages
sym <- local({ m <- read.delim(file.path(ROOT, "de_analysis/data/ens2sym.tsv"),
                               header = FALSE, col.names = c("ens", "sym"))
               setNames(m$sym, m$ens) })
GENES <- dimnames(A)[[1]]

design_for <- function(contrast) {
  s <- p$sel[[contrast]]
  mm <- p$meta[s$keep, , drop = FALSE]; mm$grp <- s$grp
  mm$plate <- droplevels(factor(mm$plate)); mm$sex <- droplevels(factor(mm$sex))
  if ("group" %in% s$cov) mm$group <- droplevels(factor(mm$group))
  cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]])) > 1, logical(1))]
  X <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse = " + "))), data = mm)
  list(X = X, samples = mm$title, grp = mm$grp) }

mkw <- function(S) { W <- 1/(S^2); f <- is.finite(W) & W > 0
  if (any(!f)) W[!f] <- min(W[f]); W/mean(W) }

ndeg <- function(Y, X, W) {
  f <- eBayes(lmFit(Y, X, weights = W), trend = TRUE)
  tt <- topTable(f, coef = "grpcase", number = Inf, sort.by = "none")
  list(n = sum(p.adjust(tt$P.Value, "BH") < FDR),
       up = sum(p.adjust(tt$P.Value, "BH") < FDR & tt$logFC > 0),
       dn = sum(p.adjust(tt$P.Value, "BH") < FDR & tt$logFC < 0),
       p = tt$P.Value, lfc = tt$logFC) }

cat(sprintf("== followup_offset_covariate ==\n%d permutations | %d-gene test subsample | BH %.2f\n",
            NPERM, NSUB, FDR))

SUMM <- list(); VAL <- list(); GENEOUT <- list()
for (CN in c("LI", "BPD")) {
  d <- design_for(CN)
  g01 <- as.integer(d$grp == "case")
  sub <- sort(sample(GENES, NSUB))
  cat(sprintf("\n%s\n== %s  (n = %d) ==\n%s\n", strrep("=", 118), CN, length(d$samples), strrep("=", 118)))

  ## ---- the offset, and how entangled it is with group ---------------------
  cat("\n--- the per-sample offset, and what adjusting for it will cost ---\n")
  OFF <- sapply(LINS, function(ct) apply(A[, ct, d$samples], 2, median))
  cat("  lineage | cor(offset, group) | offset SD | shared variance with group\n")
  for (ct in LINS) {
    r <- cor(OFF[, ct], g01)
    cat(sprintf("  %-7s |      %+6.3f        |  %.4f   |  %5.1f%%  %s\n", ct, r, sd(OFF[, ct]), 100*r^2,
        if (abs(r) > 0.7) "<- heavily entangled: adjustment removes much of any real global effect"
        else if (abs(r) > 0.4) "<- partly entangled" else "")) }

  ## ---- four treatments ------------------------------------------------------
  cat("\n--- four treatments of the same data ---\n")
  ROWS <- list()
  for (ct in LINS) {
    Y  <- A[, ct, d$samples]; W <- mkw(SE[, ct, d$samples])
    Ys <- Y[sub, ]; Ws <- W[sub, ]
    Yc <- sweep(Y, 2, OFF[, ct], "-")[sub, ]
    Xo <- cbind(d$X, offset_lin = OFF[, ct])
    ok <- qr(Xo)$rank == ncol(Xo)

    r_raw <- ndeg(Ys, d$X, Ws)
    sh    <- median(r_raw$lfc)
    f0    <- eBayes(lmFit(Ys, d$X, weights = Ws), trend = TRUE)
    se0   <- f0$stdev.unscaled[, "grpcase"] * sqrt(f0$s2.post)
    q_ph  <- p.adjust(2*pt(-abs((r_raw$lfc - sh)/se0), df = f0$df.total), "BH")
    r_cen <- ndeg(Yc, d$X, Ws)
    r_cov <- if (ok) ndeg(Ys, Xo, Ws) else list(n = NA, up = NA, dn = NA)

    cat(sprintf("  %-5s raw %4d (%3du/%3dd) | post-hoc %4d | centred %4d (%3du/%3dd) | covariate %4s (%3su/%3sd)%s\n",
        ct, r_raw$n, r_raw$up, r_raw$dn, sum(q_ph < FDR), r_cen$n, r_cen$up, r_cen$dn,
        r_cov$n, r_cov$up, r_cov$dn, if (ok) "" else "  <- DESIGN RANK-DEFICIENT"))
    ROWS[[ct]] <- data.frame(contrast = CN, lineage = ct,
      cor_offset_group = round(cor(OFF[, ct], g01), 4),
      raw = r_raw$n, raw_up = r_raw$up, raw_dn = r_raw$dn,
      posthoc = sum(q_ph < FDR), centred = r_cen$n,
      covariate = r_cov$n, covariate_up = r_cov$up, covariate_dn = r_cov$dn,
      full_rank = ok, stringsAsFactors = FALSE) }
  SUMM[[CN]] <- do.call(rbind, ROWS)

  ## ---- validate raw and covariate against the same null --------------------
  cat(sprintf("\n--- validating against %d label permutations ---\n", NPERM))
  PRE <- lapply(LINS, function(ct) list(
    Y = A[sub, ct, d$samples], W = mkw(SE[, ct, d$samples])[sub, ], off = OFF[, ct]))
  names(PRE) <- LINS
  NRAW <- NCOV <- matrix(NA_integer_, NPERM, length(LINS), dimnames = list(NULL, LINS))
  t0 <- Sys.time()
  for (b in seq_len(NPERM)) {
    pr <- sample(nrow(d$X))
    Xp <- d$X; Xp[, "grpcase"] <- Xp[pr, "grpcase"]
    for (ct in LINS) {
      NRAW[b, ct] <- ndeg(PRE[[ct]]$Y, Xp, PRE[[ct]]$W)$n
      NCOV[b, ct] <- ndeg(PRE[[ct]]$Y, cbind(Xp, offset_lin = PRE[[ct]]$off), PRE[[ct]]$W)$n }
    if (b %% 25 == 0) { el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      cat(sprintf("    %3d/%d  %.1f min, ETA %.1f\n", b, NPERM, el, el/b*(NPERM-b))); flush.console() } }

  V <- do.call(rbind, lapply(LINS, function(ct) {
    o_raw <- SUMM[[CN]]$raw[SUMM[[CN]]$lineage == ct]
    o_cov <- SUMM[[CN]]$covariate[SUMM[[CN]]$lineage == ct]
    data.frame(contrast = CN, lineage = ct,
      obs_raw = o_raw, null_raw_q95 = quantile(NRAW[, ct], .95, names = FALSE),
      null_raw_max = max(NRAW[, ct]),
      p_raw = (sum(NRAW[, ct] >= o_raw) + 1)/(NPERM + 1),
      obs_cov = o_cov, null_cov_q95 = quantile(NCOV[, ct], .95, names = FALSE),
      null_cov_max = max(NCOV[, ct]),
      p_cov = (sum(NCOV[, ct] >= o_cov) + 1)/(NPERM + 1), stringsAsFactors = FALSE) }))
  print(V, row.names = FALSE)

  cat("\n  is the covariate model calibrated? (a well-behaved null has a small 95th pct and max)\n")
  for (i in seq_len(nrow(V))) {
    r <- V[i, ]
    cat(sprintf("    %-5s null q95 %5.1f / max %4d  ->  %s\n", r$lineage, r$null_cov_q95, r$null_cov_max,
      if (r$null_cov_max > 20) "NOT CALIBRATED -- scrambled labels still yield genes"
      else if (r$obs_cov > 0 && r$p_cov < 0.05) "calibrated, and the observed result beats it"
      else if (r$obs_cov > 0) "calibrated, but observed does not beat the null"
      else "calibrated; nothing found")) }
  VAL[[CN]] <- V
}
S <- do.call(rbind, SUMM); V <- do.call(rbind, VAL)
write.csv(S, file.path(RES, "followup_D_offset_covariate.csv"), row.names = FALSE)
write.csv(V, file.path(RES, "followup_D_validation.csv"), row.names = FALSE)

cat(sprintf("\n%s\nCALIBRATION COMPARISON -- worst-case genes from scrambled labels\n%s\n",
            strrep("=", 118), strrep("=", 118)))
cmp <- rbind(
  data.frame(method = "raw (no correction)", max_null = max(V$null_raw_max)),
  data.frame(method = "offset as covariate", max_null = max(V$null_cov_max)),
  data.frame(method = "per-sample centring", max_null = 298))
print(cmp, row.names = FALSE)
cat("  (the centring figure is from followup_validate_centring.R, same 100-permutation design)\n")
cat("\nwrote results/followup_D_*.csv\n")
