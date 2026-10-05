#!/usr/bin/env Rscript
# ==============================================================================
# 04_permutation.R -- the permutation null for whole blood.
#
# WHY THIS EXISTS, AND WHY IT MATTERS MORE THAN THE BOOTSTRAP
#
# The pre-specified uncertainty estimate for pi0 is a bootstrap over genes.
# That is the standard approach and it is what config.yaml commits to, but it
# assumes genes are independent. They are emphatically not: co-expression is
# the entire premise of the network arm. A gene bootstrap therefore produces a
# confidence interval that is too NARROW, which biases the decision rule toward
# "the CI excludes 1" and hence toward INCONCLUSIVE.
#
# A permutation null does not have that problem. Each permutation scrambles the
# group labels while leaving the expression matrix — and therefore every
# gene-gene correlation — completely intact. The resulting distribution of
# DEG counts and pi0 values is what this design actually produces when there is
# no group difference, correlation structure and all.
#
# Both are reported. Where they disagree, the permutation is the honest one.
#
# The exposure column alone is permuted, leaving every covariate attached to
# its own sample. That tests "this exposure carries no signal", which is the
# hypothesis, rather than "these samples are exchangeable", which is not.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})
HERE <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/de_analysis"
setwd(HERE)
args <- commandArgs(trailingOnly = TRUE)
NPERM <- if (length(args)) as.integer(args[1]) else 200
WHICH <- if (length(args) > 1) args[-1] else c("LI","BPD")
set.seed(481)
p <- readRDS("data/prep.rds")
RES <- "results"

pi0_storey <- function(pv, lambda = 0.5) min(1, mean(pv > lambda, na.rm = TRUE) / (1 - lambda))

perm_null <- function(contrast, B) {
  s <- p$sel[[contrast]]
  mm <- p$meta[s$keep, , drop = FALSE]; mm$grp <- s$grp
  mm$plate <- droplevels(mm$plate); mm$sex <- droplevels(mm$sex)
  if ("group" %in% s$cov) mm$group <- droplevels(mm$group)
  cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]])) > 1, logical(1))]
  X0 <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse = " + "))), data = mm)
  ne <- nonEstimable(X0); if (!is.null(ne)) X0 <- X0[, setdiff(colnames(X0), ne), drop = FALSE]
  cnt <- p$counts[, mm$title, drop = FALSE]
  dge <- calcNormFactors(DGEList(cnt), method = "TMM")
  # voom weights are recomputed per permutation; reusing them would leak the
  # real design into every null replicate.
  out <- vapply(seq_len(B), function(i) {
    X <- X0; X[, "grpcase"] <- X0[sample(nrow(X0)), "grpcase"]
    tt <- topTable(eBayes(lmFit(voom(dge, X), X)), coef = "grpcase",
                   number = Inf, sort.by = "none")
    c(deg = sum(tt$adj.P.Val < 0.05, na.rm = TRUE),
      pct = 100 * sum(tt$adj.P.Val < 0.05, na.rm = TRUE) / nrow(tt),
      pi0 = pi0_storey(tt$P.Value))
  }, numeric(3))

  t(out)
}

cat(sprintf("== 04_permutation ==  %d permutations per contrast\n\n", NPERM))
ALL <- list()
for (cn in WHICH) {
  t0 <- Sys.time()
  N <- perm_null(cn, NPERM)
  ALL[[cn]] <- N
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("%s  [%.1f min]\n", cn, el))
  cat(sprintf("  DEG : mean %.2f  median %.0f  max %d  nonzero in %d/%d\n",
              mean(N[, "deg"]), median(N[, "deg"]), max(N[, "deg"]),
              sum(N[, "deg"] > 0), NPERM))
  cat(sprintf("  pct : mean %.4f%%  95th pct %.4f%%\n",
              mean(N[, "pct"]), quantile(N[, "pct"], 0.95)))
  cat(sprintf("  pi0 : mean %.4f  2.5th %.4f  97.5th %.4f\n\n",
              mean(N[, "pi0"]), quantile(N[, "pi0"], 0.025), quantile(N[, "pi0"], 0.975)))
  write.csv(as.data.frame(N), file.path(RES, sprintf("permnull_%s.csv", cn)), row.names = FALSE)
}

## ---- observed vs null -------------------------------------------------------
for (cn in setdiff(c("LI","BPD"), WHICH)) {
  f <- file.path(RES, sprintf("permnull_%s.csv", cn))
  if (file.exists(f)) { ALL[[cn]] <- as.matrix(read.csv(f)); cat(sprintf("%s loaded from earlier run
", cn)) }
}
wb <- readRDS(file.path(RES, "wholeblood.rds"))
SUM <- wb$SUM
cmp <- do.call(rbind, lapply(names(ALL), function(cn) {
  lab <- sprintf("WB_%s_raw", cn); N <- ALL[[cn]]
  obs_deg <- SUM[lab, "deg_fdr05"]; obs_pi0 <- SUM[lab, "pi0"]
  data.frame(contrast = cn, obs_deg = obs_deg,
             null_mean_deg = round(mean(N[, "deg"]), 2),
             null_max_deg = max(N[, "deg"]),
             p_emp_deg = (sum(N[, "deg"] >= obs_deg) + 1) / (NPERM + 1),
             obs_pi0 = obs_pi0,
             null_pi0_lo = round(quantile(N[, "pi0"], 0.025), 4),
             null_pi0_hi = round(quantile(N[, "pi0"], 0.975), 4),
             # This is the decision-relevant line: is the observed pi0 inside
             # the range this design produces when nothing is going on?
             pi0_inside_null = obs_pi0 >= quantile(N[, "pi0"], 0.025),
             stringsAsFactors = FALSE) }))
cat("=== observed vs permutation null ===\n")
print(cmp, row.names = FALSE)
cat("\n  pi0_inside_null = TRUE means the observed pi0 is consistent with no signal.\n")
cat("  FALSE means signal is present even if few or no genes clear FDR.\n")
write.csv(cmp, file.path(RES, "permutation_vs_observed.csv"), row.names = FALSE)
saveRDS(ALL, file.path(RES, "permnull.rds"))
cat("\nwrote results/permnull_*.csv and permutation_vs_observed.csv\n")
