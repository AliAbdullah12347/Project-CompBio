#!/usr/bin/env Rscript
# ==============================================================================
# diag_null_mechanism.R -- the heavy tail is gene-gene dependence, not a bug.
#
# diag_null_tail.R ruled out the numerical explanations: in the permutation
# that returned 1354 of 2000 genes the design was full rank, the conditioning
# was ordinary, the variance prior and residual df matched every other
# permutation, and max|t| was 4.6. What distinguished it was that 44% of genes
# had p < 0.01 simultaneously, against 0.3% normally.
#
# That is what happens when the genes are not independent. bMIND's estimates
# carry a per-sample additive offset shared by every gene in a lineage, so the
# genes move together across samples. A permutation whose labels happen to
# align with that shared axis shifts every gene at once, and BH -- which
# assumes independence or positive regression dependence -- converts a single
# sample-level coincidence into a thousand gene-level discoveries.
#
# PREDICTION, stated before the test: the number of genes a permutation returns
# should be a function of how well its labels align with the dominant
# sample-level axis of the data, and of nothing else. If so, the heavy tail is
# a property of the data, not a defect in the code, and the same mechanism
# applies to raw and centred alike.
# ==============================================================================
suppressPackageStartupMessages({library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "de_final")); set.seed(481); options(width = 175)
FDR <- 0.05; NPERM <- 150; NSUB <- 2000

p  <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
bp <- readRDS(file.path(ROOT, "de_analysis/results/bmind_profiles.rds"))
A <- bp$A; SE <- bp$SE
s <- p$sel[["LI"]]
mm <- p$meta[s$keep, , drop = FALSE]; mm$grp <- s$grp
mm$plate <- droplevels(factor(mm$plate)); mm$sex <- droplevels(factor(mm$sex))
cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]])) > 1, logical(1))]
X <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse = " + "))), data = mm)
mkw <- function(S) { W <- 1/(S^2); f <- is.finite(W) & W > 0
  if (any(!f)) W[!f] <- min(W[f]); W/mean(W) }
set.seed(481); sub <- sort(sample(dimnames(A)[[1]], NSUB))

cat(sprintf("== diag_null_mechanism ==  %d permutations, %d genes\n\n", NPERM, NSUB))

OUT <- list()
for (CT in c("gran", "T")) {
  Y <- A[sub, CT, mm$title]; W <- mkw(SE[, CT, mm$title])[sub, ]

  ## How dependent are the genes? Variance explained by the leading sample axis
  ## of the gene-centred matrix, and the mean pairwise gene correlation.
  Yc <- Y - rowMeans(Y)
  sv <- svd(Yc, nu = 1, nv = 1)
  pc1 <- sv$v[, 1]
  vexp <- sv$d[1]^2 / sum(sv$d^2)
  set.seed(1); gi <- sample(nrow(Y), 400)
  cm <- cor(t(Y[gi, ])); mean_r <- mean(cm[upper.tri(cm)])
  ## residualise PC1 on the covariates so the alignment measure is comparable
  ## to what the model actually tests
  Xn <- X[, setdiff(colnames(X), "grpcase"), drop = FALSE]
  pc1r <- residuals(lm(pc1 ~ Xn - 1))

  cat(sprintf("--- %s ---\n  leading sample axis explains %.1f%% of variance | mean pairwise gene r = %.3f\n",
              CT, 100*vexp, mean_r))

  R <- data.frame(b = seq_len(NPERM), n = NA_integer_, align = NA_real_)
  for (b in seq_len(NPERM)) {
    pr <- sample(nrow(X)); Xp <- X; Xp[, "grpcase"] <- Xp[pr, "grpcase"]
    g <- Xp[, "grpcase"]
    R$align[b] <- abs(cor(residuals(lm(g ~ Xn - 1)), pc1r))
    fit <- eBayes(lmFit(Y, Xp, weights = W), trend = TRUE)
    tt <- topTable(fit, coef = "grpcase", number = Inf, sort.by = "none")
    R$n[b] <- sum(p.adjust(tt$P.Value, "BH") < FDR) }

  rho <- cor(R$align, R$n, method = "spearman")
  cat(sprintf("  Spearman(alignment with that axis, genes returned) = %+.3f\n", rho))
  br <- cut(R$align, quantile(R$align, 0:4/4), include.lowest = TRUE,
            labels = c("Q1 lowest","Q2","Q3","Q4 highest"))
  ag <- aggregate(R$n, list(alignment = br), function(z)
          c(median = median(z), max = max(z), pct_over_50 = 100*mean(z > 50)))
  print(data.frame(alignment = ag$alignment, median = ag$x[,1], max = ag$x[,2],
                   pct_over_50 = ag$x[,3]), row.names = FALSE)
  OUT[[CT]] <- data.frame(lineage = CT, var_explained_pc1 = vexp, mean_gene_r = mean_r,
                          spearman_align_n = rho,
                          max_n = max(R$n), q95_n = quantile(R$n, .95, names = FALSE),
                          median_n = median(R$n), pct_perm_over_50 = 100*mean(R$n > 50))
  cat("\n") }
D <- do.call(rbind, OUT)
write.csv(D, "results/diag_null_mechanism.csv", row.names = FALSE)
print(D, row.names = FALSE, digits = 4)
cat("\nwrote results/diag_null_mechanism.csv\n")
