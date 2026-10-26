#!/usr/bin/env Rscript
# ==============================================================================
# diag_null_tail.R -- why does an occasional label permutation return >1000
# "significant" genes out of 2000?
#
# followup_validate_centring.R reported, for the UNCENTRED pipeline, a null
# distribution that is 0 in >=95% of permutations but whose maximum reaches
# 1354 (LI granulocytes) and 1927 (LI B cells). A 68-96% significant rate on
# scrambled labels is not a false-positive rate, it is a broken fit. Until the
# cause is known, neither the raw nor the centred calibration claim can be
# stated, because both are measured against this same null.
#
# This records, for each permutation, the quantities that could produce it:
# the moderated-variance prior, the design rank and conditioning, the smallest
# group size, and the resulting p-value distribution.
# ==============================================================================
suppressPackageStartupMessages({library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "de_final")); set.seed(481); options(width = 175)
FDR <- 0.05; NPERM <- 60; NSUB <- 2000

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
Y <- A[sub, "gran", mm$title]; W <- mkw(SE[, "gran", mm$title])[sub, ]

cat(sprintf("== diag_null_tail ==  LI granulocytes | %d genes x %d samples | %d permutations\n",
            nrow(Y), ncol(Y), NPERM))
cat(sprintf("design: %d columns, rank %d\n\n", ncol(X), qr(X)$rank))

R <- vector("list", NPERM)
for (b in seq_len(NPERM)) {
  pr <- sample(nrow(X)); Xp <- X; Xp[, "grpcase"] <- Xp[pr, "grpcase"]
  qrx <- qr(Xp)
  f <- eBayes(lmFit(Y, Xp, weights = W), trend = TRUE)
  tt <- topTable(f, coef = "grpcase", number = Inf, sort.by = "none")
  q <- p.adjust(tt$P.Value, "BH")
  R[[b]] <- data.frame(b = b, n = sum(q < FDR),
    rank = qrx$rank, ncol = ncol(Xp),
    kappa = kappa(Xp),
    n_case = sum(Xp[, "grpcase"] == 1),
    df_prior = f$df.prior[1], s2_prior = mean(f$s2.prior),
    df_total = mean(f$df.total),
    min_p = min(tt$P.Value), frac_p01 = mean(tt$P.Value < 0.01),
    max_absT = max(abs(f$t[, "grpcase"])),
    max_absLFC = max(abs(tt$logFC)),
    min_s2post = min(f$s2.post), max_absCoef = max(abs(f$coefficients[, "grpcase"])))
}
D <- do.call(rbind, R)
D$blown <- D$n > 50
cat(sprintf("permutations returning >50 genes: %d of %d\n\n", sum(D$blown), NPERM))
print(D[order(-D$n), ][1:12, ], row.names = FALSE, digits = 4)
cat("\n--- group means of every diagnostic, blown vs clean ---\n")
num <- sapply(D, is.numeric)
agg <- aggregate(D[, num & names(D) != "b"], by = list(blown = D$blown), FUN = median)
print(agg, row.names = FALSE, digits = 4)
saveRDS(D, "results/diag_null_tail.rds")
cat("\nwrote results/diag_null_tail.rds\n")
