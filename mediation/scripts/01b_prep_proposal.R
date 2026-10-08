#!/usr/bin/env Rscript
# ==============================================================================
# 01b_prep_proposal.R -- build the SECOND mediator: the proposal's own tree.
#
# The written proposal specifies a six-part composition and five balances:
#
#   "five balances follow the hematopoietic lineage: myeloid {neutrophils,
#    monocytes} vs lymphoid {CD4 naive, CD4 memory, NK, B}; then neutrophils vs
#    monocytes; T {CD4 naive, CD4 memory} vs {NK, B}; CD4 naive vs CD4 memory;
#    and NK vs B cells"
#
# The implemented mediator (01_prep.R) instead aggregates all 22 LM22 types into
# five lineages, giving four balances, because the finer parts are heavily
# zero-inflated. Rather than choose between them, we run BOTH end to end and
# report them side by side. If they agree, the choice did not matter. If they
# disagree, the reader sees it.
#
# WHAT THIS COSTS, MEASURED NOT ASSUMED
#
#   * the six parts hold 92.4% of the cell mass; the other 7.6% (macrophages,
#     dendritic, mast, eosinophils, CD8, Tfh, Treg, gamma-delta, plasma) is
#     dropped and the six re-closed to sum to one.
#   * 15.0% of samples contain at least one zero among the six (CD4 naive 3.0%,
#     B 12.2%), so zero replacement actually fires here, unlike the lineage
#     version where it never does. That makes the delta sensitivity (S6) a real
#     analysis rather than a formality.
#
# Everything else -- subjects, covariates, gene split, voom -- is identical to
# 01_prep.R, so the two runs differ in the mediator and nothing else.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
HERE <- file.path(ROOT, "mediation")
setwd(HERE); set.seed(481)
options(width = 150)
say <- function(...) cat(sprintf(...))

p <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
m <- p$meta
say("== 01b_prep_proposal ==\n")

## ---- the six parts, as the proposal names them ------------------------------
fr <- read.csv(file.path(ROOT, "data/cohort_474/fractions_bmode_474_CIBERSORTx.csv"),
               check.names = FALSE, colClasses = c(Mixture = "character"))
rownames(fr) <- fr$Mixture
fr <- as.matrix(fr[m$title, setdiff(names(fr), c("Mixture", "P-value", "Correlation", "RMSE"))])

PARTS <- list(
  neut = "Neutrophils",
  mono = "Monocytes",
  cd4n = "T cells CD4 naive",
  cd4m = c("T cells CD4 memory resting", "T cells CD4 memory activated"),
  NK   = c("NK cells resting", "NK cells activated"),
  B    = c("B cells naive", "B cells memory"))
stopifnot(all(unlist(PARTS) %in% colnames(fr)))

P <- sapply(PARTS, function(cs) rowSums(fr[, cs, drop = FALSE]))
rownames(P) <- m$title
mass <- rowSums(P)
say("six parts: %s\n", paste(names(PARTS), collapse = ", "))
say("  share of total cell mass retained: mean %.3f (range %.3f-%.3f)\n",
    mean(mass), min(mass), max(mass))
say("  dropped types (7.6%% of mass): %s\n",
    paste(setdiff(colnames(fr), unlist(PARTS)), collapse = ", "))
P <- P / mass                                   # re-close the six to sum to 1
stopifnot(max(abs(rowSums(P) - 1)) < 1e-12)

nz <- colSums(P == 0)
say("\nzeros per part (of %d samples): %s\n", nrow(P),
    paste(sprintf("%s=%d", names(nz), nz), collapse = "  "))
say("  samples with at least one zero: %d (%.1f%%)  -> zero replacement DOES fire here\n",
    sum(apply(P == 0, 1, any)), 100 * mean(apply(P == 0, 1, any)))

## ---- zero replacement, identical rule to the lineage version ----------------
# Multiplicative simple replacement (Martin-Fernandez et al. 2003): zeros become
# delta x the smallest observed value in that part, and the non-zero parts are
# rescaled down so the composition still sums to one and the ratios among the
# observed parts are preserved exactly. Additive replacement would distort those
# ratios, which is the whole quantity we are about to take logs of.
zrepl <- function(X, delta = 0.65) {
  lim <- apply(X, 2, function(v) { v <- v[v > 0]; if (!length(v)) 1e-8 else min(v) }) * delta
  z <- X <= 0
  ins <- rowSums(sweep(z, 2, lim, "*"))
  X[z] <- matrix(rep(lim, each = nrow(X)), nrow(X))[z]
  X[!z] <- (X * (1 - ins))[!z]
  X / rowSums(X)
}
Pz <- zrepl(P)
say("  replacement delta = 0.65 | max |change| %.3g | min value after %.3g (positive: %s)\n",
    max(abs(Pz - P)), min(Pz), min(Pz) > 0)

## ---- the proposal's sequential binary partition -----------------------------
SBP <- list(
  b1 = list(n = c("neut", "mono"),  d = c("cd4n", "cd4m", "NK", "B"), lab = "myeloid vs lymphoid"),
  b2 = list(n = "neut",             d = "mono",                        lab = "neutrophils vs monocytes"),
  b3 = list(n = c("cd4n", "cd4m"),  d = c("NK", "B"),                  lab = "T vs NK+B"),
  b4 = list(n = "cd4n",             d = "cd4m",                        lab = "CD4 naive vs CD4 memory"),
  b5 = list(n = "NK",               d = "B",                           lab = "NK vs B"))

# A valid SBP must split a previously-formed group each time and end with every
# part isolated. Check it rather than trust it.
groups <- list(colnames(Pz))
for (b in names(SBP)) {
  s <- SBP[[b]]; both <- c(s$n, s$d)
  hit <- which(vapply(groups, function(g) setequal(g, both), logical(1)))
  if (!length(hit)) stop("SBP step ", b, " does not split an existing group")
  groups <- c(groups[-hit], list(s$n), list(s$d))
}
stopifnot(all(lengths(groups) == 1))
say("\nsequential binary partition is valid: every split divides an existing group,\n")
say("  and the %d parts end fully separated -> exactly %d balances\n", ncol(Pz), length(SBP))
for (b in names(SBP))
  say("    %s  %-12s vs %-16s  %s\n", b, paste(SBP[[b]]$n, collapse = "+"),
      paste(SBP[[b]]$d, collapse = "+"), SBP[[b]]$lab)

psi <- t(sapply(SBP, function(b) {
  v <- setNames(numeric(ncol(Pz)), colnames(Pz))
  r <- length(b$n); s <- length(b$d)
  v[b$n] <- sqrt(s / (r * (r + s))); v[b$d] <- -sqrt(r / (s * (r + s))); v }))
orth <- max(abs(unname(psi %*% t(psi)) - diag(nrow(psi))))
say("\n  contrast matrix: rows sum to zero (max %.1g), orthonormal (max |psi psi' - I| = %.2g) -> %s\n",
    max(abs(rowSums(psi))), orth, if (orth < 1e-12) "PASS" else "FAIL")
stopifnot(orth < 1e-12, max(abs(rowSums(psi))) < 1e-12)

ILR <- log(Pz) %*% t(psi)
colnames(ILR) <- names(SBP); rownames(ILR) <- m$title

# losslessness: 5 balances must reconstruct the 6 proportions exactly
rec <- exp(ILR %*% psi); rec <- rec / rowSums(rec)
say("  reversibility check: max |reconstructed - original| = %.2g -> %s\n",
    max(abs(rec - Pz)), if (max(abs(rec - Pz)) < 1e-10) "lossless" else "FAIL")
stopifnot(max(abs(rec - Pz)) < 1e-10)

S <- data.frame(balance = colnames(ILR),
                contrast = vapply(SBP, function(b) sprintf("%s vs %s",
                  paste(b$n, collapse = "+"), paste(b$d, collapse = "+")), character(1)),
                meaning = vapply(SBP, `[[`, character(1), "lab"),
                mean = round(colMeans(ILR), 4), sd = round(apply(ILR, 2, sd), 4))
cat("\n"); print(S, row.names = FALSE)

## ---- subjects, covariates, genes: identical to 01_prep ----------------------
s  <- p$sel$LI
mm <- p$meta[s$keep, , drop = FALSE]
mm$lithium01 <- as.integer(s$grp == "case")
stopifnot(nrow(mm) == 226, sum(mm$lithium01) == 152)
stopifnot(length(unique(mm$group)) == 1)          # assessment group constant

COV <- c("age", "sex", "tobacco", "rin", "plate", "seqpc1", "seqpc2", "seqpc3")
mm$sex <- droplevels(factor(mm$sex)); mm$plate <- droplevels(factor(mm$plate))
stopifnot(!any(is.na(mm[, COV])))
Ccols <- model.matrix(as.formula(paste("~", paste(COV, collapse = " + "))), data = mm)[, -1, drop = FALSE]

M <- ILR[mm$title, , drop = FALSE]
stopifnot(identical(rownames(M), mm$title), ncol(M) == 5, !any(is.na(M)))
say("\nmediator for the 226 subjects: %d x %d | max |r| between balances %.3f\n",
    nrow(M), ncol(M), max(abs(cor(M)[upper.tri(cor(M))])))

lm22 <- unique(toupper(trimws(read.delim(file.path(ROOT, "data/LM22.txt"),
                                         check.names = FALSE, stringsAsFactors = FALSE)[[1]])))
e2s <- read.delim(file.path(ROOT, "de_analysis/data/ens2sym.tsv"), header = FALSE,
                  col.names = c("ens", "sym"), stringsAsFactors = FALSE)
sym <- setNames(toupper(e2s$sym), e2s$ens)
allg <- rownames(p$counts)
is_marker <- sym[allg] %in% lm22
gene_main <- allg[!is_marker]; gene_mark <- allg[is_marker]
say("genes: %d primary (non-LM22) + %d markers = %d\n",
    length(gene_main), length(gene_mark), length(allg))
stopifnot(length(gene_main) == 12018, length(gene_mark) == 350)

## ---- voom on the centred design that will actually be fitted ----------------
source("scripts/med.R")
cM <- centre(M); cC <- centre(Ccols)
Xd <- med_design_full(mm$lithium01, cM$X, cC$X)
say("outcome design: %d x %d | full rank: %s\n", nrow(Xd), ncol(Xd),
    qr(Xd)$rank == ncol(Xd))
stopifnot(qr(Xd)$rank == ncol(Xd))

mkvoom <- function(genes) {
  dge <- normLibSizes(DGEList(p$counts[genes, mm$title, drop = FALSE]), method = "TMM")
  voom(dge, Xd)
}
v_main <- mkvoom(gene_main); v_mark <- mkvoom(gene_mark)
say("voom: primary %s | markers %s | weights %.3f-%.3f | all finite: %s\n",
    paste(dim(v_main$E), collapse = "x"), paste(dim(v_mark$E), collapse = "x"),
    min(v_main$weights), max(v_main$weights),
    all(is.finite(v_main$weights)) && all(is.finite(v_main$E)))
stopifnot(all(is.finite(v_main$E)), all(is.finite(v_main$weights)))

saveRDS(list(meta = mm, lithium = mm$lithium01, M = M, C = Ccols,
             v_main = v_main, v_mark = v_mark,
             gene_main = gene_main, gene_mark = gene_mark,
             covariates = COV, lm22 = lm22, psi = psi, parts = P, parts_z = Pz,
             sbp = SBP, summary = S),
        "data/med_input_proposal.rds")
write.csv(S, "results/ilr_proposal_summary.csv", row.names = FALSE)
say("\nwrote data/med_input_proposal.rds and results/ilr_proposal_summary.csv\n")
