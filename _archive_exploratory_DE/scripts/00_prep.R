#!/usr/bin/env Rscript
# ==============================================================================
# 00_prep.R -- build every shared object the experiments start from.
#
# Runs once. Everything downstream loads experimentation/data/prep.rds and does
# no preparation of its own, so that a difference between two experiments is a
# difference in the experiment and never a difference in how the data was made.
#
# Reads (read-only, never written to):
#   ../data/cohort_474/counts_474.tsv.gz
#   ../data/cohort_474/metadata_474_imputed.csv
#   ../data/cohort_474/fractions_bmode_474_CIBERSORTx.csv
# Writes only inside experimentation/.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})

HERE   <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation"
DATA   <- file.path(dirname(HERE), "data", "cohort_474")
OUT    <- file.path(HERE, "data")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cat("== 00_prep ==\n")

## ---- counts ---------------------------------------------------------------
# check.names = FALSE: every sample title starts with a digit, so make.names
# would silently prefix each with "X" and break every join downstream.
counts <- as.matrix(read.delim(file.path(DATA, "counts_474.tsv.gz"),
                               row.names = 1, check.names = FALSE))
cat(sprintf("counts as read            %d x %d\n", nrow(counts), ncol(counts)))

# ENSGR rows are pseudoautosomal-region Y duplicates of ENSG0* genes. They
# double-count 1,518 reads. Dropped before anything is normalised.
ensgr  <- grepl("^ENSGR", rownames(counts))
stopifnot(sum(ensgr) == 47)
counts <- counts[!ensgr, ]
rownames(counts) <- sub("\\.\\d+$", "", rownames(counts))
stopifnot(!any(duplicated(rownames(counts))))
cat(sprintf("after ENSGR drop          %d x %d\n", nrow(counts), ncol(counts)))

## ---- metadata -------------------------------------------------------------
meta <- read.csv(file.path(DATA, "metadata_474_imputed.csv"),
                 check.names = FALSE, colClasses = c(title = "character"))
stopifnot(identical(meta$title, colnames(counts)))   # order, not just membership

m <- data.frame(
  title    = meta$title,
  dx       = factor(meta[["bipolar disorder diagnosis"]], c("Control", "BP1", "BP2")),
  lithium  = as.integer(meta[["lithium use (non-user=0, user = 1)"]]),
  age      = as.numeric(meta$age),
  sex      = factor(meta$Sex),
  group    = factor(meta[["assessment group"]]),
  rin      = as.numeric(meta$rin),
  plate    = factor(meta[["sequencing plate"]]),
  seqpc1   = as.numeric(meta[["sequencing metric pc1"]]),
  seqpc2   = as.numeric(meta[["sequencing metric pc2"]]),
  seqpc3   = as.numeric(meta[["sequencing metric pc3"]]),
  qc_pass  = as.logical(meta[["included in final analysis"]]),
  depth_ok = as.logical(meta$depth_ok),
  rin_ok   = as.logical(meta$rin_ok),
  tob_obs  = suppressWarnings(as.numeric(meta$tobacco_observed)),
  tob_imp  = as.logical(meta$tobacco_was_imputed),
  stringsAsFactors = FALSE
)
# All 20 completed tobacco columns, so the sensitivity analysis can use them.
tobcols <- grep("^tobacco_imp_\\d+$", names(meta), value = TRUE)
stopifnot(length(tobcols) == 20)
TOB <- as.matrix(meta[, tobcols]); rownames(TOB) <- meta$title
m$case <- factor(ifelse(m$dx == "Control", "Control", "BD"), c("Control", "BD"))
m$lib  <- colSums(counts)

cat(sprintf("metadata                  %d samples | dx %s\n", nrow(m),
            paste(sprintf("%s=%d", levels(m$dx), table(m$dx)), collapse = " ")))
cat(sprintf("tobacco still missing     %d (imputed in 20 completed columns)\n",
            sum(is.na(m$tob_obs))))

## ---- cell composition -----------------------------------------------------
fr <- read.csv(file.path(DATA, "fractions_bmode_474_CIBERSORTx.csv"),
               check.names = FALSE, colClasses = c(Mixture = "character"))
rownames(fr) <- fr$Mixture
qc_cs <- fr[m$title, c("P-value", "Correlation", "RMSE")]
fr    <- as.matrix(fr[m$title, setdiff(names(fr), c("Mixture", "P-value", "Correlation", "RMSE"))])
stopifnot(nrow(fr) == nrow(m), max(abs(rowSums(fr) - 1)) < 1e-6)

# Five lineages. Fine-grained LM22 output is not reliable on this dataset:
# between-method ICC is 0.389 and the Aitchison distance between two
# deconvolution methods is 1.56x the spread between subjects. Aggregating to
# lineage raises median ICC to 0.693 and the myeloid/lymphoid balance to 0.907,
# because LM22's failure mode is confusing neighbouring subtypes and that error
# cancels in the sum.
LIN <- list(
  gran = c("Neutrophils", "Eosinophils", "Mast cells resting", "Mast cells activated"),
  mono = c("Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2",
           "Dendritic cells resting", "Dendritic cells activated"),
  T    = c("T cells CD8", "T cells CD4 naive", "T cells CD4 memory resting",
           "T cells CD4 memory activated", "T cells follicular helper",
           "T cells regulatory (Tregs)", "T cells gamma delta"),
  NK   = c("NK cells resting", "NK cells activated"),
  B    = c("B cells naive", "B cells memory", "Plasma cells"))
stopifnot(setequal(unlist(LIN), colnames(fr)))
L <- sapply(LIN, function(cs) rowSums(fr[, cs, drop = FALSE]))

# Multiplicative simple replacement at 65% of each part's detection limit
# (Martin-Fernandez 2003): preserves the ratios among the non-zero parts, which
# additive replacement would distort, and every balance below is a ratio.
zrepl <- function(X, delta = 0.65) {
  X <- X / rowSums(X)
  lim <- apply(X, 2, function(v) min(v[v > 0]) * delta)
  z <- X <= 0
  ins <- rowSums(sweep(z, 2, lim, "*"))
  X[z] <- matrix(rep(lim, each = nrow(X)), nrow(X))[z]
  X[!z] <- (X * (1 - ins))[!z]
  X / rowSums(X)
}
Lz <- zrepl(L)

# ILR on a sequential binary partition of the five lineages. 5 parts -> 4
# balances. Orthonormality is checked, not assumed: a partition that is not
# sequential-and-binary still produces four numbers, just not coordinates.
sbp <- list(
  b1_myeloid_vs_lymphoid = list(num = c("gran", "mono"), den = c("T", "NK", "B")),
  b2_gran_vs_mono        = list(num = "gran",            den = "mono"),
  b3_T_vs_NKB            = list(num = "T",               den = c("NK", "B")),
  b4_NK_vs_B             = list(num = "NK",              den = "B"))
psi <- t(sapply(sbp, function(b) {
  v <- setNames(numeric(ncol(Lz)), colnames(Lz))
  r <- length(b$num); s <- length(b$den)
  v[b$num] <-  sqrt(s / (r * (r + s)))
  v[b$den] <- -sqrt(r / (s * (r + s)))
  v
}))
gram <- unname(psi %*% t(psi))
stopifnot(max(abs(gram - diag(nrow(psi)))) < 1e-12,   # orthonormal
          max(abs(rowSums(psi)))          < 1e-12)    # contrasts sum to zero
ILR <- log(Lz) %*% t(psi)
colnames(ILR) <- names(sbp)

cat(sprintf("lineages (mean)           %s\n",
            paste(sprintf("%s=%.3f", colnames(L), colMeans(L)), collapse = " ")))
cat(sprintf("ILR balances              %d, SBP verified orthonormal\n", ncol(ILR)))

## ---- contrasts ------------------------------------------------------------
# Each entry: which samples, what the exposure is, and which covariates the
# design may contain. Assessment group is CONSTANT within bipolar I (all 226
# are group A), so it must be dropped there or the design matrix is singular.
CONTRASTS <- list(
  lithium = list(
    keep = m$dx == "BP1", exposure = "lithium",
    covars = c("age", "sex", "tob", "rin", "plate", "seqpc1", "seqpc2", "seqpc3"),
    note = "bipolar I only; assessment group dropped (constant)"),
  casecon = list(
    keep = rep(TRUE, nrow(m)), exposure = "case",
    covars = c("age", "sex", "tob", "group", "rin", "plate", "seqpc1", "seqpc2", "seqpc3"),
    note = "all 474; BD (BP1+BP2) vs Control"),
  lithium_krebs = list(
    keep = m$dx != "Control", exposure = "lithium",
    covars = c("dx", "age", "sex", "tob", "group", "rin", "plate", "seqpc1", "seqpc2", "seqpc3"),
    note = "all cases, diagnosis retained -- reproduces the published contrast")
)
for (nm in names(CONTRASTS)) {
  k <- CONTRASTS[[nm]]$keep
  cat(sprintf("contrast %-14s n=%3d  %s\n", nm, sum(k), CONTRASTS[[nm]]$note))
}

saveRDS(list(counts = counts, meta = m, TOB = TOB, frac = fr, lineage = L,
             lineage_z = Lz, ILR = ILR, cs_qc = qc_cs, CONTRASTS = CONTRASTS,
             psi = psi, LIN = LIN),
        file.path(OUT, "prep.rds"))
cat(sprintf("\nwrote data/prep.rds  (%.1f MB)\n",
            file.size(file.path(OUT, "prep.rds")) / 1e6))
