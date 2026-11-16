#!/usr/bin/env Rscript
# Inputs 4 and 5 of the proposal: cell-type-specific expression estimates,
# computed per split so they can enter the prediction contest without leaking.
#
# WHY THIS IS A SEPARATE STAGE
#
# bMIND is the expensive step in the whole arm -- the DE arm's single run over
# 474 samples took 187 minutes. Permuting labels does not change a
# deconvolution, so running it inside 11_run_bpd_split.R would recompute the
# same profiles for every permutation item. This stage runs once per split and
# caches the result; the prediction runner reads it.
#
# LEAKAGE
#
# bMIND estimates its prior from the samples it is given, so it is a
# cross-sample operation and cannot be run once on the cohort. Train and test
# folds are therefore deconvolved SEPARATELY: each fold's profiles depend only
# on that fold. The labels never enter either call -- bMIND is run without the
# phenotype argument -- so no outcome information can cross the boundary.
#
# This is the same compromise the 75-split CIBERSORTx fractions already use.
# It is weaker than fitting a prior on train and applying it to test, which is
# what the leakage rule literally asks for; bMIND's interface does not expose a
# fitted prior that can be applied to new samples, and improvising one would
# mean reporting a method I cannot state exactly. The honest description is:
# each fold is deconvolved using only its own expression, and no labels.
#
# GENE UNIVERSE
#
# --ngenes restricts to the highest-variance genes IN THE TRAINING FOLD, which
# is leak-free but is a real difference from inputs 1-3, which see every gene.
# It exists because the full set does not fit the cluster's 24-hour job
# deadline. Whatever is used is recorded in the output.
#
# Usage:
#   Rscript 19_bmind_split.R --split 7 --cohort lithium [--ngenes 4000]
#                            [--impl ..] [--outdir ../prediction/cache/ct]

suppressPackageStartupMessages({ library(MIND); library(edgeR) })

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(nm, d = NULL) { i <- match(paste0("--", nm), args)
  if (is.na(i) || i == length(args)) d else args[i + 1] }
SPLIT  <- as.integer(getarg("split"))
COHORT <- getarg("cohort", "lithium")
NGENES <- as.integer(getarg("ngenes", "4000"))
impl   <- getarg("impl", "..")
OUTDIR <- getarg("outdir", file.path(impl, "prediction/cache/ct"))
NCORE  <- as.integer(Sys.getenv("BMIND_CORES", "2"))
stopifnot(!is.na(SPLIT), COHORT %in% c("lithium", "groupA", "full", "groupA_plate"))

CACHE <- if (COHORT == "lithium") "lithium_inputs.rds" else "bpd_inputs.rds"
cache <- readRDS(file.path(impl, "prediction/cache", CACHE))
sp <- read.csv(file.path(impl, sprintf("prediction/splits/splits_bpd_%s.csv", COHORT)),
               stringsAsFactors = FALSE)
s <- sp[sp$split_id == SPLIT, ]
stopifnot(nrow(s) > 0)
ord <- match(s$sample_id, cache$meta$title); stopifnot(!any(is.na(ord)))
meta <- cache$meta[ord, ]
cnt  <- cache$counts[, s$sample_id, drop = FALSE]
tr   <- s$assignment == "train"; te <- !tr

## ---- gene branch, training fold only -----------------------------------------
keep  <- rowSums(cnt[, tr, drop = FALSE] > 10) >= ceiling(0.90 * sum(tr))
cnt_k <- cnt[keep, , drop = FALSE]
lib   <- colSums(cnt_k)
f75   <- apply(cnt_k, 2, function(v) quantile(v / sum(v), 0.75))
refid <- colnames(cnt_k)[tr][which.min(abs(f75[tr] - mean(f75[tr])))]
nf    <- calcNormFactors(cnt_k, method = "TMM", refColumn = match(refid, colnames(cnt_k)))
nf    <- nf / exp(mean(log(nf[tr])))
E     <- log2(t(t(cnt_k + 0.5) / (lib * nf + 1)) * 1e6)
dimnames(E) <- dimnames(cnt_k)

if (NGENES > 0 && nrow(E) > NGENES) {
  v <- apply(E[, tr, drop = FALSE], 1, var)          # train fold only
  E <- E[order(v, decreasing = TRUE)[seq_len(NGENES)], , drop = FALSE]
}

## ---- lineage fractions --------------------------------------------------------
# Five lineages, not the proposal's six-part tree: bMIND needs every part to
# carry enough signal to estimate a profile, and splitting CD4 into naive and
# memory leaves parts too small to deconvolve stably. The composition question
# is input 2's job; this is about which lineage carries the expression.
LIN <- list(gran = c("neutrophils", "eosinophils", "mast.cells.resting"),
            mono = c("monocytes", "macrophages.m0", "macrophages.m2",
                     "dendritic.cells.resting", "dendritic.cells.activated"),
            T    = c("t.cells.cd8", "t.cells.cd4.naive", "t.cells.cd4.memory.resting",
                     "t.cells.cd4.memory.activated", "t.cells.gamma.delta"),
            NK   = c("nk.cells.resting", "nk.cells.activated"),
            B    = c("b.cells.naive", "b.cells.memory", "plasma.cells"))
L <- sapply(LIN, function(cs) rowSums(as.matrix(meta[, cs, drop = FALSE])))
L <- L / rowSums(L)
rownames(L) <- meta$title

## ---- deconvolve each fold on its own ------------------------------------------
run_fold <- function(idx, label) {
  t0 <- Sys.time()
  r <- bMIND(E[, idx, drop = FALSE], frac = L[idx, , drop = FALSE], ncore = NCORE)
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("  %-5s %d samples x %d genes -> %s  [%.1f min]\n",
              label, sum(idx), nrow(E), paste(dim(r$A), collapse = "x"), el))
  r$A                                   # genes x celltypes x samples
}
A_tr <- run_fold(tr, "train")
A_te <- run_fold(te, "test")

# Flatten to samples x (gene x celltype), the shape a classifier wants.
flat <- function(A) {
  g <- dim(A)[1]; k <- dim(A)[2]; n <- dim(A)[3]
  M <- matrix(aperm(A, c(3, 1, 2)), nrow = n)
  colnames(M) <- paste(rep(dimnames(A)[[1]], times = k),
                       rep(dimnames(A)[[2]], each = g), sep = "__")
  rownames(M) <- dimnames(A)[[3]]
  M
}
X <- rbind(flat(A_tr), flat(A_te))
X <- X[s$sample_id, , drop = FALSE]           # restore split order
stopifnot(all(is.finite(X)), identical(rownames(X), s$sample_id))

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
f <- file.path(OUTDIR, sprintf("ct_%s_g%d_split%03d.rds", COHORT, NGENES, SPLIT))
saveRDS(list(X = X, lineages = colnames(L), genes = rownames(E),
             ngenes = nrow(E), split = SPLIT, cohort = COHORT),
        f, compress = "gzip")
cat(sprintf("split %d (%s): %d samples x %d features -> %s\n",
            SPLIT, COHORT, nrow(X), ncol(X), basename(f)))
