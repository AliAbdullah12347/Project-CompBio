#!/usr/bin/env Rscript
# One split of the bipolar-status contest: preprocess inside the training
# fold, tune, fit, score the test fold once.
#
# Preprocessing happens here rather than in a cached matrix on purpose. The
# gene filter, TMM factors, log-CPM scaling and residualisation coefficients
# are all cross-sample quantities; computing them once and reusing them across
# folds is the most common way a prediction result ends up overstated. Doing
# them in this script means the test fold cannot touch any of them.
#
# Usage:
#   Rscript 11_run_bpd_split.R --split 7 --cohort groupA [--perm 0]
#                              [--impl ..] [--outdir ../results/bpd]
#                              [--models enet,rf] [--features all]
#
# --perm k (k > 0) permutes the labels across the whole cohort with seed
# 481+k before the split is applied, then runs the identical pipeline. That is
# the null: same cohort, same CV, same feature count, no real association.

suppressPackageStartupMessages({ library(edgeR); library(glmnet); library(randomForest) })

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(nm, default = NULL) {
  i <- match(paste0("--", nm), args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}
SPLIT   <- as.integer(getarg("split"))
COHORT  <- getarg("cohort", "groupA")
PERM    <- as.integer(getarg("perm", "0"))
impl    <- getarg("impl", "..")
OUTDIR  <- getarg("outdir", file.path(impl, "prediction/results/bpd"))
MODELS  <- strsplit(getarg("models", "enet,rf"), ",")[[1]]
FEATSEL <- getarg("features", "all")
COVSET  <- getarg("covset", "full")
stopifnot(!is.na(SPLIT), COHORT %in% c("groupA", "full", "groupA_plate", "lithium"))
# The lithium contest uses its own cache (226 BP1, is_case = lithium user) but
# the identical pipeline otherwise, so any difference in the answer is about
# the contrast rather than about how it was analysed.
CACHE <- if (COHORT == "lithium") "lithium_inputs.rds" else "bpd_inputs.rds"

MIN_COUNT <- 10
MIN_PROP  <- 0.90
NTREE     <- 500
RF_KEEP   <- 2000   # top-variance genes given to RF, selected within the training fold
ALPHAS    <- c(0.1, 0.5, 0.9, 1.0)
set.seed(481 + SPLIT + 1000 * PERM)

## ---- data ------------------------------------------------------------------
cache <- readRDS(file.path(impl, "prediction/cache", CACHE))
sp <- read.csv(file.path(impl, sprintf("prediction/splits/splits_bpd_%s.csv", COHORT)),
               stringsAsFactors = FALSE)
s <- sp[sp$split_id == SPLIT, ]
stopifnot(nrow(s) > 0)

ord <- match(s$sample_id, cache$meta$title)
stopifnot(!any(is.na(ord)))
meta <- cache$meta[ord, ]
cnt  <- cache$counts[, s$sample_id, drop = FALSE]
y    <- s$y
stopifnot(identical(as.integer(meta$is_case), as.integer(y)))

# The null. Permuting the whole label vector and then applying the unchanged
# split preserves class balance and fold structure, so the only thing removed
# is the association itself.
if (PERM > 0) {
  set.seed(481 + PERM)
  y <- sample(y)
}

tr <- s$assignment == "train"
te <- !tr
stopifnot(sum(tr) > 0, sum(te) > 0, length(unique(y[tr])) == 2)

## ---- gene branch, fitted on train only --------------------------------------
keep  <- rowSums(cnt[, tr, drop = FALSE] > MIN_COUNT) >= ceiling(MIN_PROP * sum(tr))
cnt_k <- cnt[keep, , drop = FALSE]

lib    <- colSums(cnt_k)
f75    <- apply(cnt_k, 2, function(v) quantile(v / sum(v), 0.75))
ref_id <- colnames(cnt_k)[tr][which.min(abs(f75[tr] - mean(f75[tr])))]
nf     <- calcNormFactors(cnt_k, method = "TMM", refColumn = match(ref_id, colnames(cnt_k)))
nf     <- nf / exp(mean(log(nf[tr])))     # recentre on train, not on all columns
E      <- log2(t(t(cnt_k + 0.5) / (lib * nf + 1)) * 1e6)
dimnames(E) <- dimnames(cnt_k)

## ---- covariate design --------------------------------------------------------
# Which covariates belong here is a decision, not a default. A confounder has
# to be causally prior to the exposure, and lithium was prescribed years
# before the blood was drawn -- so RIN, plate, library size and the sequencing
# PCs, all created by the assay, cannot confound it. They can only add noise.
# Of the variables that are prior, only age and sex are imbalanced with
# lithium (|SMD| 0.46 and 0.24; tobacco is 0.08). See 18_covariate_audit.R,
# which also shows the full block scoring 0.571 against age alone at 0.613.
#
#   confounder       age + sex           -- primary for the lithium contrast
#   confounder_tech  + RIN, plate, seqPC -- sensitivity: noise reduction only
#   full             everything          -- the original block, kept so the
#                                          earlier results stay reproducible
#
# The case-control contrast is a different matter: there the assay variables
# ARE imbalanced (plate chi-sq p = 1.6e-09), because group B was sequenced
# apart, so `full` remains its default.
COVSETS <- list(
  confounder      = c("age", "sexM"),
  confounder_tech = c("age", "sexM", "rin", "seqpc1", "seqpc2", "seqpc3", "PLATE"),
  full            = c("age", "sexM", "tobacco", "rin", "seqpc1", "seqpc2", "seqpc3", "PLATE"))
stopifnot(COVSET %in% names(COVSETS))

# Factor levels are taken from the training fold so the test matrix has exactly
# the same columns; a plate present only in test would otherwise shift them.
mk_cov <- function(m, lv, want) {
  d <- data.frame(age = as.numeric(m$age), tobacco = as.integer(m$tobacco),
                  rin = as.numeric(m$rin), seqpc1 = m$seqpc1,
                  seqpc2 = m$seqpc2, seqpc3 = m$seqpc3,
                  sexM = as.integer(m$sex == "M"))
  d <- d[, intersect(want, names(d)), drop = FALSE]
  if ("PLATE" %in% want)
    for (p in lv$plate[-1]) d[[paste0("plate_", p)]] <- as.integer(m$plate == p)
  as.matrix(d)
}
lv <- list(plate = sort(unique(as.character(meta$plate[tr]))))
Cov <- mk_cov(meta, lv, COVSETS[[COVSET]])
stopifnot(!any(is.na(Cov)), ncol(Cov) > 0)

## ---- cell composition --------------------------------------------------------
# The deposited classic CIBERSORT fractions, summed onto a haematopoietic
# lineage tree and turned into four ILR balances. Fractions sum to 1, so they
# cannot all enter a model; the balances are the compositional way round that.
#
# This is NOT a covariate set. Lithium raises neutrophils (0.378 -> 0.441 here,
# p = 2.3e-04), so composition sits on the causal path from the drug to the
# expression. Adjusting a confounder removes bias; adjusting a mediator
# removes the effect. It is carried as its own feature set so the question it
# answers -- how much of the signal is cell mixture rather than transcription
# -- is asked out loud rather than by accident.
LIN <- list(gran = c("neutrophils", "eosinophils", "mast.cells.resting"),
            mono = c("monocytes", "macrophages.m0", "macrophages.m2",
                     "dendritic.cells.resting", "dendritic.cells.activated"),
            T    = c("t.cells.cd8", "t.cells.cd4.naive", "t.cells.cd4.memory.resting",
                     "t.cells.cd4.memory.activated", "t.cells.gamma.delta"),
            NK   = c("nk.cells.resting", "nk.cells.activated"),
            B    = c("b.cells.naive", "b.cells.memory", "plasma.cells"))
L <- sapply(LIN, function(cs) rowSums(as.matrix(meta[, cs, drop = FALSE])))
L <- L / rowSums(L)

# Multiplicative zero replacement, delta = 0.65. The detection limit is the
# smallest non-zero value per lineage, and that minimum is a cross-sample
# quantity -- so it is taken over the TRAINING fold only and then applied to
# both. One BP1 subject has a zero B-cell fraction; this is the rule that
# keeps its log finite without the test fold contributing to it.
lim <- apply(L[tr, , drop = FALSE], 2, function(v) { nz <- v[v > 0]
  if (length(nz)) min(nz) * 0.65 else 1e-6 })
z   <- L <= 0
ins <- rowSums(sweep(z, 2, lim, "*"))
L[z] <- matrix(rep(lim, each = nrow(L)), nrow(L))[z]
L[!z] <- (L * (1 - ins))[!z]
L <- L / rowSums(L)

# Sequential binary partition. b1 splits myeloid from lymphoid, which is where
# lithium acts, so it is the coordinate to read first.
SBP <- list(b1 = list(n = c("gran", "mono"), d = c("T", "NK", "B")),
            b2 = list(n = "gran",            d = "mono"),
            b3 = list(n = "T",               d = c("NK", "B")),
            b4 = list(n = "NK",              d = "B"))
psi <- t(sapply(SBP, function(b) {
  v <- setNames(numeric(ncol(L)), colnames(L))
  r <- length(b$n); s <- length(b$d)
  v[b$n] <- sqrt(s / (r * (r + s))); v[b$d] <- -sqrt(r / (s * (r + s))); v }))
stopifnot(max(abs(unname(psi %*% t(psi)) - diag(nrow(psi)))) < 1e-12)
Cells <- log(L) %*% t(psi)
colnames(Cells) <- names(SBP)
stopifnot(all(is.finite(Cells)))

## ---- feature sets ------------------------------------------------------------
Xg <- t(E)                                   # samples x genes
rsd <- function(X) {                         # residualise on train-fitted coefficients
  B <- qr.solve(cbind(1, Cov[tr, , drop = FALSE]), X[tr, , drop = FALSE])
  X - cbind(1, Cov) %*% B
}
Xr <- rsd(Xg)
stopifnot(max(abs(colMeans(Xr[tr, , drop = FALSE]))) < 1e-8)

FEAT <- list(covariates   = Cov,
             genes        = Xg,
             genes_cov    = cbind(Xg, Cov),
             genes_resid  = Xr,
             cells        = Cells,
             cells_cov    = cbind(Cells, Cov),
             genes_cells  = cbind(Xg, Cells),
             # genes with the cell mixture projected out: what is left is a
             # change in transcription that the mixture cannot explain.
             genes_decell = Xg - cbind(1, Cells) %*%
               qr.solve(cbind(1, Cells[tr, , drop = FALSE]), Xg[tr, , drop = FALSE]))
if (FEATSEL != "all") FEAT <- FEAT[strsplit(FEATSEL, ",")[[1]]]
stopifnot(length(FEAT) > 0, !any(sapply(FEAT, is.null)))

## ---- helpers -----------------------------------------------------------------
auc <- function(score, lab) {
  r <- rank(score); n1 <- sum(lab == 1); n0 <- sum(lab == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[lab == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

fit_enet <- function(Xtr, ytr, Xte) {
  best <- NULL
  for (a in ALPHAS) {
    cv <- cv.glmnet(Xtr, ytr, family = "binomial", alpha = a, nfolds = 5,
                    type.measure = "auc", standardize = TRUE)
    sc <- max(cv$cvm, na.rm = TRUE)
    if (is.null(best) || sc > best$score) best <- list(score = sc, cv = cv, alpha = a)
  }
  list(pred = as.numeric(predict(best$cv, Xte, s = "lambda.min", type = "response")),
       tune = sprintf("alpha=%.1f lambda=%.4g", best$alpha, best$cv$lambda.min))
}

fit_rf <- function(Xtr, ytr, Xte) {
  # Random forest cannot be run on 12,000 features at this sample size within
  # any sane budget, so it receives the RF_KEEP highest-variance columns.
  # Variance is computed on the TRAINING fold only, so this is leak-free -- but
  # it is a real difference between the two models, not only a speed-up:
  # elastic net sees every gene and does its own selection through the L1
  # penalty, random forest sees a pre-filtered subset.
  if (ncol(Xtr) > RF_KEEP) {
    v <- apply(Xtr, 2, var)
    sel <- order(v, decreasing = TRUE)[seq_len(RF_KEEP)]
    Xtr <- Xtr[, sel, drop = FALSE]; Xte <- Xte[, sel, drop = FALSE]
  }
  # mtry is chosen on the training fold's out-of-bag AUC. OOB is computed from
  # the bootstrap samples each tree did not see, so it is an honest held-out
  # estimate and needs no inner CV.
  p <- ncol(Xtr)
  grid <- unique(pmax(1, pmin(p, round(c(sqrt(p), 2 * sqrt(p), 4 * sqrt(p))))))
  w <- c(sum(ytr == 1), sum(ytr == 0)) / length(ytr)   # inverse-prevalence weights
  best <- NULL
  f <- factor(ytr, levels = c(0, 1))
  for (mt in grid) {
    rf <- randomForest(Xtr, f, ntree = NTREE, mtry = mt, classwt = w)
    sc <- auc(rf$votes[, "1"], ytr)
    if (is.null(best) || sc > best$score) best <- list(score = sc, mtry = mt)
  }
  rf <- randomForest(Xtr, f, ntree = NTREE, mtry = best$mtry, classwt = w)
  list(pred = as.numeric(predict(rf, Xte, type = "prob")[, "1"]),
       tune = sprintf("mtry=%d oobAUC=%.3f", best$mtry, best$score))
}

## ---- run ---------------------------------------------------------------------
rows <- list(); n <- 0L
for (fs in names(FEAT)) {
  X <- FEAT[[fs]]
  Xtr <- X[tr, , drop = FALSE]; Xte <- X[te, , drop = FALSE]
  for (md in MODELS) {
    t0 <- Sys.time()
    out <- tryCatch(
      if (md == "enet") fit_enet(Xtr, y[tr], Xte) else fit_rf(Xtr, y[tr], Xte),
      error = function(e) { message("  ", fs, "/", md, " failed: ", conditionMessage(e)); NULL })
    if (is.null(out)) next
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    n <- n + 1L
    rows[[n]] <- data.frame(
      cohort = COHORT, perm = PERM, split_id = SPLIT,
      repeat_no = s$repeat_no[1], fold_no = s$fold_no[1],
      feature_set = fs, model = md, covset = COVSET,
      n_features = ncol(X), rf_keep = RF_KEEP,
      n_train = sum(tr), n_test = sum(te), n_genes_kept = nrow(cnt_k),
      sample_id = s$sample_id[te], y = y[te], pred = out$pred,
      tune = out$tune, fit_seconds = round(secs, 1),
      stringsAsFactors = FALSE)
    cat(sprintf("  %-12s %-5s  %6d feat  AUC %.3f  [%s]  %.0fs\n",
                fs, md, ncol(X), auc(out$pred, y[te]), out$tune, secs))
  }
}

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
res <- do.call(rbind, rows)
# The covariate set is in the filename, not only in a column: two runs that
# differ only by covset must never land on the same object and silently
# overwrite each other in the bucket.
f <- file.path(OUTDIR, sprintf("%s_%s_perm%02d_split%03d.csv", COHORT, COVSET, PERM, SPLIT))
write.csv(res, f, row.names = FALSE, quote = FALSE)
cat(sprintf("split %d (%s, perm %d): %d genes kept, %d train / %d test -> %s\n",
            SPLIT, COHORT, PERM, nrow(cnt_k), sum(tr), sum(te), basename(f)))
