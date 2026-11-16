#!/usr/bin/env Rscript
# Proposal inputs 4 and 5, scored through the same CV as inputs 1-3.
#
#   input 4  ct_neut  cell-type expression in the granulocyte lineage only
#   input 5  ct_all   cell-type expression in all five lineages together
#
# Reads the per-split profiles that 19_bmind_split.R cached. Those were built
# with train and test folds deconvolved separately and no labels, so the fold
# boundary is already respected before this script sees anything.
#
# "Neutrophils" in the proposal, "gran" here. bMIND was run on five lineages,
# and the granulocyte lineage is neutrophils plus eosinophils and resting mast
# cells -- both near zero in this deposit, so the lineage is neutrophils to
# three decimal places. Splitting them would mean deconvolving parts that
# carry no signal, which is worse than naming the aggregate honestly.
#
# Everything cross-sample left in this script -- feature scaling, the variance
# filter for the forest, hyperparameters -- is fitted on the training fold.
#
# Usage:
#   Rscript 21_ct_predict.R --split 7 --cohort lithium --ngenes 4000
#                           [--perm 0] [--models enet,rf] [--impl ..]

suppressPackageStartupMessages({ library(glmnet); library(randomForest) })

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(nm, d = NULL) { i <- match(paste0("--", nm), args)
  if (is.na(i) || i == length(args)) d else args[i + 1] }
SPLIT  <- as.integer(getarg("split"))
COHORT <- getarg("cohort", "lithium")
NGENES <- as.integer(getarg("ngenes", "4000"))
PERM   <- as.integer(getarg("perm", "0"))
MODELS <- strsplit(getarg("models", "enet,rf"), ",")[[1]]
impl   <- getarg("impl", "..")
CTDIR  <- getarg("ctdir", file.path(impl, "prediction/cache/ct"))
OUTDIR <- getarg("outdir", file.path(impl, "prediction/results/ct"))

NTREE   <- 500
RF_KEEP <- 2000
ALPHAS  <- c(0.1, 0.5, 0.9, 1.0)
set.seed(481 + SPLIT + 1000 * PERM)

ct <- readRDS(file.path(CTDIR, sprintf("ct_%s_g%d_split%03d.rds", COHORT, NGENES, SPLIT)))
X  <- ct$X
sp <- read.csv(file.path(impl, sprintf("prediction/splits/splits_bpd_%s.csv", COHORT)),
               stringsAsFactors = FALSE)
s <- sp[sp$split_id == SPLIT, ]
stopifnot(nrow(s) > 0, identical(rownames(X), s$sample_id))

y <- s$y
if (PERM > 0) { set.seed(481 + PERM); y <- sample(y) }   # same null as inputs 1-3
tr <- s$assignment == "train"; te <- !tr
stopifnot(sum(tr) > 0, sum(te) > 0, length(unique(y[tr])) == 2)

# Column names are gene__lineage, so input 4 is a suffix match.
lin <- sub("^.*__", "", colnames(X))
stopifnot(all(ct$lineages %in% unique(lin)))
FEAT <- list(ct_neut = X[, lin == "gran", drop = FALSE],
             ct_all  = X)

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
  if (ncol(Xtr) > RF_KEEP) {
    v <- apply(Xtr, 2, var)                         # training fold only
    sel <- order(v, decreasing = TRUE)[seq_len(RF_KEEP)]
    Xtr <- Xtr[, sel, drop = FALSE]; Xte <- Xte[, sel, drop = FALSE]
  }
  p <- ncol(Xtr)
  grid <- unique(pmax(1, pmin(p, round(c(sqrt(p), 2 * sqrt(p), 4 * sqrt(p))))))
  w <- c(sum(ytr == 1), sum(ytr == 0)) / length(ytr)
  f <- factor(ytr, levels = c(0, 1)); best <- NULL
  for (mt in grid) {
    rf <- randomForest(Xtr, f, ntree = NTREE, mtry = mt, classwt = w)
    sc <- auc(rf$votes[, "1"], ytr)
    if (is.null(best) || sc > best$score) best <- list(score = sc, mtry = mt)
  }
  rf <- randomForest(Xtr, f, ntree = NTREE, mtry = best$mtry, classwt = w)
  list(pred = as.numeric(predict(rf, Xte, type = "prob")[, "1"]),
       tune = sprintf("mtry=%d oobAUC=%.3f", best$mtry, best$score))
}

rows <- list(); n <- 0L
for (fs in names(FEAT)) {
  Xf <- FEAT[[fs]]
  Xtr <- Xf[tr, , drop = FALSE]; Xte <- Xf[te, , drop = FALSE]
  # Constant columns break glmnet's standardisation; drop on the train fold.
  keep <- apply(Xtr, 2, function(v) var(v) > 0)
  Xtr <- Xtr[, keep, drop = FALSE]; Xte <- Xte[, keep, drop = FALSE]
  for (md in MODELS) {
    t0 <- Sys.time()
    out <- tryCatch(if (md == "enet") fit_enet(Xtr, y[tr], Xte) else fit_rf(Xtr, y[tr], Xte),
                    error = function(e) { message("  ", fs, "/", md, ": ", conditionMessage(e)); NULL })
    if (is.null(out)) next
    n <- n + 1L
    rows[[n]] <- data.frame(
      cohort = COHORT, perm = PERM, split_id = SPLIT,
      repeat_no = s$repeat_no[1], fold_no = s$fold_no[1],
      feature_set = fs, model = md, covset = "none", ngenes = NGENES,
      n_features = ncol(Xtr), n_train = sum(tr), n_test = sum(te),
      sample_id = s$sample_id[te], y = y[te], pred = out$pred,
      tune = out$tune, fit_seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
      stringsAsFactors = FALSE)
    cat(sprintf("  %-8s %-5s %6d feat  AUC %.3f  [%s]\n",
                fs, md, ncol(Xtr), auc(out$pred, y[te]), out$tune))
  }
}
stopifnot(n > 0)
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
f <- file.path(OUTDIR, sprintf("ct_%s_g%d_perm%02d_split%03d.csv", COHORT, NGENES, PERM, SPLIT))
write.csv(do.call(rbind, rows), f, row.names = FALSE, quote = FALSE)
cat(sprintf("split %d (%s, perm %d) -> %s\n", SPLIT, COHORT, PERM, basename(f)))
