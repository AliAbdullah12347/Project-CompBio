#!/usr/bin/env Rscript
# ==============================================================================
# common.R -- the single differential-expression engine every experiment uses.
#
# The point of putting it here is that when two experiments disagree, the
# disagreement is the thing being varied and never an accident of how the model
# was assembled. Each experiment changes ONE argument to de_fit().
#
# source() this, then call de_fit(). Nothing here writes to disk.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma); library(statmod)})

EXP_HOME <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation"

load_prep <- function() readRDS(file.path(EXP_HOME, "data", "prep.rds"))

# ------------------------------------------------------------------ filtering
# More than `min_count` reads in at least `min_prop` of the samples IN THIS
# SUBSET. It is a cross-sample operation, so a different sample set gives a
# different gene count -- that is correct, not a bug, and freezing the list
# across subsets would leak information between them.
filter_genes <- function(counts, min_count = 10, min_prop = 0.90) {
  keep <- rowSums(counts > min_count) >= min_prop * ncol(counts)
  rownames(counts)[keep]
}

# --------------------------------------------------------------------- design
# Builds the model matrix for a contrast. `tobcol` picks which of the 20
# completed tobacco columns to use (or "obs" for observed-only, which leaves
# 30 NA samples that get dropped). `adjust` is a numeric matrix of extra
# covariates -- ILR balances, surrogate variables, whatever the experiment is
# testing. Aliased columns are detected and removed rather than left to make
# the fit silently rank-deficient.
build_design <- function(prep, contrast, tobcol = "tobacco_imp_01",
                         adjust = NULL, drop_covars = character(0),
                         samples = NULL) {
  cs <- prep$CONTRASTS[[contrast]]
  keep <- cs$keep
  if (!is.null(samples)) keep <- keep & samples
  m <- prep$meta[keep, , drop = FALSE]

  tob <- if (identical(tobcol, "obs")) m$tob_obs else prep$TOB[m$title, tobcol]
  m$tob <- as.numeric(tob)

  covars <- setdiff(cs$covars, drop_covars)
  # A covariate with one level in this subset contributes nothing and makes the
  # matrix singular. Report it rather than letting the fit fail obscurely.
  const <- covars[vapply(covars, function(v)
    length(unique(m[[v]][!is.na(m[[v]])])) < 2, logical(1))]
  covars <- setdiff(covars, const)

  ok <- stats::complete.cases(m[, c(cs$exposure, covars), drop = FALSE])
  m <- m[ok, , drop = FALSE]
  keep_idx <- which(keep)[ok]

  form <- stats::as.formula(paste("~", cs$exposure, "+",
                                  paste(c("1", covars), collapse = " + ")))
  X <- stats::model.matrix(form, data = m)

  if (!is.null(adjust)) {
    A <- as.matrix(adjust)[m$title, , drop = FALSE]
    colnames(A) <- make.names(colnames(A))
    X <- cbind(X, A)
  }
  ne <- limma::nonEstimable(X)
  if (!is.null(ne)) X <- X[, setdiff(colnames(X), ne), drop = FALSE]

  coef <- grep(paste0("^", cs$exposure), colnames(X), value = TRUE)[1]
  list(X = X, coef = coef, samples = keep_idx, meta = m,
       dropped_constant = const, dropped_aliased = ne, n = nrow(m))
}

# ------------------------------------------------------------------------ fit
# `method` selects the estimator. voom is the reference; the others exist so
# that "does this result depend on the method?" can be answered rather than
# assumed.
de_fit <- function(prep, contrast, tobcol = "tobacco_imp_01", adjust = NULL,
                   drop_covars = character(0), samples = NULL,
                   min_count = 10, min_prop = 0.90,
                   norm = "TMM", method = "voom", coef_override = NULL,
                   permute_exposure = NULL, robust = FALSE) {

  d <- build_design(prep, contrast, tobcol, adjust, drop_covars, samples)
  X <- d$X
  if (!is.null(permute_exposure)) {
    # Permute the exposure column only, leaving every covariate attached to its
    # own sample. That is the null of "this exposure carries no signal", not
    # the null of "these samples are exchangeable".
    X[, d$coef] <- X[permute_exposure, d$coef]
  }
  cnt <- prep$counts[, d$samples, drop = FALSE]
  genes <- filter_genes(cnt, min_count, min_prop)
  cnt <- cnt[genes, , drop = FALSE]

  dge <- edgeR::DGEList(counts = cnt)
  dge <- if (identical(norm, "none")) dge else edgeR::calcNormFactors(dge, method = norm)
  coef <- coef_override %||% d$coef

  if (method %in% c("voom", "voomWQW")) {
    v <- if (method == "voom") limma::voom(dge, X)
         else limma::voomWithQualityWeights(dge, X)
    fit <- limma::eBayes(limma::lmFit(v, X), robust = robust)
    tt <- limma::topTable(fit, coef = coef, number = Inf, sort.by = "none")
  } else if (method == "trend") {
    lcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
    fit <- limma::eBayes(limma::lmFit(lcpm, X), trend = TRUE, robust = robust)
    tt <- limma::topTable(fit, coef = coef, number = Inf, sort.by = "none")
  } else if (method %in% c("QLF", "LRT")) {
    dge <- edgeR::estimateDisp(dge, X, robust = robust)
    if (method == "QLF") {
      f <- edgeR::glmQLFit(dge, X, robust = robust)
      r <- edgeR::glmQLFTest(f, coef = coef)
    } else {
      f <- edgeR::glmFit(dge, X)
      r <- edgeR::glmLRT(f, coef = coef)
    }
    tt <- edgeR::topTags(r, n = Inf, sort.by = "none")$table
    names(tt)[names(tt) == "PValue"] <- "P.Value"
    names(tt)[names(tt) == "FDR"]    <- "adj.P.Val"
  } else stop("unknown method: ", method)

  tt$gene <- rownames(tt)
  list(tt = tt, design = d, n_genes = length(genes), n_samples = d$n,
       coef = coef, method = method, norm = norm)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# --------------------------------------------------------------------- counts
n_deg <- function(tt, fdr = 0.05, lfc = 0)
  sum(tt$adj.P.Val < fdr & abs(tt$logFC) >= lfc, na.rm = TRUE)

deg_ids <- function(tt, fdr = 0.05, lfc = 0)
  tt$gene[which(tt$adj.P.Val < fdr & abs(tt$logFC) >= lfc)]

# Proportion of true nulls (Storey). pi0 near 1 means almost nothing is real;
# pi0 well below 1 means signal is present even when few genes clear FDR.
pi0_storey <- function(p, lambda = 0.5) {
  p <- p[!is.na(p)]
  min(1, mean(p > lambda) / (1 - lambda))
}

# Write a run's summary so every experiment leaves the same machine-readable
# trace behind, regardless of what it varied.
write_run <- function(name, df, notes = NULL) {
  dir <- file.path(EXP_HOME, "runs", name)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(df, file.path(dir, "summary.csv"), row.names = FALSE)
  if (!is.null(notes)) writeLines(notes, file.path(dir, "NOTES.md"))
  cat(sprintf("\n[%s] wrote runs/%s/summary.csv  (%d rows)\n", name, name, nrow(df)))
  invisible(dir)
}
