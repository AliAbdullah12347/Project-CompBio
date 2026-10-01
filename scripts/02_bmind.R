#!/usr/bin/env Rscript
# ==============================================================================
# 02_bmind.R -- cell-type-specific differential expression via bMIND.
#
# DECISIONS MADE HERE, AND WHY
#
# 1. bmind_de(), not bMIND() followed by limma. bmind_de fits the phenotype
#    INSIDE the Bayesian model, so the cell-type association is tested with the
#    posterior uncertainty of the deconvolution propagated into it. Running
#    limma on bMIND's point estimates would treat shrunken posterior means as
#    if they were measured data, understate the standard errors, and inflate
#    the cell-type DEG count -- which would bias CT_BPD toward "supported"
#    exactly where we least want a thumb on the scale. It is also the function
#    Boltz et al. used, which keeps us comparable to the only published work
#    of this kind.
#
# 2. np = TRUE (non-informative prior). We have no external sorted-cell
#    reference profile for these five lineages. Supplying a wrong prior is
#    worse than supplying none: it would pull every sample toward a profile
#    that does not describe this cohort.
#
# 3. One run per contrast per adjustment -- four runs. The phenotype enters
#    the model, so runs cannot be shared across contrasts.
#
# 4. Genes are processed in chunks and written incrementally. This is a
#    multi-hour job; a crash in hour eight must not cost hours one to seven.
#
# 5. Covariates are passed as `covariate`, which bMIND applies at the bulk
#    level. Factors are expanded to dummies first because bMIND takes a
#    numeric matrix.
#
# Usage:  Rscript scripts/02_bmind.R <contrast> <adjust>
#         contrast in {LI, BPD}   adjust in {raw, ilr}
# ==============================================================================

suppressPackageStartupMessages({library(MIND)})
args <- commandArgs(trailingOnly = TRUE)
CONTRAST <- args[1]; ADJ <- args[2]
stopifnot(CONTRAST %in% c("LI", "BPD"), ADJ %in% c("raw", "ilr"))

HERE <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/de_analysis"
setwd(HERE)
p <- readRDS("data/prep.rds")
OUT <- file.path("results", sprintf("bmind_%s_%s", CONTRAST, ADJ))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

CHUNK <- 400
NCORE <- 6

s   <- p$sel[[CONTRAST]]
mm  <- p$meta[s$keep, , drop = FALSE]
grp <- s$grp
bulk <- p$logtpm[, mm$title, drop = FALSE]
frac <- p$lineage[mm$title, , drop = FALSE]

# Covariate matrix: numeric, no intercept (bMIND adds its own).
cov_terms <- s$cov
mm$plate <- droplevels(mm$plate); mm$sex <- droplevels(mm$sex)
if ("group" %in% cov_terms) mm$group <- droplevels(mm$group)
const <- cov_terms[vapply(cov_terms, function(v) length(unique(mm[[v]])) < 2, logical(1))]
if (length(const)) {
  cat(sprintf("dropping constant covariate(s): %s\n", paste(const, collapse = ", ")))
  cov_terms <- setdiff(cov_terms, const)
}
CV <- model.matrix(as.formula(paste("~", paste(cov_terms, collapse = " + "))), data = mm)[, -1, drop = FALSE]
if (ADJ == "ilr") CV <- cbind(CV, p$ILR[mm$title, , drop = FALSE])
storage.mode(CV) <- "double"

y <- as.numeric(grp == "case")       # 0 = control/non-user, 1 = case/user
cat(sprintf("== bmind_de  contrast=%s adjust=%s ==\n", CONTRAST, ADJ))
cat(sprintf("  samples %d (%d vs %d) | genes %d | lineages %d | covariates %d\n",
            ncol(bulk), sum(y == 0), sum(y == 1), nrow(bulk), ncol(frac), ncol(CV)))
cat(sprintf("  covariates: %s\n", paste(colnames(CV), collapse = ", ")))

genes <- rownames(bulk)
chunks <- split(genes, ceiling(seq_along(genes) / CHUNK))
cat(sprintf("  %d chunks of up to %d genes, ncore=%d\n\n", length(chunks), CHUNK, NCORE))

t_start <- Sys.time()
for (i in seq_along(chunks)) {
  f <- file.path(OUT, sprintf("chunk_%03d.rds", i))
  if (file.exists(f)) { cat(sprintf("  chunk %3d/%d cached\n", i, length(chunks))); next }
  g <- chunks[[i]]
  t0 <- Sys.time()
  r <- tryCatch(
    bmind_de(bulk = bulk[g, , drop = FALSE], frac = frac, y = y,
             covariate = CV, np = TRUE, ncore = NCORE),
    error = function(e) { cat(sprintf("  chunk %3d FAILED: %s\n", i, conditionMessage(e))); NULL })
  if (is.null(r)) next
  # Keep only what the DE needs: the phenotype effect and its p-value per
  # cell type. The full posterior arrays are large and are not reused.
  saveRDS(list(genes = g, pval = r$pval, coef = r$coef, qval = r$qval), f)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  done <- i; tot <- length(chunks)
  eta <- as.numeric(difftime(Sys.time(), t_start, units = "mins")) / done * (tot - done)
  cat(sprintf("  chunk %3d/%d  %5.1fs  ETA %5.1f min\n", i, tot, el, eta))
  flush.console()
}

## ---- assemble --------------------------------------------------------------
fs <- sort(list.files(OUT, "^chunk_.*\\.rds$", full.names = TRUE))
cat(sprintf("\nassembling %d chunks\n", length(fs)))
parts <- lapply(fs, readRDS)
parts <- Filter(function(x) !is.null(x$pval), parts)
P <- do.call(rbind, lapply(parts, function(x) {
  m <- as.matrix(x$pval); rownames(m) <- x$genes; m }))
cat(sprintf("  p-value matrix %d genes x %d lineages\n", nrow(P), ncol(P)))
saveRDS(list(pval = P, contrast = CONTRAST, adjust = ADJ,
             n = ncol(bulk), n_case = sum(y == 1), n_ctrl = sum(y == 0),
             covariates = colnames(CV)),
        file.path("results", sprintf("bmind_%s_%s.rds", CONTRAST, ADJ)))
cat(sprintf("wrote results/bmind_%s_%s.rds  [total %.1f min]\n", CONTRAST, ADJ,
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
