#!/usr/bin/env Rscript
# ==============================================================================
# run_de.R -- the reported differential expression analysis, one branch only.
#
# Earlier work explored many variants: three estimators, 31 gene filters, seven
# multiple-testing corrections, adjusted and unadjusted fits. Those runs are
# the EVIDENCE for the choices made here and are kept in ../de_v2/. This script
# is the single pipeline that gets written up.
#
# THE FIVE CHOICES, AND WHERE EACH IS JUSTIFIED
#
#   estimator        limma-voom                 METHODS.md section 4
#   gene filter      >10 counts in >=90%        METHODS.md section 5
#   adjustment       none (unadjusted only)     METHODS.md section 6
#   correction       Benjamini-Hochberg, 0.05   METHODS.md section 7
#   deconvolution    CIBERSORTx only            METHODS.md section 3
#
# FOUR COMPARISONS
#   WB_LI    whole blood, lithium users vs non-users within bipolar I    n=226
#   WB_BPD   whole blood, bipolar I off lithium vs healthy controls      n=308
#   CT_LI    per cell lineage, same contrast as WB_LI                    n=226
#   CT_BPD   per cell lineage, same contrast as WB_BPD                   n=308
#
# WHAT IS REUSED RATHER THAN RECOMPUTED. Two cached objects, both inputs and
# neither a result: ../de_analysis/data/prep.rds (the filtered count matrix,
# covariates and lineage fractions) and ../de_analysis/results/bmind_profiles.rds
# (3.1 hours of deconvolution). Neither depends on a threshold or on the group
# labels, so recomputing them would reproduce the same numbers exactly.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "de_final")); set.seed(481); options(width = 180)
RES <- "results"; dir.create(RES, showWarnings = FALSE, recursive = TRUE)

FDR   <- 0.05                      # the one significance threshold
LINS  <- c("gran", "mono", "T", "NK", "B")

p  <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
bp <- readRDS(file.path(ROOT, "de_analysis/results/bmind_profiles.rds"))
sym <- local({ m <- read.delim(file.path(ROOT, "de_analysis/data/ens2sym.tsv"),
                               header = FALSE, col.names = c("ens", "sym"))
               setNames(m$sym, m$ens) })
A <- bp$A; SE <- bp$SE
stopifnot(identical(dimnames(A)[[1]], rownames(p$counts)))   # same gene set, same order

cat(sprintf("== run_de ==\n%d genes (>10 counts in >=90%% of samples) | %d samples\n",
            nrow(p$counts), ncol(p$counts)))
cat(sprintf("estimator limma-voom | correction BH at FDR %.2f | unadjusted fits only\n\n", FDR))

## ---- design -----------------------------------------------------------------
# Covariates are the same for both levels so that whole blood and cell type
# differ only in what is being measured. Assessment group is dropped from the
# lithium contrast because every bipolar I subject is in group A, making it
# inestimable there; it varies once controls enter.
design_for <- function(contrast) {
  s <- p$sel[[contrast]]
  mm <- p$meta[s$keep, , drop = FALSE]; mm$grp <- s$grp
  mm$plate <- droplevels(factor(mm$plate)); mm$sex <- droplevels(factor(mm$sex))
  if ("group" %in% s$cov) mm$group <- droplevels(factor(mm$group))
  cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]])) > 1, logical(1))]
  X <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse = " + "))), data = mm)
  ne <- nonEstimable(X); if (!is.null(ne)) X <- X[, setdiff(colnames(X), ne), drop = FALSE]
  # A rank-deficient design silently drops a covariate. Fail loudly instead.
  stopifnot(qr(X)$rank == ncol(X), "grpcase" %in% colnames(X))
  list(X = X, samples = mm$title, covars = cv, dropped = setdiff(s$cov, cv),
       n_ctrl = sum(mm$grp == "ctrl"), n_case = sum(mm$grp == "case")) }

# bMIND reports a posterior SE per gene x lineage x sample. Using 1/SE^2 as a
# precision weight carries the deconvolution's own uncertainty into the test;
# treating posterior means as measured data would understate uncertainty and
# inflate the cell-type gene count.
make_weights <- function(S) {
  W <- 1 / (S^2); fin <- is.finite(W) & W > 0
  if (any(!fin)) W[!fin] <- min(W[fin]); W / mean(W) }

## ---- the two levels ---------------------------------------------------------
fit_wb <- function(contrast) {
  d <- design_for(contrast)
  dge <- calcNormFactors(DGEList(p$counts[, d$samples, drop = FALSE]), method = "TMM")
  v <- voom(dge, d$X)
  tt <- topTable(eBayes(lmFit(v, d$X)), coef = "grpcase", number = Inf, sort.by = "none")
  list(d = d, tab = data.frame(
    gene = rownames(tt), lineage = NA_character_, logFC = tt$logFC,
    AveExpr = tt$AveExpr, t = tt$t, P = tt$P.Value, BH = tt$adj.P.Val,
    stringsAsFactors = FALSE)) }

fit_ct <- function(contrast) {
  d <- design_for(contrast)
  tab <- do.call(rbind, lapply(LINS, function(ct) {
    f <- eBayes(lmFit(A[, ct, d$samples], d$X, weights = make_weights(SE[, ct, d$samples])),
                trend = TRUE)
    tt <- topTable(f, coef = "grpcase", number = Inf, sort.by = "none")
    # BH WITHIN each lineage: each lineage is its own family of 12,368 tests,
    # which answers "is there signal in granulocytes?".
    data.frame(gene = rownames(tt), lineage = ct, logFC = tt$logFC,
               AveExpr = tt$AveExpr, t = tt$t, P = tt$P.Value,
               BH = p.adjust(tt$P.Value, "BH"), stringsAsFactors = FALSE) }))
  list(d = d, tab = tab) }

## ---- run --------------------------------------------------------------------
JOBS <- data.frame(name = c("WB_LI", "WB_BPD", "CT_LI", "CT_BPD"),
                   contrast = c("LI", "BPD", "LI", "BPD"),
                   level = c("WB", "WB", "CT", "CT"), stringsAsFactors = FALSE)
OUT <- list(); SUM <- list()
for (i in seq_len(nrow(JOBS))) {
  j <- JOBS[i, ]
  r <- if (j$level == "WB") fit_wb(j$contrast) else fit_ct(j$contrast)
  tb <- r$tab
  tb$symbol <- unname(ifelse(is.na(sym[tb$gene]), tb$gene, sym[tb$gene]))
  tb$significant <- tb$BH < FDR
  # rank by p-value, within lineage for the cell-type level
  tb$rank <- if (j$level == "WB") rank(tb$P, ties.method = "min") else
    unlist(lapply(split(tb$P, tb$lineage), rank, ties.method = "min"))[order(order(tb$lineage))]
  tb <- tb[, c("symbol", "gene", "lineage", "rank", "logFC", "AveExpr", "t", "P", "BH", "significant")]
  OUT[[j$name]] <- tb

  # One extra number for the cell-type level: BH pooled over all 5 x 12,368
  # tests, collapsed to unique genes. This is the only fair basis for comparing
  # the cell-type level against whole blood, because per-lineage BH gives the
  # cell-type side five independent chances at significance.
  pooled <- NA_integer_
  if (j$level == "CT") {
    q <- p.adjust(tb$P, "BH"); pooled <- length(unique(tb$gene[q < FDR])) }

  SUM[[j$name]] <- data.frame(
    comparison = j$name, level = j$level, contrast = j$contrast,
    n = length(r$d$samples), n_group1 = r$d$n_ctrl, n_group2 = r$d$n_case,
    n_tests = nrow(tb), n_genes = length(unique(tb$gene)),
    significant_BH05 = sum(tb$significant),
    pooled_unique_genes = pooled,
    min_P = signif(min(tb$P), 3), min_BH = signif(min(tb$BH), 3),
    max_abs_logFC = round(max(abs(tb$logFC)), 3),
    covariates = paste(r$d$covars, collapse = ";"),
    design_rank = qr(r$d$X)$rank, design_cols = ncol(r$d$X), stringsAsFactors = FALSE)
  cat(sprintf("  %-7s n=%3d (%3d vs %3d) | %6d tests | %5d significant | min p %.3g\n",
              j$name, length(r$d$samples), r$d$n_ctrl, r$d$n_case, nrow(tb),
              sum(tb$significant), min(tb$P)))
}
S <- do.call(rbind, SUM)

cat("\n=== summary ===\n")
print(S[, c("comparison", "n", "n_group1", "n_group2", "n_tests", "significant_BH05",
            "pooled_unique_genes", "min_BH", "max_abs_logFC")], row.names = FALSE)

cat("\n=== cell type, significant genes per lineage ===\n")
for (k in c("CT_LI", "CT_BPD")) {
  z <- OUT[[k]]
  cat(sprintf("  %-7s %s\n", k, paste(sprintf("%s=%d", LINS,
      vapply(LINS, function(l) sum(z$significant & z$lineage == l), integer(1))), collapse = "  ")))
}

cat("\n=== covariate check ===\n")
for (i in seq_len(nrow(S)))
  cat(sprintf("  %-7s design %d x %d, rank %d (full rank: %s) | %s\n", S$comparison[i],
              S$n[i], S$design_cols[i], S$design_rank[i],
              S$design_rank[i] == S$design_cols[i], gsub(";", ", ", S$covariates[i])))

## ---- write ------------------------------------------------------------------
for (k in names(OUT))
  write.csv(OUT[[k]][order(OUT[[k]]$P), ], gzfile(file.path(RES, sprintf("DE_%s.csv.gz", k))),
            row.names = FALSE)
SIG <- do.call(rbind, lapply(names(OUT), function(k) {
  z <- OUT[[k]][OUT[[k]]$significant, ]; if (!nrow(z)) return(NULL)
  cbind(comparison = k, z[order(z$P), ]) }))
write.csv(SIG, file.path(RES, "significant_genes.csv"), row.names = FALSE)
write.csv(S, file.path(RES, "summary.csv"), row.names = FALSE)
cat(sprintf("\nwrote %d result files to results/ (%d significant rows in total)\n",
            length(OUT) + 2, nrow(SIG)))
