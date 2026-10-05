#!/usr/bin/env Rscript
# ==============================================================================
# filter_sweep_xcheck.R -- two things the main sweep cannot rule out.
#
# 1. METHOD. voom fits one mean-variance trend across whatever genes are in the
#    matrix. A loose filter puts ~25k low-count genes into that lowess, which
#    could move the weights of the well-expressed genes. If so, the filter
#    effect is partly an artefact of voom rather than of multiple testing, and
#    a method whose dispersion estimation is less globally coupled (edgeR's
#    quasi-likelihood F test, or limma-trend) should show a smaller effect.
#
# 2. CONTRAST. The lithium contrast is one 226-sample comparison. If the same
#    filter behaviour appears in the 474-sample case-control contrast, it is a
#    property of the filter; if it does not, it is a property of this contrast.
#
# Caches per fit, so this is resumable.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})
source("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation/scripts/common.R")

RUN   <- file.path(EXP_HOME, "runs", "filter_sweep")
CACHE <- file.path(RUN, "cache_x")
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
p <- load_prep(); FDR <- 0.05

run1 <- function(tag, contrast, mc, mp, method) {
  f <- file.path(CACHE, paste0(tag, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  t0 <- Sys.time()
  r <- tryCatch(de_fit(p, contrast, min_count = mc, min_prop = mp, method = method),
                error = function(e) list(error = conditionMessage(e)))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  out <- if (is.null(r$error))
    list(tag = tag, ok = TRUE, contrast = contrast, min_count = mc, min_prop = mp,
         method = method, n_genes = r$n_genes, n_samples = r$n_samples, secs = secs,
         tt = r$tt[, intersect(c("gene", "logFC", "AveExpr", "logCPM", "t", "F",
                                 "P.Value", "adj.P.Val"), names(r$tt))])
  else list(tag = tag, ok = FALSE, err = r$error, contrast = contrast,
            min_count = mc, min_prop = mp, method = method, secs = secs)
  saveRDS(out, f); out
}

jacc <- function(a, b) length(intersect(a, b)) / length(union(a, b))
rows <- list(); keep <- list()

## ---- 1. method x filter, lithium ------------------------------------------
CELLS <- list(loose = c(0, 0.10), canon = c(10, 0.90))
for (meth in c("voom", "trend", "QLF")) for (cn in names(CELLS)) {
  a <- CELLS[[cn]]; tg <- sprintf("lith_%s_%s", meth, cn)
  r <- run1(tg, "lithium", a[1], a[2], meth)
  cat(sprintf("%-20s genes=%6d DEG=%5s %.0fs\n", tg, r$n_genes %||% NA,
              if (isTRUE(r$ok)) n_deg(r$tt, FDR) else "ERR", r$secs))
  if (isTRUE(r$ok)) keep[[tg]] <- r
}
for (meth in c("voom", "trend", "QLF")) {
  lo <- keep[[sprintf("lith_%s_loose", meth)]]; cn <- keep[[sprintf("lith_%s_canon", meth)]]
  if (is.null(lo) || is.null(cn)) next
  cm <- intersect(lo$tt$gene, cn$tt$gene)
  rows[[length(rows) + 1]] <- data.frame(
    block = "method x filter (lithium)", key = meth,
    n_genes_canon = cn$n_genes, deg_canon = n_deg(cn$tt, FDR),
    n_genes_loose = lo$n_genes, deg_loose = n_deg(lo$tt, FDR),
    ratio_deg = n_deg(lo$tt, FDR) / max(1, n_deg(cn$tt, FDR)),
    pi0_canon = pi0_storey(cn$tt$P.Value), pi0_loose = pi0_storey(lo$tt$P.Value),
    # same genes, BH re-applied inside them: isolates estimation from burden
    deg_loose_fixedU = sum(p.adjust(lo$tt$P.Value[match(cm, lo$tt$gene)], "BH") < FDR),
    deg_canon_fixedU = sum(p.adjust(cn$tt$P.Value[match(cm, cn$tt$gene)], "BH") < FDR),
    spearman_p_common = cor(lo$tt$P.Value[match(cm, lo$tt$gene)],
                            cn$tt$P.Value[match(cm, cn$tt$gene)], method = "spearman"),
    jaccard = jacc(deg_ids(lo$tt, FDR), deg_ids(cn$tt, FDR)),
    stringsAsFactors = FALSE)
}

## ---- 2. same filter grid, case-control ------------------------------------
CC <- list(loose = c(0, 0.10), mid = c(10, 0.50), canon = c(10, 0.90), strict = c(50, 0.90))
ccfits <- list()
for (cn in names(CC)) {
  a <- CC[[cn]]; tg <- sprintf("cc_voom_%s", cn)
  r <- run1(tg, "casecon", a[1], a[2], "voom")
  cat(sprintf("%-20s genes=%6d DEG=%5s %.0fs\n", tg, r$n_genes %||% NA,
              if (isTRUE(r$ok)) n_deg(r$tt, FDR) else "ERR", r$secs))
  if (isTRUE(r$ok)) ccfits[[cn]] <- r
}
base <- ccfits[["canon"]]
for (cn in names(ccfits)) {
  f <- ccfits[[cn]]; cm <- intersect(f$tt$gene, base$tt$gene)
  rows[[length(rows) + 1]] <- data.frame(
    block = "filter (case-control)", key = cn,
    n_genes_canon = base$n_genes, deg_canon = n_deg(base$tt, FDR),
    n_genes_loose = f$n_genes, deg_loose = n_deg(f$tt, FDR),
    ratio_deg = n_deg(f$tt, FDR) / max(1, n_deg(base$tt, FDR)),
    pi0_canon = pi0_storey(base$tt$P.Value), pi0_loose = pi0_storey(f$tt$P.Value),
    deg_loose_fixedU = sum(p.adjust(f$tt$P.Value[match(cm, f$tt$gene)], "BH") < FDR),
    deg_canon_fixedU = sum(p.adjust(base$tt$P.Value[match(cm, base$tt$gene)], "BH") < FDR),
    spearman_p_common = cor(f$tt$P.Value[match(cm, f$tt$gene)],
                            base$tt$P.Value[match(cm, base$tt$gene)], method = "spearman"),
    jaccard = jacc(deg_ids(f$tt, FDR), deg_ids(base$tt, FDR)),
    stringsAsFactors = FALSE)
}

X <- do.call(rbind, rows)
cat("\n== cross-checks ==\n"); print(X, row.names = FALSE, digits = 4)
write.csv(X, file.path(RUN, "xcheck.csv"), row.names = FALSE)
cat("\nwrote runs/filter_sweep/xcheck.csv\n")
