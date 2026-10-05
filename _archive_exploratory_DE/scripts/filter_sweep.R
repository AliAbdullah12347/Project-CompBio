#!/usr/bin/env Rscript
# ==============================================================================
# filter_sweep.R -- how much of the lithium DE result is a property of the
# expression filter rather than of the data?
#
# Varies exactly one thing: the gene-level expression filter.
#   grid  min_count in {0,1,5,10,20,50} x min_prop in {0.10,0.25,0.50,0.75,0.90,1.00}
#   plus  edgeR::filterByExpr (three parameterisations) and no filter at all.
# Everything else -- contrast, covariates, TMM, voom, eBayes -- is the shared
# engine in common.R, untouched.
#
# The headline question is "is the DEG count driven by real signal or by how
# many genes are tested?". A DEG count on its own cannot answer that, because
# it confounds three things that all move when the filter moves:
#   (1) the size of the testing universe  -> Benjamini-Hochberg burden
#   (2) which genes are in the universe   -> different genes can be found
#   (3) the estimates themselves          -> TMM factors and the voom
#                                            mean-variance trend are both
#                                            fitted on whatever genes survive
# So every cell is also scored on a FIXED universe (the canonical 12,173
# genes), with BH re-applied inside it. That holds (1) and (2) constant, and
# whatever still moves is (3).
#
# Writes only inside experimentation/runs/filter_sweep/.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})
source("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation/scripts/common.R")

RUN   <- file.path(EXP_HOME, "runs", "filter_sweep")
CACHE <- file.path(RUN, "cache")
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

p <- load_prep()
FDR <- 0.05

# ---------------------------------------------------------------------------
# Run one cell, caching the full topTable. Caching matters because the loose
# cells take ~40 s each and the analysis below is re-run while being developed;
# a cached fit is the same fit, so no result can drift between sections.
# ---------------------------------------------------------------------------
fit_cell <- function(tag, ..., genes = NULL) {
  f <- file.path(CACHE, paste0(tag, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  t0 <- Sys.time()
  res <- tryCatch({
    if (is.null(genes)) {
      de_fit(p, "lithium", ...)
    } else {
      # Reach into the engine's lookup environment rather than duplicating
      # de_fit's body: filterByExpr picks a gene set by its own rule, and the
      # ONLY thing that should differ from a grid cell is that gene set.
      old <- get("filter_genes", envir = globalenv())
      assign("filter_genes",
             function(counts, min_count, min_prop) intersect(genes, rownames(counts)),
             envir = globalenv())
      on.exit(assign("filter_genes", old, envir = globalenv()), add = TRUE)
      de_fit(p, "lithium", ...)
    }
  }, error = function(e) list(error = conditionMessage(e)))
  res$secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  # Only the columns needed downstream; a full fit object is large and the
  # design is identical in every cell anyway.
  out <- if (is.null(res$error))
    list(tag = tag, ok = TRUE, n_genes = res$n_genes, n_samples = res$n_samples,
         secs = res$secs,
         tt = res$tt[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")])
  else list(tag = tag, ok = FALSE, err = res$error, secs = res$secs)
  saveRDS(out, f)
  out
}

# ---------------------------------------------------------------------------
# 1. Gene counts for the whole grid, before spending any time on fits.
# ---------------------------------------------------------------------------
MC <- c(0, 1, 5, 10, 20, 50)
MP <- c(0.10, 0.25, 0.50, 0.75, 0.90, 1.00)

d0  <- build_design(p, "lithium")
cnt <- p$counts[, d0$samples, drop = FALSE]
cat(sprintf("lithium subset: %d samples, %d genes before filtering\n",
            ncol(cnt), nrow(cnt)))

grid <- expand.grid(min_count = MC, min_prop = MP)
grid$n_genes_pre <- mapply(function(a, b) length(filter_genes(cnt, a, b)),
                           grid$min_count, grid$min_prop)
cat("\n-- genes retained by each grid cell --\n")
print(matrix(grid$n_genes_pre, length(MC), length(MP),
             dimnames = list(paste0("mc", MC), paste0("mp", MP))))

# ---------------------------------------------------------------------------
# 2. Fit every cell.
# ---------------------------------------------------------------------------
tag_of <- function(mc, mp) sprintf("mc%g_mp%03.0f", mc, mp * 100)
fits <- list()
for (i in seq_len(nrow(grid))) {
  mc <- grid$min_count[i]; mp <- grid$min_prop[i]
  tg <- tag_of(mc, mp)
  if (grid$n_genes_pre[i] < 20) {                 # eBayes needs genes to borrow
    cat(sprintf("  %-12s SKIPPED (%d genes)\n", tg, grid$n_genes_pre[i]))
    fits[[tg]] <- list(tag = tg, ok = FALSE, err = "too few genes",
                       n_genes = grid$n_genes_pre[i]); next
  }
  fits[[tg]] <- fit_cell(tg, min_count = mc, min_prop = mp)
  cat(sprintf("  %-12s genes=%6d  DEG=%5s  %.0fs\n", tg,
              fits[[tg]]$n_genes,
              if (fits[[tg]]$ok) n_deg(fits[[tg]]$tt, FDR) else "NA",
              fits[[tg]]$secs))
}

# ---------------------------------------------------------------------------
# 3. filterByExpr and the no-filter extreme.
#
# filterByExpr's rule is different in kind from the grid: it works on CPM, not
# raw counts, so its threshold adapts to library size, and its "in how many
# samples" requirement comes from the smallest group size in the design rather
# than from a fixed proportion. Three parameterisations because the answer it
# gives depends on what it is told about the groups.
# ---------------------------------------------------------------------------
dge0  <- DGEList(counts = cnt)
dge0  <- calcNormFactors(dge0)
lith  <- d0$meta$lithium

fbe <- list(
  fbe_design  = filterByExpr(dge0, design = d0$X),
  fbe_group   = filterByExpr(dge0, group = lith),
  fbe_nogroup = filterByExpr(dge0))
for (nm in names(fbe)) {
  g <- rownames(cnt)[fbe[[nm]]]
  cat(sprintf("  %-12s genes=%6d\n", nm, length(g)))
  fits[[nm]] <- fit_cell(nm, genes = g)
}
# No filter at all: min_count = -1 makes `counts > -1` true everywhere, so the
# engine's own filter passes all 57,773 genes through, including ~20k that are
# zero in every sample. Kept as a reference point, not as a serious option.
# dge0 holds a second copy of the 57,773 x 226 matrix and is finished with.
# On a machine with 0.4 GB free that copy is the difference between scoring
# in seconds and paging to disk.
rm(dge0); gc(verbose = FALSE)

fits[["nofilter"]] <- fit_cell("nofilter", min_count = -1, min_prop = 0)
cat(sprintf("  %-12s genes=%6d  DEG=%d\n", "nofilter",
            fits[["nofilter"]]$n_genes, n_deg(fits[["nofilter"]]$tt, FDR)))

# ---------------------------------------------------------------------------
# 4. Reference cell.
# ---------------------------------------------------------------------------
CANON <- tag_of(10, 0.90)
ref   <- fits[[CANON]]
stopifnot(ref$ok, ref$n_genes == 12173)
refU   <- ref$tt$gene                              # the 12,173-gene universe
refDEG <- deg_ids(ref$tt, FDR)
cat(sprintf("\ncanonical %s: %d genes, %d DEG, pi0=%.3f\n",
            CANON, ref$n_genes, length(refDEG), pi0_storey(ref$tt$P.Value)))

# ---------------------------------------------------------------------------
# 5. Score every cell.
# ---------------------------------------------------------------------------
jacc <- function(a, b) if (length(union(a, b)) == 0) NA_real_ else
  length(intersect(a, b)) / length(union(a, b))
ovlp <- function(a, b) if (min(length(a), length(b)) == 0) NA_real_ else
  length(intersect(a, b)) / min(length(a), length(b))

score <- function(f, min_count = NA, min_prop = NA) {
  if (!isTRUE(f$ok))
    return(data.frame(tag = f$tag, min_count = min_count, min_prop = min_prop,
                      n_genes = f$n_genes %||% NA, ok = FALSE))
  tt   <- f$tt
  degs <- deg_ids(tt, FDR)

  # --- fixed-universe rescoring: same genes, BH re-applied inside them ---
  common <- intersect(tt$gene, refU)
  sub    <- tt[match(common, tt$gene), ]
  fdr_c  <- p.adjust(sub$P.Value, "BH")
  deg_c  <- common[fdr_c < FDR]
  # the reference scored on that identical subset, for a like-for-like number
  rsub   <- ref$tt[match(common, ref$tt$gene), ]
  degr_c <- common[p.adjust(rsub$P.Value, "BH") < FDR]

  # --- genes this cell tests that the canonical universe does not ---
  extra  <- setdiff(tt$gene, refU)
  te     <- tt[match(extra, tt$gene), ]

  data.frame(
    tag = f$tag, min_count = min_count, min_prop = min_prop, ok = TRUE,
    n_genes = f$n_genes, n_samples = f$n_samples,
    n_deg_05 = length(degs),
    n_deg_01 = n_deg(tt, 0.01),
    n_deg_10 = n_deg(tt, 0.10),
    pct_deg_05 = 100 * length(degs) / f$n_genes,
    pi0 = pi0_storey(tt$P.Value),
    n_p05 = sum(tt$P.Value < 0.05, na.rm = TRUE),
    pct_p05 = 100 * mean(tt$P.Value < 0.05, na.rm = TRUE),
    # stability against the canonical DEG set
    jaccard_vs_canon = jacc(degs, refDEG),
    overlap_vs_canon = ovlp(degs, refDEG),
    n_shared_deg = length(intersect(degs, refDEG)),
    canon_deg_tested = sum(refDEG %in% tt$gene),
    canon_deg_recovered = length(intersect(degs, refDEG)),
    frac_canon_deg_recovered = length(intersect(degs, refDEG)) / length(refDEG),
    # fixed universe
    n_common = length(common),
    n_deg_fixedU = length(deg_c),
    n_deg_fixedU_canon = length(degr_c),
    jaccard_fixedU = jacc(deg_c, degr_c),
    spearman_t_common = suppressWarnings(
      cor(sub$t, rsub$t, method = "spearman", use = "complete.obs")),
    pearson_lfc_common = suppressWarnings(
      cor(sub$logFC, rsub$logFC, use = "complete.obs")),
    # what the extra genes contribute
    n_extra = length(extra),
    n_deg_extra = if (length(extra)) sum(te$adj.P.Val < FDR, na.rm = TRUE) else 0L,
    pi0_extra = if (length(extra) > 50) pi0_storey(te$P.Value) else NA_real_,
    med_aveexpr_extra = if (length(extra)) median(te$AveExpr, na.rm = TRUE) else NA_real_,
    med_aveexpr_deg = if (length(degs)) median(tt$AveExpr[tt$gene %in% degs], na.rm = TRUE) else NA_real_,
    secs = f$secs %||% NA_real_,
    stringsAsFactors = FALSE)
}

rows <- list()
for (i in seq_len(nrow(grid)))
  rows[[length(rows) + 1]] <- score(fits[[tag_of(grid$min_count[i], grid$min_prop[i])]],
                                    grid$min_count[i], grid$min_prop[i])
for (nm in c(names(fbe), "nofilter"))
  rows[[length(rows) + 1]] <- score(fits[[nm]])

S <- do.call(rbind, lapply(rows, function(r) {
  miss <- setdiff(names(rows[[which.max(sapply(rows, ncol))]]), names(r))
  for (m in miss) r[[m]] <- NA
  r[, names(rows[[which.max(sapply(rows, ncol))]])]
}))
S$is_canon <- S$tag == CANON

# ---------------------------------------------------------------------------
# 6. TMM normalisation factors -- do they move when the gene set moves?
#    This is mechanism (3) isolated from the model fit entirely.
# ---------------------------------------------------------------------------
# Only a spread of cells, not all 40. Each call copies a gene-subset of a
# 57,773 x 226 matrix, and this machine had ~0.4 GB of 7.7 GB free while the
# sweep ran (fourteen R processes). Forty copies thrash; eight do not, and
# eight spanning the full range answer the question just as well -- the
# fixed-universe columns already isolate re-estimation for every cell.
nf_of <- function(genes) {
  d <- DGEList(counts = cnt[genes, , drop = FALSE])
  nf <- calcNormFactors(d)$samples$norm.factors
  rm(d); gc(verbose = FALSE)
  nf
}
nf_ref <- nf_of(refU)
S$tmm_cor_vs_canon <- NA_real_; S$tmm_maxabsdiff_vs_canon <- NA_real_
tmm_cells <- unique(c(CANON, tag_of(0, 0.10), tag_of(0, 0.50), tag_of(5, 0.25),
                      tag_of(20, 0.50), tag_of(50, 0.90), tag_of(50, 1.00),
                      "fbe_design", "fbe_nogroup"))
for (i in which(S$tag %in% tmm_cells)) {
  f <- fits[[S$tag[i]]]
  if (!isTRUE(f$ok)) next
  nf <- nf_of(f$tt$gene)
  S$tmm_cor_vs_canon[i] <- cor(nf, nf_ref)
  S$tmm_maxabsdiff_vs_canon[i] <- max(abs(nf - nf_ref))
  cat(sprintf("  TMM %-12s cor=%.6f maxdiff=%.5f
", S$tag[i],
              S$tmm_cor_vs_canon[i], S$tmm_maxabsdiff_vs_canon[i]))
}

# ---------------------------------------------------------------------------
# 7. Is DEG count a function of universe size?
# ---------------------------------------------------------------------------
G <- S[S$ok & !is.na(S$min_count), ]
fitlm <- lm(log(n_deg_05 + 1) ~ log(n_genes), data = G[G$n_deg_05 > 0, ])
cat("\n== DEG count vs universe size (grid cells only) ==\n")
cat(sprintf("Spearman(n_genes, n_deg_05)   = %+.3f\n",
            cor(G$n_genes, G$n_deg_05, method = "spearman")))
cat(sprintf("Spearman(n_genes, pct_deg_05) = %+.3f\n",
            cor(G$n_genes, G$pct_deg_05, method = "spearman")))
cat(sprintf("Spearman(n_genes, pi0)        = %+.3f\n",
            cor(G$n_genes, G$pi0, method = "spearman")))
cat(sprintf("slope of log(DEG) on log(genes) = %+.3f  (1.0 = pure bookkeeping)\n",
            coef(fitlm)[2]))
cat(sprintf("fixed-universe DEG across cells: min %d  median %d  max %d  (canonical reference on the same subsets: %d..%d)\n",
            min(G$n_deg_fixedU, na.rm = TRUE), as.integer(median(G$n_deg_fixedU, na.rm = TRUE)),
            max(G$n_deg_fixedU, na.rm = TRUE),
            min(G$n_deg_fixedU_canon, na.rm = TRUE), max(G$n_deg_fixedU_canon, na.rm = TRUE)))

# ---------------------------------------------------------------------------
# 8. Output.
# ---------------------------------------------------------------------------
write.csv(S, file.path(RUN, "summary.csv"), row.names = FALSE)

# wide matrices, the way a reader wants to see a sweep
wide <- function(col) {
  m <- matrix(NA_real_, length(MC), length(MP),
              dimnames = list(paste0("min_count=", MC), paste0("min_prop=", MP)))
  for (i in seq_len(nrow(G))) m[match(G$min_count[i], MC), match(G$min_prop[i], MP)] <- G[[col]][i]
  m
}
for (col in c("n_genes", "n_deg_05", "pct_deg_05", "pi0", "jaccard_vs_canon",
              "n_deg_fixedU", "spearman_t_common")) {
  cat(sprintf("\n-- %s --\n", col)); print(round(wide(col), 4))
}
saveRDS(list(S = S, wide_cols = c("n_genes", "n_deg_05", "pct_deg_05", "pi0",
                                  "jaccard_vs_canon", "n_deg_fixedU"),
             MC = MC, MP = MP, CANON = CANON, refDEG = refDEG, refU = refU),
        file.path(RUN, "sweep.rds"))
cat(sprintf("\nwrote runs/filter_sweep/summary.csv (%d rows) and sweep.rds\n", nrow(S)))
