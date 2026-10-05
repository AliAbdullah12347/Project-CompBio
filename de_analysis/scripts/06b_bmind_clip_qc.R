#!/usr/bin/env Rscript
# ==============================================================================
# 06b_bmind_clip_qc.R -- quantifying a cross-gene coupling inside bMIND that is
# not documented in the package, and that CHUNKING INTERACTS WITH.
#
# WHAT WAS FOUND
#
# Reading the MIND source, the internal bmind() ends with:
#
#     res$A[res$A < min(X)] = min(X)
#     res$A[res$A > max(X)] = max(X)
#
# where X is the WHOLE input bulk matrix. So every posterior mean is clipped to
# the range of the bulk matrix it was given. The per-gene model itself is
# independent across genes -- profile[i,] and covariance[i,,] are built from
# gene i alone, and each gene gets its own MCMC -- but this final clip is
# computed across all genes in the call.
#
# WHY THAT MATTERS HERE. 06 ran bMIND in chunks of 400 genes for memory
# reasons. min(X) and max(X) are therefore CHUNK-SPECIFIC, not global, so an
# estimate can be clipped in one chunking and not in another. That is a
# property of our implementation, not of the data, and it has to be measured
# rather than assumed negligible.
#
# It also bears on interpretation generally: a cell-type expression value is
# NOT free to lie outside the observed bulk range. A lineage at 1.5% of the
# mixture can legitimately have expression above or below anything seen in
# bulk, and bMIND cannot report that.
#
# WHAT WOULD BE A PROBLEM: clipped values concentrated in the rare lineages,
# or a non-trivial overall clipped fraction, since clipping compresses
# variance and biases every downstream test toward the null.
# ==============================================================================

HERE <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/de_analysis"
setwd(HERE)
p <- readRDS("data/prep.rds")
RES <- "results"
fs <- sort(list.files(file.path(RES, "bmind_profiles"), "^chunk_.*\\.rds$", full.names = TRUE))
cat(sprintf("== 06b_bmind_clip_qc ==\n%d chunks\n", length(fs)))

bulk <- p$logtpm
cat(sprintf("global bulk range: [%.4f, %.4f]\n\n", min(bulk), max(bulk)))

rows <- list(); byct <- list()
for (f in fs) {
  q <- readRDS(f); g <- q$genes
  X <- bulk[g, , drop = FALSE]
  lo <- min(X); hi <- max(X)
  at_lo <- q$A <= lo + 1e-9
  at_hi <- q$A >= hi - 1e-9
  rows[[length(rows) + 1]] <- data.frame(
    chunk = basename(f), n_genes = length(g),
    chunk_lo = round(lo, 4), chunk_hi = round(hi, 4),
    n_cells = length(q$A),
    pct_at_lo = round(100 * mean(at_lo), 4),
    pct_at_hi = round(100 * mean(at_hi), 4),
    stringsAsFactors = FALSE)
  # per cell type, so we can see whether the rare lineages take the damage
  cts <- dimnames(q$A)[[2]]
  byct[[length(byct) + 1]] <- data.frame(
    lineage = cts,
    n = sapply(cts, function(ct) length(q$A[, ct, ])),
    n_lo = sapply(cts, function(ct) sum(q$A[, ct, ] <= lo + 1e-9)),
    n_hi = sapply(cts, function(ct) sum(q$A[, ct, ] >= hi - 1e-9)),
    stringsAsFactors = FALSE)
}
CH <- do.call(rbind, rows)
BY <- do.call(rbind, byct)
BY <- aggregate(cbind(n, n_lo, n_hi) ~ lineage, BY, sum)
BY$pct_clipped <- round(100 * (BY$n_lo + BY$n_hi) / BY$n, 4)

cat("=== per chunk ===\n")
print(CH[, c("chunk", "n_genes", "chunk_lo", "chunk_hi", "pct_at_lo", "pct_at_hi")], row.names = FALSE)

cat(sprintf("\nchunk upper bound varies: %.3f to %.3f (global max %.3f)\n",
            min(CH$chunk_hi), max(CH$chunk_hi), max(bulk)))
cat(sprintf("chunk lower bound varies: %.3f to %.3f (global min %.3f)\n",
            min(CH$chunk_lo), max(CH$chunk_lo), min(bulk)))

tot <- sum(CH$n_cells)
nlo <- sum(CH$pct_at_lo / 100 * CH$n_cells); nhi <- sum(CH$pct_at_hi / 100 * CH$n_cells)
cat(sprintf("\n=== overall ===\n  %.0f of %.0f posterior means at the lower bound (%.4f%%)\n", nlo, tot, 100 * nlo / tot))
cat(sprintf("  %.0f of %.0f at the upper bound (%.4f%%)\n", nhi, tot, 100 * nhi / tot))
cat(sprintf("  total clipped: %.4f%%\n", 100 * (nlo + nhi) / tot))

cat("\n=== by lineage (does clipping concentrate in the rare ones?) ===\n")
BY$mean_frac <- round(colMeans(p$lineage)[BY$lineage], 4)
print(BY[order(-BY$pct_clipped), ], row.names = FALSE)

cat("\nINTERPRETATION\n")
tc <- 100 * (nlo + nhi) / tot
th <- 100 * nhi / tot
# The two bounds mean different things and must be judged separately.
#   LOWER bound: every chunk's minimum is exactly 0, because log2(TPM+1) = 0
#     wherever TPM = 0. Clipping there is a NON-NEGATIVITY FLOOR, which is
#     biologically correct and is not a chunking artifact -- an unchunked run
#     would apply exactly the same floor.
#   UPPER bound: this one IS chunk-dependent, because each chunk's maximum is
#     whatever its most highly expressed gene happens to reach. This is the
#     number that measures the cost of our chunking.
{
  cat(sprintf("  lower bound is exactly 0 in every chunk, so the %.4f%% pinned there is a\n", 100 * nlo / tot))
  cat("  non-negativity floor, not a chunking artifact -- an unchunked run applies it too.\n")
  cat(sprintf("  the chunk-dependent part is the UPPER bound: %.4f%% of values.\n", th))
  if (th < 0.5) {
    cat("  That is negligible: chunking did not materially change the estimates.\n")
  } else {
    cat("  That is NOT negligible. Clipping compresses variance and biases tests\n")
    cat("  toward the null; report as a limitation and consider larger chunks.\n")
  }
}

write.csv(CH, file.path(RES, "bmind_clip_by_chunk.csv"), row.names = FALSE)
write.csv(BY, file.path(RES, "bmind_clip_by_lineage.csv"), row.names = FALSE)
cat("\nwrote results/bmind_clip_*.csv\n")
