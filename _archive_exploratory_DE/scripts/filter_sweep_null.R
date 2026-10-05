#!/usr/bin/env Rscript
# ==============================================================================
# filter_sweep_null.R -- is the filter's effect on the DEG count a change in
# real signal, or a change in false-positive behaviour?
#
# The sweep shows the DEG count moving by a factor of several across the grid.
# That is only interesting if the tests stay calibrated. Under a permuted
# lithium label there is no lithium effect by construction, so any cell that
# still returns DEGs is returning false ones, and any cell whose raw p-values
# are not uniform is mis-specified rather than merely underpowered.
#
# Only the exposure column is permuted (de_fit's `permute_exposure`), so every
# covariate stays attached to its own sample and the null being tested is
# "lithium carries no signal", not "these 226 people are exchangeable".
#
# Four cells: the loosest grid corner, the canonical cell, a strict cell, and
# the no-filter extreme. NPERM permutations each, one shared seed.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})
source("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation/scripts/common.R")

RUN <- file.path(EXP_HOME, "runs", "filter_sweep")
dir.create(RUN, showWarnings = FALSE, recursive = TRUE)
p <- load_prep()

NPERM <- 8
SEED  <- 20261218          # fixed once, here, and never changed per cell
FDR   <- 0.05

# Three cells, not four. The no-filter cell was dropped: it is a reference
# point rather than a candidate analysis (roughly 20,000 of its genes are zero
# in every sample, so voom's trend is fitted partly on structural zeros), and
# nine permuted fits on a 57,773-gene matrix would have to run on a machine
# with ~0.4 GB of 7.7 GB free. Calibrating an analysis nobody would run is not
# worth making the other three unreliable.
CELLS <- list(
  loose  = list(min_count = 0,  min_prop = 0.10),
  canon  = list(min_count = 10, min_prop = 0.90),
  strict = list(min_count = 50, min_prop = 0.90))

d0  <- build_design(p, "lithium")
n   <- d0$n
# The canonical universe, so every cell can be scored separately on the genes
# the canonical filter keeps and on the genes it throws away. That split is the
# whole point: the sweep suggests the discarded genes carry diffuse signal
# (pi0 well below 1), and the only way to tell diffuse signal from a voom
# normal-approximation failure at low counts is to look at them under the null.
REFU <- filter_genes(p$counts[, d0$samples, drop = FALSE], 10, 0.90)

set.seed(SEED)
# The same NPERM permutations are used in every cell. Reusing them means a
# difference between cells is the filter and not the draw.
PERMS <- replicate(NPERM, sample.int(n), simplify = FALSE)


# One scored row per fit. Split by whether the gene is in the canonical
# universe, because a cell's overall pi0 mixes two populations with very
# different count levels and therefore very different reliability.
row <- function(f, nm, a, k) {
  tt  <- f$tt
  lo  <- !(tt$gene %in% REFU)          # genes the canonical filter discards
  hi  <- !lo
  pv  <- function(i, q) if (sum(i) > 50) 100 * mean(tt$P.Value[i] < q, na.rm = TRUE) else NA_real_
  data.frame(
    cell = nm, min_count = a$min_count, min_prop = a$min_prop, perm = k,
    n_genes = f$n_genes, n_deg = n_deg(tt, FDR),
    pi0 = pi0_storey(tt$P.Value),
    pct_p05 = 100 * mean(tt$P.Value < 0.05, na.rm = TRUE),
    pct_p01 = 100 * mean(tt$P.Value < 0.01, na.rm = TRUE),
    min_fdr = min(tt$adj.P.Val, na.rm = TRUE),
    n_hi = sum(hi), n_lo = sum(lo),
    pi0_hi = if (sum(hi) > 50) pi0_storey(tt$P.Value[hi]) else NA_real_,
    pi0_lo = if (sum(lo) > 50) pi0_storey(tt$P.Value[lo]) else NA_real_,
    pct_p05_hi = pv(hi, 0.05), pct_p05_lo = pv(lo, 0.05),
    pct_p01_hi = pv(hi, 0.01), pct_p01_lo = pv(lo, 0.01),
    n_deg_hi = sum(tt$adj.P.Val[hi] < FDR, na.rm = TRUE),
    n_deg_lo = sum(tt$adj.P.Val[lo] < FDR, na.rm = TRUE),
    stringsAsFactors = FALSE)
}

out <- list()
for (nm in names(CELLS)) {
  a <- CELLS[[nm]]
  obs <- de_fit(p, "lithium", min_count = a$min_count, min_prop = a$min_prop)
  cat(sprintf("\n[%s] mc=%g mp=%.2f  genes=%d  observed DEG=%d  pi0=%.3f  p<.05=%.1f%%\n",
              nm, a$min_count, a$min_prop, obs$n_genes, n_deg(obs$tt, FDR),
              pi0_storey(obs$tt$P.Value), 100 * mean(obs$tt$P.Value < 0.05)))
  out[[length(out) + 1]] <- row(obs, nm, a, 0L)

  for (k in seq_len(NPERM)) {
    f <- de_fit(p, "lithium", min_count = a$min_count, min_prop = a$min_prop,
                permute_exposure = PERMS[[k]])
    out[[length(out) + 1]] <- row(f, nm, a, k)
    r <- out[[length(out)]]
    rm(f); gc(verbose = FALSE)
    cat(sprintf("   perm %2d: DEG=%5d  pi0=%.3f  p<.05=%.2f%%  |  hi: pi0=%.3f p<.05=%.2f%%  lo: pi0=%s p<.05=%s\n",
                k, r$n_deg, r$pi0, r$pct_p05, r$pi0_hi, r$pct_p05_hi,
                formatC(r$pi0_lo, format = "f", digits = 3),
                formatC(r$pct_p05_lo, format = "f", digits = 2)))
    # written every iteration so a kill under CPU contention loses one fit, not all
    write.csv(do.call(rbind, out), file.path(RUN, "null_partial.csv"), row.names = FALSE)
  }
}

N <- do.call(rbind, out)
write.csv(N, file.path(RUN, "null.csv"), row.names = FALSE)

cat("\n== permutation summary ==\n")
for (nm in names(CELLS)) {
  o <- N[N$cell == nm & N$perm == 0, ]; q <- N[N$cell == nm & N$perm > 0, ]
  cat(sprintf("%-7s genes=%6d | OBS DEG=%5d pi0=%.3f p<.05=%5.2f%% | NULL DEG med=%.0f max=%.0f  pi0 med=%.3f  p<.05 med=%5.2f%%  empirical p=%.3f\n",
              nm, o$n_genes, o$n_deg, o$pi0, o$pct_p05,
              median(q$n_deg), max(q$n_deg), median(q$pi0), median(q$pct_p05),
              (1 + sum(q$n_deg >= o$n_deg)) / (1 + nrow(q))))
  # Calibration by stratum. Under the null a well-specified test gives pi0 = 1
  # and 5.00% of p-values below 0.05 IN EVERY STRATUM. A stratum that departs
  # from that is mis-specified, and any "signal" the observed fit finds there
  # is not evidence.
  if (!is.na(o$pi0_lo))
    cat(sprintf("        canonical-universe genes (n=%5d): null pi0=%.3f p<.05=%5.2f%% p<.01=%5.2f%% | discarded genes (n=%5d): null pi0=%.3f p<.05=%5.2f%% p<.01=%5.2f%%\n",
                o$n_hi, median(q$pi0_hi), median(q$pct_p05_hi), median(q$pct_p01_hi),
                o$n_lo, median(q$pi0_lo), median(q$pct_p05_lo), median(q$pct_p01_lo)))
  if (!is.na(o$pi0_lo))
    cat(sprintf("        observed for comparison        : hi pi0=%.3f p<.05=%5.2f%%  DEG=%d | lo pi0=%.3f p<.05=%5.2f%%  DEG=%d\n",
                o$pi0_hi, o$pct_p05_hi, o$n_deg_hi, o$pi0_lo, o$pct_p05_lo, o$n_deg_lo))
}
cat(sprintf("\nwrote runs/filter_sweep/null.csv (%d rows, seed %d)\n", nrow(N), SEED))
