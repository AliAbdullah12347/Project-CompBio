#!/usr/bin/env Rscript
# ==============================================================================
# 05_celltype_de.R -- cell-type differential expression results and the
# pre-specified decisions for CT_LI and CT_BPD.
#
# DECISIONS MADE HERE, AND WHY
#
# 1. TWO multiplicity corrections, used for different questions, exactly as
#    config.yaml specifies.
#      * PER-LINEAGE BH answers "is there a signal in granulocytes?" Each
#        lineage is its own family of ~12,368 tests.
#      * POOLED BH over all 5 x 12,368 = 61,840 tests answers "is there a
#        signal anywhere?" It is the only fair basis for comparing the
#        cell-type level against whole blood, because the cell-type side
#        otherwise gets five independent chances at significance and would
#        win on multiplicity alone.
#
# 2. The cell-type percentage is computed over UNIQUE GENES out of 12,368,
#    not over the 61,840 gene x lineage tests. Whole blood reports 4 genes out
#    of 12,368; the comparable cell-type quantity is "how many distinct genes
#    are differentially expressed in at least one lineage", on the same
#    denominator. Using 61,840 as the denominator would make the two
#    percentages measure different things.
#
# 3. pi0 is computed per lineage, because the five lineages have very
#    different abundances and therefore very different power. A single pooled
#    pi0 would average a well-powered granulocyte analysis together with a
#    B-cell analysis that is nearly blind.
# ==============================================================================

suppressPackageStartupMessages({library(limma)})
HERE <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/de_analysis"
setwd(HERE); set.seed(481)
RES <- "results"
p <- readRDS("data/prep.rds")
NGENE <- nrow(p$counts)
BAND_FEW <- 0.5; BAND_MANY <- 2.0

pi0_storey <- function(pv, l = 0.5) min(1, mean(pv > l, na.rm = TRUE) / (1 - l))
pi0_ci <- function(pv, B = 1000, l = 0.5) {
  n <- length(pv)
  quantile(vapply(seq_len(B), function(i) pi0_storey(pv[sample.int(n, n, TRUE)], l),
                  numeric(1)), c(.025, .975), names = FALSE)
}

cat(sprintf("== 05_celltype_de ==\ngenes %d | bands: few <%.1f%%, many >=%.1f%%\n", NGENE, BAND_FEW, BAND_MANY))

runs <- c("LI_raw", "LI_ilr", "BPD_raw", "BPD_ilr")
have <- runs[file.exists(file.path(RES, sprintf("bmind_%s.rds", runs)))]
cat(sprintf("available runs: %s\n", paste(have, collapse = ", ")))
if (length(have) < length(runs))
  cat(sprintf("MISSING (will be skipped): %s\n", paste(setdiff(runs, have), collapse = ", ")))

PER <- list(); POOL <- list(); PV <- list()
for (r in have) {
  o <- readRDS(file.path(RES, sprintf("bmind_%s.rds", r)))
  P <- o$pval
  PV[[r]] <- P
  # per-lineage BH
  per <- do.call(rbind, lapply(colnames(P), function(ct) {
    pv <- P[, ct]
    q <- p.adjust(pv, "BH"); ci <- pi0_ci(pv)
    data.frame(run = r, lineage = ct, n_genes = length(pv),
               deg = sum(q < 0.05, na.rm = TRUE),
               pct = round(100 * sum(q < 0.05, na.rm = TRUE) / length(pv), 3),
               pi0 = round(pi0_storey(pv), 4), pi0_lo = round(ci[1], 4),
               pi0_hi = round(ci[2], 4), pi0_ci_incl_1 = ci[2] >= 1,
               mean_frac = round(mean(p$lineage[, ct]), 4),
               stringsAsFactors = FALSE) }))
  PER[[r]] <- per
  # pooled BH across lineage x gene, then collapse to unique genes
  flat <- as.vector(P); qf <- p.adjust(flat, "BH")
  Q <- matrix(qf, nrow = nrow(P), dimnames = dimnames(P))
  genes_hit <- rownames(P)[rowSums(Q < 0.05, na.rm = TRUE) > 0]
  POOL[[r]] <- data.frame(run = r, n_tests = length(flat),
                          deg_tests = sum(qf < 0.05, na.rm = TRUE),
                          deg_unique_genes = length(genes_hit),
                          pct_of_genes = round(100 * length(genes_hit) / nrow(P), 3),
                          stringsAsFactors = FALSE)
}
PER <- do.call(rbind, PER); POOL <- do.call(rbind, POOL)

cat("\n=== per lineage (BH within lineage) ===\n")
print(PER[, c("run", "lineage", "mean_frac", "deg", "pct", "pi0", "pi0_lo", "pi0_hi", "pi0_ci_incl_1")],
      row.names = FALSE)
cat("\n=== pooled across lineages (BH over all 5 x genes tests) ===\n")
print(POOL, row.names = FALSE)

## ---- decisions --------------------------------------------------------------
wb <- read.csv(file.path(RES, "wholeblood_summary.csv"), row.names = 1)
cat("\n================== PRE-SPECIFIED DECISIONS ==================\n")
DEC <- list()

if ("LI_raw" %in% have) {
  s <- PER[PER$run == "LI_raw", ]
  pooled <- POOL[POOL$run == "LI_raw", ]
  all_few <- all(s$pct < BAND_FEW)
  all_ci1 <- all(s$pi0_ci_incl_1)
  v <- if (pooled$pct_of_genes >= BAND_MANY) "REFUTED"
       else if (all_few && all_ci1) "SUPPORTED" else "INCONCLUSIVE"
  cat(sprintf("\nCT_LI  (expect FEW -- lithium should NOT act inside cells)\n"))
  cat(sprintf("  per-lineage DEG%%: %s\n", paste(sprintf("%s=%.3f%%", s$lineage, s$pct), collapse = "  ")))
  cat(sprintf("  all lineages below %.1f%%: %s | all pi0 CIs include 1: %s\n", BAND_FEW, all_few, all_ci1))
  cat(sprintf("  pooled: %d unique genes (%.3f%% of %d)\n", pooled$deg_unique_genes, pooled$pct_of_genes, NGENE))
  cat(sprintf("  VERDICT: %s\n", v))
  DEC[["CT_LI"]] <- data.frame(hypothesis = "CT_LI", expect = "few",
                               pct = pooled$pct_of_genes, verdict = v, stringsAsFactors = FALSE)
}

if ("BPD_raw" %in% have) {
  pooled <- POOL[POOL$run == "BPD_raw", ]
  wb_pct <- wb["WB_BPD_raw", "pct_of_genes"]
  v <- if (pooled$pct_of_genes >= BAND_FEW && wb_pct < BAND_FEW) "SUPPORTED"
       else if (pooled$pct_of_genes <= wb_pct) "REFUTED"
       else "INCONCLUSIVE"
  cat(sprintf("\nCT_BPD (expect MORE than whole blood -- masking in the bulk average)\n"))
  cat(sprintf("  whole blood : %.3f%% (%d genes)\n", wb_pct, wb["WB_BPD_raw", "deg_fdr05"]))
  cat(sprintf("  cell type   : %.3f%% (%d unique genes, pooled BH)\n",
              pooled$pct_of_genes, pooled$deg_unique_genes))
  cat(sprintf("  VERDICT: %s\n", v))
  DEC[["CT_BPD"]] <- data.frame(hypothesis = "CT_BPD", expect = "more than WB",
                                pct = pooled$pct_of_genes, verdict = v, stringsAsFactors = FALSE)
}

DEC <- do.call(rbind, DEC)
write.csv(PER,  file.path(RES, "celltype_per_lineage.csv"), row.names = FALSE)
write.csv(POOL, file.path(RES, "celltype_pooled.csv"), row.names = FALSE)
if (!is.null(DEC)) write.csv(DEC, file.path(RES, "celltype_decisions.csv"), row.names = FALSE)
saveRDS(PV, file.path(RES, "celltype_pvalues.rds"))
cat("\nwrote results/celltype_*.csv\n")
