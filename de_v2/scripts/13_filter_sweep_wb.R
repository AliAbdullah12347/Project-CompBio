#!/usr/bin/env Rscript
# ==============================================================================
# 13_filter_sweep_wb.R -- how much does the whole-blood answer depend on where
# we drew the expression filter?
#
# WHY THIS IS NOT JUST "RUN IT AGAIN WITH MORE GENES"
#
# Changing the filter changes the GENE UNIVERSE, and that breaks the obvious
# summary. A percentage of genes is not comparable across filters because the
# denominator moves: 1,426 of 12,368 (11.5%) and 1,426 of 24,084 (5.9%) are the
# same finding reported two ways. Worse, the multiple-testing correction itself
# depends on the number of tests, so a looser filter makes every surviving gene
# harder to detect even if nothing about that gene changed.
#
# So the headline of a filter sweep cannot be a count or a percentage. It has to
# be STABILITY:
#
#   * which genes are called differentially expressed under every filter,
#   * which appear only at particular filters,
#   * and how far a gene's RANK moves when the universe changes.
#
# A gene that sits in the top 50 under all 31 filters is a different kind of
# finding from one that crosses 0.05 at exactly one threshold. That distinction
# is invisible in a DEG count and is the actual output of this script.
#
# WHAT IS SWEPT
#
#   count-based  >C counts in >=P of samples, for C in {0,1,5,10,20,50}
#                and P in {50%, 75%, 90%, 99%, 100%}  -- 30 combinations
#   edgeR        filterByExpr(), the field-standard automatic filter, as a
#                reference point we did not choose ourselves
#
# Gene counts range from 2,087 (>50 in 100% of samples) to 24,084 (>0 in >=50%).
# The project baseline (>10 in >=90%) gives 12,368 and is marked in the output.
#
# ORDER OF OPERATIONS. The filter runs BEFORE normalisation, every time. TMM's
# scaling factors are computed from the genes that survive, so filtering after
# normalising would leave the scaling contaminated by genes we then discard.
#
# PER-GENE RECORDS. Every gene called significant under any filter is written
# out with its rank and adjusted p-value under EVERY filter, so the identity and
# ordering of findings is recoverable later without re-running anything.
#
# COST NOTE. pi0 uses the spline estimator only, not the bootstrap. The two
# agreed to within 0.02 everywhere in 01, and the bootstrap would multiply the
# runtime of 124 fits for no change in conclusion.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
HERE <- file.path(ROOT, "de_v2")
setwd(HERE); set.seed(481)
source("scripts/mtc.R")
OUT <- file.path("results", "sweep"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

p   <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
sym <- local({ m <- read.delim(file.path(ROOT, "de_analysis/data/ens2sym.tsv"),
                               header = FALSE, col.names = c("ens", "sym"))
               setNames(m$sym, m$ens) })

# Raw counts, all genes. prep.rds carries only the 12,368 that passed the
# baseline filter, so the full matrix has to be re-read here.
cnt <- as.matrix(read.delim(file.path(ROOT, "data/cohort_474/counts_474.tsv.gz"),
                            row.names = 1, check.names = FALSE))
cnt <- cnt[!grepl("^ENSGR", rownames(cnt)), ]
# counts_474.tsv.gz carries Ensembl VERSION suffixes (ENSG00000000003.10) while
# prep.rds, the bMIND output and ens2sym.tsv all use unversioned IDs. Without
# this strip every symbol lookup returns NA and every intersect with the bMIND
# gene set is empty -- which is exactly how the first run of both sweeps failed.
rownames(cnt) <- vapply(strsplit(rownames(cnt), ".", fixed = TRUE), `[`, character(1), 1)
stopifnot(!any(duplicated(rownames(cnt))))
stopifnot(identical(colnames(cnt), p$meta$title))   # refuse to proceed on a column mismatch
N <- ncol(cnt)
cat(sprintf("== 13_filter_sweep_wb ==\nraw counts %d x %d\n", nrow(cnt), N))

## ---- the filter grid --------------------------------------------------------
MC <- c(0, 1, 5, 10, 20, 50); MP <- c(0.50, 0.75, 0.90, 0.99, 1.00)
GRID <- expand.grid(min_count = MC, min_prop = MP)
GRID$fid <- sprintf("c%d_p%02d", GRID$min_count, round(GRID$min_prop * 100))
GRID$rule <- sprintf(">%d in >=%g%%", GRID$min_count, GRID$min_prop * 100)
GRID$keep <- lapply(seq_len(nrow(GRID)), function(i)
  rownames(cnt)[rowSums(cnt > GRID$min_count[i]) >= GRID$min_prop[i] * N])
# edgeR's automatic filter, included because it is what most papers use and we
# did not choose it.
fbe <- rownames(cnt)[filterByExpr(DGEList(cnt))]
GRID <- rbind(GRID, data.frame(min_count = NA, min_prop = NA, fid = "edgeR_fbe",
                               rule = "edgeR filterByExpr", keep = I(list(fbe))))
GRID$n_genes <- vapply(GRID$keep, length, integer(1))
GRID$is_baseline <- GRID$fid == "c10_p90"
GRID <- GRID[order(GRID$n_genes), ]
cat(sprintf("%d filters, %d to %d genes | baseline c10_p90 = %d\n\n",
            nrow(GRID), min(GRID$n_genes), max(GRID$n_genes),
            GRID$n_genes[GRID$is_baseline]))

## ---- design -----------------------------------------------------------------
design_for <- function(contrast, adjust) {
  s <- p$sel[[contrast]]
  mm <- p$meta[s$keep, , drop = FALSE]; mm$grp <- s$grp
  mm$plate <- droplevels(factor(mm$plate)); mm$sex <- droplevels(factor(mm$sex))
  if ("group" %in% s$cov) mm$group <- droplevels(factor(mm$group))
  cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]])) > 1, logical(1))]
  X <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse = " + "))), data = mm)
  if (adjust == "ilr") X <- cbind(X, p$ILR[mm$title, , drop = FALSE])
  ne <- nonEstimable(X); if (!is.null(ne)) X <- X[, setdiff(colnames(X), ne), drop = FALSE]
  stopifnot(qr(X)$rank == ncol(X))
  list(X = X, samples = mm$title) }

DES <- list()
for (cn in c("LI", "BPD")) for (ad in c("raw", "ilr"))
  DES[[paste0(cn, "_", ad)]] <- design_for(cn, ad)

## ---- sweep ------------------------------------------------------------------
# CHECKPOINT. The first version of this script did all 124 fits, printed the
# summary, then died in the stability block below and wrote nothing -- 45 minutes
# lost to a bug three steps downstream of the expensive part. The raw sweep is
# now saved the moment it finishes, and reloaded if it already exists.
CKPT <- file.path(OUT, "wb_sweep_raw.rds")
SUM <- list(); GENES <- list()
t0 <- Sys.time()
if (file.exists(CKPT)) {
  z <- readRDS(CKPT); S <- z$S; G <- z$G
  cat(sprintf("loaded sweep from checkpoint: %d summary rows, %d gene rows
", nrow(S), nrow(G)))
} else {
for (i in seq_len(nrow(GRID))) {
  g <- GRID$keep[[i]]; fid <- GRID$fid[i]
  for (key in names(DES)) {
    d <- DES[[key]]
    dge <- calcNormFactors(DGEList(cnt[g, d$samples, drop = FALSE]), method = "TMM")
    fit <- eBayes(lmFit(voom(dge, d$X), d$X))
    tt  <- topTable(fit, coef = "grpcase", number = Inf, sort.by = "none")
    pv  <- tt$P.Value; m <- length(pv)
    BH  <- p.adjust(pv, "BH"); BY <- p.adjust(pv, "BY"); HOLM <- p.adjust(pv, "holm")
    pi0 <- pi0_spline(pv); QS <- qvalue_storey(pv, pi0 = pi0)
    rk  <- rank(pv, ties.method = "min")
    SUM[[length(SUM) + 1]] <- data.frame(
      fid = fid, rule = GRID$rule[i], n_genes = m, is_baseline = GRID$is_baseline[i],
      analysis = key, contrast = sub("_.*", "", key), adjust = sub(".*_", "", key),
      pi0 = round(pi0, 4), min_p = signif(min(pv), 3),
      BH_05 = sum(BH < 0.05), BH_10 = sum(BH < 0.10), BY_05 = sum(BY < 0.05),
      storey_05 = sum(QS < 0.05), holm_05 = sum(HOLM < 0.05),
      pct_BH_05 = round(100 * sum(BH < 0.05) / m, 3), stringsAsFactors = FALSE)
    # per-gene record: rank and adjusted p under THIS filter, for every gene
    GENES[[length(GENES) + 1]] <- data.frame(
      gene = rownames(tt), analysis = key, fid = fid,
      rank = rk, rank_pct = round(100 * rk / m, 4), n_genes = m,
      logFC = round(tt$logFC, 5), P.Value = tt$P.Value,
      BH = BH, storey = QS, sig_BH05 = BH < 0.05, stringsAsFactors = FALSE)
  }
  cat(sprintf("  %-10s %6d genes  [%5.1f min]\n", fid, GRID$n_genes[i],
              as.numeric(difftime(Sys.time(), t0, units = "mins")))); flush.console()
}
S <- do.call(rbind, SUM)
G <- do.call(rbind, GENES)
G$symbol <- unname(ifelse(is.na(sym[G$gene]), G$gene, sym[G$gene]))
saveRDS(list(S = S, G = G), CKPT)
write.csv(S, file.path(OUT, "wb_sweep_summary.csv"), row.names = FALSE)
cat("checkpointed raw sweep before any downstream analysis
")
}

cat("\n=== DEG counts by filter (BH 0.05) ===\n")
w <- reshape(S[, c("fid", "n_genes", "analysis", "BH_05")], idvar = c("fid", "n_genes"),
             timevar = "analysis", direction = "wide")
names(w) <- sub("BH_05\\.", "", names(w))
print(w[order(w$n_genes), ], row.names = FALSE)

## ---- which genes, and how stable -------------------------------------------
# Only genes significant under at least one filter are carried forward in the
# long table; everything else would be 30 rows of "not significant" per gene.
cat("\n=== per-gene stability ===\n")
ROB <- list()
for (key in names(DES)) {
  gk <- G[G$analysis == key, ]
  nf <- length(unique(gk$fid))
  ever <- unique(gk$gene[gk$sig_BH05])
  if (!length(ever)) { cat(sprintf("  %-10s no gene significant under any filter\n", key)); next }
  sub <- gk[gk$gene %in% ever, ]
  agg <- do.call(rbind, lapply(split(sub, sub$gene), function(z) data.frame(
    gene = z$gene[1], symbol = unname(ifelse(is.na(sym[z$gene[1]]), z$gene[1], sym[z$gene[1]])),
    analysis = key,
    n_filters_tested = nrow(z), n_filters_sig = sum(z$sig_BH05),
    frac_sig = round(sum(z$sig_BH05) / nrow(z), 3),
    best_rank = min(z$rank), median_rank = median(z$rank),
    median_rank_pct = round(median(z$rank_pct), 4),
    median_logFC = round(median(z$logFC), 4),
    min_BH = signif(min(z$BH), 3), stringsAsFactors = FALSE)))
  rownames(agg) <- NULL
  agg <- agg[order(-agg$frac_sig, agg$median_rank), ]
  ROB[[key]] <- agg
  always <- sum(agg$n_filters_sig == agg$n_filters_tested)
  cat(sprintf("  %-10s %5d genes significant somewhere | %4d in EVERY filter they were tested in (%.1f%%)\n",
              key, nrow(agg), always, 100 * always / nrow(agg)))
}
R <- do.call(rbind, ROB)

## ---- overlap between filters -----------------------------------------------
cat("\n=== Jaccard overlap of DEG sets against the baseline filter ===\n")
JA <- list()
for (key in names(DES)) {
  gk <- G[G$analysis == key, ]
  base <- gk$gene[gk$fid == "c10_p90" & gk$sig_BH05]
  for (f in unique(gk$fid)) {
    s2 <- gk$gene[gk$fid == f & gk$sig_BH05]
    u <- length(union(base, s2)); it <- length(intersect(base, s2))
    JA[[length(JA) + 1]] <- data.frame(analysis = key, fid = f,
      n_genes = gk$n_genes[gk$fid == f][1], n_deg = length(s2),
      n_shared_with_baseline = it,
      jaccard = if (u == 0) NA_real_ else round(it / u, 4), stringsAsFactors = FALSE)
  }
}
J <- do.call(rbind, JA)
for (key in c("LI_raw", "BPD_raw")) {
  z <- J[J$analysis == key, ]; z <- z[order(z$n_genes), ]
  cat(sprintf("\n  %s\n", key))
  print(z[, c("fid", "n_genes", "n_deg", "n_shared_with_baseline", "jaccard")], row.names = FALSE)
}

write.csv(S, file.path(OUT, "wb_sweep_summary.csv"), row.names = FALSE)
write.csv(J, file.path(OUT, "wb_sweep_overlap.csv"), row.names = FALSE)
write.csv(R, file.path(OUT, "wb_gene_robustness.csv"), row.names = FALSE)
write.csv(G[G$sig_BH05 | G$fid == "c10_p90", ],
          gzfile(file.path(OUT, "wb_sweep_genes.csv.gz")), row.names = FALSE)
saveRDS(GRID[, c("fid", "rule", "n_genes", "is_baseline")], file.path(OUT, "filter_grid.rds"))
cat(sprintf("\nwrote results/sweep/wb_* [%.1f min]\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
