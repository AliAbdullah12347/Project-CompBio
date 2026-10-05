#!/usr/bin/env Rscript
# ==============================================================================
# 14_filter_sweep_ct.R -- the same sweep at cell-type level, plus a second
# threshold that only exists here.
#
# A HARD CONSTRAINT, STATED UP FRONT
#
# bMIND was run once, for 3.1 hours, on the 12,368 genes that passed the
# baseline filter (>10 counts in >=90% of samples). Its output exists for those
# genes and no others. So this sweep can only go in ONE direction:
#
#   STRICTER filters (subsets of the baseline) -- available, swept here
#   LOOSER  filters (supersets)                -- would need a new bMIND run
#
# Nine of the thirty grid points are nested inside the baseline: min_count >= 10
# AND min_prop >= 0.90, giving 2,087 to 12,368 genes. The looser settings that
# 13 sweeps at whole-blood level are simply not answerable here without another
# multi-hour deconvolution, and the script says so rather than quietly reporting
# the nine as if they were the whole grid.
#
# THE SECOND THRESHOLD: HOW MUCH TO TRUST AN ESTIMATE
#
# Whole blood has one filter, on expression. Cell type has a second, on
# RELIABILITY, because its values are estimates rather than measurements. Two
# are swept, because they target different failure modes:
#
#   VARIANCE FLOOR. bMIND shrinks each sample toward a per-gene prior. The
#   effect is severe and uneven -- median across-sample variance is 0.50 in
#   granulocytes and 0.0048 in B cells, a factor of ~100. A gene whose
#   cell-type estimate barely moves between subjects cannot carry a group
#   difference; it has been shrunk flat. A variance floor removes those.
#
#   POSTERIOR SE QUANTILE. bMIND reports its own uncertainty per gene x lineage
#   x sample. Keeping only the most precisely estimated fraction of gene x
#   lineage combinations asks whether the findings come from the part of the
#   output the method itself is confident about.
#
# The two sweeps are run SEPARATELY rather than crossed. Crossing them would be
# 9 x 9 x 2 contrasts x 2 adjustments x 5 lineages = 1,620 fits for a surface
# whose interesting margins are both visible without it.
#
# BH IS RECOMPUTED WITHIN WHATEVER SURVIVES. Both filters change the number of
# tests, per lineage and in the pooled set, so every correction is recomputed on
# the surviving set rather than carried over. That is the whole point of the
# sweep.
# ==============================================================================

suppressPackageStartupMessages({library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
HERE <- file.path(ROOT, "de_v2")
setwd(HERE); set.seed(481)
source("scripts/mtc.R")
OUT <- file.path("results", "sweep"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

p  <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
bp <- readRDS(file.path(ROOT, "de_analysis/results/bmind_profiles.rds"))
A <- bp$A; SE <- bp$SE; LINS <- bp$lineages
sym <- local({ m <- read.delim(file.path(ROOT, "de_analysis/data/ens2sym.tsv"),
                               header = FALSE, col.names = c("ens", "sym"))
               setNames(m$sym, m$ens) })
BASE <- dimnames(A)[[1]]
cat(sprintf("== 14_filter_sweep_ct ==\nbMIND output: %d genes x %d lineages x %d samples\n",
            dim(A)[1], dim(A)[2], dim(A)[3]))

## ---- sweep A: nested gene filters ------------------------------------------
cnt <- as.matrix(read.delim(file.path(ROOT, "data/cohort_474/counts_474.tsv.gz"),
                            row.names = 1, check.names = FALSE))
cnt <- cnt[!grepl("^ENSGR", rownames(cnt)), ]
# counts_474.tsv.gz carries Ensembl VERSION suffixes (ENSG00000000003.10) while
# prep.rds, the bMIND output and ens2sym.tsv all use unversioned IDs. Without
# this strip every symbol lookup returns NA and every intersect with the bMIND
# gene set is empty -- which is exactly how the first run of both sweeps failed.
rownames(cnt) <- vapply(strsplit(rownames(cnt), ".", fixed = TRUE), `[`, character(1), 1)
stopifnot(!any(duplicated(rownames(cnt))))
stopifnot(identical(colnames(cnt), p$meta$title))
N <- ncol(cnt)
GF <- expand.grid(min_count = c(10, 20, 50), min_prop = c(0.90, 0.99, 1.00))
GF$fid <- sprintf("c%d_p%02d", GF$min_count, round(GF$min_prop * 100))
GF$rule <- sprintf(">%d in >=%g%%", GF$min_count, GF$min_prop * 100)
GF$genes <- lapply(seq_len(nrow(GF)), function(i)
  intersect(BASE, rownames(cnt)[rowSums(cnt > GF$min_count[i]) >= GF$min_prop[i] * N]))
GF$n <- vapply(GF$genes, length, integer(1))
GF <- GF[order(GF$n), ]
cat(sprintf("\nsweep A: %d nested gene filters, %d to %d genes\n", nrow(GF), min(GF$n), max(GF$n)))
cat("  (looser filters are unavailable: bMIND has no estimates outside the baseline 12,368)\n")

## ---- per-lineage reliability statistics ------------------------------------
# Computed once over all 474 samples, so they do not depend on which contrast
# is being fitted and cannot leak a group difference into the filter.
VAR <- sapply(LINS, function(ct) apply(A[, ct, ], 1, var))
MSE <- sapply(LINS, function(ct) apply(SE[, ct, ], 1, median))
cat("\nreliability of the estimates, by lineage:\n")
for (ct in LINS)
  cat(sprintf("  %-5s median across-sample variance %.5f | median posterior SE %.4f\n",
              ct, median(VAR[, ct]), median(MSE[, ct])))

## ---- designs ----------------------------------------------------------------
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

mkw <- function(S) { W <- 1 / (S^2); fin <- is.finite(W) & W > 0
  if (any(!fin)) W[!fin] <- min(W[fin]); W / mean(W) }

# One setting = one (gene set per lineage). Returns per-lineage and pooled
# results plus the per-gene record.
run_setting <- function(sid, label, kind, genes_by_lineage) {
  per <- list(); rec <- list(); pool <- list()
  for (key in names(DES)) {
    d <- DES[[key]]
    P <- NULL
    for (ct in LINS) {
      g <- genes_by_lineage[[ct]]
      if (length(g) < 50) {                       # too few to correct meaningfully
        per[[length(per) + 1]] <- data.frame(sid = sid, label = label, kind = kind,
          analysis = key, lineage = ct, n_genes = length(g), deg_BH05 = NA_integer_,
          pi0 = NA_real_, min_p = NA_real_, stringsAsFactors = FALSE); next }
      fit <- eBayes(lmFit(A[g, ct, d$samples], d$X, weights = mkw(SE[g, ct, d$samples])),
                    trend = TRUE)
      tt <- topTable(fit, coef = "grpcase", number = Inf, sort.by = "none")
      pv <- tt$P.Value; q <- p.adjust(pv, "BH"); rk <- rank(pv, ties.method = "min")
      per[[length(per) + 1]] <- data.frame(sid = sid, label = label, kind = kind,
        analysis = key, lineage = ct, n_genes = length(pv), deg_BH05 = sum(q < 0.05),
        pi0 = round(pi0_spline(pv), 4), min_p = signif(min(pv), 3), stringsAsFactors = FALSE)
      rec[[length(rec) + 1]] <- data.frame(gene = rownames(tt), lineage = ct,
        analysis = key, sid = sid, rank = rk, rank_pct = round(100 * rk / length(pv), 4),
        n_genes = length(pv), logFC = round(tt$logFC, 5), P.Value = pv, BH = q,
        sig_BH05 = q < 0.05, stringsAsFactors = FALSE)
      P <- rbind(P, data.frame(gene = rownames(tt), p = pv, stringsAsFactors = FALSE))
    }
    if (!is.null(P)) {
      qp <- p.adjust(P$p, "BH")
      pool[[length(pool) + 1]] <- data.frame(sid = sid, label = label, kind = kind,
        analysis = key, n_tests = nrow(P), deg_tests = sum(qp < 0.05),
        deg_unique_genes = length(unique(P$gene[qp < 0.05])),
        n_unique_genes = length(unique(P$gene)), stringsAsFactors = FALSE)
    }
  }
  list(per = do.call(rbind, per), rec = do.call(rbind, rec), pool = do.call(rbind, pool)) }

PER <- list(); REC <- list(); POOL <- list(); t0 <- Sys.time()

cat("\n--- sweep A: gene filter ---\n")
for (i in seq_len(nrow(GF))) {
  g <- GF$genes[[i]]
  r <- run_setting(GF$fid[i], GF$rule[i], "gene_filter",
                   setNames(replicate(length(LINS), g, simplify = FALSE), LINS))
  PER[[length(PER) + 1]] <- r$per; REC[[length(REC) + 1]] <- r$rec; POOL[[length(POOL) + 1]] <- r$pool
  cat(sprintf("  %-10s %6d genes  [%5.1f min]\n", GF$fid[i], GF$n[i],
              as.numeric(difftime(Sys.time(), t0, units = "mins")))); flush.console()
}

cat("\n--- sweep B: variance floor on the estimates ---\n")
for (v in c(0, 0.005, 0.01, 0.05, 0.10)) {
  gl <- setNames(lapply(LINS, function(ct) BASE[VAR[, ct] > v]), LINS)
  r <- run_setting(sprintf("var%g", v), sprintf("across-sample variance > %g", v), "variance_floor", gl)
  PER[[length(PER) + 1]] <- r$per; REC[[length(REC) + 1]] <- r$rec; POOL[[length(POOL) + 1]] <- r$pool
  cat(sprintf("  var>%-6g genes kept per lineage: %s  [%5.1f min]\n", v,
              paste(sprintf("%s=%d", LINS, vapply(gl, length, integer(1))), collapse = " "),
              as.numeric(difftime(Sys.time(), t0, units = "mins")))); flush.console()
}

cat("\n--- sweep C: posterior-SE quantile ---\n")
for (qq in c(1.00, 0.90, 0.75, 0.50)) {
  gl <- setNames(lapply(LINS, function(ct)
    BASE[MSE[, ct] <= quantile(MSE[, ct], qq, names = FALSE)]), LINS)
  r <- run_setting(sprintf("se%02d", round(qq * 100)),
                   sprintf("keep most precise %g%% by posterior SE", qq * 100), "se_quantile", gl)
  PER[[length(PER) + 1]] <- r$per; REC[[length(REC) + 1]] <- r$rec; POOL[[length(POOL) + 1]] <- r$pool
  cat(sprintf("  se<=q%-5.2f genes kept per lineage: %s  [%5.1f min]\n", qq,
              paste(sprintf("%s=%d", LINS, vapply(gl, length, integer(1))), collapse = " "),
              as.numeric(difftime(Sys.time(), t0, units = "mins")))); flush.console()
}

PER <- do.call(rbind, PER); REC <- do.call(rbind, REC); POOL <- do.call(rbind, POOL)
REC$symbol <- ifelse(is.na(sym[REC$gene]), REC$gene, sym[REC$gene])

cat("\n=== per-lineage DEG counts, lithium (unadjusted) ===\n")
z <- PER[PER$analysis == "LI_raw", c("sid", "label", "lineage", "n_genes", "deg_BH05")]
w <- reshape(z[, c("sid", "lineage", "deg_BH05")], idvar = "sid", timevar = "lineage", direction = "wide")
names(w) <- sub("deg_BH05\\.", "", names(w))
print(w, row.names = FALSE)

cat("\n=== pooled, all four analyses ===\n")
print(POOL[, c("sid", "analysis", "n_tests", "deg_tests", "deg_unique_genes")], row.names = FALSE)

## ---- per-gene stability -----------------------------------------------------
cat("\n=== per-gene stability (lithium, unadjusted) ===\n")
ROB <- list()
for (key in names(DES)) {
  gk <- REC[REC$analysis == key, ]
  ever <- unique(gk[gk$sig_BH05, c("gene", "lineage")])
  if (!nrow(ever)) { cat(sprintf("  %-10s nothing significant under any setting\n", key)); next }
  gk$gl <- paste(gk$gene, gk$lineage)
  keepgl <- paste(ever$gene, ever$lineage)
  sub <- gk[gk$gl %in% keepgl, ]
  agg <- do.call(rbind, lapply(split(sub, sub$gl), function(z) data.frame(
    gene = z$gene[1], symbol = z$symbol[1], lineage = z$lineage[1], analysis = key,
    n_settings_tested = nrow(z), n_settings_sig = sum(z$sig_BH05),
    frac_sig = round(sum(z$sig_BH05) / nrow(z), 3),
    best_rank = min(z$rank), median_rank = median(z$rank),
    median_logFC = round(median(z$logFC), 4), min_BH = signif(min(z$BH), 3),
    stringsAsFactors = FALSE)))
  agg <- agg[order(-agg$frac_sig, agg$median_rank), ]
  ROB[[key]] <- agg
  cat(sprintf("  %-10s %4d gene x lineage pairs significant somewhere | %3d in every setting tested\n",
              key, nrow(agg), sum(agg$n_settings_sig == agg$n_settings_tested)))
}
R <- do.call(rbind, ROB)

write.csv(PER,  file.path(OUT, "ct_sweep_per_lineage.csv"), row.names = FALSE)
write.csv(POOL, file.path(OUT, "ct_sweep_pooled.csv"), row.names = FALSE)
write.csv(R,    file.path(OUT, "ct_gene_robustness.csv"), row.names = FALSE)
write.csv(REC[REC$sig_BH05 | REC$sid == "c10_p90", ],
          gzfile(file.path(OUT, "ct_sweep_genes.csv.gz")), row.names = FALSE)
write.csv(data.frame(gene = BASE, symbol = ifelse(is.na(sym[BASE]), BASE, sym[BASE]),
                     round(VAR, 6), check.names = FALSE),
          gzfile(file.path(OUT, "ct_estimate_variance.csv.gz")), row.names = FALSE)
cat(sprintf("\nwrote results/sweep/ct_* [%.1f min]\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
