#!/usr/bin/env Rscript
# ===========================================================================
# netmet-061  Aggregate netmet-060 and score every contrast against its null
#
# netmet-060 writes one .rds per draw holding four networks: ctrlA, ctrlB,
# nolith, lith, all at n = 74. This row turns those into answers.
#
# THE DECISION RULE, which is the whole point
#
# Each draw supplies three contrasts computed by identical code:
#
#     NULL     ctrlA  vs ctrlB     two halves of one control group
#     BIPOLAR  ctrlA  vs nolith
#     LITHIUM  nolith vs lith
#
# NULL cannot contain biology. So for every metric the reported quantity is
# not "does BIPOLAR differ" but "does BIPOLAR differ MORE THAN NULL does, on
# the same draw". The test is a paired Wilcoxon signed-rank across draws.
#
# This is what cent-040 lacked. It found 52 of 71 modules differing in degree
# between control and lithium; cent-056 found 43 from splitting the control
# group in half. Reporting the first without the second turns a measurement
# artifact into a finding. Here both numbers come out of the same loop and
# sit in the same table.
#
# WHAT IT WRITES
#   module_counts_by_contrast.csv  the cent-040 statistic, with its null beside it
#   paired_tests.csv               BIPOLAR vs NULL and LITHIUM vs NULL, per metric
#   global_similarity.csv          whole-graph agreement between groups per draw
#   gs_kme.csv                     gene significance vs module membership
#   top_genes_by_metric.csv        top 10 genes per metric per group
#   top_shifts.csv                 genes moving most, with the null shift beside
#   results_summary.csv            one row per metric: the verdict
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "netmet-061"
stopifnot(id == "netmet-061")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))
msg("=== ", id, " starting ===")

# Mirrors netmet-060's test hook: when that row was run subsetted it wrote to
# runs/netmet-060_test, and this row must read the same place or it would
# silently aggregate nothing.
TEST <- suppressWarnings(as.integer(Sys.getenv("NETMET_SUBSET", "0")))
if (is.na(TEST)) TEST <- 0L
SRC     <- if (TEST > 0L) "network/runs/netmet-060_test" else "network/runs/netmet-060"
run_dir <- file.path("network/runs", if (TEST > 0L) paste0(id, "_test") else id)
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)

ALPHA   <- 0.05
TOP_N   <- 10L
# Standardised |delta| above which a gene counts as "moved". Pre-specified so
# the threshold cannot be chosen after seeing which value favours a result.
Z_MOVE  <- 2.0

f <- sort(list.files(SRC, pattern = "^draw_\\d+\\.rds$", full.names = TRUE))
if (!length(f)) { msg("DEPENDENCY: no draws in ", SRC); quit(save = "no", status = 1L) }
msg(length(f), " draw file(s): ", paste(basename(f), collapse = ", "))

suppressPackageStartupMessages({ library(igraph) })

CONTRASTS <- list(NULL_ = c("ctrlA", "ctrlB"),
                  BIPOLAR = c("ctrlA", "nolith"),
                  LITHIUM = c("nolith", "lith"))
GENE_METRICS <- c("degree", "eigenvector", "betweenness", "mcc_log",
                  "kIN", "kME_own", "ppr_BPD", "ppr_LI")

# A metric that is -Inf (log-MCC for a gene in no clique) or NA must not turn
# a comparison into a silent NA. Both sides are filtered to finite pairs.
fin2 <- function(a, b) is.finite(a) & is.finite(b)

## ===========================================================================
## per draw
## ===========================================================================
MODCOUNT <- list(); GLOBAL <- list(); GSKME <- list(); MECOR <- list()
# Per-gene shifts are kept as genes x draws matrices, one per
# (contrast, metric), NOT as stacked data frames. The stacked form is
# 12368 genes x 8 metrics x 3 contrasts x 10 draws = 3.0M rows, and the
# aggregate()/reshape() pair that would collapse it is superlinear in the
# number of id combinations -- 4800 of them already cost seconds in a
# 600-gene test run, and the real row has 98944. rowMeans on a matrix is
# the same arithmetic in one pass.
ZMAT <- list(); GENES <- NULL; GMOD <- NULL; DRAWS <- integer(0)

for (fp in f) {
  D  <- readRDS(fp)
  k  <- D$draw
  M  <- D$metrics
  gs <- D$gs
  msg("--- draw ", k, " ---")
  mods <- names(M$ctrlA$tom_mod)

  for (cn in names(CONTRASTS)) {
    a <- CONTRASTS[[cn]][1]; b <- CONTRASTS[[cn]][2]
    GA <- M[[a]]$gene; GB <- M[[b]]$gene
    stopifnot(identical(GA$gene, GB$gene))

    ## -- the cent-040 statistic: modules differing, per metric ----------------
    for (mt in GENE_METRICS) {
      pv <- vapply(mods, function(m) {
        ix <- GA$module == m
        x  <- GA[[mt]][ix]; y <- GB[[mt]][ix]
        ok <- fin2(x, y)
        if (sum(ok) < 5L) return(NA_real_)
        tryCatch(wilcox.test(x[ok], y[ok], exact = FALSE)$p.value,
                 error = function(e) NA_real_)
      }, numeric(1))
      padj <- p.adjust(pv, "BH")
      MODCOUNT[[length(MODCOUNT) + 1L]] <- data.frame(
        draw = k, contrast = cn, metric = mt,
        n_modules = sum(!is.na(padj)),
        n_differ  = sum(padj < ALPHA, na.rm = TRUE))
    }

    ## -- 4. topological overlap: per-module and global ------------------------
    tp <- vapply(mods, function(m) M[[a]]$tom_mod[[m]] - M[[b]]$tom_mod[[m]], numeric(1))
    tom_r <- suppressWarnings(cor(M[[a]]$tom_vec, M[[b]]$tom_vec,
                                  method = "spearman", use = "complete.obs"))

    ## -- 1. module eigengene correlation --------------------------------------
    RA <- M[[a]]$me_cor; RB <- M[[b]]$me_cor
    keep <- intersect(rownames(RA), rownames(RB))
    me_r <- suppressWarnings(cor(RA[keep, keep][upper.tri(diag(length(keep)))],
                                 RB[keep, keep][upper.tri(diag(length(keep)))],
                                 method = "spearman"))

    ## -- global agreement on every gene-level metric --------------------------
    gl <- vapply(GENE_METRICS, function(mt) {
      ok <- fin2(GA[[mt]], GB[[mt]])
      if (sum(ok) < 50L) return(NA_real_)
      suppressWarnings(cor(GA[[mt]][ok], GB[[mt]][ok], method = "spearman"))
    }, numeric(1))

    GLOBAL[[length(GLOBAL) + 1L]] <- data.frame(
      draw = k, contrast = cn,
      me_mean_abs_a = M[[a]]$me_mean_abs, me_mean_abs_b = M[[b]]$me_mean_abs,
      d_me_mean_abs = M[[a]]$me_mean_abs - M[[b]]$me_mean_abs,
      me_cor_spearman = me_r,
      tom_global_spearman = tom_r,
      tom_mod_mean_abs_delta = mean(abs(tp)),
      as.list(setNames(gl, paste0("rho_", GENE_METRICS))))

    ## -- 5. differential connectivity, gene level -----------------------------
    # Standardised so metrics on different scales are comparable, and so the
    # NULL contrast's own spread sets what "moved" means. MAD rather than sd:
    # a handful of extreme genes would otherwise inflate the scale and hide
    # the very shifts this is looking for.
    if (is.null(GENES)) { GENES <<- GA$gene; GMOD <<- GA$module }
    stopifnot(identical(GA$gene, GENES))
    for (mt in GENE_METRICS) {
      dlt <- GA[[mt]] - GB[[mt]]
      ok  <- is.finite(dlt)
      if (sum(ok) < 50L) next
      sc <- stats::mad(dlt[ok]); if (!is.finite(sc) || sc == 0) sc <- sd(dlt[ok])
      if (!is.finite(sc) || sc == 0) next
      z <- rep(NA_real_, length(dlt))
      z[ok] <- (dlt[ok] - median(dlt[ok])) / sc
      key <- paste(cn, mt, sep = "|")
      if (is.null(ZMAT[[key]])) ZMAT[[key]] <- list()
      ZMAT[[key]][[as.character(k)]] <- z
    }
  }

  ## -- 3. gene significance vs module membership ------------------------------
  # The standard WGCNA question: inside a module, do the genes most central to
  # it also track the trait? Computed per contrast on that contrast's own GS.
  for (cn in c("NULL_", "BIPOLAR", "LITHIUM")) {
    gcol <- c(NULL_ = "gs_null", BIPOLAR = "gs_bpd", LITHIUM = "gs_lith")[[cn]]
    ref  <- M[[CONTRASTS[[cn]][1]]]$gene
    GSKME[[length(GSKME) + 1L]] <- do.call(rbind, lapply(mods, function(m) {
      ix <- ref$module == m
      x  <- ref$kME_own[ix]; y <- gs[[gcol]][ix]
      ok <- fin2(x, y)
      if (sum(ok) < 5L) return(NULL)
      ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
      data.frame(draw = k, contrast = cn, module = m, n = sum(ok),
                 mean_gs = mean(y[ok]), rho_gs_kme = unname(ct$estimate),
                 p = ct$p.value)
    }))
  }
  MECOR[[length(MECOR) + 1L]] <- data.frame(
    draw = k, group = names(M),
    me_mean_abs = vapply(M, function(x) x$me_mean_abs, numeric(1)),
    n_mcc_capped = vapply(M, function(x) x$n_mcc_capped, numeric(1)))
  rm(D, M); gc(verbose = FALSE)
}

MODCOUNT <- do.call(rbind, MODCOUNT)
GLOBAL   <- do.call(rbind, GLOBAL)
GSKME    <- do.call(rbind, GSKME)
MECOR    <- do.call(rbind, MECOR)
write.csv(MODCOUNT, file.path(run_dir, "module_counts_by_contrast.csv"), row.names = FALSE)
write.csv(GLOBAL,   file.path(run_dir, "global_similarity.csv"),         row.names = FALSE)
write.csv(GSKME,    file.path(run_dir, "gs_kme.csv"),                    row.names = FALSE)
write.csv(MECOR,    file.path(run_dir, "me_mean_abs_by_group.csv"),      row.names = FALSE)

## ===========================================================================
## paired tests: does a contrast exceed its own null?
## ===========================================================================
msg("paired tests over ", length(unique(MODCOUNT$draw)), " draws ...")

paired <- function(df, value, label) {
  w <- reshape(df[, c("draw", "contrast", value)], idvar = "draw",
               timevar = "contrast", direction = "wide")
  nm <- function(c_) paste0(value, ".", c_)
  out <- lapply(c("BIPOLAR", "LITHIUM"), function(c_) {
    if (!all(c(nm(c_), nm("NULL_")) %in% colnames(w))) return(NULL)
    x <- w[[nm(c_)]]; y <- w[[nm("NULL_")]]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3L) return(NULL)
    p <- tryCatch(wilcox.test(x[ok], y[ok], paired = TRUE, exact = FALSE)$p.value,
                  error = function(e) NA_real_)
    data.frame(metric = label, contrast = c_, n_draws = sum(ok),
               mean_contrast = mean(x[ok]), mean_null = mean(y[ok]),
               mean_excess = mean(x[ok] - y[ok]),
               n_draws_above_null = sum(x[ok] > y[ok]), p_paired = p)
  })
  do.call(rbind, out)
}

PT <- do.call(rbind, c(
  lapply(split(MODCOUNT, MODCOUNT$metric), function(s)
    paired(s, "n_differ", paste0("modules_differing:", s$metric[1]))),
  lapply(c("me_cor_spearman", "tom_global_spearman", "tom_mod_mean_abs_delta",
           "d_me_mean_abs", paste0("rho_", GENE_METRICS)), function(v) {
    s <- GLOBAL[, c("draw", "contrast", v)]
    names(s)[3] <- "n_differ"
    paired(s, "n_differ", v)
  })))
PT$padj <- p.adjust(PT$p_paired, "BH")
PT <- PT[order(PT$p_paired), ]
write.csv(PT, file.path(run_dir, "paired_tests.csv"), row.names = FALSE)

## ===========================================================================
## genes: the top of each metric, and the biggest shifts against the null
## ===========================================================================
msg("top genes per metric ...")
sym <- tryCatch({
  s <- read.csv("network/runs/ann-044/reference_gene_modules.csv")
  setNames(s$symbol, s$gene)
}, error = function(e) character(0))
lab <- function(g) { v <- unname(sym[g]); ifelse(is.na(v) | v == "", g, v) }

D1 <- readRDS(f[1])
TOP <- do.call(rbind, lapply(names(D1$metrics), function(grp) {
  G <- D1$metrics[[grp]]$gene
  do.call(rbind, lapply(GENE_METRICS, function(mt) {
    v  <- G[[mt]]
    ok <- is.finite(v)
    if (!any(ok)) return(NULL)
    o  <- order(v[ok], decreasing = TRUE)[seq_len(min(TOP_N, sum(ok)))]
    data.frame(group = grp, metric = mt, rank = seq_along(o),
               gene = G$gene[ok][o], symbol = lab(G$gene[ok][o]),
               module = G$module[ok][o], value = v[ok][o])
  }))
}))
write.csv(TOP, file.path(run_dir, "top_genes_by_metric.csv"), row.names = FALSE)

# A gene "moves" only if it moves further than the null contrast moves it.
# Averaged over draws first, so one draw cannot nominate a gene on its own.
msg("per-gene shifts against the null ...")
zmean <- function(cn, mt) {
  key <- paste(cn, mt, sep = "|")
  L <- ZMAT[[key]]
  if (is.null(L) || !length(L)) return(rep(NA_real_, length(GENES)))
  rowMeans(matrix(unlist(L, use.names = FALSE), nrow = length(GENES)), na.rm = TRUE)
}
METS <- unique(sub("^[^|]*[|]", "", names(ZMAT)))
WD <- do.call(rbind, lapply(METS, function(mt) {
  d <- data.frame(gene = GENES, module = GMOD, metric = mt,
                  z_NULL = zmean("NULL_", mt), z_BIPOLAR = zmean("BIPOLAR", mt),
                  z_LITHIUM = zmean("LITHIUM", mt))
  d$excess_BIPOLAR <- abs(d$z_BIPOLAR) - abs(d$z_NULL)
  d$excess_LITHIUM <- abs(d$z_LITHIUM) - abs(d$z_NULL)
  d
}))
WD$symbol <- lab(WD$gene)
write.csv(WD, file.path(run_dir, "gene_shifts.csv"), row.names = FALSE)

nmove <- do.call(rbind, lapply(METS, function(mt) {
  x <- WD[WD$metric == mt, ]
  data.frame(metric = mt,
             n_moved_NULL    = sum(abs(x$z_NULL)    > Z_MOVE, na.rm = TRUE),
             n_moved_BIPOLAR = sum(abs(x$z_BIPOLAR) > Z_MOVE, na.rm = TRUE),
             n_moved_LITHIUM = sum(abs(x$z_LITHIUM) > Z_MOVE, na.rm = TRUE))
}))
write.csv(nmove, file.path(run_dir, "genes_moved_vs_null.csv"), row.names = FALSE)

top_shift <- do.call(rbind, lapply(METS, function(mt) {
  do.call(rbind, lapply(c("BIPOLAR", "LITHIUM"), function(c_) {
    ec <- paste0("excess_", c_)
    x  <- WD[WD$metric == mt & is.finite(WD[[ec]]), ]
    if (!nrow(x)) return(NULL)
    x <- x[order(-x[[ec]]), ][seq_len(min(TOP_N, nrow(x))), ]
    data.frame(metric = mt, contrast = c_, rank = seq_len(nrow(x)),
               gene = x$gene, symbol = x$symbol, module = x$module,
               z_contrast = x[[paste0("z_", c_)]], z_null = x$z_NULL,
               excess = x[[ec]])
  }))
}))
write.csv(top_shift, file.path(run_dir, "top_shifts.csv"), row.names = FALSE)

## ===========================================================================
## verdict
## ===========================================================================
SUMM <- PT[, c("metric", "contrast", "n_draws", "mean_contrast", "mean_null",
               "mean_excess", "n_draws_above_null", "p_paired", "padj")]
SUMM$verdict <- ifelse(is.na(SUMM$padj), "untestable",
                ifelse(SUMM$padj < ALPHA, "EXCEEDS NULL", "within null"))
write.csv(SUMM, file.path(run_dir, "results_summary.csv"), row.names = FALSE)

msg("")
msg("=== modules differing, contrast vs its own null (mean over draws) ===")
mc_ <- MODCOUNT
agg2 <- aggregate(n_differ ~ metric + contrast, data = mc_, FUN = mean)
wq <- reshape(agg2, idvar = "metric", timevar = "contrast", direction = "wide")
nmod <- max(mc_$n_modules, na.rm = TRUE)
for (i in seq_len(nrow(wq)))
  msg(sprintf("  %-12s  NULL %5.1f | BIPOLAR %5.1f | LITHIUM %5.1f   (of %d modules)",
      wq$metric[i], wq$`n_differ.NULL_`[i], wq$`n_differ.BIPOLAR`[i],
      wq$`n_differ.LITHIUM`[i], nmod))

msg("")
msg("=== metrics whose contrast exceeds its null at BH ", ALPHA, " ===")
hit <- SUMM[!is.na(SUMM$padj) & SUMM$padj < ALPHA, ]
if (!nrow(hit)) {
  msg("  none of ", nrow(SUMM), " metric x contrast tests")
} else {
  for (i in seq_len(nrow(hit)))
    msg(sprintf("  %-34s %-8s excess %+.3f  p=%.3g  BH=%.3g",
        hit$metric[i], hit$contrast[i], hit$mean_excess[i],
        hit$p_paired[i], hit$padj[i]))
}

msg("")
msg("=== genes moved beyond |z| > ", Z_MOVE, " (mean over draws) ===")
for (i in seq_len(nrow(nmove)))
  msg(sprintf("  %-12s  NULL %5d | BIPOLAR %5d | LITHIUM %5d",
      nmove$metric[i], nmove$n_moved_NULL[i],
      nmove$n_moved_BIPOLAR[i], nmove$n_moved_LITHIUM[i]))

msg("")
msg("=== top 10 betweenness, control (draw ", D1$draw, ") ===")
tb <- TOP[TOP$group == "ctrlA" & TOP$metric == "betweenness", ]
for (i in seq_len(nrow(tb)))
  msg(sprintf("  %2d  %-12s %-16s %.5f", tb$rank[i], tb$symbol[i], tb$module[i], tb$value[i]))

msg("=== ", id, " COMPLETE ===")
