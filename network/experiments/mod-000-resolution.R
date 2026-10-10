#!/usr/bin/env Rscript
# ===========================================================================
# mod-000-resolution  Module resolution diagnostic
#
# Spec: base-001 found 18 modules (cutreeDynamic deepSplit=2) that
#       mergeCloseModules at 0.25 collapsed to 4. This diagnostic identifies
#       which (mergeCutHeight x deepSplit) cell yields a usable decomposition.
#
#   PRIMARY AXIS : mergeCutHeight in {0.00, 0.05, 0.10, 0.15, 0.25}
#   SECONDARY    : deepSplit in {2, 4}
#   Power        : 12 (TOM not cached — module colors only were saved)
#   Fallback     : power 16, only if no power-12 cell meets the target
#   Target       : >=10 non-grey modules AND largest module <25% of genes
#
# Design:
#   1. Build adjacency + TOM + dendrogram at power 12 (~47 min, n=234)
#      Checkpoint dissTOM as runs/mod-000-resolution/dissTOM_p12.rds
#   2. For each deepSplit: cutreeDynamic once, checkpoint labels
#   3. For each (deepSplit x mergeCutHeight): mergeCloseModules, record metrics
#   4. If no cell meets target: power-16 fallback (same steps, ~47 min more)
#
# Non-negotiable design points obeyed:
#   - enableWGCNAThreads()
#   - permutedStatisticsFile not needed (no modulePreservation here)
#   - Ensembl ids unversioned (handled by 00_load.R)
#   - Config written before any computation
#   - All parameters pre-specified in JSON before results are seen
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "mod-000-resolution"
stopifnot(id == "mod-000-resolution")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))
msg("=== ", id, " starting ===")

# ---------- 1. Directories ---------------------------------------------------
run_dir   <- file.path("network/runs",   id)
cfg_dir   <- "network/config"
cache_dir <- "network/cache"
for (d_ in c(run_dir, cfg_dir, cache_dir))
  dir.create(d_, showWarnings = FALSE, recursive = TRUE)

# ---------- 2. Parameters — ALL pre-specified BEFORE any computation ---------
POWER              <- 12L
POWER_FALLBACK     <- 16L     # only if no power-12 cell meets target
NET_TYPE           <- "signed"
COR_FNC            <- "bicor"
GENE_SET           <- "all"
REF_GROUP          <- "control"
MIN_MODULE         <- 30L
MERGE_THRESHOLDS   <- c(0.00, 0.05, 0.10, 0.15, 0.25)
DEEP_SPLITS        <- c(2L, 4L)
SEED               <- 20261009L
TARGET_N_MODS_MIN  <- 10L
TARGET_LARGEST_MAX <- 25.0    # % of total genes

# ---------- 3. Write config BEFORE any computation --------------------------
cfg_path <- file.path(cfg_dir, paste0(id, ".json"))
if (!file.exists(cfg_path)) {
  sha <- tryCatch(
    trimws(system("git rev-parse --short HEAD 2>/dev/null", intern = TRUE)[1]),
    error = function(e) "unknown"
  )
  scl <- function(v) {
    if (is.logical(v))                   return(if (isTRUE(v)) "true" else "false")
    if (is.integer(v) || is.numeric(v))  return(as.character(v[1]))
    paste0('"', gsub('"', '\\"', as.character(v[1]), fixed = TRUE), '"')
  }
  arr <- function(v) paste0("[", paste(vapply(v, scl, ""), collapse = ", "), "]")
  jsn <- function(nm, v) paste0('  "', nm, '": ', if (length(v) > 1) arr(v) else scl(v))
  entries <- list(
    id                   = id,
    family               = "moddetect",
    spec_from_queue      = "BLOCKING DIAGNOSTIC: mergeCutHeight primary axis {0.00,0.05,0.10,0.15,0.25}; deepSplit secondary {2,4}; power fixed at 12",
    power                = POWER,
    power_fallback       = POWER_FALLBACK,
    networkType          = NET_TYPE,
    corFnc               = COR_FNC,
    gene_set             = GENE_SET,
    ref_group            = REF_GROUP,
    minModuleSize        = MIN_MODULE,
    mergeCutHeights      = MERGE_THRESHOLDS,
    deepSplits           = as.integer(DEEP_SPLITS),
    target_n_mods_min    = TARGET_N_MODS_MIN,
    target_largest_max   = TARGET_LARGEST_MAX,
    seed                 = SEED,
    git_sha              = if (length(sha) && !is.na(sha)) sha else "unknown",
    r_version            = paste(R.version$major, R.version$minor, sep = "."),
    wgcna_version        = as.character(packageVersion("WGCNA"))
  )
  lines <- mapply(jsn, names(entries), entries, SIMPLIFY = TRUE)
  writeLines(c("{", paste(lines, collapse = ",\n"), "}"), cfg_path)
  msg("config written: ", cfg_path)
} else {
  msg("config exists: ", cfg_path)
}

# ---------- 4. Packages ------------------------------------------------------
suppressPackageStartupMessages({
  library(WGCNA)
  library(fastcluster)
})
enableWGCNAThreads()
options(stringsAsFactors = FALSE)
msg("WGCNA ", as.character(packageVersion("WGCNA")), " loaded; threads enabled")

# ---------- 5. Data ----------------------------------------------------------
source("network/00_load.R")
d       <- load_project()
refExpr <- expr_for(d, group = REF_GROUP)   # 234 x 12368
msg("refExpr: ", nrow(refExpr), " samples x ", ncol(refExpr), " genes")

# ---------- Helper: build_tom ------------------------------------------------
# Builds bicor signed adjacency, computes TOM, returns list(dissTOM, mean_k).
# Memory: adj (~1.22 GB) + TOM (~1.22 GB) temporarily; peak ~2.44 GB.
build_tom <- function(datExpr, power, net_type, cor_fnc, seed) {
  set.seed(seed)
  msg("Computing adjacency (power=", power, ", corFnc=", cor_fnc, ", type=", net_type, ") ...")
  adj <- adjacency(datExpr, power = power, type = net_type, corFnc = cor_fnc)
  # Mean connectivity: row-sum excluding self (diagonal = 1 before zeroing)
  diag(adj) <- 0
  mean_k <- mean(colSums(adj))
  msg("mean k = ", round(mean_k, 2),
      " (", round(mean_k / ncol(datExpr) * 100, 2), "% of genome)")
  msg("Computing TOM (~10-15 min) ...")
  TOM    <- TOMsimilarity(adj, TOMType = net_type)
  rm(adj); invisible(gc())
  dissTOM <- 1 - TOM
  rm(TOM);  invisible(gc())
  list(dissTOM = dissTOM, mean_k = mean_k)
}

# ---------- Helper: cluster_grid --------------------------------------------
# For one power's dissTOM, runs all (deepSplit x mergeCutHeight) cells.
# Checkpoints per deepSplit (labels) and per cell (result data.frame).
# Returns a data.frame with one row per cell.
cluster_grid <- function(dissTOM, refExpr_data, power, deep_splits,
                         merge_thresholds, min_module, target_n, target_pct,
                         run_d, seed) {
  msg("Building dendrogram (", ncol(dissTOM), " genes) ...")
  gene_tree <- fastcluster::hclust(as.dist(dissTOM), method = "average")
  invisible(gc())   # free the temporary dist object

  all_cells <- list()

  for (ds in deep_splits) {
    ds_ckpt <- file.path(run_d, sprintf("labels_p%d_ds%d.rds", power, ds))
    if (!file.exists(ds_ckpt)) {
      msg("cutreeDynamic deepSplit=", ds, " power=", power, " ...")
      set.seed(seed)
      raw_labels <- cutreeDynamic(
        dendro         = gene_tree,
        distM          = dissTOM,
        deepSplit      = ds,
        minClusterSize = min_module,
        method         = "hybrid",
        verbose        = 2
      )
      saveRDS(raw_labels, ds_ckpt)
      n_pre <- max(raw_labels)
      msg("deepSplit=", ds, ": ", n_pre, " modules detected (pre-merge)")
    } else {
      raw_labels <- readRDS(ds_ckpt)
      n_pre      <- max(raw_labels)
      msg("deepSplit=", ds, " checkpoint loaded: ", n_pre, " modules (pre-merge)")
    }
    raw_colors <- labels2colors(raw_labels)
    n_premerge <- sum(names(table(raw_colors)) != "grey")

    # Eigengene correlation distribution for the pre-merge module set.
    # This tells us why merge cascades when cutHeight is too large.
    ME_data <- tryCatch(
      moduleEigengenes(refExpr_data, colors = raw_colors, verbose = 0)$eigengenes,
      error = function(e) { msg("moduleEigengenes error: ", conditionMessage(e)); NULL }
    )
    ME_cols <- if (!is.null(ME_data)) setdiff(colnames(ME_data), "MEgrey") else character(0)

    if (length(ME_cols) > 1) {
      me_cor    <- cor(ME_data[, ME_cols])
      upper_tri <- me_cor[upper.tri(me_cor)]
      ec_p25 <- round(quantile(upper_tri, 0.25), 3)
      ec_p50 <- round(quantile(upper_tri, 0.50), 3)
      ec_p75 <- round(quantile(upper_tri, 0.75), 3)
      ec_max <- round(max(abs(upper_tri)), 3)
      msg("Pre-merge ME correlations (|r|): median=", ec_p50,
          " p75=", ec_p75, " max=", ec_max)
    } else {
      ec_p25 <- ec_p50 <- ec_p75 <- ec_max <- NA_real_
    }

    for (mch in merge_thresholds) {
      cell_id   <- sprintf("p%d_ds%d_mch%04.2f", power, ds, mch)
      cell_ckpt <- file.path(run_d, sprintf("cell_%s.rds", cell_id))

      if (!file.exists(cell_ckpt)) {
        msg("Cell power=", power, " deepSplit=", ds, " mergeCutHeight=", mch, " ...")
        if (mch > 0) {
          merged_result <- tryCatch(
            mergeCloseModules(
              refExpr_data, raw_colors,
              cutHeight  = mch,
              corFnc     = "cor",   # Pearson is standard for eigengene correlation
              verbose    = 0
            ),
            error = function(e) {
              msg("mergeCloseModules error: ", conditionMessage(e))
              NULL
            }
          )
          final_colors <- if (!is.null(merged_result)) merged_result$colors else raw_colors
        } else {
          final_colors <- raw_colors   # mch=0 means no merging
        }

        tbl         <- table(final_colors)
        non_grey    <- tbl[names(tbl) != "grey"]
        n_mods      <- length(non_grey)
        n_grey_g    <- as.integer(if ("grey" %in% names(tbl)) tbl["grey"] else 0L)
        tot_genes   <- sum(tbl)
        largest_pct <- if (n_mods > 0) round(max(non_grey) / tot_genes * 100, 2) else 100.0
        grey_pct    <- round(n_grey_g / tot_genes * 100, 2)
        med_size    <- if (n_mods > 0) round(median(as.integer(non_grey)), 0) else 0L
        meets       <- as.integer(n_mods >= target_n && largest_pct < target_pct)

        msg("  -> ", n_mods, " modules; largest=", largest_pct,
            "%; grey=", grey_pct, "%; meets_target=", meets)

        cell_df <- data.frame(
          power              = power,
          deepSplit          = ds,
          mergeCutHeight     = mch,
          n_modules_premerge = n_premerge,
          n_modules          = n_mods,
          n_grey             = n_grey_g,
          largest_pct        = largest_pct,
          grey_pct           = grey_pct,
          median_size        = med_size,
          eigene_cor_p25     = ec_p25,
          eigene_cor_p50     = ec_p50,
          eigene_cor_p75     = ec_p75,
          eigene_cor_max     = ec_max,
          meets_target       = meets,
          stringsAsFactors   = FALSE
        )
        saveRDS(cell_df, cell_ckpt)
      } else {
        cell_df <- readRDS(cell_ckpt)
        msg("Cell ", cell_id, " loaded from checkpoint: ",
            cell_df$n_modules, " modules, meets_target=", cell_df$meets_target)
      }
      all_cells[[cell_id]] <- cell_df
    }
  }
  do.call(rbind, all_cells)
}

# ---------- 6. Power-12 grid -------------------------------------------------
tom_ckpt_12 <- file.path(run_dir, "dissTOM_p12.rds")
if (!file.exists(tom_ckpt_12)) {
  msg("Building power-12 TOM (~47 min for n=234) ...")
  tom12 <- build_tom(refExpr, POWER, NET_TYPE, COR_FNC, SEED)
  saveRDS(tom12, tom_ckpt_12, compress = FALSE)   # ~1.2 GB; compress=FALSE for speed
  msg("Power-12 TOM checkpoint saved: ", tom_ckpt_12)
} else {
  msg("Loading power-12 TOM from checkpoint ...")
  tom12 <- readRDS(tom_ckpt_12)
}
dissTOM_12 <- tom12$dissTOM
mean_k_12  <- tom12$mean_k
rm(tom12); invisible(gc())

msg("Running resolution grid at power ", POWER, " (mean k=", round(mean_k_12, 2), ") ...")
grid12        <- cluster_grid(
  dissTOM          = dissTOM_12,
  refExpr_data     = refExpr,
  power            = POWER,
  deep_splits      = DEEP_SPLITS,
  merge_thresholds = MERGE_THRESHOLDS,
  min_module       = MIN_MODULE,
  target_n         = TARGET_N_MODS_MIN,
  target_pct       = TARGET_LARGEST_MAX,
  run_d            = run_dir,
  seed             = SEED
)
grid12$mean_k <- round(mean_k_12, 2)
rm(dissTOM_12); invisible(gc())

# ---------- 7. Power-16 fallback (only if no power-12 cell meets target) ----
n_meeting_12 <- sum(grid12$meets_target, na.rm = TRUE)
msg("Power-12: ", n_meeting_12, " of ", nrow(grid12),
    " cells meet target (n_modules>=", TARGET_N_MODS_MIN,
    " AND largest_pct<", TARGET_LARGEST_MAX, "%)")

grid16 <- NULL
if (n_meeting_12 == 0L) {
  msg("No power-12 cell meets target. Attempting power-", POWER_FALLBACK, " fallback ...")
  tom_ckpt_16 <- file.path(run_dir, "dissTOM_p16.rds")
  if (!file.exists(tom_ckpt_16)) {
    msg("Building power-16 TOM (~47 min) ...")
    tom16 <- build_tom(refExpr, POWER_FALLBACK, NET_TYPE, COR_FNC, SEED + 1L)
    saveRDS(tom16, tom_ckpt_16, compress = FALSE)
    msg("Power-16 TOM checkpoint saved")
  } else {
    msg("Loading power-16 TOM from checkpoint ...")
    tom16 <- readRDS(tom_ckpt_16)
  }
  dissTOM_16 <- tom16$dissTOM
  mean_k_16  <- tom16$mean_k
  rm(tom16); invisible(gc())

  grid16        <- cluster_grid(
    dissTOM          = dissTOM_16,
    refExpr_data     = refExpr,
    power            = POWER_FALLBACK,
    deep_splits      = DEEP_SPLITS,
    merge_thresholds = MERGE_THRESHOLDS,
    min_module       = MIN_MODULE,
    target_n         = TARGET_N_MODS_MIN,
    target_pct       = TARGET_LARGEST_MAX,
    run_d            = run_dir,
    seed             = SEED + 1L
  )
  grid16$mean_k <- round(mean_k_16, 2)
  rm(dissTOM_16); invisible(gc())

  n_meeting_16 <- sum(grid16$meets_target, na.rm = TRUE)
  msg("Power-16: ", n_meeting_16, " of ", nrow(grid16), " cells meet target")
} else {
  msg("Power-12 has meeting cells; skipping power-", POWER_FALLBACK, " fallback")
}

# ---------- 8. Combine and write results ------------------------------------
all_grid <- if (!is.null(grid16)) rbind(grid12, grid16) else grid12

results_path <- file.path(run_dir, "resolution_grid.csv")
write.csv(all_grid, results_path, row.names = FALSE)
msg("Results written: ", results_path)

# Summary to log
msg("\n=== RESOLUTION GRID (power x deepSplit x mergeCutHeight) ===")
print(all_grid[, c("power", "deepSplit", "mergeCutHeight",
                   "n_modules_premerge", "n_modules",
                   "largest_pct", "grey_pct", "median_size",
                   "eigene_cor_p75", "eigene_cor_max",
                   "meets_target")],
      row.names = FALSE)

# Cells meeting target, sorted by n_modules descending (most informative first)
target_cells <- all_grid[all_grid$meets_target == 1L, ]
if (nrow(target_cells) > 0L) {
  target_cells <- target_cells[order(-target_cells$n_modules), ]
  msg("\n=== CELLS MEETING TARGET (n_modules>=", TARGET_N_MODS_MIN,
      " AND largest_pct<", TARGET_LARGEST_MAX, "%) ===")
  print(target_cells[, c("power", "deepSplit", "mergeCutHeight",
                         "n_modules", "largest_pct", "grey_pct", "median_size")],
        row.names = FALSE)
  best <- target_cells[1L, ]
  msg("RECOMMENDED: power=", best$power,
      " deepSplit=", best$deepSplit,
      " mergeCutHeight=", best$mergeCutHeight,
      " -> ", best$n_modules, " modules, largest=", best$largest_pct, "%")
} else {
  msg("\nNO CELL MEETS THE TARGET at any tested parameter combination.")
  msg("This data does not admit a usable WGCNA decomposition at these settings.")
  msg("Report as negative result: further investigation needed (different gene sets,")
  msg("residualised expression, or higher powers beyond 16).")
}

msg("=== ", id, " COMPLETE ===")
# Exit 0 → driver marks the row DONE.
