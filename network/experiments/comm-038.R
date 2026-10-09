#!/usr/bin/env Rscript
# ===========================================================================
# comm-038  Louvain community detection on adjacency matrix
#
# Spec: "Louvain on adjacency; compare to WGCNA modules (ARI, NMI)"
#
# Design:
#   - Load base-001 WGCNA reference network from cache.
#   - For each group (control all, bp_nolith all, bp_lith subsample n=74):
#       1. Compute adjacency via WGCNA adjacency() (signed bicor, power 12).
#       2. Threshold adjacency at 0.1 for computational feasibility.
#       3. Build igraph from thresholded adjacency.
#       4. Run Louvain at resolution 0.5, 1.0, 2.0 (robustness check).
#       5. Compare Louvain communities to base-001 WGCNA moduleColors via ARI
#          and NMI (implemented with contingency tables; no external packages).
#   - bp_lith: ONE subsample draw (not 20) for speed; documented in config.
#   - Outputs:
#       network/runs/comm-038/louvain_membership_control.csv
#       network/runs/comm-038/louvain_membership_nolith.csv
#       network/runs/comm-038/louvain_membership_lith.csv
#       network/runs/comm-038/community_comparison.csv
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "comm-038"
stopifnot(id == "comm-038")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

run_dir   <- file.path("network/runs",   id)
cfg_dir   <- "network/config"
cache_dir <- "network/cache"
for (d_ in c(run_dir, cfg_dir, cache_dir))
  dir.create(d_, showWarnings = FALSE, recursive = TRUE)

# ---- parameters (ALL pre-specified) ----------------------------------------
SEED_BASE          <- 20261009L + 38000L
POWER              <- 12L
NET_TYPE           <- "signed"
COR_FNC            <- "bicor"
ADJ_THRESHOLD      <- 0.1        # threshold adjacency for igraph construction
LOUVAIN_RESOLUTIONS <- c(0.5, 1, 2)  # robustness check
SUBSAMPLE_N        <- 74L
# bp_lith: ONE fixed subsample for speed
LITH_SINGLE_SEED   <- SEED_BASE + 1L
REF_CACHE          <- "network/cache/ref_control_all_bicor_signed_p12.rds"

# ---- config ----------------------------------------------------------------
cfg_path <- file.path(cfg_dir, paste0(id, ".json"))
if (!file.exists(cfg_path)) {
  scl <- function(v) {
    if (is.logical(v)) return(if (isTRUE(v)) "true" else "false")
    if (is.integer(v) || is.numeric(v)) return(as.character(v[1]))
    paste0('"', gsub('"', '\\"', as.character(v[1]), fixed = TRUE), '"')
  }
  arr <- function(v) paste0("[", paste(vapply(v, scl, ""), collapse = ", "), "]")
  jsn <- function(nm, v) paste0('  "', nm, '": ', if (length(v) > 1) arr(v) else scl(v))
  sha <- tryCatch(trimws(system("git rev-parse --short HEAD 2>/dev/null", intern=TRUE)[1]),
                  error = function(e) "unknown")
  cfg <- list(
    id                    = id,
    family                = "community",
    spec                  = "Louvain on adjacency; compare to WGCNA modules (ARI, NMI)",
    seed_base             = SEED_BASE,
    power                 = POWER,
    networkType           = NET_TYPE,
    corFnc                = COR_FNC,
    adjacency_threshold   = ADJ_THRESHOLD,
    adjacency_threshold_note = "Adjacency values below 0.1 set to zero before igraph construction; reduces edge count for computational feasibility while retaining strong co-expression links",
    louvain_resolution    = LOUVAIN_RESOLUTIONS,
    louvain_resolution_note = "Louvain run at three resolutions (0.5, 1, 2) for robustness; resolution=1 is the igraph default and primary result",
    groups_tested         = c("control", "bp_nolith", "bp_lith"),
    subsample_n           = SUBSAMPLE_N,
    bp_lith_draws         = 1L,
    bp_lith_single_draw_note = "bp_lith uses ONE fixed subsample (not 20) for computational speed; documents variability not averaged here",
    ref_cache             = REF_CACHE,
    comparison_metrics    = c("ARI", "NMI"),
    comparison_note       = "ARI and NMI computed from contingency tables using custom R functions; no external packages required",
    git_sha               = if (length(sha) && !is.na(sha)) sha else "unknown",
    r_version             = paste(R.version$major, R.version$minor, sep = "."),
    wgcna_version         = as.character(packageVersion("WGCNA")),
    igraph_version        = as.character(packageVersion("igraph"))
  )
  lines <- mapply(jsn, names(cfg), cfg, SIMPLIFY = TRUE)
  writeLines(c("{", paste(lines, collapse = ",\n"), "}"), cfg_path)
  msg("config written: ", cfg_path)
}

# ---- packages --------------------------------------------------------------
suppressPackageStartupMessages({
  library(WGCNA)
  library(fastcluster)
  library(igraph)
})
enableWGCNAThreads()
options(stringsAsFactors = FALSE)
msg("WGCNA ", as.character(packageVersion("WGCNA")),
    " | igraph ", as.character(packageVersion("igraph")), " loaded")

# ---- data ------------------------------------------------------------------
source("network/00_load.R")
d <- load_project()

# ---- load base-001 reference modules ----------------------------------------
if (!file.exists(REF_CACHE)) {
  stop("Base-001 reference cache not found: ", REF_CACHE,
       "\nRun base-001.R first to populate the cache.")
}
ref_net      <- readRDS(REF_CACHE)
wgcna_colors <- ref_net$moduleColors   # named character vector: gene -> colour
genes        <- names(wgcna_colors)
stopifnot(length(genes) == nrow(d$logtpm))
msg("base-001 reference loaded: ", length(unique(wgcna_colors)) - 1L,
    " non-grey WGCNA modules; ", length(genes), " genes")

# ---- ARI and NMI helpers (no external packages) ----------------------------
contingency_ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  sum_a <- rowSums(tab); sum_b <- colSums(tab)
  sum_c2 <- sum(choose(tab, 2))
  a2 <- sum(choose(sum_a, 2)); b2 <- sum(choose(sum_b, 2))
  expected <- a2 * b2 / choose(n, 2)
  (sum_c2 - expected) / ((a2 + b2) / 2 - expected)
}

nmi <- function(a, b) {
  tab <- as.matrix(table(a, b)); n <- sum(tab)
  pa <- rowSums(tab) / n; pb <- colSums(tab) / n
  ha <- -sum(pa[pa > 0] * log(pa[pa > 0]))
  hb <- -sum(pb[pb > 0] * log(pb[pb > 0]))
  pab <- tab / n
  mi <- sum(pab[pab > 0] * log(pab[pab > 0] / outer(pa, pb))[pab > 0])
  mi / sqrt(ha * hb)
}

# ---- main computation function per group -----------------------------------
# Returns: list(membership_df, comparison_rows)
run_louvain_for_group <- function(group_label, expr_mat, wgcna_col_vec, seed) {
  msg("--- ", group_label, ": computing adjacency (", nrow(expr_mat), " samples) ---")

  # 1. WGCNA adjacency: genes x genes
  set.seed(seed)
  A <- adjacency(expr_mat, power = POWER, type = "signed", corFnc = "bicor")
  # A is genes x genes (12368 x 12368); diag = 1 by WGCNA convention
  diag(A) <- 0   # zero diagonal before thresholding and igraph construction

  # 2. Threshold at ADJ_THRESHOLD
  A_thresh           <- A
  A_thresh[A_thresh < ADJ_THRESHOLD] <- 0
  n_edges_before     <- sum(A > 0) / 2
  n_edges_after      <- sum(A_thresh > 0) / 2
  msg("  edges before threshold: ", n_edges_before,
      "  after: ", n_edges_after,
      "  (", round(100 * n_edges_after / max(n_edges_before, 1), 1), "% retained)")
  rm(A); gc()

  # 3. Build igraph
  g <- graph_from_adjacency_matrix(A_thresh, mode = "undirected",
                                    weighted = TRUE, diag = FALSE)
  rm(A_thresh); gc()
  msg("  igraph: ", vcount(g), " vertices, ", ecount(g), " edges")

  # 4. Louvain at multiple resolutions
  comparison_rows <- list()
  best_comm       <- NULL   # resolution=1 is primary

  for (res in LOUVAIN_RESOLUTIONS) {
    set.seed(seed + as.integer(res * 100))
    comm <- cluster_louvain(g, weights = E(g)$weight, resolution = res)
    mem_vec <- membership(comm)   # named integer: gene -> community id
    n_comm  <- length(unique(mem_vec))
    msg("  Louvain res=", res, ": ", n_comm, " communities")

    # align with wgcna_col_vec (same gene order as rownames of logtpm)
    ari_val <- contingency_ari(mem_vec, wgcna_col_vec)
    nmi_val <- nmi(mem_vec, wgcna_col_vec)
    msg("    ARI=", round(ari_val, 4), "  NMI=", round(nmi_val, 4))

    comparison_rows[[length(comparison_rows) + 1L]] <- data.frame(
      group                 = group_label,
      resolution            = res,
      n_louvain_communities = n_comm,
      ARI_vs_wgcna          = ari_val,
      NMI_vs_wgcna          = nmi_val,
      n_genes               = length(mem_vec),
      n_samples             = nrow(expr_mat),
      stringsAsFactors      = FALSE
    )

    if (res == 1) best_comm <- mem_vec
  }

  # 5. Membership data frame (using resolution=1 as primary)
  mem_df <- data.frame(
    gene             = genes,
    wgcna_module     = wgcna_col_vec,
    louvain_community = as.integer(best_comm),
    stringsAsFactors = FALSE
  )

  list(membership = mem_df, comparison = do.call(rbind, comparison_rows))
}

# ---- control (all 234) -----------------------------------------------------
ctrl_ckpt <- file.path(run_dir, "louvain_ctrl_result.rds")
if (!file.exists(ctrl_ckpt)) {
  E_ctrl <- expr_for(d, group = "control")   # 234 x 12368
  res_ctrl <- run_louvain_for_group("control", E_ctrl, wgcna_colors, SEED_BASE)
  saveRDS(res_ctrl, ctrl_ckpt)
  write.csv(res_ctrl$membership,
            file.path(run_dir, "louvain_membership_control.csv"), row.names = FALSE)
  msg("control done; membership written")
  rm(E_ctrl, res_ctrl); gc()
} else {
  msg("control checkpoint exists, skipping computation")
  res_ctrl <- readRDS(ctrl_ckpt)
  write.csv(res_ctrl$membership,
            file.path(run_dir, "louvain_membership_control.csv"), row.names = FALSE)
}

# ---- bp_nolith (all 74) ----------------------------------------------------
nolith_ckpt <- file.path(run_dir, "louvain_nolith_result.rds")
if (!file.exists(nolith_ckpt)) {
  E_nolith <- expr_for(d, group = "bp_nolith")   # 74 x 12368
  res_nolith <- run_louvain_for_group("bp_nolith", E_nolith, wgcna_colors,
                                       SEED_BASE + 10000L)
  saveRDS(res_nolith, nolith_ckpt)
  write.csv(res_nolith$membership,
            file.path(run_dir, "louvain_membership_nolith.csv"), row.names = FALSE)
  msg("bp_nolith done; membership written")
  rm(E_nolith, res_nolith); gc()
} else {
  msg("bp_nolith checkpoint exists, skipping computation")
  res_nolith <- readRDS(nolith_ckpt)
  write.csv(res_nolith$membership,
            file.path(run_dir, "louvain_membership_nolith.csv"), row.names = FALSE)
}

# ---- bp_lith (ONE subsample, n=74) -----------------------------------------
# Design note: bp_lith has 152 samples; we subsample to 74 (LITH_SINGLE_SEED)
# for one fixed draw. This is intentional (speed); draw variance not averaged here.
lith_ckpt <- file.path(run_dir, "louvain_lith_result.rds")
if (!file.exists(lith_ckpt)) {
  set.seed(LITH_SINGLE_SEED)
  lith_all_sids <- colnames(d$logtpm)[!is.na(d$groups) & d$groups == "bp_lith"]
  stopifnot(length(lith_all_sids) == 152L)
  sel_sids <- sample(lith_all_sids, SUBSAMPLE_N, replace = FALSE)
  E_lith   <- t(d$logtpm[, sel_sids])   # 74 x 12368
  msg("bp_lith subsample: ", length(sel_sids), " of 152 samples (seed=",
      LITH_SINGLE_SEED, ")")
  res_lith <- run_louvain_for_group("bp_lith", E_lith, wgcna_colors,
                                     LITH_SINGLE_SEED + 1000L)
  res_lith$lith_samples <- paste(sel_sids, collapse = ",")
  saveRDS(res_lith, lith_ckpt)
  write.csv(res_lith$membership,
            file.path(run_dir, "louvain_membership_lith.csv"), row.names = FALSE)
  msg("bp_lith done; membership written")
  rm(E_lith, res_lith); gc()
} else {
  msg("bp_lith checkpoint exists, skipping computation")
  res_lith <- readRDS(lith_ckpt)
  write.csv(res_lith$membership,
            file.path(run_dir, "louvain_membership_lith.csv"), row.names = FALSE)
}

# ---- aggregate community_comparison.csv ------------------------------------
all_done <- file.exists(ctrl_ckpt) && file.exists(nolith_ckpt) && file.exists(lith_ckpt)
if (!all_done) {
  msg("not all groups complete; exiting with status 1 for resume")
  quit(save = "no", status = 1L)
}

res_ctrl   <- readRDS(ctrl_ckpt)
res_nolith <- readRDS(nolith_ckpt)
res_lith   <- readRDS(lith_ckpt)

comparison <- rbind(res_ctrl$comparison, res_nolith$comparison, res_lith$comparison)
# flag primary resolution
comparison$is_primary_resolution <- comparison$resolution == 1

write.csv(comparison,
          file.path(run_dir, "community_comparison.csv"), row.names = FALSE)
msg("community_comparison.csv written (", nrow(comparison), " rows)")

# brief summary to log
pri <- comparison[comparison$is_primary_resolution, ]
msg("\n--- Community comparison summary (resolution=1) ---")
for (i in seq_len(nrow(pri))) {
  r <- pri[i, ]
  msg("  ", r$group, ": n_communities=", r$n_louvain_communities,
      "  ARI=", round(r$ARI_vs_wgcna, 4),
      "  NMI=", round(r$NMI_vs_wgcna, 4))
}

msg("=== ", id, " COMPLETE ===")
