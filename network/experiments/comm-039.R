#!/usr/bin/env Rscript
# ===========================================================================
# comm-039  Leiden community detection on adjacency matrix
#
# Spec: "Leiden on adjacency; compare to Louvain and to WGCNA (ARI, NMI)"
#
# Design:
#   - Same adjacency pipeline as comm-038 (WGCNA adjacency, threshold 0.1).
#   - Run Leiden (cluster_leiden) if available in igraph 2.3.4; else fall back
#     to Louvain with a documented note (tryCatch detection at script start).
#   - Same ARI/NMI helpers (contingency table; no external packages).
#   - Also compare to comm-038 Louvain results if available.
#   - bp_lith: ONE fixed subsample (n=74, same seed derivation as comm-038).
#
# Outputs:
#   network/runs/comm-039/leiden_membership_control.csv
#   network/runs/comm-039/leiden_membership_nolith.csv
#   network/runs/comm-039/leiden_membership_lith.csv
#   network/runs/comm-039/community_comparison.csv
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "comm-039"
stopifnot(id == "comm-039")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

run_dir   <- file.path("network/runs",   id)
cfg_dir   <- "network/config"
cache_dir <- "network/cache"
for (d_ in c(run_dir, cfg_dir, cache_dir))
  dir.create(d_, showWarnings = FALSE, recursive = TRUE)

# ---- packages (load before detection check) --------------------------------
suppressPackageStartupMessages({
  library(WGCNA)
  library(fastcluster)
  library(igraph)
})
enableWGCNAThreads()
options(stringsAsFactors = FALSE)
msg("WGCNA ", as.character(packageVersion("WGCNA")),
    " | igraph ", as.character(packageVersion("igraph")), " loaded")

# ---- detect Leiden availability BEFORE writing config ----------------------
has_leiden <- exists("cluster_leiden", where = "package:igraph", inherits = FALSE)
if (!has_leiden) {
  msg("cluster_leiden not found in igraph ", as.character(packageVersion("igraph")),
      " -- using cluster_louvain() as fallback")
  msg("NOTE: this is documented in config as leiden_available=false; results are",
      " Louvain-based and should be interpreted accordingly")
  ALGORITHM_USED <- "louvain_fallback"
} else {
  msg("cluster_leiden found in igraph ", as.character(packageVersion("igraph")))
  ALGORITHM_USED <- "leiden"
}

# ---- parameters (ALL pre-specified) ----------------------------------------
SEED_BASE          <- 20261009L + 39000L
POWER              <- 12L
NET_TYPE           <- "signed"
COR_FNC            <- "bicor"
ADJ_THRESHOLD      <- 0.1
RESOLUTIONS        <- c(0.5, 1, 2)
SUBSAMPLE_N        <- 74L
LITH_SINGLE_SEED   <- SEED_BASE + 1L
REF_CACHE          <- "network/cache/ref_control_all_bicor_signed_p12.rds"
COMM038_DIR        <- "network/runs/comm-038"

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
    spec                  = "Leiden on adjacency; compare to Louvain and to WGCNA (ARI, NMI)",
    seed_base             = SEED_BASE,
    power                 = POWER,
    networkType           = NET_TYPE,
    corFnc                = COR_FNC,
    adjacency_threshold   = ADJ_THRESHOLD,
    adjacency_threshold_note = "Same as comm-038: values below 0.1 set to zero before igraph construction",
    resolutions_tested    = RESOLUTIONS,
    leiden_available      = has_leiden,
    algorithm_used        = ALGORITHM_USED,
    leiden_fallback_note  = if (!has_leiden) paste0(
      "cluster_leiden not found in igraph ",
      as.character(packageVersion("igraph")),
      "; cluster_louvain() used as fallback. Results are Louvain-based. ",
      "Re-run with igraph >= 2.x that ships Leiden to get true Leiden results."
    ) else "cluster_leiden found; true Leiden algorithm used",
    groups_tested         = c("control", "bp_nolith", "bp_lith"),
    subsample_n           = SUBSAMPLE_N,
    bp_lith_draws         = 1L,
    bp_lith_single_draw_note = "bp_lith uses ONE fixed subsample (n=74) for speed; same design as comm-038",
    comm038_comparison    = "ARI/NMI vs comm-038 Louvain computed if network/runs/comm-038/ output exists; skipped with note if absent",
    ref_cache             = REF_CACHE,
    comparison_metrics    = c("ARI", "NMI"),
    git_sha               = if (length(sha) && !is.na(sha)) sha else "unknown",
    r_version             = paste(R.version$major, R.version$minor, sep = "."),
    wgcna_version         = as.character(packageVersion("WGCNA")),
    igraph_version        = as.character(packageVersion("igraph"))
  )
  lines <- mapply(jsn, names(cfg), cfg, SIMPLIFY = TRUE)
  writeLines(c("{", paste(lines, collapse = ",\n"), "}"), cfg_path)
  msg("config written: ", cfg_path)
}

# ---- data ------------------------------------------------------------------
source("network/00_load.R")
d <- load_project()

# ---- load base-001 reference modules ----------------------------------------
if (!file.exists(REF_CACHE)) {
  stop("Base-001 reference cache not found: ", REF_CACHE,
       "\nRun base-001.R first.")
}
ref_net      <- readRDS(REF_CACHE)
wgcna_colors <- ref_net$moduleColors
genes        <- names(wgcna_colors)
msg("base-001 reference loaded: ", length(unique(wgcna_colors)) - 1L,
    " non-grey WGCNA modules")

# ---- ARI and NMI helpers ---------------------------------------------------
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
  mi  <- sum(pab[pab > 0] * log(pab[pab > 0] / outer(pa, pb))[pab > 0])
  mi / sqrt(ha * hb)
}

# ---- community detection wrapper -------------------------------------------
run_cluster <- function(g, resolution, seed) {
  if (has_leiden) {
    set.seed(seed)
    cluster_leiden(g, resolution_parameter = resolution,
                   weights = E(g)$weight)
  } else {
    # Louvain fallback: resolution parameter supported in igraph >= 1.3
    set.seed(seed)
    cluster_louvain(g, weights = E(g)$weight, resolution = resolution)
  }
}

# ---- main computation function per group -----------------------------------
run_leiden_for_group <- function(group_label, expr_mat, wgcna_col_vec, seed,
                                  louvain_mem_path = NULL) {
  msg("--- ", group_label, ": computing adjacency (", nrow(expr_mat), " samples) ---")

  # 1. Adjacency
  set.seed(seed)
  A <- adjacency(expr_mat, power = POWER, type = "signed", corFnc = "bicor")
  diag(A) <- 0

  # 2. Threshold
  A_thresh <- A
  A_thresh[A_thresh < ADJ_THRESHOLD] <- 0
  n_edges_before <- sum(A > 0) / 2
  n_edges_after  <- sum(A_thresh > 0) / 2
  msg("  edges before threshold: ", n_edges_before,
      "  after: ", n_edges_after,
      "  (", round(100 * n_edges_after / max(n_edges_before, 1), 1), "% retained)")
  rm(A); gc()

  # 3. igraph
  g <- graph_from_adjacency_matrix(A_thresh, mode = "undirected",
                                    weighted = TRUE, diag = FALSE)
  rm(A_thresh); gc()
  msg("  igraph: ", vcount(g), " vertices, ", ecount(g), " edges")

  # 4. Leiden/Louvain at multiple resolutions
  comparison_rows <- list()
  best_comm       <- NULL

  for (res in RESOLUTIONS) {
    comm    <- run_cluster(g, resolution = res, seed = seed + as.integer(res * 100))
    mem_vec <- membership(comm)
    n_comm  <- length(unique(mem_vec))
    msg("  ", ALGORITHM_USED, " res=", res, ": ", n_comm, " communities")

    ari_wgcna <- contingency_ari(mem_vec, wgcna_col_vec)
    nmi_wgcna <- nmi(mem_vec, wgcna_col_vec)
    msg("    vs WGCNA: ARI=", round(ari_wgcna, 4), "  NMI=", round(nmi_wgcna, 4))

    # compare to comm-038 Louvain if available
    ari_louvain <- NA_real_; nmi_louvain <- NA_real_
    if (!is.null(louvain_mem_path) && file.exists(louvain_mem_path) && res == 1) {
      louvain_df  <- read.csv(louvain_mem_path, stringsAsFactors = FALSE)
      louvain_mem <- louvain_df$louvain_community[match(genes, louvain_df$gene)]
      if (!anyNA(louvain_mem)) {
        ari_louvain <- contingency_ari(mem_vec, louvain_mem)
        nmi_louvain <- nmi(mem_vec, louvain_mem)
        msg("    vs Louvain(comm-038, res=1): ARI=", round(ari_louvain, 4),
            "  NMI=", round(nmi_louvain, 4))
      } else {
        msg("    comm-038 Louvain membership has NAs for some genes; skipping comparison")
      }
    } else if (!is.null(louvain_mem_path) && !file.exists(louvain_mem_path) && res == 1) {
      msg("    comm-038 output not found at ", louvain_mem_path, "; skipping Louvain comparison")
    }

    comparison_rows[[length(comparison_rows) + 1L]] <- data.frame(
      group                      = group_label,
      algorithm                  = ALGORITHM_USED,
      resolution                 = res,
      n_communities              = n_comm,
      ARI_vs_wgcna               = ari_wgcna,
      NMI_vs_wgcna               = nmi_wgcna,
      ARI_vs_louvain_comm038     = ari_louvain,
      NMI_vs_louvain_comm038     = nmi_louvain,
      n_genes                    = length(mem_vec),
      n_samples                  = nrow(expr_mat),
      stringsAsFactors           = FALSE
    )

    if (res == 1) best_comm <- mem_vec
  }

  mem_df <- data.frame(
    gene              = genes,
    wgcna_module      = wgcna_col_vec,
    leiden_community  = as.integer(best_comm),
    stringsAsFactors  = FALSE
  )

  list(membership = mem_df, comparison = do.call(rbind, comparison_rows))
}

# ---- control (all 234) -----------------------------------------------------
ctrl_ckpt <- file.path(run_dir, "leiden_ctrl_result.rds")
if (!file.exists(ctrl_ckpt)) {
  E_ctrl <- expr_for(d, group = "control")
  louvain_ctrl_path <- file.path(COMM038_DIR, "louvain_membership_control.csv")
  res_ctrl <- run_leiden_for_group("control", E_ctrl, wgcna_colors,
                                    SEED_BASE, louvain_ctrl_path)
  saveRDS(res_ctrl, ctrl_ckpt)
  write.csv(res_ctrl$membership,
            file.path(run_dir, "leiden_membership_control.csv"), row.names = FALSE)
  msg("control done; membership written")
  rm(E_ctrl, res_ctrl); gc()
} else {
  msg("control checkpoint exists, skipping")
  res_ctrl <- readRDS(ctrl_ckpt)
  write.csv(res_ctrl$membership,
            file.path(run_dir, "leiden_membership_control.csv"), row.names = FALSE)
}

# ---- bp_nolith (all 74) ----------------------------------------------------
nolith_ckpt <- file.path(run_dir, "leiden_nolith_result.rds")
if (!file.exists(nolith_ckpt)) {
  E_nolith <- expr_for(d, group = "bp_nolith")
  louvain_nolith_path <- file.path(COMM038_DIR, "louvain_membership_nolith.csv")
  res_nolith <- run_leiden_for_group("bp_nolith", E_nolith, wgcna_colors,
                                      SEED_BASE + 10000L, louvain_nolith_path)
  saveRDS(res_nolith, nolith_ckpt)
  write.csv(res_nolith$membership,
            file.path(run_dir, "leiden_membership_nolith.csv"), row.names = FALSE)
  msg("bp_nolith done; membership written")
  rm(E_nolith, res_nolith); gc()
} else {
  msg("bp_nolith checkpoint exists, skipping")
  res_nolith <- readRDS(nolith_ckpt)
  write.csv(res_nolith$membership,
            file.path(run_dir, "leiden_membership_nolith.csv"), row.names = FALSE)
}

# ---- bp_lith (ONE subsample, n=74) -----------------------------------------
lith_ckpt <- file.path(run_dir, "leiden_lith_result.rds")
if (!file.exists(lith_ckpt)) {
  set.seed(LITH_SINGLE_SEED)
  lith_all_sids <- colnames(d$logtpm)[!is.na(d$groups) & d$groups == "bp_lith"]
  stopifnot(length(lith_all_sids) == 152L)
  sel_sids <- sample(lith_all_sids, SUBSAMPLE_N, replace = FALSE)
  E_lith   <- t(d$logtpm[, sel_sids])
  msg("bp_lith subsample: ", length(sel_sids), " of 152 (seed=", LITH_SINGLE_SEED, ")")
  louvain_lith_path <- file.path(COMM038_DIR, "louvain_membership_lith.csv")
  res_lith <- run_leiden_for_group("bp_lith", E_lith, wgcna_colors,
                                    LITH_SINGLE_SEED + 1000L, louvain_lith_path)
  res_lith$lith_samples <- paste(sel_sids, collapse = ",")
  saveRDS(res_lith, lith_ckpt)
  write.csv(res_lith$membership,
            file.path(run_dir, "leiden_membership_lith.csv"), row.names = FALSE)
  msg("bp_lith done; membership written")
  rm(E_lith, res_lith); gc()
} else {
  msg("bp_lith checkpoint exists, skipping")
  res_lith <- readRDS(lith_ckpt)
  write.csv(res_lith$membership,
            file.path(run_dir, "leiden_membership_lith.csv"), row.names = FALSE)
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
comparison$is_primary_resolution <- comparison$resolution == 1

write.csv(comparison,
          file.path(run_dir, "community_comparison.csv"), row.names = FALSE)
msg("community_comparison.csv written (", nrow(comparison), " rows)")

pri <- comparison[comparison$is_primary_resolution, ]
msg("\n--- Community comparison summary (resolution=1, algorithm=", ALGORITHM_USED, ") ---")
for (i in seq_len(nrow(pri))) {
  r <- pri[i, ]
  msg("  ", r$group, ": n_communities=", r$n_communities,
      "  ARI_wgcna=", round(r$ARI_vs_wgcna, 4),
      "  NMI_wgcna=", round(r$NMI_vs_wgcna, 4),
      if (!is.na(r$ARI_vs_louvain_comm038))
        paste0("  ARI_louvain=", round(r$ARI_vs_louvain_comm038, 4)) else "")
}

msg("=== ", id, " COMPLETE ===")
