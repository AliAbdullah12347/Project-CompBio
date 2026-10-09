#!/usr/bin/env Rscript
# ===========================================================================
# hub-037  Differential connectivity per gene between groups
#
# Spec (from queue):
#   differential connectivity per gene between groups, same permutation null
#
# Design (K revised to 20 — same budget constraint as base-001):
#   1. Requires base-001 reference cache (for moduleColors) AND hub-036 results
#      (kME per group). Dependency checks: quit(status=1) if either is missing.
#   2. For each gene: compare kME in control vs bp_lith.
#        - delta_kME = kME_lith_mean - kME_ctrl (per module column)
#        - Also compute delta_kME_ctrl_vs_nolith for reference
#      Connectivity summarised per gene: intramodular kME (own module column).
#   3. Permutation null (K=20 revised from K=500):
#        - For each of K permuted draws: randomly select 74 samples from all 474
#          as pseudo-"bp_lith"; compute kME vs control kME.
#        - delta_kME under null for each gene
#        - Empirical p-value per gene: fraction of K null |delta_kME| >= observed
#        - BH correction across genes within each module separately.
#
# Non-negotiable design points obeyed:
#   - K implemented = 20 (spec = 500; documented in config)
#   - SEED_BASE = 20261009L + 37000L
#   - Config written before any computation
#   - Checkpoint per permuted draw; skip if file exists
#   - Empirical p minimum = 1/K = 0.05 (noted in config)
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "hub-037"
stopifnot(id == "hub-037")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

# ---------- 1. Directories ---------------------------------------------------
run_dir <- file.path("network/runs", id)
cfg_dir <- "network/config"
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- 2. Parameters — ALL pre-specified BEFORE any computation ---------
SEED_BASE      <- 20261009L + 37000L
K              <- 20L        # spec = 500; revised for budget (see METHODS.md 2026-10-09)
SUBSAMPLE_N    <- 74L
REF_CACHE      <- "network/cache/ref_control_all_bicor_signed_p12.rds"
HUB036_CTRL    <- "network/runs/hub-036/kme_control.csv"
HUB036_NOLITH  <- "network/runs/hub-036/kme_nolith.csv"
HUB036_LITH    <- "network/runs/hub-036/kme_lith_mean.csv"

# ---------- 3. Dependency checks ---------------------------------------------
if (!file.exists(REF_CACHE)) {
  msg("DEPENDENCY: base-001 cache missing: ", REF_CACHE)
  quit(save = "no", status = 1L)
}
if (!file.exists(HUB036_CTRL)) {
  msg("DEPENDENCY: hub-036 not done — missing: ", HUB036_CTRL)
  quit(save = "no", status = 1L)
}

# ---------- 4. Write config BEFORE any computation ---------------------------
cfg_path <- file.path(cfg_dir, paste0(id, ".json"))
if (!file.exists(cfg_path)) {
  scl <- function(v) {
    if (is.logical(v)) return(if (isTRUE(v)) "true" else "false")
    if (is.integer(v) || is.numeric(v)) return(as.character(v[1]))
    paste0('"', gsub('"', '\\"', as.character(v[1]), fixed = TRUE), '"')
  }
  arr <- function(v) paste0("[", paste(vapply(v, scl, ""), collapse = ", "), "]")
  jsn <- function(nm, v) paste0('  "', nm, '": ', if (length(v) > 1) arr(v) else scl(v))
  sha <- tryCatch(
    trimws(system("git rev-parse --short HEAD 2>/dev/null", intern = TRUE)[1]),
    error = function(e) "unknown"
  )
  cfg_entries <- list(
    id                  = id,
    family              = "hub",
    spec                = "differential connectivity per gene between groups, same permutation null",
    K_spec              = 500L,
    K_implemented       = K,
    K_revision_note     = "K=500 x compute time infeasible; revised to 20 (same budget constraint as base-001, see METHODS.md 2026-10-09). Minimum empirical p = 1/K = 0.05.",
    ref_cache           = REF_CACHE,
    depends_on          = "hub-036",
    pool_for_permutation = "all 474 samples",
    subsample_n         = SUBSAMPLE_N,
    delta_kme_direction = "kME_lith_mean - kME_ctrl (positive = more connected in lith)",
    correction          = "BH within each module separately, across genes",
    min_empirical_p     = 1.0 / K,
    seed_base           = SEED_BASE,
    git_sha             = if (length(sha) && !is.na(sha)) sha else "unknown",
    r_version           = paste(R.version$major, R.version$minor, sep = "."),
    wgcna_version       = as.character(packageVersion("WGCNA"))
  )
  lines <- mapply(jsn, names(cfg_entries), cfg_entries, SIMPLIFY = TRUE)
  writeLines(c("{", paste(lines, collapse = ",\n"), "}"), cfg_path)
  msg("config written: ", cfg_path)
}

# ---------- 5. Packages ------------------------------------------------------
suppressPackageStartupMessages({
  library(WGCNA)
  library(fastcluster)
})
enableWGCNAThreads()
options(stringsAsFactors = FALSE)
msg("WGCNA ", as.character(packageVersion("WGCNA")), " loaded; threads enabled")

# ---------- 6. Data ----------------------------------------------------------
source("network/00_load.R")
d <- load_project()

ref_net      <- readRDS(REF_CACHE)
moduleColors <- ref_net$moduleColors
modules      <- sort(unique(moduleColors[moduleColors != "grey"]))
genes        <- rownames(d$logtpm)
msg(length(modules), " non-grey modules; ", length(genes), " genes")

# ---------- 7. Load hub-036 kME tables ---------------------------------------
msg("loading hub-036 kME tables ...")
kme_ctrl   <- read.csv(HUB036_CTRL,   row.names = 1, check.names = FALSE,
                        stringsAsFactors = FALSE)
kme_nolith <- read.csv(HUB036_NOLITH, row.names = 1, check.names = FALSE,
                        stringsAsFactors = FALSE)
kme_lith   <- read.csv(HUB036_LITH,   row.names = 1, check.names = FALSE,
                        stringsAsFactors = FALSE)
msg("kME tables loaded: ", nrow(kme_ctrl), " genes x ", ncol(kme_ctrl), " modules")

# Ensure row alignment: genes must be in the same order across tables.
stopifnot(identical(rownames(kme_ctrl), rownames(kme_nolith)))
stopifnot(identical(rownames(kme_ctrl), rownames(kme_lith)))

# ---------- 8. Compute observed delta_kME ------------------------------------
# Per gene, per module: delta_kME = kME_lith - kME_ctrl (and nolith - ctrl)
delta_lith   <- kme_lith   - kme_ctrl
delta_nolith <- kme_nolith - kme_ctrl

# Intramodular kME per gene: gene's kME in its OWN module
gene_module   <- moduleColors[match(genes, rownames(d$logtpm))]
names(gene_module) <- genes

# For each gene, pull kME from its own module column in each group
pull_own_kme <- function(kme_mat, gene_mod_vec) {
  sapply(seq_along(gene_mod_vec), function(i) {
    gn  <- names(gene_mod_vec)[i]
    mod <- gene_mod_vec[i]
    if (is.na(mod) || mod == "grey" || !gn %in% rownames(kme_mat) ||
        !mod %in% colnames(kme_mat)) return(NA_real_)
    kme_mat[gn, mod]
  })
}

msg("extracting intramodular kME per gene ...")
kme_own_ctrl   <- pull_own_kme(kme_ctrl,   gene_module)
kme_own_nolith <- pull_own_kme(kme_nolith, gene_module)
kme_own_lith   <- pull_own_kme(kme_lith,   gene_module)

delta_own_lith   <- kme_own_lith   - kme_own_ctrl
delta_own_nolith <- kme_own_nolith - kme_own_ctrl

# Wide-format delta table: one row per gene per module
delta_kme_rows <- do.call(rbind, lapply(modules, function(m) {
  g_in_mod <- names(gene_module)[gene_module == m]
  g_in_mod <- g_in_mod[!is.na(g_in_mod)]
  if (!length(g_in_mod) || !m %in% colnames(delta_lith)) return(NULL)
  data.frame(
    gene                     = g_in_mod,
    module                   = m,
    kME_ctrl                 = kme_ctrl[g_in_mod, m],
    kME_nolith               = kme_nolith[g_in_mod, m],
    kME_lith_mean            = kme_lith[g_in_mod, m],
    delta_kME_ctrl_vs_nolith = delta_nolith[g_in_mod, m],
    delta_kME_ctrl_vs_lith   = delta_lith[g_in_mod, m],
    stringsAsFactors         = FALSE
  )
}))

write.csv(delta_kme_rows, file.path(run_dir, "delta_kme.csv"), row.names = FALSE)
msg("delta_kme.csv written: ", nrow(delta_kme_rows), " gene-module rows")

# ---------- 9. kME helper for permuted draws ----------------------------------
refExpr <- expr_for(d, group = "control")   # 234 x 12368 (needed for MEs)
all_sids <- colnames(d$logtpm)

# Compute kME for an arbitrary expr matrix using the reference module eigengenes.
# Returns named vector of per-gene intramodular kME (own module column only).
compute_own_kme_perm <- function(E_test) {
  # E_test: samples x genes
  MEs_test <- moduleEigengenes(E_test, colors = moduleColors,
                               excludeGrey = TRUE)$eigengenes
  kme_test <- signedKME(E_test, MEs_test)
  colnames(kme_test) <- sub("^kME", "", colnames(kme_test))
  # Per gene: own module kME
  pull_own_kme(kme_test, gene_module)
}

# ---------- 10. Permuted draws (K=20 null distribution) ----------------------
for (k in seq_len(K)) {
  perm_ckpt <- file.path(run_dir, sprintf("perm_%03d.rds", k))
  if (file.exists(perm_ckpt)) {
    msg("permuted draw ", k, "/", K, " checkpoint exists, skipping")
    next
  }
  msg("permuted draw ", k, "/", K, ": sampling 74 from all 474 ...")
  set.seed(SEED_BASE + k)
  sel_sids <- sample(all_sids, SUBSAMPLE_N, replace = FALSE)
  E_perm   <- t(d$logtpm[, sel_sids, drop = FALSE])   # 74 x 12368
  kme_own_perm <- tryCatch(
    compute_own_kme_perm(E_perm),
    error = function(e) {
      msg("ERROR in permuted draw ", k, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(kme_own_perm)) {
    rm(E_perm); gc()
    next
  }
  delta_perm <- kme_own_perm - kme_own_ctrl
  saveRDS(list(delta_own = delta_perm, samples = sel_sids, k = k), perm_ckpt)
  msg("permuted draw ", k, "/", K, " done")
  rm(E_perm, kme_own_perm, delta_perm); gc()
}

done_k <- vapply(seq_len(K),
  function(i) file.exists(file.path(run_dir, sprintf("perm_%03d.rds", i))),
  logical(1))
if (sum(done_k) < K) {
  msg("not all permuted draws complete: ", sum(done_k), "/", K,
      "; exiting for resume in a later session")
  quit(save = "no", status = 1L)
}

# ---------- 11. Aggregate null distribution and compute empirical p-values ----
msg("aggregating null distribution across ", K, " permuted draws ...")
# null_delta_mat: genes x K matrix of delta_kME under null
null_list <- lapply(seq_len(K), function(i) {
  x <- readRDS(file.path(run_dir, sprintf("perm_%03d.rds", i)))
  x$delta_own   # named numeric vector, length = n_genes
})
null_mat <- do.call(cbind, null_list)   # genes x K
rownames(null_mat) <- genes

# For each gene in a non-grey module: compute empirical p-value
# p = fraction of K null |delta| >= observed |delta|
msg("computing empirical p-values and BH correction per module ...")
diff_conn_rows <- do.call(rbind, lapply(modules, function(m) {
  g_in_mod <- names(gene_module)[gene_module == m]
  g_in_mod <- g_in_mod[!is.na(g_in_mod)]
  if (!length(g_in_mod)) return(NULL)

  obs_delta  <- delta_own_lith[g_in_mod]
  null_g     <- null_mat[g_in_mod, , drop = FALSE]   # genes x K

  emp_pval <- vapply(seq_along(g_in_mod), function(i) {
    obs_abs <- abs(obs_delta[i])
    if (is.na(obs_abs)) return(NA_real_)
    null_abs <- abs(null_g[i, ])
    mean(null_abs >= obs_abs, na.rm = TRUE)
  }, numeric(1))

  padj <- p.adjust(emp_pval, method = "BH")

  data.frame(
    gene           = g_in_mod,
    module         = m,
    kME_ctrl       = kme_own_ctrl[g_in_mod],
    kME_lith       = kme_own_lith[g_in_mod],
    delta_kME      = obs_delta,
    abs_delta_kME  = abs(obs_delta),
    direction      = ifelse(obs_delta > 0, "more_connected_in_lith",
                            ifelse(obs_delta < 0, "less_connected_in_lith", "unchanged")),
    empirical_pval = emp_pval,
    padj           = padj,
    stringsAsFactors = FALSE
  )
}))

write.csv(diff_conn_rows, file.path(run_dir, "diff_connectivity.csv"), row.names = FALSE)
msg("diff_connectivity.csv written: ", nrow(diff_conn_rows), " genes")

n_sig <- sum(diff_conn_rows$padj < 0.05, na.rm = TRUE)
n_more <- sum(diff_conn_rows$direction == "more_connected_in_lith" &
              !is.na(diff_conn_rows$padj) & diff_conn_rows$padj < 0.05, na.rm = TRUE)
n_less <- sum(diff_conn_rows$direction == "less_connected_in_lith" &
              !is.na(diff_conn_rows$padj) & diff_conn_rows$padj < 0.05, na.rm = TRUE)
msg(n_sig, " genes with differential connectivity (BH 5%): ",
    n_more, " more connected in lith; ", n_less, " less connected")
msg("Note: minimum empirical p = 1/K = ", round(1 / K, 3), " with K=", K,
    " permuted draws")
msg("=== ", id, " COMPLETE ===")
