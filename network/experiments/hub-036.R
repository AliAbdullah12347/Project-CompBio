#!/usr/bin/env Rscript
# ===========================================================================
# hub-036  kME hub identity and stability across draws
#
# Spec (from queue):
#   kME hub identity; stability across 100 draws; do hubs change between groups?
#
# Design (draws revised to 20 — same budget revision as base-001):
#   1. For EACH GROUP (control, bp_nolith, bp_lith):
#        - control:   kME from all 234 control samples (full group, no subsampling)
#        - bp_nolith: kME from all 74 bp_nolith samples (full group, no subsampling)
#        - bp_lith:   20 draws of 74 samples (seeds: SEED_BASE_DRAW + draw);
#                     kME averaged and SD computed across draws
#   2. For each module: rank genes by kME within each group. Hub = top-kME gene.
#   3. Stability metric: Spearman correlation of kME ranks between control and
#      each test group (bp_nolith, bp_lith).
#   4. Check: does the top hub gene change between control and bp_lith?
#   5. No TOM computation — just eigengenes + kME (FAST).
#
# Non-negotiable design points obeyed:
#   - Draws implemented = 20 (spec = 100; documented in config)
#   - SEED_BASE = 20261009L + 36000L
#   - Config written before any computation
#   - Checkpoint per draw; skip if draw file exists
#   - Dependency check: REF_CACHE must exist; quit(status=1) if missing
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "hub-036"
stopifnot(id == "hub-036")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

# ---------- 1. Directories ---------------------------------------------------
run_dir <- file.path("network/runs", id)
cfg_dir <- "network/config"
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg_dir, showWarnings = FALSE, recursive = TRUE)

# ---------- 2. Parameters — ALL pre-specified BEFORE any computation ---------
SEED_BASE    <- 20261009L + 36000L
DRAWS        <- 20L        # spec = 100; revised for budget (see METHODS.md 2026-10-09)
SUBSAMPLE_N  <- 74L
REF_CACHE    <- "network/cache/ref_control_all_bicor_signed_p12.rds"

# ---------- 3. Dependency check ----------------------------------------------
if (!file.exists(REF_CACHE)) {
  msg("DEPENDENCY: base-001 cache missing: ", REF_CACHE)
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
    spec                = "kME hub identity; stability across 100 draws; do hubs change between groups?",
    draws_spec          = 100L,
    draws_implemented   = DRAWS,
    draws_revision_note = "100 draws x compute time infeasible; revised to 20 (same budget constraint as base-001, see METHODS.md 2026-10-09)",
    ref_cache           = REF_CACHE,
    groups              = c("control", "bp_nolith", "bp_lith"),
    subsample_n_lith    = SUBSAMPLE_N,
    subsample_lith      = TRUE,
    subsample_nolith    = FALSE,
    subsample_control   = FALSE,
    kme_method          = "signedKME",
    hub_definition      = "top kME gene per module per group",
    stability_metric    = "Spearman rho of kME ranks between control and test group",
    tom_computed        = FALSE,
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

# ---------- 7. Helper: compute signedKME for a samples-x-genes matrix --------
# Returns a data.frame: genes as rows, modules as columns.
# Only columns for non-grey modules are returned.
compute_kme <- function(E) {
  # E: samples x genes (WGCNA orientation)
  MEs <- moduleEigengenes(E, colors = moduleColors, excludeGrey = TRUE)$eigengenes
  kme <- signedKME(E, MEs)
  # signedKME returns a data.frame; column names are like "kMEblue"
  # standardise column names to module color only
  colnames(kme) <- sub("^kME", "", colnames(kme))
  kme
}

# ---------- 8. control kME ---------------------------------------------------
ctrl_ckpt <- file.path(run_dir, "kme_control.csv")
if (!file.exists(ctrl_ckpt)) {
  msg("computing kME for control ...")
  E_ctrl <- expr_for(d, group = "control")   # 234 x 12368
  kme_ctrl <- compute_kme(E_ctrl)
  write.csv(kme_ctrl, ctrl_ckpt, row.names = TRUE)
  msg("kme_control.csv written: ", nrow(kme_ctrl), " genes x ", ncol(kme_ctrl), " modules")
  rm(E_ctrl); gc()
} else {
  msg("kme_control.csv exists; loading from checkpoint")
  kme_ctrl <- read.csv(ctrl_ckpt, row.names = 1, check.names = FALSE,
                       stringsAsFactors = FALSE)
}

# ---------- 9. bp_nolith kME -------------------------------------------------
nolith_ckpt <- file.path(run_dir, "kme_nolith.csv")
if (!file.exists(nolith_ckpt)) {
  msg("computing kME for bp_nolith ...")
  E_nolith <- expr_for(d, group = "bp_nolith")   # 74 x 12368
  kme_nolith <- compute_kme(E_nolith)
  write.csv(kme_nolith, nolith_ckpt, row.names = TRUE)
  msg("kme_nolith.csv written: ", nrow(kme_nolith), " genes x ", ncol(kme_nolith), " modules")
  rm(E_nolith); gc()
} else {
  msg("kme_nolith.csv exists; loading from checkpoint")
  kme_nolith <- read.csv(nolith_ckpt, row.names = 1, check.names = FALSE,
                         stringsAsFactors = FALSE)
}

# ---------- 10. bp_lith draws (20 draws x signedKME) ------------------------
lith_sids <- colnames(d$logtpm)[!is.na(d$groups) & d$groups == "bp_lith"]
stopifnot(length(lith_sids) == 152L)

for (draw in seq_len(DRAWS)) {
  draw_ckpt <- file.path(run_dir, sprintf("draw_%03d.rds", draw))
  if (file.exists(draw_ckpt)) {
    msg("bp_lith draw ", draw, "/", DRAWS, " checkpoint exists, skipping")
    next
  }
  msg("bp_lith draw ", draw, "/", DRAWS, ": subsampling 74 from 152 ...")
  set.seed(SEED_BASE + draw)
  sel_sids <- sample(lith_sids, SUBSAMPLE_N, replace = FALSE)
  E_draw   <- t(d$logtpm[, sel_sids, drop = FALSE])   # 74 x 12368
  kme_draw <- tryCatch(
    compute_kme(E_draw),
    error = function(e) {
      msg("ERROR in draw ", draw, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(kme_draw)) {
    rm(E_draw); gc()
    next
  }
  saveRDS(list(kme = kme_draw, samples = sel_sids, draw = draw), draw_ckpt)
  msg("draw ", draw, "/", DRAWS, " done (", nrow(kme_draw), " genes x ",
      ncol(kme_draw), " modules)")
  rm(E_draw, kme_draw); gc()
}

# ---------- 11. Check all draws complete before aggregating ------------------
done_lith <- vapply(seq_len(DRAWS),
  function(i) file.exists(file.path(run_dir, sprintf("draw_%03d.rds", i))),
  logical(1))
n_lith_done <- sum(done_lith)
if (n_lith_done < DRAWS) {
  msg("not all bp_lith draws complete: ", n_lith_done, "/", DRAWS,
      "; exiting for resume in a later session")
  quit(save = "no", status = 1L)
}

# ---------- 12. Aggregate bp_lith draws → mean and SD kME matrices -----------
msg("all draws complete; aggregating bp_lith kME across ", DRAWS, " draws ...")
draw_kme_list <- lapply(seq_len(DRAWS), function(i) {
  x <- readRDS(file.path(run_dir, sprintf("draw_%03d.rds", i)))
  as.matrix(x$kme)
})
# Stack into array: genes x modules x draws
n_genes   <- nrow(draw_kme_list[[1]])
n_modules <- ncol(draw_kme_list[[1]])
kme_array <- array(
  unlist(draw_kme_list),
  dim      = c(n_genes, n_modules, DRAWS),
  dimnames = list(rownames(draw_kme_list[[1]]),
                  colnames(draw_kme_list[[1]]),
                  NULL)
)
kme_lith_mean <- apply(kme_array, c(1, 2), mean, na.rm = TRUE)
kme_lith_sd   <- apply(kme_array, c(1, 2), sd,   na.rm = TRUE)
kme_lith_mean <- as.data.frame(kme_lith_mean)
kme_lith_sd   <- as.data.frame(kme_lith_sd)

write.csv(kme_lith_mean, file.path(run_dir, "kme_lith_mean.csv"), row.names = TRUE)
write.csv(kme_lith_sd,   file.path(run_dir, "kme_lith_sd.csv"),   row.names = TRUE)
msg("kme_lith_mean.csv and kme_lith_sd.csv written")

# ---------- 13. Hub identity per module per group ----------------------------
# Hub = gene with highest kME in its own module column.
find_hubs <- function(kme_mat, group_label) {
  do.call(rbind, lapply(modules, function(m) {
    if (!m %in% colnames(kme_mat)) {
      return(data.frame(module = m, hub_gene = NA_character_,
                        hub_kME = NA_real_, group = group_label,
                        stringsAsFactors = FALSE))
    }
    v      <- kme_mat[[m]]
    idx    <- which.max(v)
    hub_gn <- rownames(kme_mat)[idx]
    data.frame(module = m, hub_gene = hub_gn, hub_kME = v[idx],
               group = group_label, stringsAsFactors = FALSE)
  }))
}

hubs_ctrl   <- find_hubs(kme_ctrl,      "control")
hubs_nolith <- find_hubs(kme_nolith,    "bp_nolith")
hubs_lith   <- find_hubs(kme_lith_mean, "bp_lith")

hub_identity           <- rbind(hubs_ctrl, hubs_nolith, hubs_lith)
hub_identity$is_ref_hub <- hub_identity$hub_gene ==
  hub_identity$hub_gene[match(hub_identity$module, hubs_ctrl$module)]

write.csv(hub_identity, file.path(run_dir, "hub_identity.csv"), row.names = FALSE)
msg("hub_identity.csv written")

# ---------- 14. Stability: Spearman rho of kME ranks -------------------------
# Per module: rank genes by kME in control vs test group.
# Spearman rho between the two rank vectors.
spearman_kme_rank <- function(kme_ref, kme_test, mod) {
  if (!mod %in% colnames(kme_ref) || !mod %in% colnames(kme_test)) return(NA_real_)
  v_ref  <- kme_ref[[mod]]
  v_test <- kme_test[[mod]]
  if (length(v_ref) != length(v_test)) return(NA_real_)
  cor(v_ref, v_test, method = "spearman", use = "pairwise.complete.obs")
}

hub_stability <- do.call(rbind, lapply(modules, function(m) {
  rho_nolith <- spearman_kme_rank(kme_ctrl, kme_nolith,    m)
  rho_lith   <- spearman_kme_rank(kme_ctrl, kme_lith_mean, m)
  # Hub change check: does the top hub gene differ?
  hub_ctrl_gene   <- hubs_ctrl$hub_gene[hubs_ctrl$module == m]
  hub_nolith_gene <- hubs_nolith$hub_gene[hubs_nolith$module == m]
  hub_lith_gene   <- hubs_lith$hub_gene[hubs_lith$module == m]
  data.frame(
    module               = m,
    spearman_ctrl_nolith = rho_nolith,
    spearman_ctrl_lith   = rho_lith,
    hub_ctrl             = if (length(hub_ctrl_gene))   hub_ctrl_gene[1]   else NA,
    hub_nolith           = if (length(hub_nolith_gene)) hub_nolith_gene[1] else NA,
    hub_lith             = if (length(hub_lith_gene))   hub_lith_gene[1]   else NA,
    hub_changed_nolith   = !identical(hub_ctrl_gene, hub_nolith_gene),
    hub_changed_lith     = !identical(hub_ctrl_gene, hub_lith_gene),
    stringsAsFactors     = FALSE
  )
}))

write.csv(hub_stability, file.path(run_dir, "hub_stability.csv"), row.names = FALSE)
msg("hub_stability.csv written")

# ---------- 15. Summary report -----------------------------------------------
n_changed_nolith <- sum(hub_stability$hub_changed_nolith, na.rm = TRUE)
n_changed_lith   <- sum(hub_stability$hub_changed_lith,   na.rm = TRUE)
msg("Hub gene changed between control and bp_nolith: ", n_changed_nolith, "/",
    length(modules), " modules")
msg("Hub gene changed between control and bp_lith:   ", n_changed_lith, "/",
    length(modules), " modules")

rho_nolith_med <- median(hub_stability$spearman_ctrl_nolith, na.rm = TRUE)
rho_lith_med   <- median(hub_stability$spearman_ctrl_lith,   na.rm = TRUE)
msg("Median Spearman rho (ctrl vs bp_nolith): ", round(rho_nolith_med, 3))
msg("Median Spearman rho (ctrl vs bp_lith):   ", round(rho_lith_med, 3))

msg("=== ", id, " COMPLETE ===")
