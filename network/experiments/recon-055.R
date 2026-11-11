#!/usr/bin/env Rscript
# ===========================================================================
# recon-055  Is it the partition? (NEXT_STEPS item 3)
#
# THE DISCREPANCY
#
# On the fixed verified 71-module reference, exactly 1 of 71 lithium-module
# associations survives adjustment for cell composition (trait-053). On
# grp-051's group-specific 69-module control_A network, 28 of 69 do. Same
# data, same subjects, wildly different answers.
#
# trait-053 already ruled out two of the three candidate explanations by
# running both models on both sample sets. What it could not rule out is the
# MODULE SET itself, because it only ever had one.
#
# THE TEST
#
# Re-score grp-051's partition through trait-053's code path, unchanged in
# every other respect: same expression matrix, same sample sets, same three
# models, same BH correction. The only thing that differs is which genes are
# grouped together. If 28 survives, a group-specific decomposition is finding
# composition-independent lithium signal the fixed reference misses, and that
# is worth pursuing. If it collapses toward 1, grp-051b's count came from
# something else in its pipeline and the fixed-reference answer stands.
#
# WHY NOT JUST EDIT trait-053
#
# trait-053 is a published row whose numbers are cited in the findings. This
# row borrows its model verbatim and swaps the partition, so the two can be
# compared directly and neither overwrites the other.
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "recon-055"
stopifnot(id == "recon-055")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

run_dir <- file.path("network/runs", id)
cfg_dir <- "network/config"
for (d_ in c(run_dir, cfg_dir)) dir.create(d_, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({ library(WGCNA) })
enableWGCNAThreads()
options(stringsAsFactors = FALSE)
source("network/00_load.R")
d <- load_project(cache = TRUE, verbose = FALSE)

## ---- the two partitions -----------------------------------------------------
E <- expr_for(d, group = NULL)                 # offset-removed, all samples
genes_E <- colnames(E)
msg("expression: ", nrow(E), " samples x ", ncol(E), " genes")

REF <- "network/cache/ref_control_all_bicor_signed_p12_ds4_mch00_offsetrm.rds"
if (!ref_cache_ok(REF)) stop("verified reference missing or unstamped: ", REF)
ref <- readRDS(REF)
cols_fixed <- ref$moduleColors
# base-001 cached its colours positionally, by design -- the preservation rows
# consume them that way. Attach ids locally, as comm-038/039 and ann-044 do.
if (is.null(names(cols_fixed))) names(cols_fixed) <- genes_E
stopifnot(length(cols_fixed) == ncol(E))

g <- read.csv("network/runs/grp-051/genes_control_A.csv")
stopifnot(all(c("gene", "module") %in% names(g)))
cols_grp <- setNames(g$module, g$gene)[genes_E]
if (anyNA(cols_grp)) stop(sum(is.na(cols_grp)), " genes in E are absent from the ",
                          "grp-051 partition; the two are not over the same gene set")
msg("fixed reference : ", length(setdiff(unique(cols_fixed), "grey")), " non-grey modules")
msg("grp-051 control_A: ", length(setdiff(unique(cols_grp), "grey")), " non-grey modules")
msg("gene sets identical: ", identical(sort(genes_E), sort(g$gene)))

## ---- trait-053's model, verbatim --------------------------------------------
md <- d$meta
md$li <- as.integer(d$groups == "bp_lith")
md$dx <- as.integer(d$groups %in% c("bp_nolith", "bp_lith"))
ilr <- d$ilr[match(md$title, d$ilr$sample), c("b1", "b2", "b3", "b4"), drop = FALSE]

sample_sets <- list(
  BP1_only = !is.na(d$groups) & d$groups %in% c("bp_lith", "bp_nolith"),
  all      = !is.na(d$groups))

fit_and_count <- function(ME, keep, model) {
  p <- vapply(colnames(ME), function(m) {
    y  <- ME[[m]][keep]
    df <- data.frame(y = y, li = md$li[keep], dx = md$dx[keep], age = md$age[keep],
                     sex = md$sex[keep], rin = md$rin[keep],
                     plate = factor(md$plate[keep]),
                     b1 = ilr$b1[keep], b2 = ilr$b2[keep],
                     b3 = ilr$b3[keep], b4 = ilr$b4[keep])
    f <- switch(model,
      marginal = lm(y ~ li, data = df),
      M_001    = lm(y ~ li + b1 + b2 + b3 + b4, data = df),
      M_051b   = if (all(df$dx == df$dx[1]))      # dx constant inside BP1: drop it
                   lm(y ~ li + age + sex + rin + plate + b1 + b2 + b3 + b4, data = df)
                 else
                   lm(y ~ li + dx + age + sex + rin + plate + b1 + b2 + b3 + b4, data = df))
    s <- summary(f)$coefficients
    if ("li" %in% rownames(s)) s["li", 4] else NA_real_
  }, numeric(1))
  list(p = p, n_bh = sum(p.adjust(p, "BH") < 0.05, na.rm = TRUE))
}

## ---- config -----------------------------------------------------------------
cfg_path <- file.path(cfg_dir, paste0(id, ".json"))
if (!file.exists(cfg_path)) {
  sha <- tryCatch(trimws(system("git rev-parse --short HEAD", intern = TRUE)[1]),
                  error = function(e) "unknown", warning = function(w) "unknown")
  writeLines(c("{",
    paste0('  "id": "', id, '",'),
    '  "family": "module-trait",',
    '  "spec": "re-score grp-051 control_A partition through trait-053 model; only the partition differs",',
    '  "models": "marginal; M_001 (li + b1..b4); M_051b (li + dx + technical + b1..b4)",',
    '  "partitions": "fixed verified reference; grp-051 control_A",',
    paste0('  "git_sha": "', if (length(sha) && !is.na(sha)) sha else "unknown", '",'),
    paste0('  "r_version": "', paste(R.version$major, R.version$minor, sep = "."), '",'),
    paste0('  "wgcna_version": "', as.character(packageVersion("WGCNA")), '"'),
    "}"), cfg_path)
  msg("config written: ", cfg_path)
}

## ---- run both partitions through it -----------------------------------------
partitions <- list(fixed_reference = cols_fixed, grp051_control_A = cols_grp)

out <- list(); per_mod <- list()
for (pn in names(partitions)) {
  msg("")
  msg("--- partition: ", pn, " ---")
  ME <- moduleEigengenes(E, colors = partitions[[pn]], excludeGrey = TRUE)$eigengenes
  msg("  eigengenes: ", nrow(ME), " x ", ncol(ME))
  for (ss in names(sample_sets)) {
    keep <- sample_sets[[ss]]
    for (mm in c("marginal", "M_001", "M_051b")) {
      r <- fit_and_count(ME, keep, mm)
      per_mod[[length(per_mod) + 1]] <- data.frame(
        partition = pn, sample_set = ss, n = sum(keep), model = mm,
        module = sub("^ME", "", names(r$p)), p = unname(r$p),
        p_bh = unname(p.adjust(r$p, "BH")))
      out[[length(out) + 1]] <- data.frame(
        partition = pn, sample_set = ss, n = sum(keep), model = mm,
        n_modules = ncol(ME), n_BH05 = r$n_bh)
      msg(sprintf("  %-9s n=%-3d  %-9s -> %d of %d modules at BH<0.05",
                  ss, sum(keep), mm, r$n_bh, ncol(ME)))
    }
  }
}

R  <- do.call(rbind, out)
PM <- do.call(rbind, per_mod)
write.csv(R,  file.path(run_dir, "reconciliation_by_partition.csv"), row.names = FALSE)
write.csv(PM, file.path(run_dir, "per_module_p.csv"), row.names = FALSE)

## ---- the comparison this row exists to make ---------------------------------
msg("")
msg("=== the same model on two partitions ===")
for (ss in unique(R$sample_set)) for (mm in c("M_001", "M_051b")) {
  a <- R$n_BH05[R$partition == "fixed_reference"  & R$sample_set == ss & R$model == mm]
  b <- R$n_BH05[R$partition == "grp051_control_A" & R$sample_set == ss & R$model == mm]
  na <- R$n_modules[R$partition == "fixed_reference"  & R$sample_set == ss & R$model == mm]
  nb <- R$n_modules[R$partition == "grp051_control_A" & R$sample_set == ss & R$model == mm]
  msg(sprintf("  %-9s %-7s: fixed %2d/%d   grp-051 %2d/%d", ss, mm, a, na, b, nb))
}

surv <- PM[PM$p_bh < 0.05 & PM$model != "marginal", ]
for (pn in unique(surv$partition)) for (ss in unique(surv$sample_set)) {
  for (mm in unique(surv$model)) {
    s <- surv[surv$partition == pn & surv$sample_set == ss & surv$model == mm, ]
    if (nrow(s)) msg(sprintf("  survivors %-17s %-9s %-7s: %s",
                             pn, ss, mm, paste(s$module, collapse = ", ")))
  }
}

msg("")
msg("Read it this way: if grp-051's partition still gives ~28 under M_051b")
msg("while the fixed reference gives 1, the difference IS the module set and")
msg("a group-specific decomposition is seeing something the fixed one misses.")
msg("If it collapses toward 1, grp-051b's 28 came from elsewhere in its")
msg("pipeline and the fixed-reference answer stands.")
msg("=== ", id, " COMPLETE ===")
