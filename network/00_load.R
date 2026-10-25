#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Canonical data loader for the network-analysis arm.
#
# Everything the arm needs is reassembled from the CSVs in
# data/analysis_matrices/, which ARE tracked in git. Nothing here depends on
# de_analysis/data/prep.rds, which is NOT tracked (74 MB binary cache) — so a
# fresh clone can load the data with no rebuild step and no GEO download.
#
#   source("network/00_load.R")
#   d <- load_project()        # cached after the first call
#
# Returns a list:
#   $logtpm   12368 x 474 numeric, log2(TPM + 1), genes x samples
#   $counts   12368 x 474 integer, filtered raw counts
#   $meta     474 rows of covariates, in the same order as the matrix columns
#   $groups   factor over the 474: "control" / "bp_nolith" / "bp_lith" / NA
#   $ilr      474 x 4 ILR balances (b1..b4)
#   $frac     474 x 5 lineage fractions (gran, mono, T, NK, B)
#   $tmm      474 x 3 library size and TMM factors
#   $tobacco  474 x 20 multiply-imputed tobacco status
#   $sym      named character vector, Ensembl gene id -> HGNC symbol
#   $lm22     character vector of LM22 signature gene symbols (circularity check)
#
# Every load is checked against the verified facts in CLAUDE.md section 2.
# If a check fails the function stops rather than returning a silently wrong
# object — an unattended run must not analyse a mangled matrix for two days.
# ---------------------------------------------------------------------------

.proj_root <- function() {
  # Walk up until we find the repo markers, so this works from any wd.
  d <- normalizePath(getwd(), winslash = "/")
  for (i in 1:6) {
    if (all(file.exists(file.path(d, c("data/analysis_matrices", "de_v2"))))) return(d)
    d <- dirname(d)
  }
  stop("Cannot locate the repository root. Run from inside the Implementation/ checkout.")
}

.read_split <- function(root, stem) {
  # Matrices over 14 MB were split into .partN.csv on export; glue them back.
  dir <- file.path(root, "data/analysis_matrices")
  parts <- sort(list.files(dir, pattern = paste0("^", stem, "\\.part[0-9]+\\.csv$"), full.names = TRUE))
  if (!length(parts)) {
    one <- file.path(dir, paste0(stem, ".csv"))
    if (!file.exists(one)) stop("Missing matrix: ", stem, " (looked in ", dir, ")")
    parts <- one
  }
  chunks <- lapply(parts, function(p) read.csv(p, check.names = FALSE, stringsAsFactors = FALSE))
  df <- do.call(rbind, chunks)
  m <- as.matrix(df[, -1, drop = FALSE])
  rownames(m) <- df[[1]]
  storage.mode(m) <- "double"
  m
}

load_project <- function(cache = TRUE, verbose = TRUE) {
  root <- .proj_root()
  cf <- file.path(root, "network", "cache", "project.rds")
  if (cache && file.exists(cf)) {
    if (verbose) message("load_project: reusing cache ", cf)
    return(readRDS(cf))
  }
  say <- function(...) if (verbose) message(...)

  say("load_project: reading matrices from data/analysis_matrices/ ...")
  logtpm <- .read_split(root, "log2TPM_12368x474")
  counts <- .read_split(root, "counts_filtered_12368x474")

  rd <- function(f) read.csv(file.path(root, "data/analysis_matrices", f),
                             check.names = FALSE, stringsAsFactors = FALSE)
  meta <- rd("covariates_474.csv")
  mem  <- rd("contrast_membership.csv")
  ilr  <- rd("ILR_coordinates.csv")
  frac <- rd("fractions_5lineages.csv")
  tmm  <- rd("TMM_normalisation_factors.csv")
  tob  <- rd("tobacco_20_imputations.csv")

  # --- align everything to the column order of the expression matrix --------
  sid <- colnames(logtpm)
  ord <- function(df, key) {
    i <- match(sid, df[[key]])
    if (anyNA(i)) stop("Sample ids do not line up for one of the covariate tables.")
    df[i, , drop = FALSE]
  }
  meta <- ord(meta, "title"); mem <- ord(mem, "sample"); ilr <- ord(ilr, "sample")
  frac <- ord(frac, "sample"); tmm <- ord(tmm, "sample"); tob <- ord(tob, "sample")
  rownames(meta) <- rownames(mem) <- rownames(ilr) <- NULL
  rownames(frac) <- rownames(tmm) <- rownames(tob) <- NULL

  # --- the three network groups --------------------------------------------
  # LI  contrast: case = BP1 on lithium,     ctrl = BP1 off lithium
  # BPD contrast: case = BP1 off lithium,    ctrl = healthy control
  g <- rep(NA_character_, length(sid))
  g[mem$BPD == "ctrl"] <- "control"
  g[mem$LI  == "ctrl"] <- "bp_nolith"
  g[mem$LI  == "case"] <- "bp_lith"
  groups <- factor(g, levels = c("control", "bp_nolith", "bp_lith"))

  # --- gene symbols ---------------------------------------------------------
  s2 <- file.path(root, "de_analysis/data/ens2sym.tsv")
  sym <- character(0)
  if (file.exists(s2)) {
    m <- read.delim(s2, header = FALSE, col.names = c("ens", "sym"),
                    stringsAsFactors = FALSE)
    sym <- setNames(m$sym, m$ens)
  } else {
    warning("ens2sym.tsv not found; module genes will be Ensembl ids only.")
  }

  # --- LM22 signature genes, for the circularity sensitivity check ----------
  lm <- file.path(root, "data/LM22.txt")
  lm22 <- if (file.exists(lm)) read.delim(lm, check.names = FALSE)[[1]] else character(0)

  d <- list(logtpm = logtpm, counts = counts, meta = meta, groups = groups,
            ilr = ilr, frac = frac, tmm = tmm, tobacco = tob, sym = sym,
            lm22 = lm22, root = root)

  # --- verified facts: CLAUDE.md section 2. Fail loudly, never silently. ----
  chk <- function(cond, msg) if (!isTRUE(cond)) stop("DATA CHECK FAILED: ", msg)
  chk(nrow(logtpm) == 12368, sprintf("expected 12,368 genes, got %d", nrow(logtpm)))
  chk(ncol(logtpm) == 474,   sprintf("expected 474 samples, got %d", ncol(logtpm)))
  chk(identical(dim(counts), dim(logtpm)), "counts and log2TPM disagree in shape")
  chk(identical(rownames(counts), rownames(logtpm)), "counts and log2TPM gene order differs")
  chk(identical(colnames(counts), colnames(logtpm)), "counts and log2TPM sample order differs")
  chk(sum(groups == "control",   na.rm = TRUE) == 234, "expected 234 healthy controls")
  chk(sum(groups == "bp_nolith", na.rm = TRUE) ==  74, "expected 74 BP1 off lithium")
  chk(sum(groups == "bp_lith",   na.rm = TRUE) == 152, "expected 152 BP1 on lithium")
  chk(!any(duplicated(rownames(logtpm))), "duplicate gene ids")
  chk(!anyNA(logtpm), "NA in the expression matrix")
  chk(all(grepl("^ENSG[0-9]+$", rownames(logtpm))), "gene ids are not unversioned Ensembl")
  chk(nrow(meta) == 474 && identical(meta$title, sid), "covariates are not aligned to the matrix")

  say(sprintf("load_project: OK  %d genes x %d samples | control %d, bp_nolith %d, bp_lith %d",
              nrow(logtpm), ncol(logtpm),
              sum(groups == "control", na.rm = TRUE),
              sum(groups == "bp_nolith", na.rm = TRUE),
              sum(groups == "bp_lith", na.rm = TRUE)))
  say(sprintf("load_project: symbols for %d/%d genes | LM22 signature genes: %d",
              sum(rownames(logtpm) %in% names(sym)), nrow(logtpm), length(lm22)))

  if (cache) {
    dir.create(dirname(cf), showWarnings = FALSE, recursive = TRUE)
    saveRDS(d, cf)
    say("load_project: cached to ", cf)
  }
  d
}

# Expression matrix as samples x genes, which is the orientation WGCNA wants.
# Optionally restricted to one group and/or residualised for covariates.
# --- normalisation mode -----------------------------------------------------
# Defined in its own file so experiment scripts can source it cheaply before
# they build CACHE_KEY. See network/00_exprmode.R for the rationale.
source("network/00_exprmode.R")

# Subset by explicit sample IDs with the SAME normalisation expr_for() applies.
# Exists because 42 of 51 experiment scripts reached for `t(d$logtpm[, sids])`
# directly to build a test set, which silently bypassed offset removal and
# scored a raw test matrix against an offset-removed reference. Always use this
# (or expr_for) - never d$logtpm - for anything that becomes expression input.
# Verify a matrix REALLY carries the normalisation it is claimed to.
#
# A cache key records intent, not content. On 2026-10-25 a reference network
# built from raw log2(TPM+1) sat for a week in a file whose key ended
# "_offsetrm", because the key names the mode the script meant to use, not the
# data it read. This reads the artifact instead of trusting its label.
#
# It reconstructs a random sample of the matrix via expr_samples() under the
# declared mode and compares. A gene-subsetted matrix cannot be checked by
# "per-sample mean is 0" - the offset is removed over all 12,368 genes, so any
# subset has a non-zero mean - which is why this compares values, not means.
# Gate a reference cache on the way IN. assert_expr_mode() stops a bad reference
# being built; this stops an already-bad one being reused.
#
# It returns FALSE rather than raising, and moves the offending file aside, so
# the caller's existing "if the cache is not usable, build it" branch simply
# rebuilds. Raising here would wedge the row instead: it would fail every
# attempt and reach BLOCKED without ever regenerating the thing it objected to.
#
# A cache with no expr_mode stamp predates 2026-10-25 and cannot be verified
# from its contents, so it is refused rather than trusted - that is precisely
# the window in which a raw network sat under an "_offsetrm" key.
ref_cache_ok <- function(path, mode = EXPR_MODE, quiet = FALSE) {
  if (!file.exists(path)) return(FALSE)
  r <- tryCatch(readRDS(path), error = function(e) NULL)
  bad <- if (is.null(r)) "unreadable"
         else if (is.null(r$expr_mode)) "no expr_mode stamp (pre-2026-10-25, unverifiable)"
         else if (!identical(r$expr_mode, mode)) paste0("built with expr_mode='", r$expr_mode,
                                                        "', this run uses '", mode, "'")
         else NA_character_
  if (is.na(bad)) return(TRUE)
  q <- file.path(dirname(path), "quarantine")
  dir.create(q, showWarnings = FALSE, recursive = TRUE)
  dest <- file.path(q, paste0(basename(path), ".", format(Sys.time(), "%Y%m%dT%H%M%SZ"), ".rds"))
  file.rename(path, dest)
  if (!quiet) message("ref_cache_ok: QUARANTINED ", basename(path), " - ", bad,
                      "\n  moved to ", dest, "; it will be rebuilt.")
  FALSE
}

assert_expr_mode <- function(M, d, mode = EXPR_MODE, what = "expression matrix",
                             n_probe = 400, tol = 1e-8) {
  if (is.null(dim(M))) stop("assert_expr_mode: ", what, " is not a matrix")
  # Orientation: these matrices are samples x genes (expr_for default) or
  # genes x samples. Decide by which margin's names match the project's samples.
  sids_all <- colnames(d$logtpm)
  by_row <- sum(rownames(M) %in% sids_all); by_col <- sum(colnames(M) %in% sids_all)
  if (max(by_row, by_col) == 0)
    stop("assert_expr_mode: ", what, " has no dimnames matching project samples")
  samples_in_rows <- by_row >= by_col
  sids <- if (samples_in_rows) rownames(M) else colnames(M)
  gids <- if (samples_in_rows) colnames(M) else rownames(M)
  gids <- intersect(gids, rownames(d$logtpm))
  if (!length(gids)) stop("assert_expr_mode: ", what, " shares no genes with the project")
  probe <- if (length(gids) > n_probe) sample(gids, n_probe) else gids

  expect <- expr_samples(d, sids, transpose = FALSE, mode = mode)[probe, , drop = FALSE]
  got <- if (samples_in_rows) t(M[, probe, drop = FALSE]) else M[probe, , drop = FALSE]
  worst <- max(abs(got - expect))

  if (!is.finite(worst) || worst > tol) {
    alt <- if (identical(mode, "raw")) "offset_removed" else "raw"
    expect2 <- expr_samples(d, sids, transpose = FALSE, mode = alt)[probe, , drop = FALSE]
    looks_like <- if (max(abs(got - expect2)) <= tol) alt else "neither"
    stop(sprintf(
      "assert_expr_mode: %s does NOT match mode '%s' (max abs deviation %.3g). It looks like: %s.\n  This is the 2026-10-25 failure mode - see network/lint_expr.sh and\n  findings/2026-10-25-raw-test-matrices-bypassed-normalisation.md",
      what, mode, worst, looks_like))
  }
  invisible(list(mode = mode, n_samples = length(sids), n_genes = length(gids),
                 probe = length(probe), max_dev = worst,
                 per_sample_offset_sd = stats::sd(colMeans(
                   expr_samples(d, sids, transpose = FALSE, mode = mode)))))
}

expr_samples <- function(d, sids, transpose = TRUE, mode = EXPR_MODE) {
  miss <- setdiff(sids, colnames(d$logtpm))
  if (length(miss)) stop("expr_samples: unknown sample ids: ", paste(head(miss, 3), collapse = ", "))
  M <- d$logtpm[, sids, drop = FALSE]
  if (!identical(mode, "raw")) {
    M <- sweep(M, 2, colMeans(d$logtpm)[sids], "-")
  }
  if (transpose) t(M) else M
}

expr_for <- function(d, group = NULL, transpose = TRUE, residualise = NULL,
                     mode = EXPR_MODE) {
  keep <- if (is.null(group)) rep(TRUE, ncol(d$logtpm)) else
    !is.na(d$groups) & d$groups %in% group
  M <- d$logtpm[, keep, drop = FALSE]
  # The offset is computed over ALL 12,368 genes and the FULL sample set, then
  # subset - so every experiment removes the identical quantity regardless of
  # which genes or which group it happens to be looking at.
  if (!identical(mode, "raw")) {
    M <- sweep(M, 2, colMeans(d$logtpm)[keep], "-")
  }
  if (!is.null(residualise)) {
    md <- d$meta[keep, , drop = FALSE]
    if (!is.null(d$ilr)) md <- cbind(md, d$ilr[keep, -1, drop = FALSE])
    X <- model.matrix(residualise, data = md)
    if (nrow(X) != ncol(M)) stop("residualise: model frame dropped rows (missing covariates?)")
    M <- t(residuals(lm.fit(X, t(M))))
  }
  if (transpose) t(M) else M
}

if (sys.nframe() == 0L) {
  d <- load_project(cache = FALSE)
  cat("\nself-test passed\n")
}
