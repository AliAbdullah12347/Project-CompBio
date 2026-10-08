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
expr_for <- function(d, group = NULL, transpose = TRUE, residualise = NULL) {
  keep <- if (is.null(group)) rep(TRUE, ncol(d$logtpm)) else
    !is.na(d$groups) & d$groups %in% group
  M <- d$logtpm[, keep, drop = FALSE]
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
