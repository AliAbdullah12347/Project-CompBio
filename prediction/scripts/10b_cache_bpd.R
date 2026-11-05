#!/usr/bin/env Rscript
# Cache the raw counts and covariates for the bipolar-status contrast.
#
# Everything here is a column subset or a row drop -- both per-sample, so
# nothing cross-sample is baked in and the cache cannot leak between folds.
# The gene filter, TMM and everything else that looks across samples happens
# inside 11_run_bpd_split.R, within each training fold.
#
# The cache covers the full 308 (74 cases + all 234 controls); the group-A
# cohort is a subset of it, so one cache serves both.

args <- commandArgs(trailingOnly = TRUE)
impl <- if (length(args)) args[1] else ".."
outd <- file.path(impl, "prediction/cache")
dir.create(outd, showWarnings = FALSE, recursive = TRUE)

meta <- read.csv(file.path(impl, "data/cohort_474/metadata_474.csv"),
                 check.names = FALSE, stringsAsFactors = FALSE)
meta$li <- as.integer(meta[["lithium use (non-user=0, user = 1)"]])
meta$dx <- meta[["bipolar disorder diagnosis"]]

keep <- (meta$dx == "BP1" & meta$li == 0) | meta$dx == "Control"
b <- meta[keep, ]
stopifnot(nrow(b) == 308, sum(b$dx == "BP1") == 74, sum(b$dx == "Control") == 234)

# 29 subjects have no recorded tobacco status and every one of them is a
# control -- 0 of 74 cases, 29 of 234 controls. So "tobacco is missing" names a
# control with certainty, and any encoding of it leaks: a missingness flag
# obviously, but also an imputed value, because a tobacco of 0.24 sitting among
# 0s and 1s is the same flag by another name and a tree splits on it instantly.
# Dropping them is the only handling a model cannot exploit. It costs 29
# controls and no cases, and shrinks the primary cohort's controls from 125 to
# 112; that is the price of a number we can defend.
miss <- is.na(b[["tobacco use"]]) | b[["tobacco use"]] == ""
cat(sprintf("dropping %d subjects with no tobacco status (%d cases, %d controls)\n",
            sum(miss), sum(miss & b$dx == "BP1"), sum(miss & b$dx == "Control")))
b <- b[!miss, ]
stopifnot(nrow(b) == 279, sum(b$dx == "BP1") == 74)

m <- data.frame(
  title   = b$title,
  is_case = as.integer(b$dx == "BP1"),
  group   = b[["assessment group"]],
  age     = as.numeric(b$age),
  sex     = b$Sex,
  tobacco = as.integer(b[["tobacco use"]]),
  rin     = as.numeric(b$rin),
  plate   = as.character(b[["sequencing plate"]]),
  seqpc1  = as.numeric(b[["sequencing metric pc1"]]),
  seqpc2  = as.numeric(b[["sequencing metric pc2"]]),
  seqpc3  = as.numeric(b[["sequencing metric pc3"]]),
  depth_ok = b$depth_ok %in% c("True", "TRUE"),
  rin_ok   = b$rin_ok   %in% c("True", "TRUE"),
  stringsAsFactors = FALSE)
stopifnot(!any(is.na(m)))

cnt <- read.delim(gzfile(file.path(impl, "data/cohort_474/counts_474.tsv.gz")),
                  row.names = 1, check.names = FALSE)
stopifnot(nrow(cnt) == 57820)
cnt <- cnt[!grepl("^ENSGR", rownames(cnt)), , drop = FALSE]   # 47 PAR-Y duplicates
stopifnot(nrow(cnt) == 57773)
cnt <- as.matrix(cnt[, m$title, drop = FALSE])
storage.mode(cnt) <- "integer"
stopifnot(ncol(cnt) == 279, identical(colnames(cnt), m$title))

saveRDS(list(counts = cnt, meta = m), file.path(outd, "bpd_inputs.rds"), compress = "gzip")
cat(sprintf("cached %d genes x %d samples (%d cases, %d controls)\n",
            nrow(cnt), ncol(cnt), sum(m$is_case), sum(1 - m$is_case)))
cat(sprintf("  group A: %d  group B: %d (controls only)\n",
            sum(m$group == "A"), sum(m$group == "B")))
cat(sprintf("  library size %.2fM - %.2fM\n", min(colSums(cnt))/1e6, max(colSums(cnt))/1e6))
cat(sprintf("  -> %s\n", file.path(outd, "bpd_inputs.rds")))
