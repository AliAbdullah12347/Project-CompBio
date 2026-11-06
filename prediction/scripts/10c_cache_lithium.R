#!/usr/bin/env Rscript
# Cache for the lithium contest: 226 BP1 subjects, 152 users vs 74 non-users.
#
# Same structure as the bipolar cache so 11_run_bpd_split.R can run either
# without changes -- `is_case` carries lithium status here rather than
# diagnosis. Keeping one runner for both contrasts is the point: the whole
# comparison is that the same pipeline, the same models and the same metric
# are applied to each, so any difference in the answer is about the contrast
# and not about how it was analysed.
#
# This contrast is far cleaner than the case-control one. Every BP1 subject is
# in assessment group A, so that confound is constant; sequencing plate is
# balanced (p = 0.887) where it is near-deterministic for case-control
# (p = 1.6e-09); tobacco is balanced (p = 0.70) where it is not (p = 4.7e-04);
# and no BP1 subject is missing a covariate. Age is the one shared confounder
# (p = 0.00073).
#
# Column subsets and row drops only, so nothing cross-sample is baked in.

args <- commandArgs(trailingOnly = TRUE)
impl <- if (length(args)) args[1] else ".."
outd <- file.path(impl, "prediction/cache")
dir.create(outd, showWarnings = FALSE, recursive = TRUE)

meta <- read.csv(file.path(impl, "data/cohort_474/metadata_474.csv"),
                 check.names = FALSE, stringsAsFactors = FALSE)
meta$li <- as.integer(meta[["lithium use (non-user=0, user = 1)"]])
b <- meta[meta[["bipolar disorder diagnosis"]] == "BP1", ]
stopifnot(nrow(b) == 226, sum(b$li == 1) == 152, sum(b$li == 0) == 74)

m <- data.frame(
  title   = b$title,
  is_case = b$li,                       # lithium user = 1
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
stopifnot(!any(is.na(m)), length(unique(m$group)) == 1)   # group A only, by construction

cnt <- read.delim(gzfile(file.path(impl, "data/cohort_474/counts_474.tsv.gz")),
                  row.names = 1, check.names = FALSE)
stopifnot(nrow(cnt) == 57820)
cnt <- cnt[!grepl("^ENSGR", rownames(cnt)), , drop = FALSE]   # 47 PAR-Y duplicates
stopifnot(nrow(cnt) == 57773)
cnt <- as.matrix(cnt[, m$title, drop = FALSE])
storage.mode(cnt) <- "integer"
stopifnot(ncol(cnt) == 226, identical(colnames(cnt), m$title))

saveRDS(list(counts = cnt, meta = m), file.path(outd, "lithium_inputs.rds"), compress = "gzip")
cat(sprintf("cached %d genes x %d samples (%d lithium users, %d non-users)\n",
            nrow(cnt), ncol(cnt), sum(m$is_case), sum(1 - m$is_case)))
cat(sprintf("  assessment group: all %s (constant -- no confound here)\n", unique(m$group)))
cat(sprintf("  -> %s\n", file.path(outd, "lithium_inputs.rds")))
