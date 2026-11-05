#!/usr/bin/env Rscript
# Train/test partitions for the bipolar-status contrast.
#
# Two cohorts, both written, both run:
#
#   groupA  74 cases vs 125 controls   PRIMARY
#           Controls restricted to assessment group A. Group B holds 109
#           controls and zero cases, and assessment group is predictable from
#           expression at AUC 0.939 among controls -- so including group B
#           hands a classifier 47% of the controls for free. See BPD_DESIGN.md.
#
#   full    74 cases vs 234 controls   SECONDARY
#           Reported beside the primary. The gap between the two is the
#           measure of how much that confound was worth.
#
# 25 repeats of stratified 3-fold CV each, so every subject is out-of-fold
# exactly once per repeat and the pooled evaluation set is the whole cohort.
#
# Usage: Rscript 10_make_splits_bpd.R [impl_root]

args <- commandArgs(trailingOnly = TRUE)
impl <- if (length(args)) args[1] else ".."

SEED    <- 481
REPEATS <- 25
FOLDS   <- 3

meta <- read.csv(file.path(impl, "data/cohort_474/metadata_474.csv"),
                 check.names = FALSE, stringsAsFactors = FALSE)
meta$li  <- as.integer(meta[["lithium use (non-user=0, user = 1)"]])
meta$dx  <- meta[["bipolar disorder diagnosis"]]
meta$grp <- meta[["assessment group"]]

is_case <- meta$dx == "BP1" & meta$li == 0          # BP2 excluded, see design §1
tobacco_known <- !(is.na(meta[["tobacco use"]]) | meta[["tobacco use"]] == "")
is_ctrl <- meta$dx == "Control"

# Tobacco is missing for 29 subjects, all of them controls, so its absence
# identifies a control. Those subjects are excluded rather than imputed; see
# 10b_cache_bpd.R and BPD_DESIGN.md.
# Sequencing plate is near-deterministic for this contrast: plate 1 holds 15
# cases and 0 controls, plate 2 holds 15 cases and 3 controls (chi-sq
# p = 1.6e-09). Plate is a batch variable and is therefore also written into
# the expression, so a model can read diagnosis off it without any biology.
# CLAUDE.md records plate as balanced, but that is for the lithium contrast
# inside BP1; it does not hold here.
#
# groupA_plate keeps only the plates where both classes are represented in
# reasonable numbers (3, 4, 5). It costs 30 of 74 cases, which is why it is a
# sensitivity analysis rather than the primary.
plate_ok <- meta[["sequencing plate"]] %in% c("3", "4", "5")
cohorts <- list(
  groupA       = (is_case | (is_ctrl & meta$grp == "A")) & tobacco_known,
  full         = (is_case | is_ctrl) & tobacco_known,
  groupA_plate = (is_case | (is_ctrl & meta$grp == "A")) & tobacco_known & plate_ok)

# Round-robin within each class after shuffling: fold sizes differ by at most
# one per class, which is the tightest balance 3 folds allow.
assign_folds <- function(idx, k) {
  shuffled <- sample(idx)
  ((seq_along(shuffled) - 1L) %% k)[order(order(shuffled))] + 1L
}

dir.create(file.path(impl, "prediction/splits"), showWarnings = FALSE, recursive = TRUE)

for (nm in names(cohorts)) {
  keep <- cohorts[[nm]]
  ids  <- meta$title[keep]
  y    <- as.integer(is_case[keep])
  stopifnot(!any(is.na(y)), length(ids) == length(y))

  set.seed(SEED)                     # same seed per cohort; the cohorts differ,
  rows <- vector("list", REPEATS * FOLDS)   # so the partitions do too
  n <- 0L
  for (r in seq_len(REPEATS)) {
    fold <- integer(length(ids))
    for (lv in c(0L, 1L)) {
      w <- which(y == lv)
      fold[w] <- assign_folds(w, FOLDS)
    }
    for (k in seq_len(FOLDS)) {
      n <- n + 1L
      rows[[n]] <- data.frame(split_id = n, repeat_no = r, fold_no = k,
                              sample_id = ids, y = y,
                              assignment = ifelse(fold == k, "test", "train"),
                              stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)

  # Every subject must be tested exactly once per repeat, or the pooled
  # out-of-fold set is not the cohort the design claims.
  chk <- aggregate(assignment ~ sample_id + repeat_no, out, function(v) sum(v == "test"))
  stopifnot(all(chk$assignment == 1L))

  f <- file.path(impl, sprintf("prediction/splits/splits_bpd_%s.csv", nm))
  write.csv(out, f, row.names = FALSE, quote = FALSE)

  te <- out[out$assignment == "test", ]
  ncase <- aggregate(y ~ split_id, te, sum)$y
  nall  <- aggregate(y ~ split_id, te, length)$y
  cat(sprintf("%-7s %3d cases / %3d controls (%.1f%%)  |  %d splits, test folds %d-%d with %d-%d cases\n",
              nm, sum(y), sum(1 - y), 100 * mean(y), n,
              min(nall), max(nall), min(ncase), max(ncase)))
  cat(sprintf("        -> %s\n", f))
}
