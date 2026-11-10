#!/usr/bin/env Rscript
# ===========================================================================
# prep_bmind.R  --  per-lineage expression matrices for bmind-047, from the
#                   bMIND profiles that already exist.
#
# WHY THIS EXISTS RATHER THAN RE-RUNNING bMIND
#
# bmind-046's spec is "REBUILD bMIND profiles ... (~3.1h)". That rebuild was
# already carried out on 2026-10-01 by de_analysis/scripts/06_bmind_profiles.R
# -- all 31 chunks, 12368 genes x 5 lineages x 474 samples, no NAs, SE in
# 0.163-6.401 -- and the script has not changed since (last touched
# 2026-09-30 22:34, profiles written 2026-10-01 01:42).
#
# Re-running it would not reproduce that file. bMIND is a Gibbs sampler, so a
# second run returns a different posterior draw, and every published DE-arm
# result (celltype_limma, celltype_mde, celltype_null) was computed from the
# 1 Oct draw. Rebuilding would silently desynchronise the network arm from the
# DE arm for no gain. The profiles are the input to bmind-047; they are not
# themselves a result anyone is waiting on.
#
# WHAT THIS WRITES
#
# bmind-047 reads one lineage at a time and uses only the posterior mean A.
# The full object is 441 MB because it also carries SE, which 047 never opens
# and which would double the S3 transfer. This writes A alone, reshaped to the
# samples x genes layout WGCNA wants, as a named list of five matrices --
# the one shape bmind-047's normalise_bmind() handles without ambiguity.
#
#   Rscript network/cluster/prep_bmind.R
# ===========================================================================

CHUNK_DIR <- "de_analysis/results/bmind_profiles"
OUT       <- "de_analysis/results/bmind_lineage_expr.rds"
LINEAGES  <- c("gran", "mono", "T", "NK", "B")

fs <- sort(list.files(CHUNK_DIR, "^chunk_[0-9]+[.]rds$", full.names = TRUE))
if (!length(fs)) stop("no bMIND chunks under ", CHUNK_DIR)
cat(sprintf("== prep_bmind == %d chunks\n", length(fs)))

# Read the first chunk only to fix the axes, then stream the rest into a
# preallocated result. Reading the assembled 441 MB object instead would need
# ~900 MB of peak memory to produce a 234 MB answer.
q1 <- readRDS(fs[1])
dn <- dimnames(q1$A)
stopifnot(length(dim(q1$A)) == 3L, identical(dn[[2]], LINEAGES))
samples <- dn[[3]]
genes   <- unlist(lapply(fs, function(f) readRDS(f)$genes), use.names = FALSE)
stopifnot(!anyDuplicated(genes))
cat(sprintf("   %d genes x %d lineages x %d samples\n",
            length(genes), length(LINEAGES), length(samples)))

out <- setNames(lapply(LINEAGES, function(l)
  matrix(NA_real_, length(samples), length(genes),
         dimnames = list(samples, genes))), LINEAGES)

k <- 0L
for (f in fs) {
  q <- readRDS(f)
  stopifnot(identical(dimnames(q$A)[[3]], samples))
  idx <- (k + 1L):(k + length(q$genes)); k <- k + length(q$genes)
  # A is genes x lineages x samples; t() of a gene x sample slice is the
  # samples x genes orientation blockwiseModules() and modulePreservation()
  # both expect.
  for (l in LINEAGES) out[[l]][, idx] <- t(q$A[, l, ])
  rm(q)
}
stopifnot(k == length(genes))

na <- vapply(out, function(m) sum(is.na(m)), numeric(1))
if (any(na > 0)) stop("NAs left after fill: ", paste(names(na)[na > 0], collapse = ", "))
cat(sprintf("   ranges  %s\n",
            paste(sprintf("%s [%.2f, %.2f]", LINEAGES,
                          vapply(out, min, 0), vapply(out, max, 0)), collapse = "  ")))

saveRDS(out, OUT)
cat(sprintf("wrote %s (%.0f MB)\n", OUT, file.size(OUT) / 1048576))
