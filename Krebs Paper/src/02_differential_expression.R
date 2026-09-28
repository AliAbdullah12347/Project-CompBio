#!/usr/bin/env Rscript
# Differential expression, reproducing Krebs et al. 2020 Psychological Medicine.
#
# Models are taken verbatim from the Supplementary Methods ("Differential
# expression analysis"):
#
#   BD:      expression ~ age + sex + lithium + tobacco + group + RIN
#                         + plate + seqPC1-3 + diagnosis
#   Lithium: expression ~ diagnosis + age + sex + tobacco + group + RIN
#                         + plate + seqPC1-3 + lithium
#
# Note the lithium model retains `ascertainment group`. That variable is
# constant within cases (all cases are group A), so the lithium contrast cannot
# have been run on cases only -- the design would be rank-deficient. We
# therefore fit both sample sets and report which one reproduces the deposited
# supplementary table, rather than assuming.

suppressPackageStartupMessages({
  library(edgeR); library(limma)
})

# Resolve the project root from the script's own location so the pipeline can
# be run from any working directory.
script_path <- sub("^--file=", "",
                   grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
ROOT <- if (!is.na(script_path)) dirname(dirname(normalizePath(script_path))) else "."
if (!dir.exists(file.path(ROOT, "data"))) ROOT <- "."
PROC <- file.path(ROOT, "data", "processed")
REF  <- file.path(ROOT, "data", "reference")
RES  <- file.path(ROOT, "results")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------- load
cat("Loading data ...\n")
counts <- as.matrix(read.delim(file.path(PROC, "counts_filtered.tsv"),
                               row.names = 1, check.names = FALSE))
meta <- read.csv(file.path(PROC, "sample_metadata.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
meta <- meta[meta$qc_pass %in% c(TRUE, "True", "TRUE"), ]

# Join by key, never by position (titles all start with a digit, so R would
# otherwise mangle them via make.names).
stopifnot(!anyDuplicated(meta$title))
meta <- meta[match(colnames(counts), meta$title), ]
stopifnot(identical(as.character(meta$title), colnames(counts)))
cat(sprintf("  counts: %d genes x %d samples\n", nrow(counts), ncol(counts)))

meta$case      <- factor(ifelse(meta$diagnosis == "Control", "control", "case"),
                         levels = c("control", "case"))
meta$sex       <- factor(meta$sex)
meta$tobacco   <- factor(meta$tobacco)
meta$group     <- factor(meta$group)
meta$plate     <- factor(meta$plate)
meta$lithium   <- as.numeric(meta$lithium)

cat(sprintf("  cases %d / controls %d ; lithium users %d / non-users %d\n",
            sum(meta$case == "case"), sum(meta$case == "control"),
            sum(meta$lithium == 1), sum(meta$lithium == 0)))

# ---------------------------------------------------------------- pipeline
#' Fit one limma-voom model and return the full topTable.
#'
#' Filtering happened upstream on the full 444; TMM is computed on whatever
#' sample set is passed in, matching a per-analysis normalisation.
run_de <- function(keep, formula, coef, label) {
  y <- DGEList(counts = counts[, keep, drop = FALSE])
  y <- calcNormFactors(y, method = "TMM")
  design <- model.matrix(formula, data = meta[keep, , drop = FALSE])

  # Drop aliased columns so a constant covariate cannot silently break the fit.
  qr_rank <- qr(design)$rank
  if (qr_rank < ncol(design)) {
    keepcol <- qr(design)$pivot[seq_len(qr_rank)]
    dropped <- setdiff(colnames(design), colnames(design)[keepcol])
    cat(sprintf("    [%s] dropping rank-deficient column(s): %s\n",
                label, paste(dropped, collapse = ", ")))
    design <- design[, sort(keepcol), drop = FALSE]
  }
  stopifnot(coef %in% colnames(design))

  v   <- voom(y, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  tt  <- topTable(fit, coef = coef, number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  cat(sprintf("    [%s] n=%d  genes=%d  FDR<0.05: %d\n",
              label, sum(keep), nrow(tt), sum(tt$adj.P.Val < 0.05)))
  tt[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
}

all444 <- rep(TRUE, nrow(meta))
cases  <- meta$case == "case"

F_BD <- ~ age + sex + lithium + tobacco + group + rin + plate +
          seqpc1 + seqpc2 + seqpc3 + case
F_LI <- ~ case + age + sex + tobacco + group + rin + plate +
          seqpc1 + seqpc2 + seqpc3 + lithium

cat("\n=== BD case/control DE (all 444) ===\n")
bd <- run_de(all444, F_BD, "casecase", "BD")

cat("\n=== Lithium DE ===\n")
li_all   <- run_de(all444, F_LI, "lithium", "lithium | all 444")
li_cases <- run_de(cases,  F_LI, "lithium", "lithium | cases only")

# ---------------------------------------------------------------- validate
cat("\n=== Agreement with deposited Supplementary File S1 ===\n")
ref_pre <- read.csv(file.path(REF, "File_S1_DEGs__Li_DEGs_PreCellTypeCorrection.csv"),
                    stringsAsFactors = FALSE)
ref_bd  <- read.csv(file.path(REF, "File_S1_DEGs__BD_DEGs.csv"),
                    stringsAsFactors = FALSE)

compare <- function(ours, theirs, label) {
  m <- merge(ours, theirs, by = "gene", suffixes = c(".ours", ".theirs"))
  r_lfc <- cor(m$logFC.ours, m$logFC.theirs)
  r_t   <- cor(m$t.ours, m$t.theirs)
  conc  <- mean(sign(m$logFC.ours) == sign(m$logFC.theirs))
  n_ours   <- sum(ours$adj.P.Val < 0.05)
  n_theirs <- sum(theirs$adj.P.Val < 0.05)
  ov <- length(intersect(ours$gene[ours$adj.P.Val < 0.05],
                         theirs$gene[theirs$adj.P.Val < 0.05]))
  cat(sprintf(paste0("  %-24s genes matched %5d | r(logFC)=%.4f r(t)=%.4f ",
                     "| sign concordance %.1f%%\n%26s FDR<0.05  ours %d  theirs %d  overlap %d\n"),
              label, nrow(m), r_lfc, r_t, 100 * conc, "", n_ours, n_theirs, ov))
  invisible(c(r = r_lfc, n = nrow(m)))
}

compare(li_all,   ref_pre, "lithium | all 444")
compare(li_cases, ref_pre, "lithium | cases only")

cat("\n  BD model vs their 6 published BD DEGs:\n")
bd_hit <- bd[bd$gene %in% ref_bd$gene, c("gene", "logFC", "adj.P.Val")]
bd_cmp <- merge(bd_hit, ref_bd[, c("gene", "logFC", "adj.P.Val")],
                by = "gene", suffixes = c(".ours", ".theirs"))
print(bd_cmp, row.names = FALSE, digits = 4)
cat(sprintf("  our BD DEGs at FDR<0.05: %d (paper reports 6)\n",
            sum(bd$adj.P.Val < 0.05)))

# ---------------------------------------------------------------- write
write.csv(bd,       file.path(RES, "DE_BD_all444.csv"), row.names = FALSE)
write.csv(li_all,   file.path(RES, "DE_lithium_all444.csv"), row.names = FALSE)
write.csv(li_cases, file.path(RES, "DE_lithium_casesonly.csv"), row.names = FALSE)
cat("\nWrote 3 result tables to results/\n")
