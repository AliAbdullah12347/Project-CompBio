#!/usr/bin/env Rscript
# Which lithium model produced the deposited Supplementary File S1?
#
# The BD model from the Supplementary Methods reproduces the paper's six DEGs
# almost exactly, so the pipeline is sound. The lithium sheet does not match the
# documented lithium model, so this script fits a grid of plausible variants and
# reports which one agrees with the deposited table. Diagnostic, not part of the
# main pipeline.

suppressPackageStartupMessages({ library(edgeR); library(limma) })

script_path <- sub("^--file=", "",
                   grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
ROOT <- if (!is.na(script_path)) dirname(dirname(normalizePath(script_path))) else "."
PROC <- file.path(ROOT, "data", "processed")
REF  <- file.path(ROOT, "data", "reference")

counts <- as.matrix(read.delim(file.path(PROC, "counts_filtered.tsv"),
                               row.names = 1, check.names = FALSE))
meta <- read.csv(file.path(PROC, "sample_metadata.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
meta <- meta[meta$qc_pass %in% c(TRUE, "True", "TRUE"), ]
meta <- meta[match(colnames(counts), meta$title), ]

meta$case    <- factor(ifelse(meta$diagnosis == "Control", "control", "case"),
                       levels = c("control", "case"))
meta$sex     <- factor(meta$sex);   meta$tobacco <- factor(meta$tobacco)
meta$group   <- factor(meta$group); meta$plate   <- factor(meta$plate)
meta$lithium <- as.numeric(meta$lithium)

ref_pre  <- read.csv(file.path(REF, "File_S1_DEGs__Li_DEGs_PreCellTypeCorrection.csv"))
ref_post <- read.csv(file.path(REF, "File_S1_DEGs__Li_DEGs_PostCellTypeCorrection.csv"))

CELLS <- c("b.cells.naive","b.cells.memory","plasma.cells","t.cells.cd8",
           "t.cells.cd4.naive","t.cells.cd4.memory.resting",
           "t.cells.cd4.memory.activated","t.cells.gamma.delta",
           "nk.cells.resting","nk.cells.activated","monocytes",
           "macrophages.m0","macrophages.m2","dendritic.cells.resting",
           "dendritic.cells.activated","mast.cells.resting","eosinophils",
           "neutrophils")   # 18 non-structurally-zero LM22 types

fit_one <- function(keep, rhs, coef) {
  y <- DGEList(counts = counts[, keep, drop = FALSE])
  y <- suppressMessages(calcNormFactors(y, method = "TMM"))
  design <- model.matrix(as.formula(paste("~", rhs)), data = meta[keep, ])
  r <- qr(design)$rank
  if (r < ncol(design)) design <- design[, sort(qr(design)$pivot[seq_len(r)]), drop = FALSE]
  if (!coef %in% colnames(design)) return(NULL)
  v <- voom(y, design, plot = FALSE)
  tt <- topTable(eBayes(lmFit(v, design)), coef = coef, number = Inf, sort.by = "none")
  data.frame(gene = rownames(tt), logFC = tt$logFC, t = tt$t,
             adj.P.Val = tt$adj.P.Val, stringsAsFactors = FALSE)
}

score <- function(ours, theirs) {
  if (is.null(ours)) return(NULL)
  m <- merge(ours, theirs, by = "gene", suffixes = c(".o", ".t"))
  c(r = cor(m$logFC.o, m$logFC.t),
    conc = mean(sign(m$logFC.o) == sign(m$logFC.t)),
    n_sig = sum(ours$adj.P.Val < 0.05),
    overlap = length(intersect(ours$gene[ours$adj.P.Val < 0.05],
                               theirs$gene[theirs$adj.P.Val < 0.05])))
}

BASE   <- "age + sex + tobacco + group + rin + plate + seqpc1 + seqpc2 + seqpc3"
NOGRP  <- "age + sex + tobacco + rin + plate + seqpc1 + seqpc2 + seqpc3"
CELLRHS <- paste(CELLS, collapse = " + ")

all444 <- rep(TRUE, nrow(meta)); cases <- meta$case == "case"

variants <- list(
  list("documented: cases+controls, with diagnosis", all444,
       paste("case +", BASE, "+ lithium"), "lithium"),
  list("all 444, WITHOUT diagnosis",                 all444,
       paste(BASE, "+ lithium"), "lithium"),
  list("all 444, no diagnosis, no group",            all444,
       paste(NOGRP, "+ lithium"), "lithium"),
  list("cases only (group dropped)",                 cases,
       paste(NOGRP, "+ lithium"), "lithium"),
  list("all 444, diagnosis, no tobacco",             all444,
       paste("case + age + sex + group + rin + plate + seqpc1 + seqpc2 + seqpc3 + lithium"), "lithium"),
  list("all 444, minimal (age+sex+lithium)",         all444,
       "age + sex + lithium", "lithium")
)

cat("=== vs deposited PRE-cell-type sheet (3,031 sig) ===\n")
cat(sprintf("%-44s %8s %8s %7s %8s\n", "model", "r(logFC)", "sign%", "n_sig", "overlap"))
for (v in variants) {
  s <- score(fit_one(v[[2]], v[[3]], v[[4]]), ref_pre)
  if (!is.null(s))
    cat(sprintf("%-44s %8.4f %7.1f%% %7d %8d\n", v[[1]], s["r"],
                100 * s["conc"], s["n_sig"], s["overlap"]))
}

cat("\n=== vs deposited POST-cell-type sheet (233 sig) ===\n")
cat(sprintf("%-44s %8s %8s %7s %8s\n", "model", "r(logFC)", "sign%", "n_sig", "overlap"))
post_variants <- list(
  list("documented + 18 cell types, all 444", all444,
       paste("case +", BASE, "+", CELLRHS, "+ lithium"), "lithium"),
  list("no diagnosis + 18 cell types, all 444", all444,
       paste(BASE, "+", CELLRHS, "+ lithium"), "lithium"),
  list("cases only + 18 cell types", cases,
       paste(NOGRP, "+", CELLRHS, "+ lithium"), "lithium")
)
for (v in post_variants) {
  s <- score(fit_one(v[[2]], v[[3]], v[[4]]), ref_post)
  if (!is.null(s))
    cat(sprintf("%-44s %8.4f %7.1f%% %7d %8d\n", v[[1]], s["r"],
                100 * s["conc"], s["n_sig"], s["overlap"]))
}
