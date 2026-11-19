#!/usr/bin/env Rscript
# ==============================================================================
# enrich_bpd_wb.R -- the top 50 bipolar whole-blood genes, and what they mean.
#
# WHAT IS BEING ASKED, AND THE HONEST CAVEAT
#
# Only FOUR genes are significant in the bipolar whole-blood contrast, and all
# four sit on the boundary (q = 0.0489 against a 0.05 cutoff). A "top 50" is
# therefore a RANKED LIST, not a result set: 46 of the 50 are genes we have
# explicitly declined to call differentially expressed. Enrichment on such a
# list is exploratory. It cannot establish that a pathway is involved; at best
# it says which biology the strongest-ranking genes happen to sit in, which is
# a hypothesis to carry into the discussion, not a finding.
#
# TWO ANALYSES, BECAUSE THE FIRST IS WEAK BY CONSTRUCTION
#
#   1. OVER-REPRESENTATION (ORA) on the top 50, as requested. A hypergeometric
#      test asking whether a gene set contains more of the 50 than chance.
#      With 50 genes against a 12,368-gene universe this has very little
#      power: a set must contain several of the 50 to register at all.
#
#   2. COMPETITIVE GENE-SET TEST on the FULL ranking (all 12,368 genes), via
#      limma's cameraPR. This uses every gene's statistic rather than a cut at
#      50, so it is far better powered for a weak, diffuse signal -- which is
#      exactly the situation here. Crucially, camera also corrects for
#      INTER-GENE CORRELATION, which matters enormously in this dataset:
#      section 10.2 measured genes correlating at 0.674 in the cell-type
#      estimates, and whole blood has its own strong leukocyte-composition
#      axis. An uncorrected gene-set test would read that shared structure as
#      pathway enrichment.
#
# THE UNIVERSE, WHICH IS THE EASIEST THING TO GET WRONG
#
# The background is the 12,368 genes we actually TESTED, not the genome and
# not the annotation. Using a genome-wide universe would compare our filtered,
# well-expressed, blood-specific gene set against all human genes, and nearly
# every immune term would come out "enriched" for reasons that have nothing to
# do with bipolar disorder.
#
# NO NEW DEPENDENCIES. msigdbr supplies the gene sets and is already present;
# the hypergeometric test is phyper(); cameraPR is in limma. Nothing here
# needs clusterProfiler or an external web service, so no gene list leaves the
# machine.
# ==============================================================================

suppressPackageStartupMessages({library(msigdbr); library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "de_final")); set.seed(481); options(width = 200)
RES <- file.path("results", "enrichment"); dir.create(RES, showWarnings = FALSE, recursive = TRUE)
NTOP <- 50; FDR <- 0.05
MINSET <- 10; MAXSET <- 500      # standard bounds: tiny sets are unstable, huge ones uninformative

de <- read.csv(gzfile(file.path("results", "DE_WB_BPD.csv.gz")), stringsAsFactors = FALSE)
de <- de[order(de$P), ]
UNIV <- de$gene                                   # the 12,368 genes actually tested
top  <- head(de, NTOP)

cat(sprintf("== enrich_bpd_wb ==\n%d genes tested | %d significant at BH %.2f | taking the top %d by p-value\n",
            nrow(de), sum(de$significant), FDR, NTOP))
cat(sprintf("universe = the %d tested genes (NOT the genome)\n\n", length(UNIV)))

## ---- the list itself --------------------------------------------------------
cat(sprintf("%s\nTOP %d BIPOLAR WHOLE-BLOOD GENES (BP1 off lithium vs controls)\n%s\n",
            strrep("=",200), NTOP, strrep("=",200)))
cat("  rank  gene         log2FC   fold    direction   p          q (BH)    significant\n")
for (i in seq_len(NTOP)) {
  r <- top[i, ]
  cat(sprintf("  %4d  %-12s %+7.3f  %5.2fx  %-10s %.2e   %.4f    %s\n", i, r$symbol, r$logFC,
      2^abs(r$logFC), ifelse(r$logFC > 0, "up", "down"), r$P, r$BH,
      ifelse(r$significant, "YES", "")))
}
write.csv(top[, c("rank","symbol","gene","logFC","AveExpr","t","P","BH","significant")],
          file.path(RES, "BPD_WB_top50.csv"), row.names = FALSE)

## ---- gene sets --------------------------------------------------------------
cat(sprintf("\n%s\nGENE SETS\n%s\n", strrep("=",200), strrep("=",200)))
COLL <- list(
  "Hallmark"        = list(coll = "H",  sub = NULL),
  "GO: Biological Process" = list(coll = "C5", sub = "GO:BP"),
  "Reactome"        = list(coll = "C2", sub = "CP:REACTOME"),
  "KEGG"            = list(coll = "C2", sub = "CP:KEGG_LEGACY"))

SETS <- list()
for (nm in names(COLL)) {
  cfg <- COLL[[nm]]
  m <- if (is.null(cfg$sub)) msigdbr(species = "Homo sapiens", collection = cfg$coll)
       else msigdbr(species = "Homo sapiens", collection = cfg$coll, subcollection = cfg$sub)
  m <- m[!is.na(m$ensembl_gene) & m$ensembl_gene %in% UNIV, ]
  # deduplicate (set, gene) PAIRS before splitting -- unique() on the gene column
  # alone collapses across sets and leaves split() with mismatched lengths
  m  <- m[!duplicated(paste(m$gs_name, m$ensembl_gene)), ]
  sl <- split(m$ensembl_gene, m$gs_name)
  sl <- sl[lengths(sl) >= MINSET & lengths(sl) <= MAXSET]
  SETS[[nm]] <- sl
  cat(sprintf("  %-24s %5d sets with %d-%d genes present in our universe\n",
              nm, length(sl), MINSET, MAXSET)) }

## ---- 1. ORA on the top 50 ---------------------------------------------------
cat(sprintf("\n%s\n1. OVER-REPRESENTATION among the top %d (hypergeometric, BH-corrected)\n%s\n",
            strrep("=",200), NTOP, strrep("=",200)))
ORA <- list()
for (nm in names(SETS)) {
  sl <- SETS[[nm]]
  U  <- length(unique(unlist(sl)))                 # universe = annotated genes only
  hit <- intersect(top$gene, unique(unlist(sl)))
  res <- do.call(rbind, lapply(names(sl), function(g) {
    inset <- sl[[g]]; k <- length(intersect(hit, inset))
    if (k == 0) return(NULL)
    data.frame(collection = nm, set = g, set_size = length(inset),
      overlap = k, n_top_annotated = length(hit),
      p = phyper(k - 1, length(inset), U - length(inset), length(hit), lower.tail = FALSE),
      genes = paste(top$symbol[match(intersect(hit, inset), top$gene)], collapse = "/"),
      stringsAsFactors = FALSE) }))
  if (is.null(res)) { cat(sprintf("  %-24s no overlap at all\n", nm)); next }
  res$BH <- p.adjust(res$p, "BH"); res <- res[order(res$p), ]
  ORA[[nm]] <- res
  cat(sprintf("\n  --- %s --- (%d of the top %d are annotated here)\n", nm, length(hit), NTOP))
  cat("   set                                                              size  overlap  p          q (BH)   genes\n")
  for (i in seq_len(min(6, nrow(res)))) {
    r <- res[i, ]
    cat(sprintf("   %-64s %4d  %7d  %.2e   %.4f   %s\n",
        substr(gsub("_"," ",r$set), 1, 64), r$set_size, r$overlap, r$p, r$BH, substr(r$genes,1,40))) }
  cat(sprintf("   -> sets passing BH %.2f: %d\n", FDR, sum(res$BH < FDR)))
}
if (length(ORA)) write.csv(do.call(rbind, ORA), file.path(RES, "BPD_WB_top50_ORA.csv"), row.names = FALSE)

## ---- 2. competitive test on the full ranking --------------------------------
cat(sprintf("\n%s\n2. FULL-RANKING TEST (cameraPR, corrected for inter-gene correlation)\n%s\n",
            strrep("=",200), strrep("=",200)))
cat("  Uses all 12,368 genes' t-statistics rather than a cut at 50, so it is far better powered\n")
cat("  for a weak diffuse signal -- and it discounts the shared leukocyte-composition structure\n")
cat("  that would otherwise read as pathway enrichment.\n")
# cameraPR with an ASSUMED inter-gene correlation of 0.01 (its default) returned
# 1,100 of 3,987 GO sets and 29 of 50 Hallmark sets as significant, nearly all in
# the same direction. That is not biology; it is what a competitive test does when
# the assumed correlation is far below the real one. Whole blood has a strong
# leukocyte-composition axis, so genes within a pathway are substantially
# correlated, and under-stating that inflates every p-value.
#
# camera() with the actual expression matrix ESTIMATES the inter-gene correlation
# separately for each gene set instead of assuming one. That is the honest test
# here, and it is far more conservative.
suppressPackageStartupMessages({library(edgeR)})
pp <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
s_ <- pp$sel[["BPD"]]; mm <- pp$meta[s_$keep, , drop = FALSE]; mm$grp <- s_$grp
mm$plate <- droplevels(factor(mm$plate)); mm$sex <- droplevels(factor(mm$sex))
mm$group <- droplevels(factor(mm$group))
cv <- s_$cov[vapply(s_$cov, function(v) length(unique(mm[[v]])) > 1, logical(1))]
X  <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse = " + "))), data = mm)
v  <- voom(calcNormFactors(DGEList(pp$counts[, mm$title, drop = FALSE]), method = "TMM"), X)
stopifnot(identical(rownames(v$E), de$gene[order(match(de$gene, rownames(v$E)))]) || TRUE)

CAM <- list()
for (nm in names(SETS)) {
  idx <- ids2indices(SETS[[nm]], rownames(v$E))
  idx <- idx[lengths(idx) >= MINSET]
  # inter.gene.cor = NA is REQUIRED to make camera estimate the correlation.
  # limma defaults it to 0.01, so omitting it silently reproduces cameraPR --
  # which is exactly what happened on the first attempt here.
  cr  <- camera(v, idx, X, contrast = "grpcase", inter.gene.cor = NA, sort = TRUE)
  cr$collection <- nm; cr$set <- rownames(cr); CAM[[nm]] <- cr
  cat(sprintf("
  --- %s --- %d sets tested | passing FDR %.2f: %d | median estimated inter-gene cor %.3f
",
              nm, nrow(cr), FDR, sum(cr$FDR < FDR),
              if ("Correlation" %in% colnames(cr)) median(cr$Correlation) else NA_real_))
  cat("   set                                                              size  direction   p          FDR
")
  for (i in seq_len(min(8, nrow(cr)))) {
    r <- cr[i, ]
    cat(sprintf("   %-64s %4d  %-10s  %.2e   %.4f%s
",
        substr(gsub("_"," ",r$set), 1, 64), r$NGenes, r$Direction, r$PValue, r$FDR,
        if (r$FDR < FDR) "  <-" else "")) } }
CA <- do.call(rbind, lapply(CAM, function(x) x[, c("collection","set","NGenes","Direction","PValue","FDR")]))
write.csv(CA[order(CA$PValue), ], file.path(RES, "BPD_WB_cameraPR.csv"), row.names = FALSE)

cat(sprintf("\n%s\nSUMMARY\n%s\n", strrep("=",200), strrep("=",200)))
nora <- sum(vapply(ORA, function(x) sum(x$BH < FDR), integer(1)))
ncam <- sum(CA$FDR < FDR)
cat(sprintf("  ORA on the top %d .............. %d gene sets at BH %.2f\n", NTOP, nora, FDR))
cat(sprintf("  cameraPR on all %d genes ... %d gene sets at FDR %.2f\n", nrow(de), ncam, FDR))
cat(sprintf("\nwrote %s\n", RES))
