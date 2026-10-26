#!/usr/bin/env Rscript
# ==============================================================================
# sdfilter_leakage.R -- are the SD-filter cell-type hits real, or composition?
#
# WHY THIS IS NEEDED
#
# The 200-gene SD cut produced the first non-zero bipolar cell-type result in
# the project: 11 genes in T cells, permutation p = 0.018. Before that can be
# reported, the gene names have to be taken seriously:
#
#   CPA3, FCER1A    mast cell / basophil
#   TUBB1, TLN1     platelet
#   GZMB, PRF1      cytotoxic granule
#   MME (CD10)      neutrophil / progenitor
#
# Those are cell-identity markers, and finding them inside a *T cell* profile
# is the signature of deconvolution leakage rather than T-cell biology. bMIND
# estimates each lineage's profile per sample; if a sample's true mixture is
# not perfectly captured, the expression of another cell type bleeds into the
# T-cell estimate. A group difference in the amount of that other cell type
# then appears as a "T cell" differential expression result.
#
# This is the same circularity hazard that CLAUDE.md section 4.2 flags for
# mediation, arriving by a different route: the SD filter selects the highest-
# variance genes, and cell-identity markers are exactly the genes whose
# variance across a cohort is dominated by how much of that cell type each
# person has.
#
# THREE TESTS, IN ORDER OF SEVERITY
#
#   1. Are the genes LM22 signature markers? If the deconvolution used a gene
#      to decide the mixture, that gene's lineage profile is not independent
#      evidence about the lineage.
#   2. Does the gene track a cell FRACTION rather than the diagnosis? Reported
#      as the strongest correlation with any of the 22 fractions.
#   3. Does the group effect survive adjusting for cell composition? This is
#      decisive. If adding the fractions to the design kills the effect, it was
#      composition.
# ==============================================================================

suppressPackageStartupMessages({library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "de_final")); set.seed(481); options(width = 205)
RES <- "results"; FDR <- 0.05

p  <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
bp <- readRDS(file.path(ROOT, "de_analysis/results/bmind_profiles.rds"))
A <- bp$A; SE <- bp$SE
sym <- local({ m <- read.delim(file.path(ROOT, "de_analysis/data/ens2sym.tsv"),
                               header = FALSE, col.names = c("ens", "sym"))
               setNames(m$sym, m$ens) })

# LM22 signature genes, by HGNC symbol
lm22 <- read.delim(file.path(ROOT, "data/LM22.txt"), check.names = FALSE)
LM22G <- toupper(trimws(lm22[[1]]))
cat(sprintf("LM22 signature: %d marker genes\n", length(unique(LM22G))))

FR <- p$frac                                   # 474 x 22 cell fractions
hits <- read.csv(gzfile(file.path(RES, "sdfilter200_genes.csv.gz")), stringsAsFactors = FALSE)
hits <- hits[hits$significant & hits$level == "CT", ]
cat(sprintf("cell-type hits at the 200-gene cut: %d\n\n", nrow(hits)))

design_for <- function(contrast) {
  s <- p$sel[[contrast]]
  mm <- p$meta[s$keep, , drop = FALSE]; mm$grp <- s$grp
  mm$plate <- droplevels(factor(mm$plate)); mm$sex <- droplevels(factor(mm$sex))
  if ("group" %in% s$cov) mm$group <- droplevels(factor(mm$group))
  cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]])) > 1, logical(1))]
  X <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse = " + "))), data = mm)
  ne <- nonEstimable(X); if (!is.null(ne)) X <- X[, setdiff(colnames(X), ne), drop = FALSE]
  list(X = X, samples = mm$title) }
mkw <- function(S) { W <- 1/(S^2); f <- is.finite(W) & W > 0
  if (any(!f)) W[!f] <- min(W[f]); W/mean(W) }

# The five lineages as fraction columns, for the composition adjustment. Using
# four of five (one dropped) because they sum to 1 and would be collinear.
LINMAP <- list(
  gran = c("Neutrophils","Eosinophils","Mast cells resting","Mast cells activated"),
  mono = c("Monocytes","Macrophages M0","Macrophages M1","Macrophages M2",
           "Dendritic cells resting","Dendritic cells activated"),
  T    = c("T cells CD8","T cells CD4 naive","T cells CD4 memory resting",
           "T cells CD4 memory activated","T cells follicular helper",
           "T cells regulatory (Tregs)","T cells gamma delta"),
  NK   = c("NK cells resting","NK cells activated"),
  B    = c("B cells naive","B cells memory","Plasma cells"))

OUT <- list()
for (CN in unique(hits$contrast)) {
  d <- design_for(CN)
  LF <- sapply(names(LINMAP), function(L) rowSums(FR[d$samples, LINMAP[[L]], drop = FALSE]))
  Xadj <- cbind(d$X, LF[, c("gran","mono","NK","B")])   # drop T to avoid the simplex
  ok <- qr(Xadj)$rank == ncol(Xadj)
  cat(sprintf("\n%s\n== %s == composition-adjusted design full rank: %s\n%s\n",
              strrep("=",205), CN, ok, strrep("=",205)))
  cat("  lineage gene       | LM22 marker? | strongest fraction correlation        | logFC raw | logFC adj | %% lost | q adj  | verdict\n")

  for (ct in unique(hits$lineage[hits$contrast == CN])) {
    gs <- hits$gene[hits$contrast == CN & hits$lineage == ct]
    # A[genes, one lineage, samples] drops to a vector for one gene and keeps a
    # 3-D shape with drop=FALSE, so build the genes x samples matrix explicitly.
    Y <- matrix(A[gs, ct, d$samples], nrow = length(gs),
                dimnames = list(gs, d$samples))
    W <- mkw(SE[, ct, d$samples])[gs, , drop = FALSE]

    f1 <- eBayes(lmFit(Y, d$X, weights = W), trend = FALSE)
    t1 <- topTable(f1, coef = "grpcase", number = Inf, sort.by = "none")
    f2 <- eBayes(lmFit(Y, Xadj, weights = W), trend = FALSE)
    t2 <- topTable(f2, coef = "grpcase", number = Inf, sort.by = "none")
    q2 <- p.adjust(t2$P.Value, "BH")

    for (i in seq_along(gs)) {
      g <- gs[i]; sy <- toupper(ifelse(is.na(sym[g]), g, sym[g]))
      cors <- cor(A[g, ct, d$samples], FR[d$samples, ], use = "complete.obs")
      j <- which.max(abs(cors))
      lost <- 100 * (1 - abs(t2$logFC[i]) / abs(t1$logFC[i]))
      vd <- if (q2[i] < FDR && lost < 50) "survives adjustment" else
            if (q2[i] < FDR) "survives, but effect much reduced" else "LOST -- composition"
      cat(sprintf("  %-7s %-10s | %-12s | %-28s %+.2f | %+8.3f | %+9.3f | %5.1f | %.4f | %s\n",
                  ct, sy, if (sy %in% LM22G) "YES" else "no",
                  colnames(FR)[j], cors[j], t1$logFC[i], t2$logFC[i], lost, q2[i], vd))
      OUT[[length(OUT)+1]] <- data.frame(contrast = CN, lineage = ct, gene = g, symbol = sy,
        lm22_marker = sy %in% LM22G, top_fraction = colnames(FR)[j], top_fraction_cor = cors[j],
        logFC_raw = t1$logFC[i], logFC_adj = t2$logFC[i], pct_lost = lost,
        q_adj = q2[i], survives = q2[i] < FDR, stringsAsFactors = FALSE) } }
}
R <- do.call(rbind, OUT)
write.csv(R, file.path(RES, "sdfilter200_leakage.csv"), row.names = FALSE)

cat(sprintf("\n%s\nSUMMARY\n%s\n", strrep("=",205), strrep("=",205)))
for (CN in unique(R$contrast)) {
  z <- R[R$contrast == CN, ]
  cat(sprintf("  %-4s %d hits | %d are LM22 markers (%.0f%%) | %d survive composition adjustment (%.0f%%)\n",
              CN, nrow(z), sum(z$lm22_marker), 100*mean(z$lm22_marker),
              sum(z$survives), 100*mean(z$survives))) }
cat("\n  For reference, LM22 covers %0.1f%% of the tested transcriptome, so a marker\n")
cat(sprintf("  share far above %.1f%% indicates the filter is selecting signature genes.\n",
            100*mean(toupper(sym[dimnames(A)[[1]]]) %in% LM22G)))
cat("\nwrote results/sdfilter200_leakage.csv\n")
