#!/usr/bin/env Rscript
# ==============================================================================
# 18_estimator_and_covariates.R -- two questions, both of which must be settled
# before the results can be narrowed to a single reported pipeline.
#
# PART A. IF WE REPORT ONE ESTIMATOR, HOW MUCH IS THAT CHOICE DOING?
#
# 03 reported only DEG counts per estimator (voom 1,426 / trend 1,305 /
# QLF 1,341). A count comparison is nearly worthless here: three methods could
# each return ~1,400 genes and share almost none of them. The question is
# whether they agree GENE BY GENE, which needs correlation of the effect
# estimates, correlation of the rankings, and overlap of the called sets.
#
# PART B. ARE THE COVARIATES ACTUALLY DOING ANYTHING?
#
# Naming covariates in a formula is not evidence that they were fitted. This
# has already bitten this project once: bmind_de() accepts a covariate argument
# and silently discards it, which was only caught because two runs differing by
# four variables returned byte-identical p-values.
#
# So three checks, in increasing strength:
#   1. the design matrix contains the columns  -- necessary, nearly worthless
#   2. the design is full rank and grpcase is the tested coefficient
#   3. REMOVING the covariates changes the answer -- the only one that proves
#      they were used. If adjusted and unadjusted coefficients agree to machine
#      precision, the covariates are decorative.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "de_v2")); set.seed(481); options(width = 175)
OUT <- "results"
p <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
sym <- local({ m <- read.delim(file.path(ROOT, "de_analysis/data/ens2sym.tsv"),
                               header = FALSE, col.names = c("ens", "sym"))
               setNames(m$sym, m$ens) })

design_for <- function(contrast, covars = TRUE) {
  s <- p$sel[[contrast]]
  mm <- p$meta[s$keep, , drop = FALSE]; mm$grp <- s$grp
  mm$plate <- droplevels(factor(mm$plate)); mm$sex <- droplevels(factor(mm$sex))
  if ("group" %in% s$cov) mm$group <- droplevels(factor(mm$group))
  cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]])) > 1, logical(1))]
  f <- if (covars) paste("~ grp +", paste(cv, collapse = " + ")) else "~ grp"
  X <- model.matrix(as.formula(f), data = mm)
  ne <- nonEstimable(X); if (!is.null(ne)) X <- X[, setdiff(colnames(X), ne), drop = FALSE]
  list(X = X, samples = mm$title, cov_used = if (covars) cv else character(0),
       cov_dropped = setdiff(s$cov, cv)) }

fit_one <- function(contrast, method, covars = TRUE) {
  d <- design_for(contrast, covars)
  dge <- calcNormFactors(DGEList(p$counts[, d$samples, drop = FALSE]), method = "TMM")
  if (method == "voom") {
    tt <- topTable(eBayes(lmFit(voom(dge, d$X), d$X)), coef = "grpcase", number = Inf, sort.by = "none")
    data.frame(gene = rownames(tt), logFC = tt$logFC, t = tt$t, P = tt$P.Value,
               BH = tt$adj.P.Val, stringsAsFactors = FALSE)
  } else if (method == "trend") {
    f <- eBayes(lmFit(cpm(dge, log = TRUE, prior.count = 3), d$X), trend = TRUE)
    tt <- topTable(f, coef = "grpcase", number = Inf, sort.by = "none")
    data.frame(gene = rownames(tt), logFC = tt$logFC, t = tt$t, P = tt$P.Value,
               BH = tt$adj.P.Val, stringsAsFactors = FALSE)
  } else {
    dge <- estimateDisp(dge, d$X)
    r <- topTags(glmQLFTest(glmQLFit(dge, d$X), coef = "grpcase"), n = Inf, sort.by = "none")$table
    data.frame(gene = rownames(r), logFC = r$logFC, t = r$F, P = r$PValue,
               BH = r$FDR, stringsAsFactors = FALSE)
  } }

## ============================ PART A =========================================
cat("================ PART A: do the three estimators agree gene by gene? ================\n")
MET <- c("voom", "trend", "QLF")
AGREE <- list(); TOPOV <- list()
for (cn in c("LI", "BPD")) {
  F <- lapply(MET, function(m) fit_one(cn, m)); names(F) <- MET
  g <- Reduce(intersect, lapply(F, `[[`, "gene"))
  L <- sapply(MET, function(m) F[[m]]$logFC[match(g, F[[m]]$gene)])
  R <- sapply(MET, function(m) rank(F[[m]]$P[match(g, F[[m]]$gene)]))
  DEG <- lapply(MET, function(m) F[[m]]$gene[F[[m]]$BH < 0.05]); names(DEG) <- MET

  cat(sprintf("\n--- %s (n=%d, %d genes) ---\n", cn, length(design_for(cn)$samples), length(g)))
  cat("log2FC, Pearson:\n"); print(round(cor(L), 4))
  cat("p-value rank, Spearman:\n"); print(round(cor(R), 4))
  cat(sprintf("DEG at BH 0.05: %s\n", paste(sprintf("%s=%d", MET, lengths(DEG)), collapse = "  ")))
  for (i in 1:2) for (j in (i+1):3) {
    a <- DEG[[i]]; b <- DEG[[j]]
    u <- length(union(a, b)); it <- length(intersect(a, b))
    AGREE[[length(AGREE)+1]] <- data.frame(contrast = cn, m1 = MET[i], m2 = MET[j],
      n1 = length(a), n2 = length(b), shared = it,
      jaccard = if (u == 0) NA else round(it/u, 4),
      pct_of_smaller = if (min(length(a),length(b)) == 0) NA else round(100*it/min(length(a),length(b)), 1),
      logFC_pearson = round(cor(L[,i], L[,j]), 4),
      rank_spearman = round(cor(R[,i], R[,j], method = "spearman"), 4),
      stringsAsFactors = FALSE) }
  # concordance among the genes each method ranks highest
  for (N in c(50, 100, 500)) {
    tops <- lapply(MET, function(m) g[order(F[[m]]$P[match(g, F[[m]]$gene)])][1:N])
    TOPOV[[length(TOPOV)+1]] <- data.frame(contrast = cn, topN = N,
      shared_all_three = length(Reduce(intersect, tops)),
      pct = round(100*length(Reduce(intersect, tops))/N, 1), stringsAsFactors = FALSE) }
}
A <- do.call(rbind, AGREE); TO <- do.call(rbind, TOPOV)
cat("\n--- pairwise DEG-set overlap ---\n"); print(A, row.names = FALSE)
cat("\n--- genes shared by ALL THREE methods among each one's top N ---\n"); print(TO, row.names = FALSE)
write.csv(A,  file.path(OUT, "estimator_agreement.csv"), row.names = FALSE)
write.csv(TO, file.path(OUT, "estimator_top_overlap.csv"), row.names = FALSE)

## ============================ PART B =========================================
cat("\n\n================ PART B: are the covariates actually fitted? ================\n")
VER <- list()
for (cn in c("LI", "BPD")) {
  d <- design_for(cn, TRUE); d0 <- design_for(cn, FALSE)
  cat(sprintf("\n--- %s ---\n", cn))
  cat(sprintf("  design: %d x %d | rank %d | full rank: %s\n", nrow(d$X), ncol(d$X),
              qr(d$X)$rank, qr(d$X)$rank == ncol(d$X)))
  cat(sprintf("  columns: %s\n", paste(colnames(d$X), collapse = ", ")))
  cat(sprintf("  covariates requested and used : %s\n", paste(d$cov_used, collapse = ", ")))
  cat(sprintf("  dropped as constant           : %s\n",
              if (length(d$cov_dropped)) paste(d$cov_dropped, collapse = ", ") else "none"))
  cat(sprintf("  tested coefficient present    : %s\n", "grpcase" %in% colnames(d$X)))

  fa <- fit_one(cn, "voom", TRUE); fu <- fit_one(cn, "voom", FALSE)
  g <- intersect(fa$gene, fu$gene)
  la <- fa$logFC[match(g, fa$gene)]; lu <- fu$logFC[match(g, fu$gene)]
  identical_coef <- isTRUE(all.equal(la, lu, tolerance = 1e-12))
  cat(sprintf("\n  PROOF TEST -- does removing the covariates change the answer?\n"))
  cat(sprintf("    adjusted vs unadjusted log2FC identical : %s   <-- must be FALSE\n", identical_coef))
  cat(sprintf("    Pearson r between them                  : %.5f\n", cor(la, lu)))
  cat(sprintf("    median |difference| in log2FC           : %.5f\n", median(abs(la - lu))))
  cat(sprintf("    max |difference| in log2FC              : %.5f\n", max(abs(la - lu))))
  cat(sprintf("    DEG adjusted %d vs unadjusted %d\n",
              sum(fa$BH < 0.05), sum(fu$BH < 0.05)))
  VER[[length(VER)+1]] <- data.frame(contrast = cn, n = nrow(d$X), n_coef = ncol(d$X),
    rank = qr(d$X)$rank, full_rank = qr(d$X)$rank == ncol(d$X),
    covariates = paste(d$cov_used, collapse = ";"),
    dropped = paste(d$cov_dropped, collapse = ";"),
    coef_identical_without_covariates = identical_coef,
    r_adj_vs_unadj = round(cor(la, lu), 5),
    median_abs_diff_logFC = round(median(abs(la - lu)), 5),
    deg_adjusted = sum(fa$BH < 0.05), deg_unadjusted = sum(fu$BH < 0.05),
    stringsAsFactors = FALSE)
}
V <- do.call(rbind, VER)
write.csv(V, file.path(OUT, "covariate_verification.csv"), row.names = FALSE)
cat("\n--- verdict ---\n")
for (i in seq_len(nrow(V)))
  cat(sprintf("  %-4s %s\n", V$contrast[i],
      if (!V$coef_identical_without_covariates[i] && V$full_rank[i])
        "covariates ARE fitted: full rank, and dropping them changes every coefficient"
      else "PROBLEM -- covariates may not be entering the model"))
cat("\nwrote results/estimator_agreement.csv, estimator_top_overlap.csv, covariate_verification.csv\n")
