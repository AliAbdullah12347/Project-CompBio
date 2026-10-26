#!/usr/bin/env Rscript
# ==============================================================================
# sdfilter_leakage2.R -- the sharper composition control.
#
# sdfilter_leakage.R adjusted for the five LINEAGE fractions and every hit
# survived. That control is too coarse for what the correlations showed: the
# BPD T-cell hits track 22-type fractions that the five-lineage sums hide --
# CPA3 and FCER1A with resting mast cells (r = +0.47, +0.46), PRF1 with
# resting NK (+0.48), GZMB with activated CD4 memory (+0.58). Mast cells are
# folded inside the granulocyte sum, so adjusting for that sum barely touches
# them.
#
# This re-runs the adjustment against all 22 CIBERSORTx fractions (21 columns,
# one dropped for the simplex), which is the strongest composition control the
# data supports, and tests the LM22-marker enrichment formally.
# ==============================================================================
suppressPackageStartupMessages({library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "de_final")); set.seed(481); options(width = 190)
RES <- "results"; FDR <- 0.05
p  <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
bp <- readRDS(file.path(ROOT, "de_analysis/results/bmind_profiles.rds"))
A <- bp$A; SE <- bp$SE; FR <- p$frac
sym <- local({ m <- read.delim(file.path(ROOT, "de_analysis/data/ens2sym.tsv"),
                 header = FALSE, col.names = c("ens","sym")); setNames(m$sym, m$ens) })
LM22G <- toupper(trimws(read.delim(file.path(ROOT,"data/LM22.txt"), check.names=FALSE)[[1]]))
design_for <- function(cn) {
  s <- p$sel[[cn]]; mm <- p$meta[s$keep, , drop=FALSE]; mm$grp <- s$grp
  mm$plate <- droplevels(factor(mm$plate)); mm$sex <- droplevels(factor(mm$sex))
  if ("group" %in% s$cov) mm$group <- droplevels(factor(mm$group))
  cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]]))>1, logical(1))]
  X <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse=" + "))), data=mm)
  ne <- nonEstimable(X); if (!is.null(ne)) X <- X[, setdiff(colnames(X), ne), drop=FALSE]
  list(X=X, samples=mm$title) }
mkw <- function(S){W<-1/(S^2);f<-is.finite(W)&W>0;if(any(!f))W[!f]<-min(W[f]);W/mean(W)}

## --- formal LM22 enrichment among the hits ---
allsym <- toupper(ifelse(is.na(sym[dimnames(A)[[1]]]), dimnames(A)[[1]], sym[dimnames(A)[[1]]]))
bg <- mean(allsym %in% LM22G)
h <- read.csv(gzfile(file.path(RES,"sdfilter200_genes.csv.gz")), stringsAsFactors=FALSE)
h <- h[h$significant & h$level=="CT", ]
hs <- toupper(h$symbol); k <- sum(hs %in% LM22G)
ft <- fisher.test(matrix(c(k, length(hs)-k,
       sum(allsym %in% LM22G)-k, sum(!(allsym %in% LM22G))-(length(hs)-k)), 2))
cat(sprintf("LM22 markers among the %d cell-type hits: %d (%.1f%%) vs %.1f%% background\n",
            length(hs), k, 100*k/length(hs), 100*bg))
cat(sprintf("  Fisher exact: odds ratio %.1f, p = %.2e\n\n", ft$estimate, ft$p.value))
# and among the top-200 SD sets themselves, which is the real question
for (CN in c("LI","BPD")) { d <- design_for(CN)
  for (ct in c("gran","mono","T","NK","B")) {
    Y <- A[, ct, d$samples]; keep <- order(apply(Y,1,sd), decreasing=TRUE)[1:200]
    ks <- toupper(ifelse(is.na(sym[rownames(Y)[keep]]), rownames(Y)[keep], sym[rownames(Y)[keep]]))
    cat(sprintf("  top-200 SD set, %s %-5s: %3d of 200 are LM22 markers (%.1f%%, %.1fx background)\n",
                CN, ct, sum(ks %in% LM22G), 100*mean(ks %in% LM22G), mean(ks %in% LM22G)/bg) ) } }

## --- adjustment against all 22 fractions ---
cat(sprintf("\n%s\nadjusting for ALL 22 cell fractions (21 columns + the base design)\n%s\n",
            strrep("=",190), strrep("=",190)))
OUT <- list()
for (CN in unique(h$contrast)) {
  d <- design_for(CN)
  F22 <- as.matrix(FR[d$samples, , drop=FALSE])
  F22 <- F22[, apply(F22, 2, sd) > 0, drop=FALSE]          # drop structural zeros
  F22 <- F22[, -ncol(F22), drop=FALSE]                      # drop one for the simplex
  Xa <- cbind(d$X, F22)
  qq <- qr(Xa); if (qq$rank < ncol(Xa)) {
    keepc <- qq$pivot[seq_len(qq$rank)]; Xa <- Xa[, keepc, drop=FALSE] }
  cat(sprintf("\n  %s: design %d x %d (rank %d), %d fraction columns\n",
              CN, nrow(Xa), ncol(Xa), qr(Xa)$rank, ncol(Xa)-ncol(d$X)))
  cat("   lineage gene          | LM22 | logFC raw | logFC adj22 |  %lost | q adj22 | verdict\n")
  for (ct in unique(h$lineage[h$contrast==CN])) {
    gs <- h$gene[h$contrast==CN & h$lineage==ct]
    Y <- matrix(A[gs, ct, d$samples], nrow=length(gs), dimnames=list(gs, d$samples))
    W <- mkw(SE[, ct, d$samples])[gs, , drop=FALSE]
    t1 <- topTable(eBayes(lmFit(Y, d$X, weights=W)), coef="grpcase", number=Inf, sort.by="none")
    t2 <- topTable(eBayes(lmFit(Y, Xa,  weights=W)), coef="grpcase", number=Inf, sort.by="none")
    q2 <- p.adjust(t2$P.Value, "BH")
    for (i in seq_along(gs)) {
      sy <- toupper(ifelse(is.na(sym[gs[i]]), gs[i], sym[gs[i]]))
      lost <- 100*(1 - abs(t2$logFC[i])/abs(t1$logFC[i]))
      cat(sprintf("   %-7s %-13s | %-4s | %+9.3f | %+11.3f | %6.1f | %7.4f | %s\n",
          ct, sy, if (sy %in% LM22G) "YES" else "no", t1$logFC[i], t2$logFC[i], lost, q2[i],
          if (q2[i] < FDR) "survives" else "LOST -- composition"))
      OUT[[length(OUT)+1]] <- data.frame(contrast=CN, lineage=ct, symbol=sy,
        lm22=sy %in% LM22G, logFC_raw=t1$logFC[i], logFC_adj22=t2$logFC[i],
        pct_lost=lost, q_adj22=q2[i], survives=q2[i] < FDR, stringsAsFactors=FALSE) } } }
R <- do.call(rbind, OUT)
write.csv(R, file.path(RES,"sdfilter200_leakage22.csv"), row.names=FALSE)
cat(sprintf("\n%s\n", strrep("=",190)))
for (CN in unique(R$contrast)) { z <- R[R$contrast==CN, ]
  cat(sprintf("  %-4s %2d hits | %2d survive adjustment for all 22 fractions (%.0f%%) | median effect lost %.1f%%\n",
              CN, nrow(z), sum(z$survives), 100*mean(z$survives), median(z$pct_lost))) }
cat("\nwrote results/sdfilter200_leakage22.csv\n")
