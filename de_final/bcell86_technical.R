#!/usr/bin/env Rscript
# ==============================================================================
# bcell86_technical.R -- are the 86 B-cell genes biology or RNA quality?
#
# The mean filter at 1,000 genes gave 86 significant genes in B cells, 97% of
# which survive adjustment for all 22 cell fractions -- so, unlike every other
# filtering result, not composition. But the gene classes are a warning:
# SCARNA10, SNORA53, RPL39P3, RPS23P8, RPL32, RPS21, HIST1H4C. Small nucleolar
# RNAs, ribosomal protein genes, their pseudogenes, and a histone.
#
# That is the classic signature of RNA degradation and library composition,
# not of B-cell biology. snoRNA and ribosomal-protein content shifts with RNA
# integrity and with how much of the library is taken up by rRNA-adjacent
# species -- and cell-fraction adjustment cannot touch it, which is exactly
# why these survived when the T-cell set did not.
#
# Three checks: what fraction of the 86 are snoRNA/scaRNA/ribosomal/pseudogene
# by symbol; whether their expression tracks RIN; and whether the group
# difference survives when RIN is allowed a more flexible fit than the single
# linear term the main model gives it.
# ==============================================================================
suppressPackageStartupMessages({library(limma); library(splines)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "de_final")); set.seed(481); options(width = 190)
RES <- "results"; FDR <- 0.05
p  <- readRDS(file.path(ROOT,"de_analysis/data/prep.rds"))
bp <- readRDS(file.path(ROOT,"de_analysis/results/bmind_profiles.rds"))
A <- bp$A; SE <- bp$SE
sym <- local({m<-read.delim(file.path(ROOT,"de_analysis/data/ens2sym.tsv"),header=FALSE,
                col.names=c("ens","sym"));setNames(m$sym,m$ens)})
s <- p$sel[["BPD"]]; mm <- p$meta[s$keep,]; mm$grp <- s$grp
mm$plate<-droplevels(factor(mm$plate)); mm$sex<-droplevels(factor(mm$sex))
if ("group" %in% s$cov) mm$group<-droplevels(factor(mm$group))
cv <- s$cov[vapply(s$cov,function(v)length(unique(mm[[v]]))>1,logical(1))]
X <- model.matrix(as.formula(paste("~ grp +",paste(cv,collapse=" + "))),data=mm)
ne<-nonEstimable(X); if(!is.null(ne)) X<-X[,setdiff(colnames(X),ne),drop=FALSE]
mkw <- function(S){W<-1/(S^2);f<-is.finite(W)&W>0;if(any(!f))W[!f]<-min(W[f]);W/mean(W)}
Y <- A[, "B", mm$title]; W <- mkw(SE[, "B", mm$title])
keep <- order(rowMeans(Y), decreasing=TRUE)[seq_len(1000)]
Ys <- Y[keep,,drop=FALSE]; Ws <- W[keep,,drop=FALSE]
t1 <- topTable(eBayes(lmFit(Ys,X,weights=Ws),trend=TRUE),coef="grpcase",number=Inf,sort.by="none")
q1 <- p.adjust(t1$P.Value,"BH"); sg <- which(q1 < FDR)
gs <- rownames(Ys)[sg]; sy <- toupper(ifelse(is.na(sym[gs]), gs, sym[gs]))

## 1. gene classes
pat <- "^(SNOR|SCARNA|RNU|RNY|RPL|RPS|MRPL|MRPS|HIST)|P[0-9]+$|^RP[0-9]+-"
cls <- grepl(pat, sy)
cat(sprintf("== the %d B-cell genes ==\n", length(gs)))
cat(sprintf("  snoRNA / scaRNA / ribosomal-protein / histone / pseudogene by symbol: %d of %d (%.0f%%)\n",
            sum(cls), length(sy), 100*mean(cls)))
# the same share among all 1,000 tested, for comparison
sy_all <- toupper(ifelse(is.na(sym[rownames(Ys)]), rownames(Ys), sym[rownames(Ys)]))
cat(sprintf("  same share among all 1,000 genes tested:                              %.0f%%  (%.1fx enrichment)\n",
            100*mean(grepl(pat, sy_all)), mean(cls)/mean(grepl(pat, sy_all))))

## 2. do they track RIN?
rin <- mm$rin
r_hit <- apply(Ys[sg,,drop=FALSE], 1, function(y) cor(y, rin, method="spearman"))
r_all <- apply(Ys, 1, function(y) cor(y, rin, method="spearman"))
cat(sprintf("\n  Spearman with RIN: hits median %+.3f (|r| > 0.2 in %d of %d)\n",
            median(r_hit), sum(abs(r_hit) > 0.2), length(r_hit)))
cat(sprintf("                     all 1,000 median %+.3f\n", median(r_all)))
cat(sprintf("  Wilcoxon, |r| of hits vs all: p = %.2e\n",
            wilcox.test(abs(r_hit), abs(r_all))$p.value))
cat(sprintf("  RIN itself: ctrl %.2f vs case %.2f, t-test p = %.4f\n",
            mean(rin[mm$grp=="ctrl"]), mean(rin[mm$grp=="case"]),
            t.test(rin ~ mm$grp)$p.value))

## 3. give RIN a flexible fit instead of one linear term
Xs <- cbind(X[, setdiff(colnames(X), "rin"), drop=FALSE], ns(rin, df=4))
qs <- qr(Xs); if (qs$rank < ncol(Xs)) Xs <- Xs[, qs$pivot[seq_len(qs$rank)], drop=FALSE]
t3 <- topTable(eBayes(lmFit(Ys,Xs,weights=Ws),trend=TRUE),coef="grpcase",number=Inf,sort.by="none")
q3 <- p.adjust(t3$P.Value,"BH")
cat(sprintf("\n  with RIN as a 4-df spline instead of a linear term: %d of %d remain significant\n",
            sum(q3[sg] < FDR), length(sg)))
cat(sprintf("  median %% of effect lost: %.1f%%\n",
            median(100*(1-abs(t3$logFC[sg])/abs(t1$logFC[sg])))))
write.csv(data.frame(gene=gs, symbol=sy, technical_class=cls, logFC=t1$logFC[sg],
  q=q1[sg], cor_rin=r_hit, q_rin_spline=q3[sg], survives_rin_spline=q3[sg] < FDR,
  stringsAsFactors=FALSE), file.path(RES,"bcell86_technical.csv"), row.names=FALSE)
cat("wrote results/bcell86_technical.csv\n")
