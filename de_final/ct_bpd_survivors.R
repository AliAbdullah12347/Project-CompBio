#!/usr/bin/env Rscript
# ==============================================================================
# ct_bpd_survivors.R -- composition test for the cell-type cells that survived
# the selection-aware null.
#
# meanfilter_validate.R left four survivors. Two are trivial or already
# understood; two are large and need the same composition test that section
# 10.5 applied to the 200-gene hits:
#
#   bipolar T cells,  SD filter at 2,000 genes :   0 -> 192  (p = 0.0348)
#   bipolar B cells,  MEAN filter at 1,000     :   0 ->  86  (p = 0.0149)
#
# "Survives the selection-aware null" only means the count is unlikely under
# permuted labels. It says nothing about whether the genes are cell-identity
# markers tracking a composition difference, which is what the 200-gene
# version of this result turned out to be (6 of 11 lost on adjustment, LM22
# enrichment 9.2x among hits). The same two questions therefore get asked of
# these larger sets.
# ==============================================================================
suppressPackageStartupMessages({library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "de_final")); set.seed(481); options(width = 190)
RES <- "results"; FDR <- 0.05
p  <- readRDS(file.path(ROOT,"de_analysis/data/prep.rds"))
bp <- readRDS(file.path(ROOT,"de_analysis/results/bmind_profiles.rds"))
A <- bp$A; SE <- bp$SE; FR <- p$frac
sym <- local({m<-read.delim(file.path(ROOT,"de_analysis/data/ens2sym.tsv"),header=FALSE,
                col.names=c("ens","sym"));setNames(m$sym,m$ens)})
LM22G <- toupper(trimws(read.delim(file.path(ROOT,"data/LM22.txt"),check.names=FALSE)[[1]]))
ALLSYM <- toupper(ifelse(is.na(sym[dimnames(A)[[1]]]), dimnames(A)[[1]], sym[dimnames(A)[[1]]]))
BG <- mean(ALLSYM %in% LM22G)

s <- p$sel[["BPD"]]; mm <- p$meta[s$keep,]; mm$grp <- s$grp
mm$plate<-droplevels(factor(mm$plate)); mm$sex<-droplevels(factor(mm$sex))
if ("group" %in% s$cov) mm$group<-droplevels(factor(mm$group))
cv <- s$cov[vapply(s$cov,function(v)length(unique(mm[[v]]))>1,logical(1))]
X <- model.matrix(as.formula(paste("~ grp +",paste(cv,collapse=" + "))),data=mm)
ne<-nonEstimable(X); if(!is.null(ne)) X<-X[,setdiff(colnames(X),ne),drop=FALSE]
F22 <- as.matrix(FR[mm$title,]); F22 <- F22[,apply(F22,2,sd)>0,drop=FALSE]
Xa <- cbind(X, F22[,-ncol(F22),drop=FALSE])
qa <- qr(Xa); if (qa$rank < ncol(Xa)) Xa <- Xa[, qa$pivot[seq_len(qa$rank)], drop=FALSE]
mkw <- function(S){W<-1/(S^2);f<-is.finite(W)&W>0;if(any(!f))W[!f]<-min(W[f]);W/mean(W)}

JOBS <- list(list(ct="T", stat="SD",   n=2000), list(ct="B", stat="MEAN", n=1000))
OUT <- NULL
for (J in JOBS) {
  Y <- A[, J$ct, mm$title]; W <- mkw(SE[, J$ct, mm$title])
  v <- if (J$stat=="SD") apply(Y,1,sd) else rowMeans(Y)
  keep <- order(v, decreasing=TRUE)[seq_len(J$n)]
  Ys <- Y[keep,,drop=FALSE]; Ws <- W[keep,,drop=FALSE]
  t1 <- topTable(eBayes(lmFit(Ys,X,weights=Ws),trend=TRUE),coef="grpcase",number=Inf,sort.by="none")
  q1 <- p.adjust(t1$P.Value,"BH"); sg <- which(q1 < FDR)
  t2 <- topTable(eBayes(lmFit(Ys,Xa,weights=Ws),trend=TRUE),coef="grpcase",number=Inf,sort.by="none")
  q2 <- p.adjust(t2$P.Value,"BH")
  ks <- toupper(ifelse(is.na(sym[rownames(Ys)[sg]]), rownames(Ys)[sg], sym[rownames(Ys)[sg]]))
  lost <- 100*(1-abs(t2$logFC[sg])/abs(t1$logFC[sg]))
  surv <- q2[sg] < FDR
  cat(sprintf("\n== bipolar %s cells, %s filter at %d genes: %d significant ==\n",
              J$ct, J$stat, J$n, length(sg)))
  cat(sprintf("  LM22 markers among hits : %d of %d (%.1f%%) vs %.1f%% background = %.1fx\n",
              sum(ks %in% LM22G), length(ks), 100*mean(ks %in% LM22G), 100*BG, mean(ks %in% LM22G)/BG))
  cat(sprintf("  direction               : %d down, %d up\n", sum(t1$logFC[sg]<0), sum(t1$logFC[sg]>0)))
  cat(sprintf("  survive adjustment for all 22 fractions: %d of %d (%.0f%%)\n",
              sum(surv), length(sg), 100*mean(surv)))
  cat(sprintf("  median %% of effect lost : %.1f%%\n", median(lost)))
  cat(sprintf("  strongest 12            : %s\n", paste(ks[order(q1[sg])][1:min(12,length(ks))], collapse=", ")))
  OUT <- rbind(OUT, data.frame(lineage=J$ct, filter=J$stat, n_kept=J$n, n_sig=length(sg),
    lm22_frac=mean(ks %in% LM22G), lm22_enrichment=mean(ks %in% LM22G)/BG,
    n_down=sum(t1$logFC[sg]<0), n_survive=sum(surv), pct_survive=100*mean(surv),
    median_pct_lost=median(lost), stringsAsFactors=FALSE)) }
write.csv(OUT, file.path(RES,"ct_bpd_survivors.csv"), row.names=FALSE)
cat("\nwrote results/ct_bpd_survivors.csv\n")
