#!/usr/bin/env Rscript
# ==============================================================================
# meanfilter_bpd30.R -- the 30 bipolar whole-blood genes from the mean filter.
#
# meanfilter_validate.R found the one unambiguous win of the whole filtering
# exercise: filtering whole blood to the 200 most-expressed genes takes the
# bipolar contrast from 4 significant genes to 30, and that survives the
# selection-aware null (max over 9 cutoffs x 2 statistics, p = 0.0050).
#
# Unlike the SD results this is not suspect by construction -- mean expression
# does not preferentially select cell-identity markers, and the independence
# precondition holds. But three things still have to be checked before these
# are called findings, and they are the same three that the rest of the
# project applies to everything:
#
#   1. Are they composition? Re-fit with all 22 cell fractions in the design.
#   2. Are they LM22 markers, i.e. genes the deconvolution leaned on?
#   3. Does any single sample carry them? Leave-one-out on the effect size.
#
# The 4 genes found without any filter are also reported, to show whether the
# filtered set contains them or replaces them.
# ==============================================================================
suppressPackageStartupMessages({library(edgeR); library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "de_final")); set.seed(481); options(width = 195)
RES <- "results"; FDR <- 0.05; NTOP <- 200
p <- readRDS(file.path(ROOT,"de_analysis/data/prep.rds")); FR <- p$frac
sym <- local({m<-read.delim(file.path(ROOT,"de_analysis/data/ens2sym.tsv"),header=FALSE,
                col.names=c("ens","sym"));setNames(m$sym,m$ens)})
LM22G <- toupper(trimws(read.delim(file.path(ROOT,"data/LM22.txt"),check.names=FALSE)[[1]]))

s <- p$sel[["BPD"]]; mm <- p$meta[s$keep,]; mm$grp <- s$grp
mm$plate<-droplevels(factor(mm$plate)); mm$sex<-droplevels(factor(mm$sex))
if ("group" %in% s$cov) mm$group<-droplevels(factor(mm$group))
cv <- s$cov[vapply(s$cov,function(v)length(unique(mm[[v]]))>1,logical(1))]
X <- model.matrix(as.formula(paste("~ grp +",paste(cv,collapse=" + "))),data=mm)
ne<-nonEstimable(X); if(!is.null(ne)) X<-X[,setdiff(colnames(X),ne),drop=FALSE]
v <- voom(calcNormFactors(DGEList(p$counts[,mm$title]),method="TMM"), X)
E <- v$E; VW <- v$weights; dimnames(VW)<-dimnames(E)

unf <- topTable(eBayes(lmFit(E,X,weights=VW)),coef="grpcase",number=Inf,sort.by="none")
unf_sig <- rownames(E)[p.adjust(unf$P.Value,"BH") < FDR]

keep <- order(rowMeans(E), decreasing=TRUE)[seq_len(NTOP)]
Y <- E[keep,,drop=FALSE]; W <- VW[keep,,drop=FALSE]
t1 <- topTable(eBayes(lmFit(Y,X,weights=W)),coef="grpcase",number=Inf,sort.by="none")
q1 <- p.adjust(t1$P.Value,"BH"); gs <- rownames(Y)[q1 < FDR]
cat(sprintf("mean-expression cut at %d genes (min mean log-CPM %.2f)\n", NTOP, min(rowMeans(E)[keep])))
cat(sprintf("significant: %d  | unfiltered found %d (%s)\n",
    length(gs), length(unf_sig),
    paste(ifelse(is.na(sym[unf_sig]), unf_sig, sym[unf_sig]), collapse=", ")))
cat(sprintf("overlap between the two sets: %d\n\n", sum(unf_sig %in% gs)))

F22 <- as.matrix(FR[mm$title,]); F22 <- F22[,apply(F22,2,sd)>0,drop=FALSE]
Xa <- cbind(X, F22[,-ncol(F22),drop=FALSE])
qa <- qr(Xa); if (qa$rank < ncol(Xa)) Xa <- Xa[, qa$pivot[seq_len(qa$rank)], drop=FALSE]
t2 <- topTable(eBayes(lmFit(Y,Xa,weights=W)),coef="grpcase",number=Inf,sort.by="none")
q2 <- p.adjust(t2$P.Value,"BH")

idx <- match(gs, rownames(Y))
cat("  gene      | LM22 | logFC  |   q    | logFC adj22 | % lost | q adj22 | worst LOO (%) | sign | verdict\n")
OUT <- NULL
for (i in idx) {
  g <- rownames(Y)[i]; sy <- toupper(ifelse(is.na(sym[g]), g, sym[g]))
  y <- Y[i,]; wt <- W[i,]
  loo <- vapply(seq_along(y), function(k)
    unname(coef(lm(y[-k] ~ X[-k,-1,drop=FALSE], weights=wt[-k]))[2]), numeric(1))
  worst <- loo[which.max(abs(loo - t1$logFC[i]))]
  lost <- 100*(1-abs(t2$logFC[i])/abs(t1$logFC[i]))
  cat(sprintf("  %-9s | %-4s | %+.3f | %.4f | %+11.3f | %6.1f | %7.4f | %+13.0f | %-4s | %s\n",
      sy, if (sy %in% LM22G) "YES" else "no", t1$logFC[i], q1[i], t2$logFC[i], lost, q2[i],
      100*(worst-t1$logFC[i])/abs(t1$logFC[i]),
      if (sign(worst)==sign(t1$logFC[i])) "ok" else "FLIP",
      if (q2[i] < FDR) "survives" else "LOST -- composition"))
  OUT <- rbind(OUT, data.frame(gene=g, symbol=sy, lm22=sy %in% LM22G,
    logFC=t1$logFC[i], q=q1[i], logFC_adj22=t2$logFC[i], pct_lost=lost, q_adj22=q2[i],
    loo_pct=100*(worst-t1$logFC[i])/abs(t1$logFC[i]),
    sign_stable=sign(worst)==sign(t1$logFC[i]),
    survives=q2[i] < FDR, in_unfiltered=g %in% unf_sig, stringsAsFactors=FALSE)) }
OUT <- OUT[order(OUT$q), ]
write.csv(OUT, file.path(RES,"meanfilter_bpd30_genes.csv"), row.names=FALSE)
cat(sprintf("\n  %d of %d survive full composition adjustment | %d are LM22 markers (background 2.8%%) | %d sign-stable\n",
    sum(OUT$survives), nrow(OUT), sum(OUT$lm22), sum(OUT$sign_stable)))
cat(sprintf("  direction: %d up, %d down in unmedicated bipolar I\n",
    sum(OUT$logFC>0), sum(OUT$logFC<0)))
cat("wrote results/meanfilter_bpd30_genes.csv\n")
