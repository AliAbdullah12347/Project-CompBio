#!/usr/bin/env Rscript
# ==============================================================================
# bestcut5_composition.R -- what are the 5 bipolar whole-blood genes?
#
# The 50-gene SD cut returned TMTC1, HDC, ALOX15, CPA3 and GATA2, and that set
# survived a selection-aware null (p = 0.015) and a leave-one-out check (no
# gene moves more than 15% when its most influential sample is removed).
#
# Four of the five are the same compartment. HDC makes histamine; CPA3 is the
# mast-cell carboxypeptidase; GATA2 is the transcription factor that specifies
# the basophil/mast/eosinophil lineage; ALOX15 is the basophil and eosinophil
# lipoxygenase. All four move DOWN in unmedicated bipolar I.
#
# A coordinated drop in four markers of one lineage is what a smaller lineage
# looks like, not what altered transcription looks like. That distinction is
# the entire point of this project, so it gets tested rather than asserted:
# are these LM22 markers, do they track a measured cell fraction, and does the
# effect survive adjusting for composition?
# ==============================================================================
suppressPackageStartupMessages({library(edgeR); library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "de_final")); set.seed(481); options(width = 190)
RES <- "results"; FDR <- 0.05
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

G <- read.csv(file.path(RES,"sdfilter_bestcut_genes.csv"), stringsAsFactors=FALSE)
gs <- G$gene; cat(sprintf("the %d genes: %s\n\n", length(gs), paste(G$symbol, collapse=", ")))

F22 <- as.matrix(FR[mm$title,]); F22 <- F22[,apply(F22,2,sd)>0,drop=FALSE]
Xa <- cbind(X, F22[,-ncol(F22),drop=FALSE])
q <- qr(Xa); if (q$rank < ncol(Xa)) Xa <- Xa[, q$pivot[seq_len(q$rank)], drop=FALSE]

Y <- E[gs,,drop=FALSE]; W <- VW[gs,,drop=FALSE]
t1 <- topTable(eBayes(lmFit(Y,X,weights=W)), coef="grpcase", number=Inf, sort.by="none")
t2 <- topTable(eBayes(lmFit(Y,Xa,weights=W)),coef="grpcase", number=Inf, sort.by="none")
q2 <- p.adjust(t2$P.Value,"BH")

cat("  gene     | LM22 | strongest fraction correlation          | logFC raw | logFC adj22 | % lost | q adj22 | verdict\n")
OUT <- NULL
for (i in seq_along(gs)) {
  sy <- toupper(G$symbol[i]); cors <- cor(E[gs[i],], F22); j <- which.max(abs(cors))
  lost <- 100*(1-abs(t2$logFC[i])/abs(t1$logFC[i]))
  cat(sprintf("  %-8s | %-4s | %-30s %+.2f | %+9.3f | %+11.3f | %6.1f | %7.4f | %s\n",
      sy, if (sy %in% LM22G) "YES" else "no", colnames(F22)[j], cors[j],
      t1$logFC[i], t2$logFC[i], lost, q2[i],
      if (q2[i] < FDR) "survives" else "LOST -- composition"))
  OUT <- rbind(OUT, data.frame(symbol=sy, lm22=sy %in% LM22G, top_fraction=colnames(F22)[j],
    cor=cors[j], logFC_raw=t1$logFC[i], logFC_adj22=t2$logFC[i], pct_lost=lost,
    q_adj22=q2[i], survives=q2[i]<FDR, stringsAsFactors=FALSE)) }
write.csv(OUT, file.path(RES,"bestcut5_composition.csv"), row.names=FALSE)

# Do the group differ in the basophil/mast compartment itself?
cat("\n--- the compartment these genes mark ---\n")
for (ctp in c("Mast cells resting","Mast cells activated","Eosinophils","Neutrophils")) {
  if (!ctp %in% colnames(FR)) next
  x <- FR[mm$title, ctp]; if (sd(x)==0) { cat(sprintf("  %-22s structurally zero\n", ctp)); next }
  tt <- t.test(x ~ mm$grp)
  cat(sprintf("  %-22s ctrl %.4f vs case %.4f | t-test p = %.4f\n",
      ctp, mean(x[mm$grp=="ctrl"]), mean(x[mm$grp=="case"]), tt$p.value)) }
cat(sprintf("\n  %d of %d survive full composition adjustment\n", sum(OUT$survives), nrow(OUT)))
cat("wrote results/bestcut5_composition.csv\n")
