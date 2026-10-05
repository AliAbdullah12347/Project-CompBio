# Why did my COMPOS survival rate disagree with the agent's?
# The spike I added, c * b1_centred, lies exactly in the column space of the
# ADJUSTED design (which contains b1). So it must be absorbed exactly, and the
# adjusted lithium coefficient must be identical with and without the spike.
# If that is true, my "survivors" are simply genes that were already DE after
# adjustment -- a contaminated denominator, not a failure of the adjustment.
source("scripts/common.R"); p <- load_prep(); set.seed(20260929)
m <- p$meta[p$meta$dx=="BP1",]; m$tob <- p$TOB[m$title,"tobacco_imp_01"]
cv <- c("age","sex","tob","rin","plate","seqpc1","seqpc2","seqpc3")
m <- m[complete.cases(m[,c("lithium",cv)]),]; m$plate<-droplevels(m$plate); m$sex<-droplevels(m$sex)
b1 <- p$ILR[m$title,"b1_myeloid_vs_lymphoid"]; b1c <- b1-mean(b1)
delta <- coef(lm(as.formula(paste("b1 ~ lithium +",paste(cv,collapse=" + "))),data=cbind(m,b1=b1)))["lithium"]
cn <- p$counts[,m$title,drop=FALSE]; cn <- cn[filter_genes(cn,10,0.90),,drop=FALSE]
E <- cpm(calcNormFactors(DGEList(cn),"TMM"),log=TRUE,prior.count=3)
Xm <- model.matrix(as.formula(paste("~ lithium +",paste(cv,collapse=" + "))),data=m)
Xa <- cbind(Xm, b1=b1)
spk <- sample(rownames(E),1200)
Y <- E; Y[spk,] <- Y[spk,] + matrix(rep((0.20/delta)*b1c, each=length(spk)), nrow=length(spk))
adj_clean  <- topTable(eBayes(lmFit(E,Xa),trend=TRUE),coef="lithium",number=Inf,sort.by="none")
adj_spiked <- topTable(eBayes(lmFit(Y,Xa),trend=TRUE),coef="lithium",number=Inf,sort.by="none")
d <- max(abs(adj_spiked[spk,"logFC"] - adj_clean[spk,"logFC"]))
cat(sprintf("max |adjusted logFC spiked - unspiked| over the 1200 spiked genes: %.3e\n", d))
cat(sprintf("  -> the compositional spike is absorbed %s by the adjustment\n",
            if (d < 1e-8) "EXACTLY" else "only partly"))
base_adj_deg <- intersect(spk, rownames(adj_clean)[adj_clean$adj.P.Val<0.05])
cat(sprintf("spiked genes already DE after adjustment BEFORE any spike: %d\n", length(base_adj_deg)))
marg <- topTable(eBayes(lmFit(Y,Xm),trend=TRUE),coef="lithium",number=Inf,sort.by="none")
det <- intersect(spk, rownames(marg)[marg$adj.P.Val<0.05])
sur <- intersect(det, rownames(adj_spiked)[adj_spiked$adj.P.Val<0.05])
cat(sprintf("detected marginally %d | survivors %d | of those already-DE-at-baseline %d (%.0f%%)\n",
            length(det), length(sur), length(intersect(sur,base_adj_deg)),
            100*length(intersect(sur,base_adj_deg))/max(length(sur),1)))
clean_spk <- setdiff(spk, base_adj_deg)
det2 <- intersect(clean_spk, rownames(marg)[marg$adj.P.Val<0.05])
sur2 <- intersect(det2, rownames(adj_spiked)[adj_spiked$adj.P.Val<0.05])
cat(sprintf("\nCORRECTED (spiking only genes NOT already DE after adjustment):\n"))
cat(sprintf("  detected %d | survive %d | survival rate %.3f  <- compare agent's 0.000\n",
            length(det2), length(sur2), length(sur2)/max(length(det2),1)))
