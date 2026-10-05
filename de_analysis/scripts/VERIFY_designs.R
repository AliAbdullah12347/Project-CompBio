# Verify the covariate adjustment actually happened: print the exact design
# matrix used by every fit, confirm it is full rank, and confirm the tested
# coefficient is the group term and not something else.
suppressPackageStartupMessages({library(limma); library(edgeR)})
setwd("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/de_analysis")
p <- readRDS("data/prep.rds")
for (cn in c("LI","BPD")) for (ad in c("raw","ilr")) {
  s <- p$sel[[cn]]; mm <- p$meta[s$keep,,drop=FALSE]; mm$grp <- s$grp
  mm$plate <- droplevels(mm$plate); mm$sex <- droplevels(mm$sex)
  if ("group" %in% s$cov) mm$group <- droplevels(mm$group)
  cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]]))>1, logical(1))]
  dropped <- setdiff(s$cov, cv)
  X <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse=" + "))), data=mm)
  if (ad=="ilr") X <- cbind(X, p$ILR[mm$title,,drop=FALSE])
  ne <- nonEstimable(X); if (!is.null(ne)) X <- X[, setdiff(colnames(X), ne), drop=FALSE]
  cat(sprintf("\n=== WB_%s_%s ===\n", cn, ad))
  cat(sprintf("  samples %d | design %d x %d | rank %d | FULL RANK: %s\n",
      nrow(mm), nrow(X), ncol(X), qr(X)$rank, qr(X)$rank==ncol(X)))
  cat(sprintf("  tested coefficient: 'grpcase'  present: %s\n", "grpcase" %in% colnames(X)))
  cat(sprintf("  covariates requested : %s\n", paste(s$cov, collapse=", ")))
  cat(sprintf("  dropped as constant  : %s\n", if(length(dropped)) paste(dropped,collapse=", ") else "none"))
  cat(sprintf("  dropped as aliased   : %s\n", if(!is.null(ne)) paste(ne,collapse=", ") else "none"))
  cat(sprintf("  DESIGN COLUMNS (%d): %s\n", ncol(X), paste(colnames(X), collapse=" | ")))
}
# Does adjustment actually change anything? Fit LI with and without covariates.
cat("\n\n=== does the adjustment do any work? (contrast LI) ===\n")
s <- p$sel$LI; mm <- p$meta[s$keep,]; mm$grp <- s$grp
mm$plate<-droplevels(mm$plate); mm$sex<-droplevels(mm$sex)
dge <- calcNormFactors(DGEList(p$counts[,mm$title]), method="TMM")
for (lab in c("grp ONLY (no covariates)","grp + all covariates")) {
  f <- if (grepl("ONLY", lab)) ~ grp else ~ grp+age+sex+tobacco+rin+plate+seqpc1+seqpc2+seqpc3
  X <- model.matrix(f, data=mm)
  tt <- topTable(eBayes(lmFit(voom(dge,X),X)), coef="grpcase", number=Inf, sort.by="none")
  cat(sprintf("  %-26s cols=%2d  DEG=%5d  pi0=%.3f\n", lab, ncol(X),
      sum(tt$adj.P.Val<0.05), min(1,mean(tt$P.Value>.5)/.5)))
}
