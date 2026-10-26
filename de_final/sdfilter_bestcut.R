#!/usr/bin/env Rscript
# ==============================================================================
# sdfilter_bestcut.R -- are the 5 bipolar whole-blood genes at the 50-gene cut
# real, or the product of having tried several cutoffs?
#
# The SD sweep found, for bipolar whole blood: 5 significant genes at a 50-gene
# cutoff, against 1, 0, 0, 0, 0 and 2 at the 100, 200, 500, 1,000, 2,000 and
# 5,000 cutoffs, and 4 with no filter at all. Reporting the 5 means reporting
# the best of seven tries, and a per-cutoff p-value does not price that in.
#
# The honest null is therefore the distribution of the MAXIMUM count across the
# same seven filtered cutoffs under permuted labels -- each permutation gets
# exactly the same seven chances we gave ourselves.
#
# The unfiltered 12,368-gene fit is excluded from the grid: it is the baseline
# being compared against, not one of the filter choices, and including it would
# make every permutation 60x more expensive for no change in the comparison.
#
# Anything that survives is then checked for dependence on single samples,
# because a 50-gene cut keeps the 50 most variable genes in the cohort, and
# high variance is what outliers look like.
# ==============================================================================
suppressPackageStartupMessages({library(edgeR); library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "de_final")); set.seed(481); options(width = 190)
RES <- "results"; FDR <- 0.05
CUTS <- c(50, 100, 200, 500, 1000, 2000, 5000); NPERM <- 1000

p <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
sym <- local({ m <- read.delim(file.path(ROOT,"de_analysis/data/ens2sym.tsv"),
                 header=FALSE, col.names=c("ens","sym")); setNames(m$sym, m$ens) })
design_for <- function(cn) {
  s <- p$sel[[cn]]; mm <- p$meta[s$keep, , drop=FALSE]; mm$grp <- s$grp
  mm$plate <- droplevels(factor(mm$plate)); mm$sex <- droplevels(factor(mm$sex))
  if ("group" %in% s$cov) mm$group <- droplevels(factor(mm$group))
  cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]]))>1, logical(1))]
  X <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse=" + "))), data=mm)
  ne <- nonEstimable(X); if (!is.null(ne)) X <- X[, setdiff(colnames(X), ne), drop=FALSE]
  list(X=X, samples=mm$title, grp=mm$grp) }
fitq <- function(Y, X, W) {
  tt <- topTable(eBayes(lmFit(Y, X, weights=W)), coef="grpcase", number=Inf, sort.by="none")
  list(q=p.adjust(tt$P.Value,"BH"), lfc=tt$logFC, g=rownames(tt)) }

CN <- "BPD"
d <- design_for(CN)
v <- voom(calcNormFactors(DGEList(p$counts[, d$samples]), method="TMM"), d$X)
E <- v$E; VW <- v$weights; dimnames(VW) <- dimnames(E)
sdg <- apply(E,1,sd); ord <- order(sdg, decreasing=TRUE)
KEEP <- lapply(CUTS, function(m) ord[seq_len(m)])

unf <- sum(fitq(E, d$X, VW)$q < FDR)
obs <- vapply(KEEP, function(k) sum(fitq(E[k,,drop=FALSE], d$X, VW[k,,drop=FALSE])$q < FDR), integer(1))
bi <- which.max(obs)
cat(sprintf("== bipolar whole blood ==\n  unfiltered: %d\n  per cutoff: %s\n  best: %d genes at the %d-gene cutoff\n",
            unf, paste(sprintf("%d:%d", CUTS, obs), collapse="  "), obs[bi], CUTS[bi]))

NMAX <- NAT <- integer(NPERM); t0 <- Sys.time()
for (b in seq_len(NPERM)) {
  Xp <- d$X; Xp[,"grpcase"] <- Xp[sample(nrow(Xp)),"grpcase"]
  nb <- vapply(KEEP, function(k) sum(fitq(E[k,,drop=FALSE], Xp, VW[k,,drop=FALSE])$q < FDR), integer(1))
  NMAX[b] <- max(nb); NAT[b] <- nb[bi]
  if (b %% 250 == 0) cat(sprintf("    %4d/%d  %.1f min\n", b, NPERM,
                                 as.numeric(difftime(Sys.time(), t0, units="mins")))) }
p_at  <- (sum(NAT  >= obs[bi]) + 1)/(NPERM+1)
p_max <- (sum(NMAX >= obs[bi]) + 1)/(NPERM+1)
cat(sprintf("\n  null AT the %d-gene cutoff : median %2.0f, q95 %5.1f, max %4d -> p = %.4f\n",
            CUTS[bi], median(NAT), quantile(NAT,.95,names=FALSE), max(NAT), p_at))
cat(sprintf("  null MAX over all %d cutoffs: median %2.0f, q95 %5.1f, max %4d -> p = %.4f  <-- the honest one\n",
            length(CUTS), median(NMAX), quantile(NMAX,.95,names=FALSE), max(NMAX), p_max))

k <- KEEP[[bi]]; r <- fitq(E[k,,drop=FALSE], d$X, VW[k,,drop=FALSE])
sig <- r$g[r$q < FDR]
GEN <- NULL
if (length(sig)) {
  cat(sprintf("\n  the %d genes, and how far one sample can move each estimate:\n", length(sig)))
  cat("   gene        logFC     q    | max |z| | worst leave-one-out logFC | sign stable\n")
  for (gg in sig) {
    y <- E[gg,]; wt <- VW[gg,]; lf <- r$lfc[match(gg, r$g)]
    loo <- vapply(seq_along(y), function(i)
      unname(coef(lm(y[-i] ~ d$X[-i,-1,drop=FALSE], weights=wt[-i]))[2]), numeric(1))
    worst <- loo[which.max(abs(loo - lf))]
    cat(sprintf("   %-11s %+.3f  %.4f |  %5.2f  | %+.3f (%+.0f%%)            | %s\n",
      ifelse(is.na(sym[gg]), gg, sym[gg]), lf, r$q[match(gg,r$g)], max(abs(scale(y))),
      worst, 100*(worst-lf)/abs(lf), if (sign(worst)==sign(lf)) "yes" else "NO"))
    GEN <- rbind(GEN, data.frame(gene=gg, symbol=unname(ifelse(is.na(sym[gg]),gg,sym[gg])),
      logFC=lf, q=r$q[match(gg,r$g)], sd=sdg[gg], max_abs_z=max(abs(scale(y))),
      loo_worst=worst, sign_stable=sign(worst)==sign(lf), stringsAsFactors=FALSE)) } }

write.csv(data.frame(contrast=CN, unfiltered=unf, best_cut=CUTS[bi], observed=obs[bi],
  null_at_q95=quantile(NAT,.95,names=FALSE), null_at_max=max(NAT), p_at_cut=p_at,
  null_max_q95=quantile(NMAX,.95,names=FALSE), null_max_max=max(NMAX), p_max_over_cuts=p_max),
  file.path(RES,"sdfilter_bestcut_validation.csv"), row.names=FALSE)
if (!is.null(GEN)) write.csv(GEN, file.path(RES,"sdfilter_bestcut_genes.csv"), row.names=FALSE)
cat(sprintf("\n=== verdict: %s ===\n",
  if (p_max < 0.05) "survives the selection-aware null" else
  "DOES NOT survive once the choice of cutoff is priced in"))
cat("wrote results/sdfilter_bestcut_*.csv\n")
