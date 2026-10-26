#!/usr/bin/env Rscript
# ==============================================================================
# meanfilter_validate.R -- do any of the apparent filtering gains survive?
#
# WHAT NEEDS DEFENDING
#
# meanfilter.R swept nine cutoffs x two filter statistics across twelve
# analyses. Twenty of those cells returned more significant genes than no
# filtering at all, and some by a lot:
#
#   bipolar, T cells,  SD   at 1,000 genes :   0 -> 167
#   bipolar, T cells,  SD   at 2,000 genes :   0 -> 154
#   bipolar, B cells,  MEAN at   500 genes :   0 ->  86
#   bipolar, whole blood, MEAN at 200      :   4 ->  30
#
# Taken at face value these would be the biggest findings in the project. They
# are almost certainly not real, for two reasons that have both already been
# measured in this dataset:
#
#   1. SELECTION. Each number is the best of nine cutoffs x two statistics.
#      Eighteen chances per analysis, and a per-cell p-value prices in none of
#      it. The honest null is the distribution of the MAXIMUM count across the
#      same grid under permuted labels.
#
#   2. DEPENDENCE. diag_null_mechanism.R showed one sample-level axis explains
#      71% of variance in the cell-type estimates and any two genes correlate
#      at 0.674. Under that much dependence the number of BH discoveries is
#      wildly unstable: the permutation nulls already measured reached 105 and
#      168 genes in these very lineages. A count of 167 may be an ordinary
#      draw from that null.
#
# THE TEST
#
# One fit per permutation on the full gene set, then BH applied within each of
# the eighteen subsets. This is "filter after fitting" (what DESeq2 does), and
# meanfilter.R showed it tracks the refit version to within a gene or two,
# while costing eighteen times less -- which is what makes 500 permutations
# per analysis affordable.
#
# The independence diagnostic is also redone properly here. meanfilter.R
# estimated it from a single permutation, which is far too noisy to support a
# claim; this averages over 50.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "de_final")); set.seed(481); options(width = 200)
RES <- "results"
FDR   <- 0.05
CUTS  <- c(200, 500, 1000, 2000, 4000, 6000, 8000, 10000, 12368)
NPERM <- 200
NRHO  <- 50

p  <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
bp <- readRDS(file.path(ROOT, "de_analysis/results/bmind_profiles.rds"))
A <- bp$A; SE <- bp$SE
sym <- local({ m <- read.delim(file.path(ROOT, "de_analysis/data/ens2sym.tsv"),
                 header = FALSE, col.names = c("ens","sym")); setNames(m$sym, m$ens) })

design_for <- function(cn) {
  s <- p$sel[[cn]]; mm <- p$meta[s$keep, , drop = FALSE]; mm$grp <- s$grp
  mm$plate <- droplevels(factor(mm$plate)); mm$sex <- droplevels(factor(mm$sex))
  if ("group" %in% s$cov) mm$group <- droplevels(factor(mm$group))
  cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]])) > 1, logical(1))]
  X <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse = " + "))), data = mm)
  ne <- nonEstimable(X); if (!is.null(ne)) X <- X[, setdiff(colnames(X), ne), drop = FALSE]
  list(X = X, samples = mm$title) }
mkw <- function(S) { W <- 1/(S^2); f <- is.finite(W) & W > 0
  if (any(!f)) W[!f] <- min(W[f]); W/mean(W) }
pv <- function(Y, X, W, tr) {
  tt <- topTable(eBayes(lmFit(Y, X, weights = W), trend = tr),
                 coef = "grpcase", number = Inf, sort.by = "none")
  setNames(tt$P.Value, rownames(tt)) }

# counts at every (filter, cutoff) from one vector of p-values
grid_counts <- function(pvec, ORD) {
  vapply(names(ORD), function(FS)
    vapply(CUTS, function(m) sum(p.adjust(pvec[ORD[[FS]][seq_len(min(m, length(pvec)))]], "BH") < FDR),
           integer(1)), integer(length(CUTS))) }

# Only the analyses where some cutoff actually beat no filtering need a null --
# there is nothing to defend where filtering already lost. WB_LI is kept as the
# reference case because it is the headline analysis and its independence
# diagnostic is worth having. The six others are every unit from
# meanfilter_curve.csv with at least one winning cell.
UNITS <- rbind(
  data.frame(level="WB", contrast="LI",   lineage=NA_character_, stringsAsFactors=FALSE),
  data.frame(level="WB", contrast="BPD",  lineage=NA_character_, stringsAsFactors=FALSE),
  data.frame(level="CT", contrast="LI",   lineage="mono",  stringsAsFactors=FALSE),
  data.frame(level="CT", contrast="LI",   lineage="T",     stringsAsFactors=FALSE),
  data.frame(level="CT", contrast="BPD",  lineage="mono",  stringsAsFactors=FALSE),
  data.frame(level="CT", contrast="BPD",  lineage="T",     stringsAsFactors=FALSE),
  data.frame(level="CT", contrast="BPD",  lineage="B",     stringsAsFactors=FALSE))

cat(sprintf("== meanfilter_validate ==\n%d permutations | %d cutoffs x 2 filters = %d chances per analysis\n",
            NPERM, length(CUTS), 2*length(CUTS)))

OUT <- list(); RHO <- list()
for (u in seq_len(nrow(UNITS))) {
  U <- UNITS[u, ]; d <- design_for(U$contrast)
  if (U$level == "WB") {
    v <- voom(calcNormFactors(DGEList(p$counts[, d$samples]), method = "TMM"), d$X)
    Y <- v$E; W <- v$weights; dimnames(W) <- dimnames(Y); tr <- FALSE
  } else {
    Y <- A[, U$lineage, d$samples]; W <- mkw(SE[, U$lineage, d$samples]); tr <- TRUE }
  ORD <- list(SD = rownames(Y)[order(apply(Y,1,sd), decreasing=TRUE)],
              MEAN = rownames(Y)[order(rowMeans(Y), decreasing=TRUE)])
  stat <- list(SD = apply(Y,1,sd), MEAN = rowMeans(Y))

  p0  <- pv(Y, d$X, W, tr)
  oc  <- grid_counts(p0, ORD)                      # 9 cutoffs x 2 filters
  unf <- sum(p.adjust(p0, "BH") < FDR)
  obs_max <- apply(oc, 2, max)
  lab <- sprintf("%s %s%s", U$level, U$contrast, if (is.na(U$lineage)) "" else paste0(" ", U$lineage))

  NM <- matrix(0L, NPERM, 2, dimnames = list(NULL, c("SD","MEAN")))
  RH <- matrix(NA_real_, NRHO, 2, dimnames = list(NULL, c("SD","MEAN")))
  t0 <- Sys.time()
  for (b in seq_len(NPERM)) {
    Xp <- d$X; Xp[, "grpcase"] <- Xp[sample(nrow(Xp)), "grpcase"]
    pb <- pv(Y, Xp, W, tr)
    NM[b, ] <- apply(grid_counts(pb, ORD), 2, max)
    if (b <= NRHO) for (FS in c("SD","MEAN"))
      RH[b, FS] <- cor(stat[[FS]], -log10(pb), method = "spearman") }

  for (FS in c("SD","MEAN")) {
    pmax <- (sum(NM[, FS] >= obs_max[FS]) + 1)/(NPERM + 1)
    best <- CUTS[which.max(oc[, FS])]
    OUT[[length(OUT)+1]] <- data.frame(level=U$level, contrast=U$contrast,
      lineage=ifelse(is.na(U$lineage), "", U$lineage), filter=FS, unfiltered=unf,
      best_cut=best, observed_max=obs_max[FS],
      null_median=median(NM[,FS]), null_q95=quantile(NM[,FS],.95,names=FALSE),
      null_max=max(NM[,FS]), p_max_over_grid=pmax,
      rho_null_mean=mean(RH[,FS]), rho_null_sd=sd(RH[,FS]),
      beats_unfiltered=obs_max[FS] > unf, stringsAsFactors=FALSE) }
  cat(sprintf("  %-14s unfiltered %4d | SD best %4d @%5d (p %.4f) | MEAN best %4d @%5d (p %.4f) | %.1f min\n",
      lab, unf, obs_max["SD"], CUTS[which.max(oc[,"SD"])],
      (sum(NM[,"SD"] >= obs_max["SD"])+1)/(NPERM+1),
      obs_max["MEAN"], CUTS[which.max(oc[,"MEAN"])],
      (sum(NM[,"MEAN"] >= obs_max["MEAN"])+1)/(NPERM+1),
      as.numeric(difftime(Sys.time(), t0, units="mins"))))
  flush.console() }

R <- do.call(rbind, OUT)
write.csv(R, file.path(RES, "meanfilter_validation.csv"), row.names = FALSE)

cat(sprintf("\n%s\nSELECTION-AWARE RESULT (null = max over the same 9 cutoffs)\n%s\n",
            strrep("=",200), strrep("=",200)))
R$verdict <- ifelse(!R$beats_unfiltered, "filtering does not beat no filter",
             ifelse(R$p_max_over_grid < 0.05, "SURVIVES -- a real gain",
                    "does NOT survive: within the null once the cutoff choice is priced in"))
print(R[, c("level","contrast","lineage","filter","unfiltered","best_cut","observed_max",
            "null_q95","null_max","p_max_over_grid","verdict")], row.names = FALSE)

cat(sprintf("\n%s\nINDEPENDENCE OF THE FILTER STATISTIC (mean over %d permutations)\n%s\n",
            strrep("=",200), NRHO, strrep("=",200)))
ag <- aggregate(cbind(rho_null_mean, rho_null_sd) ~ level + filter, data = R, FUN = mean)
print(ag, row.names = FALSE, digits = 3)
cat("  both statistics should sit at ~0; a non-zero value means the filter peeks at the answer\n")
cat("\nwrote results/meanfilter_validation.csv\n")
