#!/usr/bin/env Rscript
# ==============================================================================
# method_sweep_timing.R -- how the five estimators scale with the number of
# genes, on the real design.
#
# Written because one voomWQW fit on the full 12,173 genes ran for well over
# ten CPU minutes and it was not obvious whether that was the honest cost of
# the estimator or a pathology. Measuring the cost at 500/1000/2000/4000 genes
# answers that: a linear trend extrapolating to the observed full-matrix time
# means the cost is real and is simply the price of
# `arrayWeights(method = "genebygene", maxiter = 50)` looping over genes in R.
#
# The gene subsets are random draws from the genes that survive the real
# filter, so each subset has the same count distribution as the real matrix and
# the timing is not flattered by using only high-expression genes. This script
# measures runtime only -- none of its p-values are used anywhere.
#
# Writes only inside experimentation/.
# ==============================================================================

source("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation/scripts/common.R")
OUT <- file.path(EXP_HOME, "runs", "method_sweep")
log <- function(...) { cat(sprintf(...), file = stderr()); flush(stderr()) }

p   <- load_prep()
d   <- build_design(p, "lithium")
cnt <- p$counts[, d$samples, drop = FALSE]
g   <- filter_genes(cnt, 10, 0.90)
log("full filtered matrix: %d genes x %d samples\n", length(g), ncol(cnt))

set.seed(20261218)
SIZES   <- c(500, 1000, 2000, 4000)
METHODS <- c("voom", "voomWQW", "trend", "QLF", "LRT")

rows <- list(); k <- 0
for (nsub in SIZES) {
  sub <- sample(g, nsub)
  dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = cnt[sub, , drop = FALSE]),
                                method = "TMM")
  for (mth in METHODS) {
    t0 <- proc.time()
    if (mth == "voom") {
      limma::eBayes(limma::lmFit(limma::voom(dge, d$X), d$X))
    } else if (mth == "voomWQW") {
      limma::eBayes(limma::lmFit(limma::voomWithQualityWeights(dge, d$X, plot = FALSE), d$X))
    } else if (mth == "trend") {
      limma::eBayes(limma::lmFit(edgeR::cpm(dge, log = TRUE, prior.count = 3), d$X),
                    trend = TRUE)
    } else if (mth == "QLF") {
      dd <- edgeR::estimateDisp(dge, d$X)
      edgeR::glmQLFTest(edgeR::glmQLFit(dd, d$X), coef = d$coef)
    } else {
      dd <- edgeR::estimateDisp(dge, d$X)
      edgeR::glmLRT(edgeR::glmFit(dd, d$X), coef = d$coef)
    }
    el <- proc.time() - t0
    k <- k + 1
    rows[[k]] <- data.frame(n_genes = nsub, method = mth,
                            cpu = round(el[["user.self"]] + el[["sys.self"]], 2),
                            wall = round(el[["elapsed"]], 2), stringsAsFactors = FALSE)
    log("  %5d genes  %-8s cpu %7.2fs\n", nsub, mth, rows[[k]]$cpu)
  }
}
tm <- do.call(rbind, rows)
write.csv(tm, file.path(OUT, "scaling_raw.csv"), row.names = FALSE)

# Fit cpu = a * n_genes through the origin for each method and extrapolate to
# the full 12,173. A ratio near 1 against the measured full-matrix time means
# the cost is linear in genes and therefore not a pathology.
ext <- do.call(rbind, lapply(split(tm, tm$method), function(x) {
  a <- coef(lm(cpu ~ 0 + n_genes, data = x))[[1]]
  data.frame(method = x$method[1],
             cpu_per_1000_genes = round(a * 1000, 3),
             predicted_cpu_at_12173 = round(a * length(g), 1),
             r2_through_origin = round(1 - sum(resid(lm(cpu ~ 0 + n_genes, data = x))^2) /
                                         sum(x$cpu^2), 4),
             stringsAsFactors = FALSE)
}))
ext <- ext[order(-ext$predicted_cpu_at_12173), ]
log("\n== scaling, extrapolated to the full %d genes ==\n", length(g))
print(ext, row.names = FALSE)
write.csv(ext, file.path(OUT, "scaling_extrapolated.csv"), row.names = FALSE)
log("\n[timing] done\n")
