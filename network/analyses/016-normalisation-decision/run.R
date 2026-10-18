# DECISION TEST: does removing the per-sample library offset change the network?
# Bounded on purpose - 4000 most-variable genes, 1 thread - so it does not starve
# the base-001 preservation run sharing these 4 cores.
suppressMessages(library(WGCNA))
WGCNAnThreads <- 1; options(stringsAsFactors = FALSE)
source("network/00_load.R")
OUT <- "network/analyses/016-normalisation-decision/results"
say <- function(...) cat(sprintf(...), "\n")
set.seed(42)

d <- load_project(cache = TRUE)
E <- expr_for(d, group = "control", transpose = FALSE)      # genes x samples
v <- apply(E, 1, var); keep <- order(v, decreasing = TRUE)[1:4000]
raw <- E[keep, , drop = FALSE]
adj <- sweep(raw, 2, colMeans(E), "-")                      # remove per-sample offset
say("genes=%d samples=%d", nrow(raw), ncol(raw))

summarise <- function(M, tag) {
  X <- t(M)                                                  # samples x genes
  C <- cor(t(M[sample(nrow(M), 1200), ])); diag(C) <- NA   # t(): gene-gene, not sample-sample
  dens  <- mean(C > 0.5, na.rm = TRUE)
  meanr <- mean(C, na.rm = TRUE)
  sft <- pickSoftThreshold(X, powerVector = c(4,6,8,10,12,14,16,18,20),
                           networkType = "signed", verbose = 0)
  r2 <- -sign(sft$fitIndices[,3]) * sft$fitIndices[,2]
  best <- max(r2, na.rm = TRUE); bestp <- sft$fitIndices[which.max(r2), 1]
  net <- blockwiseModules(X, power = 12, networkType = "signed", corType = "bicor",
                          deepSplit = 4, minModuleSize = 30, mergeCutHeight = 0,
                          numericLabels = TRUE, maxBlockSize = 4500, verbose = 0)
  nmod <- length(setdiff(unique(net$colors), 0))
  grey <- mean(net$colors == 0)
  me <- moduleEigengenes(X, colors = labels2colors(net$colors))$eigengenes
  me <- me[, setdiff(names(me), "MEgrey"), drop = FALSE]
  R  <- abs(cor(me)); diag(R) <- NA
  medr <- median(R, na.rm = TRUE)
  ev <- eigen(cor(me), only.values = TRUE)$values
  meff <- (sum(ev))^2 / sum(ev^2)                            # Li & Ji style
  lg <- max(table(net$colors[net$colors != 0])) / sum(net$colors != 0)
  say("")
  say("### %s", tag)
  say("  mean gene-gene r            %+.3f", meanr)
  say("  %% random pairs r>0.5        %.1f%%", 100*dens)
  say("  best scale-free R2          %.3f  (at power %d)", best, bestp)
  say("  modules detected            %d", nmod)
  say("  %% genes grey                %.1f%%", 100*grey)
  say("  largest module share        %.1f%%", 100*lg)
  say("  median |r| between MEs      %.3f", medr)
  say("  effective dimensions        %.2f of %d", meff, ncol(me))
  data.frame(set=tag, mean_r=meanr, pct_gt_half=dens, sft_r2=best, sft_power=bestp,
             n_modules=nmod, pct_grey=grey, largest_share=lg, median_me_r=medr, meff=meff)
}

a <- summarise(raw, "RAW log2(TPM+1) -- what the arm used")
b <- summarise(adj, "PER-SAMPLE OFFSET REMOVED")
res <- rbind(a, b)
write.csv(res, file.path(OUT, "normalisation_comparison.csv"), row.names = FALSE)
say(""); say("=== VERDICT INPUTS ==="); print(res, row.names = FALSE)
say(""); say("ANALYSIS_016_DONE")
