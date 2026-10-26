#!/usr/bin/env Rscript
# ==============================================================================
# 10_circularity.R -- the circularity sensitivity analysis (PRESPEC.md section 7).
#
# THE QUESTION
#
# The cell fractions (the mediator) were estimated FROM the expression matrix,
# using the LM22 signature genes. For those genes the mediator is close to a
# function of the outcome, so mediation through cell mix is partly built in.
# How much does that matter?
#
# THE DESIGN -- ONE COMMON FIT
#
# TMM, library sizes, voom weights and the empirical-Bayes variance prior are
# all computed ONCE, on the full 12,368-gene filtered set. The analyses then
# differ only in which genes are summarised:
#
#   C1  non-marker genes                         12,018
#   C2  LM22 signature genes that are expressed     350
#   C3  everything, markers included             12,368
#   C4  every LM22 gene with any counts             539   (point estimates only)
#
# Holding the fit fixed is what makes this a clean sensitivity analysis. The
# earlier marker check (S1) normalised the markers on their own counts, so its
# library size was the sum of marker counts -- which itself tracks cell
# composition. That put composition into the denominator. S1 is superseded.
#
# Both mediator trees are run.
# ==============================================================================

ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "mediation")); set.seed(481)
suppressPackageStartupMessages({library(edgeR); library(limma)})
source("scripts/med.R")
options(width = 170)
say  <- function(...) cat(sprintf(...))
rule <- function(t) cat(sprintf("\n%s\n%s\n", t, strrep("=", nchar(t))))
dir.create("results/circularity", showWarnings = FALSE, recursive = TRUE)

prep <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
e2s  <- read.delim(file.path(ROOT, "de_analysis/data/ens2sym.tsv"), header = FALSE,
                   col.names = c("ens", "sym"), stringsAsFactors = FALSE)
sym  <- setNames(toupper(e2s$sym), e2s$ens)
LM   <- read.delim(file.path(ROOT, "data/LM22.txt"), check.names = FALSE, stringsAsFactors = FALSE)
lm22 <- unique(toupper(trimws(LM[[1]])))
LMm  <- as.matrix(LM[, -1]); rownames(LMm) <- toupper(LM[[1]])

## =============================================================================
## 1. Account for all 547 signature genes (the table in PRESPEC 7.1, regenerated)
## =============================================================================
rule("1. Where the 547 LM22 signature genes are")
raw <- read.delim(gzfile(file.path(ROOT, "data/raw_as_deposited/GSE124326_count_matrix.txt.gz")),
                  check.names = FALSE, row.names = 1)
raw <- raw[!grepl("^ENSGR", rownames(raw)), ]
rownames(raw) <- sub("[.].*$", "", rownames(raw))
colnames(raw) <- sub("[.]counts$", "", colnames(raw))
rsym <- sym[rownames(raw)]
filt <- rownames(prep$counts)
acc <- data.frame(symbol = lm22, stringsAsFactors = FALSE)
acc$in_annotation <- acc$symbol %in% rsym
acc$ensembl <- names(rsym)[match(acc$symbol, rsym)]
acc$passes_filter <- acc$ensembl %in% filt
acc$top_celltype <- colnames(LMm)[apply(LMm[acc$symbol, , drop = FALSE], 1, which.max)]
acc$specificity <- apply(LMm[acc$symbol, , drop = FALSE], 1, function(v) max(v) / sum(v))
rc <- as.matrix(raw[match(acc$ensembl, rownames(raw)), prep$meta$title])
acc$share_samples_gt10 <- ifelse(acc$in_annotation, rowMeans(rc > 10), NA)
acc$all_zero <- ifelse(acc$in_annotation, rowSums(rc) == 0, NA)
acc$status <- ifelse(!acc$in_annotation, "absent from annotation",
              ifelse(acc$passes_filter, "expressed (C2)",
              ifelse(acc$all_zero, "all-zero", "present, fails filter")))
print(as.data.frame(table(status = acc$status)), row.names = FALSE)
write.csv(acc, "results/circularity/lm22_accounting.csv", row.names = FALSE)

## =============================================================================
## 2. The common fit, per tree
## =============================================================================
VAR <- list(lineage = "data/med_input.rds", proposal = "data/med_input_proposal.rds")
ALL <- list()

wls_rss <- function(Y, W, X) {
  r <- numeric(nrow(Y))
  for (g in seq_len(nrow(Y))) {
    w <- W[g, ]; ch <- chol(crossprod(X, X * w))
    b <- backsolve(ch, backsolve(ch, crossprod(X, Y[g, ] * w), transpose = TRUE))
    e <- Y[g, ] - X %*% b; r[g] <- sum(w * e^2)
  }
  r
}

summarise_set <- function(tb, idx, label, it_p, pr2) {
  s <- tb[idx, ]
  qTE <- p.adjust(s$TE_p, "BH"); qNDE <- p.adjust(s$NDE_p, "BH")
  qNIE <- p.adjust(s$NIE_p, "BH"); qINT <- p.adjust(it_p[idx], "BH")
  sig <- qTE < 0.05
  fw <- c(CDE = mean(abs(s$CDE[sig])), INTref = mean(abs(s$INTref[sig])),
          INTmed = mean(abs(s$INTmed[sig])), PIE = mean(abs(s$PIE[sig])))
  data.frame(set = label, n_genes = length(idx),
             n_TE = sum(sig), n_NDE = sum(qNDE < 0.05), n_NIE_delta = sum(qNIE < 0.05),
             n_interaction = sum(qINT < 0.05),
             median_PM = if (any(sig)) median(s$PM[sig]) else NA_real_,
             PM_q25 = if (any(sig)) quantile(s$PM[sig], .25, names = FALSE) else NA_real_,
             PM_q75 = if (any(sig)) quantile(s$PM[sig], .75, names = FALSE) else NA_real_,
             share_PM_gt_half = if (any(sig)) mean(s$PM[sig] > 0.5) else NA_real_,
             median_abs_NIE = median(abs(s$NIE)), median_abs_TE = median(abs(s$TE)),
             pct_mediated_fourway = 100 * (fw["PIE"] + fw["INTmed"]) / sum(fw),
             median_partialR2_cellmix = median(pr2[idx]),
             stringsAsFactors = FALSE)
}

for (vn in names(VAR)) {
  rule(sprintf("2. COMMON FIT -- %s tree", vn))
  d <- readRDS(VAR[[vn]]); K <- ncol(d$M)
  stopifnot(setequal(c(d$gene_main, d$gene_mark), filt))
  CNT <- prep$counts[filt, d$meta$title, drop = FALSE]
  x <- d$lithium; cM <- centre(d$M); cC <- centre(d$C)
  Xd <- med_design_full(x, cM$X, cC$X)
  dge <- normLibSizes(DGEList(CNT), method = "TMM")          # ONE normalisation
  v   <- voom(dge, Xd)                                         # ONE set of weights
  mf  <- mediator_fit(cM$X, x, cC$X)
  of  <- outcome_fit(v$E, v$weights, Xd, se = TRUE, moderate = TRUE)  # ONE prior
  tb  <- med_table(med_effects(of, mf), rownames(v$E))
  tb$symbol <- sym[tb$gene]
  tb$set <- ifelse(tb$gene %in% d$gene_mark, "marker", "non-marker")
  it  <- interaction_test(v$E, v$weights, Xd, K)
  tb$INT_p <- it$p
  say("  %d genes | voom weights %.2f-%.2f | eBayes prior df %.1f\n",
      nrow(tb), min(v$weights), max(v$weights), of$df_prior)

  ## ---- circularity signature: how much of each gene does the cell mix explain?
  # partial R^2 of the balance block in a model without the interaction:
  #   full    Y ~ lithium + balances + covariates
  #   reduced Y ~ lithium + covariates
  Xf <- cbind(1, x, cM$X, cC$X); Xr <- cbind(1, x, cC$X)
  rss_f <- wls_rss(v$E, v$weights, Xf); rss_r <- wls_rss(v$E, v$weights, Xr)
  tb$partialR2_cellmix <- (rss_r - rss_f) / rss_r

  ## ---- the three sets, from one fit
  i1 <- which(tb$set == "non-marker"); i2 <- which(tb$set == "marker"); i3 <- seq_len(nrow(tb))
  S <- rbind(summarise_set(tb, i1, "C1 non-markers", it$p, tb$partialR2_cellmix),
             summarise_set(tb, i2, "C2 LM22 markers (expressed)", it$p, tb$partialR2_cellmix),
             summarise_set(tb, i3, "C3 all filtered genes", it$p, tb$partialR2_cellmix))
  S$tree <- vn
  cat("\n"); print(transform(S[, c("set", "n_genes", "n_TE", "n_NDE", "n_NIE_delta", "n_interaction",
                                   "median_PM", "share_PM_gt_half", "median_abs_NIE",
                                   "median_partialR2_cellmix")],
                             median_PM = round(median_PM, 4), share_PM_gt_half = round(share_PM_gt_half, 3),
                             median_abs_NIE = signif(median_abs_NIE, 3),
                             median_partialR2_cellmix = round(median_partialR2_cellmix, 4)),
                   row.names = FALSE)

  ## ---- test 1 of PRESPEC 7.4: is the mediator-outcome link stronger for markers?
  w <- wilcox.test(tb$partialR2_cellmix[i2], tb$partialR2_cellmix[i1])
  say("\n  CIRCULARITY SIGNATURE -- partial R^2 of the cell mix (markers vs non-markers)\n")
  say("    markers     median %.4f  IQR %.4f-%.4f\n", median(tb$partialR2_cellmix[i2]),
      quantile(tb$partialR2_cellmix[i2], .25), quantile(tb$partialR2_cellmix[i2], .75))
  say("    non-markers median %.4f  IQR %.4f-%.4f\n", median(tb$partialR2_cellmix[i1]),
      quantile(tb$partialR2_cellmix[i1], .25), quantile(tb$partialR2_cellmix[i1], .75))
  say("    ratio of medians %.2fx | Wilcoxon p = %.3g\n",
      median(tb$partialR2_cellmix[i2]) / median(tb$partialR2_cellmix[i1]), w$p.value)
  # among markers: are the most cell-type-specific signature genes the most explained?
  spec <- acc$specificity[match(tb$symbol[i2], acc$symbol)]
  ct <- cor.test(spec, tb$partialR2_cellmix[i2], method = "spearman", exact = FALSE)
  say("    among markers, LM22 specificity vs partial R^2: Spearman rho %.3f, p = %.3g\n",
      ct$estimate, ct$p.value)

  ## ---- markers judged within the FULL family as well as their own
  qfull <- p.adjust(tb$TE_p, "BH")
  say("\n  LM22 markers with a significant total effect: %d within their own family (350 tests),\n",
      sum(p.adjust(tb$TE_p[i2], "BH") < 0.05))
  say("    %d within the full 12,368-test family\n", sum(qfull[i2] < 0.05))

  ## ---- normalisation-basis check against the ORIGINAL primary (TMM on 12,018)
  orig <- readRDS(sprintf("results/primary/%s_fit.rds", vn))$tb
  m <- match(orig$gene, tb$gene)
  nb <- data.frame(tree = vn,
    corr_NIE = cor(orig$NIE, tb$NIE[m]), corr_TE = cor(orig$TE, tb$TE[m]),
    corr_NDE = cor(orig$NDE, tb$NDE[m]),
    max_abs_diff_TE = max(abs(orig$TE - tb$TE[m])),
    median_PM_original = median(orig$PM[orig$TE_q_BH < 0.05]),
    median_PM_commonfit = S$median_PM[1],
    n_TE_original = sum(orig$TE_q_BH < 0.05), n_TE_commonfit = S$n_TE[1],
    stringsAsFactors = FALSE)
  say("\n  NORMALISATION BASIS -- original primary (TMM on 12,018) vs common fit (TMM on 12,368)\n")
  say("    correlation of per-gene NIE %.5f | TE %.5f | NDE %.5f\n", nb$corr_NIE, nb$corr_TE, nb$corr_NDE)
  say("    median PM %.4f vs %.4f | genes with significant TE %d vs %d\n",
      nb$median_PM_original, nb$median_PM_commonfit, nb$n_TE_original, nb$n_TE_commonfit)

  ## ---- C4: every LM22 gene with any counts (adds the 189 unexpressed ones)
  ext <- acc$ensembl[acc$status == "present, fails filter"]
  CNT4 <- rbind(CNT, as.matrix(raw[ext, d$meta$title]))
  storage.mode(CNT4) <- "integer"
  dge4 <- DGEList(CNT4, lib.size = dge$samples$lib.size, norm.factors = dge$samples$norm.factors)
  v4  <- voom(dge4, Xd)
  of4 <- outcome_fit(v4$E, v4$weights, Xd, se = TRUE, moderate = TRUE)
  tb4 <- med_table(med_effects(of4, mf), rownames(v4$E))
  i4  <- which(tb4$gene %in% c(d$gene_mark, ext))
  it4 <- interaction_test(v4$E, v4$weights, Xd, K)
  pr4 <- rep(NA_real_, nrow(tb4))
  S4 <- summarise_set(tb4, i4, "C4 every LM22 gene with counts", it4$p, pr4)
  S4$median_partialR2_cellmix <- NA_real_; S4$tree <- vn
  say("\n  C4 -- %d LM22 genes (350 expressed + %d unexpressed): significant TE %d | median PM %.4f\n",
      length(i4), length(ext), S4$n_TE, S4$median_PM)
  unexpr_sig <- sum(p.adjust(tb4$TE_p[i4], "BH")[tb4$gene[i4] %in% ext] < 0.05)
  say("    of the %d unexpressed genes, %d reach a significant total effect\n", length(ext), unexpr_sig)
  S <- rbind(S, S4)

  ## ---- save
  write.csv(tb, sprintf("results/circularity/%s_common_fit_gene_effects.csv", vn), row.names = FALSE)
  write.csv(S,  sprintf("results/circularity/%s_summary.csv", vn), row.names = FALSE)
  saveRDS(list(tb = tb, S = S, nb = nb, mf = mf, sets = list(C1 = i1, C2 = i2, C3 = i3),
               sig = list(C1 = i1[p.adjust(tb$TE_p[i1], "BH") < 0.05],
                          C2 = i2[p.adjust(tb$TE_p[i2], "BH") < 0.05],
                          C3 = i3[p.adjust(tb$TE_p[i3], "BH") < 0.05]),
               wilcox_p = w$p.value, spec_rho = unname(ct$estimate), spec_p = ct$p.value,
               c4 = tb4[i4, ]),
          sprintf("results/circularity/%s_circularity.rds", vn))
  ALL[[vn]] <- list(S = S, nb = nb)
}

write.csv(do.call(rbind, lapply(ALL, `[[`, "nb")), "results/circularity/normalisation_basis.csv",
          row.names = FALSE)
write.csv(do.call(rbind, lapply(ALL, `[[`, "S")), "results/circularity/summary_both_trees.csv",
          row.names = FALSE)
cat("\n10_circularity complete.\n")
