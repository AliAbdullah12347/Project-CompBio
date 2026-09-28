#!/usr/bin/env Rscript
# Evaluate bMIND output the way Boltz et al. (2024) did -- Figure 1B, Figure S2, Table S3, and the
# bmind_de counts (p. 333) -- plus the controls their evaluation lacks.
#
# Paper: Fig 1B "correlations (measured by R2 of mean expression across samples) between the eight
# main cell types"; Fig S2 "Principal component analysis of cell type expression ... each individual
# is represented nine times"; Table S3 / Fig S3 "we compare the median TPM values for protein coding
# genes using both scRNA-Seq and computationally deconvoluted bulk RNA-Seq" (Schmiedel = DICE).
#
# Added controls (not in the paper): a high R2 between "bMIND monocytes" and DICE monocytes is only
# evidence of deconvolution if it beats (a) the raw bulk against the same DICE profile and (b) a
# *different* bMIND cell type against it. Every expressed gene is broadly similar across blood cell
# types, so both baselines are expected to be high; the paper reports neither.

root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
A <- readRDS(file.path(root, "data/processed/bmind_cts.rds"))$A   # gene x cell x sample, log2(TPM + 1)
s <- list(mean = apply(A, 1:2, mean), median = apply(A, 1:2, median))
x <- readRDS(file.path(root, "data/processed/log2tpm_pcs.rds"))
stopifnot(identical(dimnames(A)[[3]], colnames(x$logx)))
ann <- read.csv(file.path(root, "data/processed/gene_annotation.csv"), row.names = 1)
tgt <- with(read.csv(file.path(root, "data/reference/paper_targets.csv")), setNames(value, quantity))
strip <- function(v) sub("\\..*$", "", v)
genes <- rownames(s$mean)
cells <- colnames(s$mean)
r2 <- function(a, b) cor(a, b, use = "complete.obs")^2

# ---- Figure 1B: R2 of per-gene mean expression between contexts --------------------------------
ctx_mean <- cbind(s$mean, bulk = rowMeans(x$logx[genes, ]))
fig1b <- cor(ctx_mean)^2
cat("=== Figure 1B: R2 of per-gene mean expression between contexts ===\n")
print(round(fig1b, 3))
write.csv(fig1b, file.path(root, "results/fig1B_celltype_R2.csv"))
png(file.path(root, "results/fig1B_celltype_R2.png"), width = 900, height = 800, res = 120)
par(mar = c(10, 10, 2, 2))
image(1:ncol(fig1b), 1:ncol(fig1b), fig1b, col = hcl.colors(50, "YlOrRd", rev = TRUE), axes = FALSE,
      xlab = "", ylab = "", main = "R2 of mean expression across samples")
axis(1, 1:ncol(fig1b), colnames(fig1b), las = 2, cex.axis = 0.8)
axis(2, 1:ncol(fig1b), colnames(fig1b), las = 2, cex.axis = 0.8)
text(rep(1:ncol(fig1b), ncol(fig1b)), rep(1:ncol(fig1b), each = ncol(fig1b)), sprintf("%.2f", fig1b), cex = 0.6)
invisible(dev.off())

# ---- Figure S2: PCA of sample x context profiles -------------------------------------------------
g <- dimnames(A)[[1]]
prof <- do.call(rbind, c(lapply(cells, function(k) t(A[, k, ])), list(t(x$logx[g, ]))))
ctx <- factor(rep(c(cells, "bulk"), each = dim(A)[3]))
keep <- apply(prof, 2, var) > 0
pc <- prcomp(prof[, keep], center = TRUE, scale. = FALSE)
ve <- summary(pc)$importance[2, 1:2]
ctx_r2 <- sapply(1:2, function(j) summary(lm(pc$x[, j] ~ ctx))$r.squared)
cat(sprintf("\n=== Figure S2: PCA of %d profiles over %d genes ===\n", nrow(prof), sum(keep)))
cat(sprintf("PC1 %.1f%% / PC2 %.1f%% of variance; share of each explained by context: %.2f / %.2f\n",
            100 * ve[1], 100 * ve[2], ctx_r2[1], ctx_r2[2]))
png(file.path(root, "results/figS2_celltype_expression_PCA.png"), width = 1000, height = 800, res = 120)
pal <- setNames(c(hcl.colors(length(cells), "Dark 3"), "black"), levels(factor(c(cells, "bulk"))))
plot(pc$x[, 1], pc$x[, 2], pch = 16, cex = 0.4, col = pal[as.character(ctx)],
     xlab = sprintf("PC1 (%.1f%%)", 100 * ve[1]), ylab = sprintf("PC2 (%.1f%%)", 100 * ve[2]),
     main = "PCA of bMIND cell-type expression (each sample x 9 contexts)")
legend("topright", names(pal), col = pal, pch = 16, cex = 0.7, bty = "n")
invisible(dev.off())

# ---- Table S3: bMIND vs DICE (Schmiedel 2018) median TPM, protein-coding genes -------------------
dice_files <- c(monocytes = "MONOCYTES_TPM.csv.gz", t.cells.cd4.naive = "CD4_NAIVE_TPM.csv.gz",
                b.cells.naive = "B_CELL_NAIVE_TPM.csv.gz")
pcg <- genes[ann[genes, "biotype"] == "protein_coding"]
med_tpm <- 2^s$median - 1                         # median commutes with the monotone log transform
bulk_med <- setNames(2^apply(x$logx[pcg, ], 1, median) - 1, strip(pcg))
rownames(med_tpm) <- strip(rownames(med_tpm))
rows <- list()
for (k in names(dice_files)) {
  d <- read.csv(file.path(root, "data/external/DICE_Schmiedel2018", dice_files[[k]]), check.names = FALSE)
  dmed <- setNames(apply(as.matrix(d[, -(1:3)]), 1, median), strip(d[[1]]))
  common <- intersect(strip(pcg), names(dmed))
  lg <- function(v) log2(v + 1)
  add <- function(label, est) rows[[length(rows) + 1]] <<- data.frame(
    dice_cell_type = k, estimate = label, n_genes = length(common),
    R2_log = r2(lg(est[common]), lg(dmed[common])), R2_linear = r2(est[common], dmed[common]))
  if (k %in% cells) add(paste("bMIND", k), med_tpm[, k]) else add(paste("bMIND", k, "(not estimated)"), NA * med_tpm[, 1])
  add("bulk (no deconvolution)", bulk_med)
  for (o in setdiff(cells, k)) add(paste("bMIND", o, "(mismatched)"), med_tpm[, o])
}
s3 <- do.call(rbind, rows)
s3$boltz_R2 <- ifelse(grepl("^bMIND", s3$estimate) & !grepl("mismatched|not est", s3$estimate),
                      tgt[paste0("r2_bmind_schmiedel_", s3$dice_cell_type)], NA)
cat("\n=== Table S3: R2 with DICE median TPM (protein-coding) ===\n")
for (k in names(dice_files)) {
  z <- s3[s3$dice_cell_type == k, ]
  mm <- z[grepl("mismatched", z$estimate), ]
  cat(sprintf("%-26s matched %.3f | bulk %.3f | best mismatched %.3f (%s) | Boltz %.2f   [log scale, %d genes]\n",
              k, z$R2_log[1], z$R2_log[2], max(mm$R2_log), sub("bMIND (.*) \\(mismatched\\)", "\\1", mm$estimate[which.max(mm$R2_log)]),
              z$boltz_R2[1], z$n_genes[1]))
}
write.csv(s3, file.path(root, "results/tableS3_bmind_vs_DICE.csv"), row.names = FALSE)

# ---- bmind_de: counts vs Boltz, and do Boltz's cell-type DEGs replicate? -------------------------
# Boltz's cell-type names in Table S15 ('Bcellsnaive') vs LM22 names here.
map <- c(b.cells.naive = "Bcellsnaive", b.cells.memory = "Bcellsmemory", t.cells.cd8 = "TcellsCD8",
         t.cells.cd4.naive = "TcellsCD4naive", t.cells.cd4.memory.resting = "TcellsCD4memoryresting",
         nk.cells.resting = "NKcellsresting", monocytes = "Monocytes", neutrophils = "Neutrophils")
de_rows <- list()
for (tag in c("casecontrol", "lithium")) {
  q <- read.csv(file.path(root, sprintf("results/bmind_de_%s_qvals.csv", tag)), check.names = FALSE, row.names = 1)
  p <- read.csv(file.path(root, sprintf("results/bmind_de_%s_pvals.csv", tag)), check.names = FALSE, row.names = 1)
  bq <- read.csv(file.path(root, sprintf("data/reference/S15_bmind_diff_expr_qvals_%s.csv",
                                        if (tag == "lithium") "Li" else "case_cont")), check.names = FALSE)
  rownames(p) <- strip(rownames(p))
  for (k in colnames(q)) {
    bk <- map[k]
    b_sig <- if (!is.na(bk)) strip(bq$gene[bq[[bk]] < 0.05]) else character(0)
    rep_n <- if (length(b_sig)) sum(p[intersect(b_sig, rownames(p)), k] < 0.05) else NA
    de_rows[[length(de_rows) + 1]] <- data.frame(contrast = tag, cell_type = k, degs_here = sum(q[[k]] < 0.05),
      boltz_degs = if (!is.na(bk)) sum(bq[[bk]] < 0.05) else NA, boltz_degs_tested_here = length(intersect(b_sig, rownames(p))),
      boltz_degs_nominal_p05_here = rep_n)
  }
}
de <- do.call(rbind, de_rows)
cat("\n=== bmind_de: DEGs (q < 0.05) here vs Boltz ===\n")
print(de, row.names = FALSE)
write.csv(de, file.path(root, "results/bmind_de_summary.csv"), row.names = FALSE)
cat("\nWrote results/fig1B_*, figS2_*, tableS3_bmind_vs_DICE.csv, bmind_de_summary.csv\n")
