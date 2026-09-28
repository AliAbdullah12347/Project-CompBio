#!/usr/bin/env Rscript
# Bulk differential expression, Boltz et al. (2024) p. 326 and p. 332-333 / Figure 4C / Figure S8.
#
# Paper: "We used the limma eBayes function with trend = true ... We include only those genes with at
# least 1 TPM in at least 436 individuals (about 25% ...), leaving 17,194 genes to be tested. We then
# log2-transformed this matrix and computed the first 50 expression principal components to be
# included as covariates. In the lithium user vs. non-user analysis, we included only diagnosed
# individuals ..., whereas in the case-control analysis, all individuals diagnosed with BP or SCZ
# were included as cases." And (p. 333): "we tested the inclusion of the cell-type proportions as
# covariates (in addition to the 50 expression PCs) in the bulk lithium differential-expression test."
#
# Unstated details, fixed here and recorded in REPRODUCIBILITY.md: log2(TPM + 1) pseudocount;
# PCs from prcomp on genes centred (not scaled) across all 444 samples, computed once, as the
# sentence order implies; no other covariates (the paper lists none for DE).

suppressPackageStartupMessages(library(limma))
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
krebs <- file.path(root, "..", "Krebs Paper")
N_PC <- 50

meta <- read.csv(file.path(root, "data/processed/sample_metadata.csv"), check.names = FALSE,
                 colClasses = c(title = "character"))
tpm <- read.delim(file.path(root, "data/processed/tpm_filtered.tsv.gz"), row.names = 1, check.names = FALSE)
stopifnot(identical(colnames(tpm), meta$title))  # joined by key in script 00; confirm, don't assume
cells <- readLines(file.path(root, "data/processed/celltypes_selected.txt"))
ann <- read.csv(file.path(root, "data/processed/gene_annotation.csv"), row.names = 1)
tgt <- with(read.csv(file.path(root, "data/reference/paper_targets.csv")), setNames(value, quantity))
strip <- function(x) sub("\\..*$", "", x)

logx <- log2(as.matrix(tpm) + 1)
cat(sprintf("Genes tested: %d (Boltz: %d)\n", nrow(logx), tgt["n_genes_tested"]))
pcs <- prcomp(t(logx), center = TRUE, scale. = FALSE)$x[, 1:N_PC]
colnames(pcs) <- paste0("PC", 1:N_PC)
saveRDS(list(logx = logx, pcs = pcs), file.path(root, "data/processed/log2tpm_pcs.rds"))

run_de <- function(rows, coef_name, y, extra = NULL) {
  d <- data.frame(y = y[rows], pcs[rows, ], check.names = FALSE)
  if (!is.null(extra)) d <- cbind(d, meta[rows, extra, drop = FALSE])
  design <- model.matrix(~ ., data = d)
  colnames(design)[2] <- coef_name
  stopifnot(qr(design)$rank == ncol(design))
  fit <- eBayes(lmFit(logx[, rows], design), trend = TRUE)
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  tt$symbol <- ann[tt$gene, "symbol"]
  tt[, c("gene", "symbol", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
}
summ <- function(label, tt, n_boltz) {
  s <- tt[tt$adj.P.Val < 0.05, ]
  cat(sprintf("%-34s DEGs (FDR<0.05): %4d   Boltz: %4d   logFC range [%.3f, %.3f]\n", label, nrow(s),
              n_boltz, if (nrow(s)) min(s$logFC) else NA, if (nrow(s)) max(s$logFC) else NA))
}

cases <- meta$case == 1
cat(sprintf("Lithium contrast: %d users vs %d non-users (cases only)\n", sum(meta$lithium[cases] == 1),
            sum(meta$lithium[cases] == 0)))
de_li  <- run_de(cases, "lithium", meta$lithium)
de_cc  <- run_de(rep(TRUE, nrow(meta)), "case", meta$case)
de_adj <- run_de(cases, "lithium", meta$lithium, extra = cells)

cat("\n=== Bulk DE ===\n")
summ("Lithium users vs non-users", de_li, tgt["n_deg_lithium"])
summ("Case vs control", de_cc, tgt["n_deg_casecontrol"])
summ("Lithium, + cell proportions", de_adj, tgt["n_deg_lithium_celladj"])
sig <- function(tt) tt$gene[tt$adj.P.Val < 0.05]
cat(sprintf("Lithium & case/control overlap: %d (Boltz: %d)\n", length(intersect(sig(de_li), sig(de_cc))),
            tgt["n_deg_overlap_cc_li"]))
cat(sprintf("Cell-adjusted lithium DEGs also in unadjusted: %d of %d (Boltz: %d of %d)\n",
            length(intersect(sig(de_adj), sig(de_li))), length(sig(de_adj)),
            tgt["n_deg_lithium_celladj_shared"], tgt["n_deg_lithium_celladj"]))

# Sensitivity (not in Boltz): in GSE124326 every case is in assessment group A while controls span
# A and B, so group is a collection batch confounded with diagnosis. Boltz's PC-only model cannot
# separate the two; adjusting for group, or restricting to group A, shows how much of the
# case/control signal is batch. (Boltz's own cohort may not share this structure.)
de_cc_grp <- run_de(rep(TRUE, nrow(meta)), "case", meta$case, extra = "group")
de_cc_A <- run_de(meta$group == "A", "case", meta$case)
cat("\n--- Sensitivity: case/control vs assessment-group batch ---\n")
summ("Case vs control, + group", de_cc_grp, tgt["n_deg_casecontrol"])
summ("Case vs control, group A only", de_cc_A, tgt["n_deg_casecontrol"])
write.csv(de_cc_grp, file.path(root, "results/DE_casecontrol_bulk_groupadjusted.csv"), row.names = FALSE)

write.csv(de_li, file.path(root, "results/DE_lithium_bulk.csv"), row.names = FALSE)
write.csv(de_cc, file.path(root, "results/DE_casecontrol_bulk.csv"), row.names = FALSE)
write.csv(de_adj, file.path(root, "results/DE_lithium_bulk_celladjusted.csv"), row.names = FALSE)

# ---- Replication of Boltz's own DEGs (Table S15) in this cohort ------------------------------------
# Their gene IDs carry GENCODE v33lift37 versions ("ENSG...10_3"); match on the unversioned ID.
cat("\n=== Do Boltz's published DEGs replicate here? ===\n")
rep_rows <- list()
for (x in list(list("lithium", "S15_lithium_DEGs_limma_bulk.csv", de_li),
               list("case_control", "S15_case_control_DEGs_limma_bulk.csv", de_cc))) {
  b <- read.csv(file.path(root, "data/reference", x[[2]]))
  ours <- x[[3]]
  m <- match(strip(b$gene), strip(ours$gene))
  ok <- !is.na(m)
  conc <- mean(sign(b$logFC[ok]) == sign(ours$logFC[m[ok]]))
  r <- cor(b$logFC[ok], ours$logFC[m[ok]])
  rep_rows[[x[[1]]]] <- data.frame(contrast = x[[1]], boltz_degs = nrow(b), tested_here = sum(ok),
    sign_concordant = conc, logFC_r = r, nominal_p05_same_sign = mean(ours$P.Value[m[ok]] < 0.05 &
    sign(b$logFC[ok]) == sign(ours$logFC[m[ok]])), fdr05_here = sum(ours$adj.P.Val[m[ok]] < 0.05))
  cat(sprintf("%-13s %3d Boltz DEGs, %3d tested here: sign concordance %.0f%%, logFC r = %.2f, nominal p<0.05 same sign %.0f%%, FDR<0.05 here %d\n",
              x[[1]], nrow(b), sum(ok), 100 * conc, r, 100 * rep_rows[[x[[1]]]]$nominal_p05_same_sign,
              rep_rows[[x[[1]]]]$fdr05_here))
}
write.csv(do.call(rbind, rep_rows), file.path(root, "results/replication_of_boltz_degs.csv"), row.names = FALSE)

# ---- Overlap with Krebs et al. (2020) -----------------------------------------------------------
# Boltz report 33 of their 100 lithium DEGs "previously reported in Krebs et al.", Fisher OR 6.43,
# p 4.74e-14, without saying which Krebs list. Krebs deposited two (pre- and post-cell-type
# correction) and their text quotes a third number (976) that matches neither. Test Boltz's own 100
# against each, to find the list that reproduces 33 -- and ours against the same lists.
fisher_overlap <- function(query, kr) {
  U <- intersect(strip(kr$gene), universe)
  q <- intersect(strip(query), U)
  s <- intersect(strip(kr$gene[kr$adj.P.Val < 0.05]), U)
  a <- length(intersect(q, s))
  ft <- fisher.test(matrix(c(a, length(q) - a, length(s) - a, length(U) - length(s) - length(q) + a), 2))
  c(universe = length(U), query_in_universe = length(q), krebs_sig = length(s), overlap = a,
    OR = unname(ft$estimate), p = ft$p.value)
}
universe <- strip(rownames(logx))
krebs_lists <- c(
  deposited_pre_celltype  = "data/reference/File_S1_DEGs__Li_DEGs_PreCellTypeCorrection.csv",
  deposited_post_celltype = "data/reference/File_S1_DEGs__Li_DEGs_PostCellTypeCorrection.csv",
  recreated_all444        = "results/DE_lithium_all444.csv",
  recreated_cases_only    = "results/DE_lithium_casesonly.csv",
  recreated_celltype      = "results/DE_lithium_celltype.csv")
boltz_li <- read.csv(file.path(root, "data/reference/S15_lithium_DEGs_limma_bulk.csv"))$gene
ov <- do.call(rbind, lapply(names(krebs_lists), function(k) {
  f <- file.path(krebs, krebs_lists[[k]])
  if (!file.exists(f)) return(NULL)  # the Krebs recreation is optional for this script
  kr <- read.csv(f)
  rbind(data.frame(query = "boltz_published_100", krebs_list = k, t(fisher_overlap(boltz_li, kr))),
        data.frame(query = "ours_lithium", krebs_list = k, t(fisher_overlap(sig(de_li), kr))))
}))
cat(sprintf("\n=== Overlap with Krebs lithium DEGs (Boltz: %d, OR %.2f, p %.3g) ===\n",
            tgt["n_deg_lithium_in_krebs"], tgt["fisher_or_krebs"], tgt["fisher_p_krebs"]))
print(format(ov, digits = 3), row.names = FALSE)
write.csv(ov, file.path(root, "results/overlap_with_krebs.csv"), row.names = FALSE)

# ---- Figure 4C (lithium) and Figure S8 (case/control): volcano + MA ----------------------------
plot_de <- function(tt, file, title) {
  png(file.path(root, "results", file), width = 1400, height = 650, res = 130)
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  s <- tt$adj.P.Val < 0.05
  plot(tt$logFC, -log10(tt$P.Value), pch = 16, cex = 0.4, col = ifelse(s, "red3", "grey55"),
       xlab = "logFC", ylab = "-log10 p", main = sprintf("%s: %d DEGs (FDR < 0.05)", title, sum(s)))
  plot(tt$AveExpr, tt$logFC, pch = 16, cex = 0.4, col = ifelse(s, "red3", "grey55"),
       xlab = "Average log2(TPM + 1)", ylab = "logFC", main = title)
  abline(h = 0, lty = 2)
  invisible(dev.off())
}
plot_de(de_li, "fig4C_lithium_DE.png", "Lithium users vs non-users")
plot_de(de_cc, "figS8_casecontrol_DE.png", "Cases vs controls")
cat("\nWrote results/DE_*.csv, replication_of_boltz_degs.csv, overlap_with_krebs.csv, fig4C/figS8 PNGs\n")
