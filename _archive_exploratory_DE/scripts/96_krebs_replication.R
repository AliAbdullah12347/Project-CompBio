#!/usr/bin/env Rscript
# ==============================================================================
# 96_krebs_replication.R -- replicate the source paper, and use their own
# cell-type correction as an external check on ours.
#
# Krebs et al. deposited three things that make this unusually testable:
#   File_S1_DEGs__Li_DEGs_PreCellTypeCorrection.csv   lithium DE, 12,353 genes
#   File_S1_DEGs__Li_DEGs_PostCellTypeCorrection.csv  the same, after adjusting
#                                                     for cell proportions
#   File_S1_DEGs__BD_DEGs.csv                         their case-control DEGs
#
# Three questions:
#   1. Do our lithium results replicate theirs? (logFC correlation, overlap)
#   2. Their Pre vs Post pair IS a published mediation estimate. How much did
#      their correction remove, and does our number agree? An independent group,
#      a different adjustment, same direction of answer would be strong.
#   3. Their case-control DEG list is tiny. Do those genes show ANY signal in
#      bipolar I subjects who are NOT on lithium? If the published bipolar
#      signature evaporates in unmedicated patients, that is the project's
#      thesis demonstrated on the original authors' own gene list.
#
# Read-only on everything outside experimentation/.
# ==============================================================================

source("scripts/common.R")
p <- load_prep()
REF <- file.path(dirname(EXP_HOME), "Krebs Paper", "data", "reference")
OUT <- file.path(EXP_HOME, "runs", "krebs_replication")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

rd <- function(f) {
  x <- read.csv(file.path(REF, f), stringsAsFactors = FALSE)
  x$gene <- sub("\\..*$", "", x$gene)          # strip Ensembl version suffixes
  x[!duplicated(x$gene), ]
}
pre  <- rd("File_S1_DEGs__Li_DEGs_PreCellTypeCorrection.csv")
post <- rd("File_S1_DEGs__Li_DEGs_PostCellTypeCorrection.csv")
bd   <- rd("File_S1_DEGs__BD_DEGs.csv")
cat(sprintf("Krebs tables: pre %d genes, post %d genes, BD DEG list %d genes\n",
            nrow(pre), nrow(post), nrow(bd)))
cat(sprintf("  their lithium DEG at FDR .05: pre = %d, post = %d  (%.1f%% lost to their own correction)\n",
            sum(pre$adj.P.Val < 0.05), sum(post$adj.P.Val < 0.05),
            100 * (1 - sum(post$adj.P.Val < 0.05) / sum(pre$adj.P.Val < 0.05))))

## ---- our fits --------------------------------------------------------------
# lithium_krebs is the closest match to their design: all cases, diagnosis kept.
ours     <- de_fit(p, "lithium_krebs")$tt
ours_ilr <- de_fit(p, "lithium_krebs", adjust = p$ILR)$tt
cat(sprintf("\nours (lithium_krebs): %d genes, DEG pre = %d, post-ILR = %d (%.1f%% lost)\n",
            nrow(ours), n_deg(ours), n_deg(ours_ilr),
            100 * (1 - n_deg(ours_ilr) / n_deg(ours))))

## ---- 1. replication --------------------------------------------------------
u <- intersect(ours$gene, pre$gene)
a <- ours[match(u, ours$gene), ]; b <- pre[match(u, pre$gene), ]
da <- intersect(deg_ids(ours), u); db <- b$gene[b$adj.P.Val < 0.05]
inter <- intersect(da, db)
cat("\n=== 1. REPLICATION of the published lithium result ===\n")
cat(sprintf("  shared gene universe        %d\n", length(u)))
cat(sprintf("  our DEG %d | their DEG %d | overlap %d\n", length(da), length(db), length(inter)))
cat(sprintf("  expected overlap if independent %.0f  -> %.1fx enrichment\n",
            length(da) * length(db) / length(u),
            length(inter) / (length(da) * length(db) / length(u))))
cat(sprintf("  logFC Pearson %+.3f | Spearman %+.3f\n",
            cor(a$logFC, b$logFC), cor(a$logFC, b$logFC, method = "spearman")))
cat(sprintf("  sign concordance on shared DEGs %.3f\n",
            mean(sign(a$logFC[match(inter, a$gene)]) == sign(b$logFC[match(inter, b$gene)]))))

## ---- 2. their mediation estimate vs ours -----------------------------------
up <- intersect(pre$gene, post$gene)
pp <- pre[match(up, pre$gene), ]; qq <- post[match(up, post$gene), ]
cat("\n=== 2. THEIR cell-type correction vs OURS ===\n")
cat(sprintf("  Krebs   : %d -> %d DEG  (%.1f%% removed)   logFC shrinkage %.1f%%\n",
            sum(pp$adj.P.Val < 0.05), sum(qq$adj.P.Val < 0.05),
            100 * (1 - sum(qq$adj.P.Val < 0.05) / sum(pp$adj.P.Val < 0.05)),
            100 * (1 - mean(abs(qq$logFC)) / mean(abs(pp$logFC)))))
ui <- intersect(ours$gene, ours_ilr$gene)
o1 <- ours[match(ui, ours$gene), ]; o2 <- ours_ilr[match(ui, ours_ilr$gene), ]
cat(sprintf("  ours    : %d -> %d DEG  (%.1f%% removed)   logFC shrinkage %.1f%%\n",
            n_deg(o1), n_deg(o2), 100 * (1 - n_deg(o2) / n_deg(o1)),
            100 * (1 - mean(abs(o2$logFC)) / mean(abs(o1$logFC)))))
cat("  Two independent groups, different adjustments, same dataset.\n")

## ---- 3. the decisive test --------------------------------------------------
# Does the PUBLISHED bipolar signature exist in unmedicated patients?
m <- p$meta; m$tob <- p$TOB[m$title, "tobacco_imp_01"]
fit_grp <- function(keep, grpfun) {
  d <- m[keep, ]; d$grp <- grpfun(d)
  cv <- c("age", "sex", "tob", "rin", "plate", "seqpc1", "seqpc2", "seqpc3")
  if (nlevels(droplevels(d$group)) > 1) cv <- c(cv, "group")
  d$group <- droplevels(d$group); d$plate <- droplevels(d$plate); d$sex <- droplevels(d$sex)
  d <- d[complete.cases(d[, c("grp", cv)]), ]
  X <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse = " + "))), data = d)
  ne <- nonEstimable(X); if (!is.null(ne)) X <- X[, setdiff(colnames(X), ne), drop = FALSE]
  cn <- p$counts[, d$title, drop = FALSE]; cn <- cn[filter_genes(cn, 10, 0.90), , drop = FALSE]
  tt <- topTable(eBayes(lmFit(voom(calcNormFactors(DGEList(cn), "TMM"), X), X)),
                 coef = "grpcase", number = Inf, sort.by = "none")
  tt$gene <- rownames(tt); tt
}
ctrl_lab <- function(d) factor(ifelse(d$dx == "Control", "ctrl", "case"), c("ctrl", "case"))
off <- fit_grp((m$dx == "BP1" & m$lithium == 0) | m$dx == "Control", ctrl_lab)
on  <- fit_grp((m$dx == "BP1" & m$lithium == 1) | m$dx == "Control", ctrl_lab)

cat("\n=== 3. does the PUBLISHED bipolar signature exist off lithium? ===\n")
cat(sprintf("  Krebs' own case-control DEG list: %d genes\n", nrow(bd)))
for (nm in c("off", "on")) {
  t <- get(nm)
  g <- intersect(bd$gene, t$gene)
  if (!length(g)) next
  s <- t[match(g, t$gene), ]
  cat(sprintf("  %-4s lithium: %d/%d of their genes testable | median p = %.3f | %d at FDR .05 | mean |logFC| %.3f\n",
              nm, length(g), nrow(bd), median(s$P.Value), sum(s$adj.P.Val < 0.05),
              mean(abs(s$logFC))))
}
# Same question with their far larger lithium DEG list, which has the power to
# show a graded answer rather than a 6-gene coin flip.
kd <- pre$gene[pre$adj.P.Val < 0.05]
cat(sprintf("\n  using their %d-gene LITHIUM DEG list as the probe set:\n", length(kd)))
for (nm in c("off", "on")) {
  t <- get(nm); g <- intersect(kd, t$gene); s <- t[match(g, t$gene), ]
  bgp <- median(t$P.Value)
  cat(sprintf("  %-4s lithium vs control: n=%d genes | median p %.4f (background %.4f) | %d at FDR .05\n",
              nm, length(g), median(s$P.Value), bgp, sum(s$adj.P.Val < 0.05)))
}

summ <- data.frame(
  metric = c("krebs_deg_pre", "krebs_deg_post", "krebs_pct_removed",
             "our_deg_pre", "our_deg_post_ilr", "our_pct_removed",
             "logFC_cor_with_krebs", "deg_overlap", "off_lithium_deg", "on_lithium_deg"),
  value = c(sum(pre$adj.P.Val < 0.05), sum(post$adj.P.Val < 0.05),
            round(100 * (1 - sum(post$adj.P.Val < 0.05) / sum(pre$adj.P.Val < 0.05)), 1),
            n_deg(ours), n_deg(ours_ilr), round(100 * (1 - n_deg(ours_ilr) / n_deg(ours)), 1),
            round(cor(a$logFC, b$logFC), 3), length(inter),
            n_deg(off), n_deg(on)), stringsAsFactors = FALSE)
write.csv(summ, file.path(OUT, "summary.csv"), row.names = FALSE)
cat("\nwrote runs/krebs_replication/summary.csv\n")
