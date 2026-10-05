#!/usr/bin/env Rscript
# ==============================================================================
# cellcomp_annot.R -- what KIND of gene does composition adjustment kill?
#
# cellcomp.R establishes that adjustment removes most of the lithium DEGs and
# cellcomp_power.R establishes that it is not mainly a power artifact. Neither
# says anything about gene identity. If the composition story is right, the
# casualties should be cell-identity genes -- the markers that separate a
# neutrophil from a lymphocyte -- while anything that survives should not be
# lineage-restricted.
#
# Gene sets are hand-curated canonical lineage markers, written down here
# before any result was inspected, not derived from this dataset. They are NOT
# the LM22 signature (that file is not in this repository), so they do not
# reproduce the circularity by construction -- a point in their favour, since a
# marker set taken from LM22 would be guaranteed to look composition-driven.
#
# Reads runs/cellcomp/gene_level.csv and survivor_genes.csv, plus the GENCODE
# v19 annotation built into experimentation/data/gene_annot.csv.
# Writes only to experimentation/runs/cellcomp/.
# ==============================================================================

source("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation/scripts/common.R")
OUT <- file.path(EXP_HOME, "runs", "cellcomp")
cat("== cellcomp_annot ==\n")

ann <- read.csv(file.path(EXP_HOME, "data", "gene_annot.csv"), stringsAsFactors = FALSE)
gl  <- read.csv(file.path(OUT, "gene_level.csv"), stringsAsFactors = FALSE)
gl  <- merge(gl, ann, by = "gene", all.x = TRUE, sort = FALSE)
cat(sprintf("%d tested genes, %d with a symbol\n", nrow(gl), sum(!is.na(gl$symbol))))

## ---------------------------------------------------------------------------
## Direction. Lithium raises the myeloid share (path a, positive). So if the
## signature is compositional, a gene's lithium logFC should track its loading
## on the myeloid/lymphoid balance. A strong positive correlation here is the
## single cleanest picture of the whole result.
## ---------------------------------------------------------------------------
dtab <- data.frame(
  set = c("all tested genes", "canonical DEGs", "lost DEGs", "surviving DEGs"),
  n = c(nrow(gl), sum(gl$status != "never" & gl$status != "new"),
        sum(gl$status == "lost"), sum(gl$status == "survives")),
  r_logFC_vs_b1 = c(
    cor(gl$logFC_canonical, gl$b1_coef),
    cor(gl$logFC_canonical[gl$status %in% c("lost", "survives")],
        gl$b1_coef[gl$status %in% c("lost", "survives")]),
    cor(gl$logFC_canonical[gl$status == "lost"], gl$b1_coef[gl$status == "lost"]),
    cor(gl$logFC_canonical[gl$status == "survives"], gl$b1_coef[gl$status == "survives"])))
dtab$r_logFC_vs_b1 <- round(dtab$r_logFC_vs_b1, 4)
cat("\n--- correlation of canonical lithium logFC with myeloid/lymphoid loading ---\n")
print(dtab, row.names = FALSE)

# Same thing as a variance share: how much of the spread in canonical logFC
# across genes is explained by the composition term alone.
r2 <- summary(lm(logFC_canonical ~ indirect_pred, data = gl))$r.squared
cat(sprintf("R2 of canonical logFC on the predicted indirect (composition) effect: %.4f\n", r2))

## ---------------------------------------------------------------------------
## HOW MUCH OF EACH GENE IS COMPOSITION? The spike-in simulations kept
## degenerating on the compositional arm: planting an effect that runs entirely
## through the balances made the simulated genes undetectable, because the
## mediator's variance lands in their residual. Yet 1,382 genes ARE detected in
## the real data and most of them die on adjustment. Those two facts are only
## compatible if the real casualties are genes with an unusually LARGE share of
## their variance explained by composition -- far larger than the randomly
## chosen genes the spike-in used.
##
## That is an empirical question, so measure it instead of simulating it: the
## incremental R^2 from adding the four balances to the covariate-only design,
## per gene, and compare across outcome classes. No simulation, no assumption.
## ---------------------------------------------------------------------------
p <- load_prep()
d <- build_design(p, "lithium"); Xa <- d$X; mt <- d$meta$title
cnt <- p$counts[, d$samples, drop = FALSE]
cnt <- cnt[filter_genes(cnt, 10, 0.90), , drop = FALSE]
E <- voom(calcNormFactors(DGEList(counts = cnt), method = "TMM"), Xa)$E
stopifnot(all(gl$gene %in% rownames(E)))
E <- E[gl$gene, , drop = FALSE]
ILR4 <- p$ILR[mt, , drop = FALSE]

# residual sum of squares for every gene against a design, vectorised
rss_all <- function(X) {
  Q <- qr.Q(qr(X))
  rowSums(E^2) - rowSums((E %*% Q)^2)
}
r_full <- rss_all(cbind(Xa, ILR4))
r_base <- rss_all(Xa)
gl$comp_R2 <- (r_base - r_full) / r_base

compR2 <- do.call(rbind, lapply(c("lost", "survives", "new", "never"), function(s) {
  x <- gl$comp_R2[gl$status == s]
  data.frame(status = s, n = length(x),
             median_comp_R2 = round(median(x), 4),
             q25 = round(quantile(x, .25), 4), q75 = round(quantile(x, .75), 4),
             pct_over_10 = round(100 * mean(x > 0.10), 1))
}))
cat("\n--- variance explained by the 4 balances, over and above covariates ---\n")
print(compR2, row.names = FALSE)
cat(sprintf("all tested genes: median %.4f\n", median(gl$comp_R2)))
cat(sprintf("lost vs surviving, Wilcoxon p = %s\n",
            format.pval(wilcox.test(gl$comp_R2[gl$status == "lost"],
                                    gl$comp_R2[gl$status == "survives"])$p.value, 3)))
write.csv(compR2, file.path(OUT, "composition_R2_by_status.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## HOW DEAD are the casualties? "Lost a DEG" is a threshold statement, and a
## gene that slipped from p=0.049 to p=0.051 has not been refuted. Split the
## lost genes by how far they actually fell, and by how much of their effect
## size remains. If most land near the threshold the headline is soft; if most
## land at p>0.5 with the effect gone, it is not.
## ---------------------------------------------------------------------------
lost <- gl[gl$status == "lost", ]
retained_frac <- lost$logFC_ILR4 / lost$logFC_canonical
depth <- data.frame(
  band = c("still nominally sig (p<0.05 raw)", "adjP 0.05-0.20", "adjP 0.20-0.50",
           "adjP > 0.50 (signal gone)"),
  n = c(NA, sum(lost$adjP_ILR4 >= 0.05 & lost$adjP_ILR4 < 0.20),
        sum(lost$adjP_ILR4 >= 0.20 & lost$adjP_ILR4 < 0.50),
        sum(lost$adjP_ILR4 >= 0.50)))
depth$pct <- round(100 * depth$n / nrow(lost), 1)
cat(sprintf("\n--- how far the %d lost DEGs fell ---\n", nrow(lost)))
print(depth[-1, ], row.names = FALSE)
cat(sprintf("median adjusted FDR among lost genes: %.3f\n", median(lost$adjP_ILR4)))
cat(sprintf("effect retained by lost genes: median %.1f%%, %.1f%% of them keep <half\n",
            100 * median(retained_frac),
            100 * mean(retained_frac < 0.5)))
cat(sprintf("lost genes with effect reversed in sign: %d (%.1f%%)\n",
            sum(retained_frac < 0), 100 * mean(retained_frac < 0)))
write.csv(depth[-1, ], file.path(OUT, "attenuation_depth.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## Marker-set enrichment. Pre-declared, curated from textbook lineage markers.
## ---------------------------------------------------------------------------
MARK <- list(
  granulocyte = c("FCGR3B","CSF3R","CXCR2","CXCR1","S100A8","S100A9","S100A12",
                  "MMP8","MMP9","DEFA4","LCN2","CEACAM8","ELANE","MPO","ARG1",
                  "SIGLEC5","CEACAM3","ALPL","TREM1","FPR1","BASP1","CDA"),
  monocyte    = c("CD14","LYZ","VCAN","FCN1","CSF1R","MNDA","CEBPD","ITGAM",
                  "CD68","CLEC4E","SERPINA1","TYROBP","AIF1","AOAH"),
  Tcell       = c("CD3D","CD3E","CD3G","CD2","CD28","LCK","IL7R","TRAC","TRBC1",
                  "CD6","THEMIS","SKAP1","ITK","ZAP70","CD5","LAT","CD27","TCF7"),
  Bcell       = c("CD19","MS4A1","CD79A","CD79B","BLNK","PAX5","BANK1",
                  "TNFRSF13C","VPREB3","FCRL1","FCRL2","CD22"),
  NK          = c("KLRD1","KLRF1","NKG7","GNLY","PRF1","NCR1","FGFBP2","SH2D1B",
                  "KLRB1","GZMB","GZMH"),
  # Lithium's documented molecular pharmacology: inositol monophosphatase
  # pathway, GSK3, and the circadian genes lithium is known to lengthen.
  # If any part of the signature is a genuine drug effect rather than a cell
  # shift, this is where it should concentrate.
  lithium_pharm = c("IMPA1","IMPA2","INPP1","INPP5A","GSK3A","GSK3B","BCL2",
                    "NR1D1","NR1D2","PER1","PER2","PER3","ARNTL","CRY1","CRY2",
                    "CLOCK","DDIT4","AQP9","MARCKS","PPP3CA","CTNNB1","MCL1")
)
MARK$lymphoid <- unique(c(MARK$Tcell, MARK$Bcell, MARK$NK))
MARK$myeloid  <- unique(c(MARK$granulocyte, MARK$monocyte))

enrich <- do.call(rbind, lapply(names(MARK), function(k) {
  inset <- gl$symbol %in% MARK[[k]]
  n_tested <- sum(inset)
  is_deg  <- gl$status %in% c("lost", "survives")
  # hypergeometric: are set members over-represented among canonical DEGs?
  p_deg <- phyper(sum(inset & is_deg) - 1, sum(is_deg), sum(!is_deg), n_tested,
                  lower.tail = FALSE)
  # and, among the canonical DEGs, are set members preferentially LOST?
  lost_in  <- sum(inset & gl$status == "lost")
  surv_in  <- sum(inset & gl$status == "survives")
  lost_out <- sum(!inset & gl$status == "lost")
  surv_out <- sum(!inset & gl$status == "survives")
  ft <- if (lost_in + surv_in > 0)
    fisher.test(matrix(c(lost_in, surv_in, lost_out, surv_out), 2)) else NULL
  data.frame(set = k, n_in_set = length(MARK[[k]]), n_tested = n_tested,
             n_canonical_deg = sum(inset & is_deg),
             pct_of_set_deg = round(100 * sum(inset & is_deg) / max(1, n_tested), 1),
             p_deg_enrich = signif(p_deg, 3),
             n_lost = lost_in, n_survive = surv_in,
             pct_lost = round(100 * lost_in / max(1, lost_in + surv_in), 1),
             OR_lost = if (is.null(ft)) NA else round(unname(ft$estimate), 3),
             p_lost  = if (is.null(ft)) NA else signif(ft$p.value, 3),
             stringsAsFactors = FALSE)
}))
cat("\n--- marker-set enrichment (background: all tested genes / all canonical DEGs) ---\n")
print(enrich, row.names = FALSE)
cat(sprintf("\nbaseline loss rate across all canonical DEGs: %.1f%%\n",
            100 * sum(gl$status == "lost") / sum(gl$status %in% c("lost", "survives"))))

## ---------------------------------------------------------------------------
## Biotype. A composition shift moves whole expression programmes, so the
## casualties should look like ordinary protein-coding immune genes rather than
## anything technical.
## ---------------------------------------------------------------------------
bt <- as.data.frame.matrix(table(gl$biotype, gl$status))
bt <- bt[order(-rowSums(bt)), , drop = FALSE]
bt$pct_lost <- round(100 * bt$lost / pmax(1, bt$lost + bt$survives), 1)
cat("\n--- status by gene biotype (top 8) ---\n"); print(head(bt, 8))

## ---------------------------------------------------------------------------
## Named survivors: the candidate genuinely transcriptional lithium genes.
## ---------------------------------------------------------------------------
sv <- read.csv(file.path(OUT, "survivor_genes.csv"), stringsAsFactors = FALSE)
sv <- merge(sv, ann[, c("gene", "symbol", "biotype")], by = "gene", all.x = TRUE)
sv <- sv[order(sv$adjP_adjusted), ]
cat("\n--- 30 strongest composition-independent lithium genes, named ---\n")
print(head(sv[, c("symbol", "gene", "logFC_canonical", "logFC_adjusted",
                  "retained", "adjP_adjusted", "biotype")], 30), row.names = FALSE)

# and the most strongly composition-driven casualties, for contrast
lo <- gl[gl$status == "lost", ]
lo <- lo[order(-abs(lo$indirect_pred)), ]
cat("\n--- 20 DEGs whose entire lithium effect was composition ---\n")
print(head(lo[, c("symbol", "gene", "logFC_canonical", "logFC_ILR4",
                  "indirect_pred", "adjP_ILR4", "biotype")], 20), row.names = FALSE)

write.csv(dtab,   file.path(OUT, "direction_check.csv"),   row.names = FALSE)
write.csv(enrich, file.path(OUT, "marker_enrichment.csv"), row.names = FALSE)
write.csv(bt,     file.path(OUT, "biotype_by_status.csv"))
write.csv(sv,     file.path(OUT, "survivor_genes_named.csv"), row.names = FALSE)
write.csv(head(lo, 200), file.path(OUT, "top_composition_driven.csv"), row.names = FALSE)
cat("\ndone.\n")
