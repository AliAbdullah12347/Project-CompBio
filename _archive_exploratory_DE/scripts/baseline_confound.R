#!/usr/bin/env Rscript
# ==============================================================================
# baseline_confound.R -- stage three of the baseline.
#
# The canonical fits turned up something the summary tables hide: the top genes
# of the case/control contrast ARE the top genes of the lithium contrast
# (TSPAN2, RFX2, ARMC2, LINC00877, FAR2 in both top tens). This script asks how
# far that goes, because it is the premise the whole project rests on.
#
# The sharp version is a prediction with a number attached. Suppose bipolar
# diagnosis has no transcriptional effect of its own, and cases differ from
# controls only because a fraction f of them take lithium. Then the unadjusted
# case-vs-control log fold change should be f times the within-BP1 lithium log
# fold change, gene for gene -- and f is known, not fitted: 152/239 = 0.636.
#
# So: regress one on the other and read off the slope. A slope near 0.64 with a
# high R-squared says the case/control signature in this dataset is a lithium
# signature with almost nothing left over. Anything else bounds how much is
# genuinely illness.
#
# Uses only fits already cached by baseline.R and baseline_gap.R.
# Writes only inside experimentation/runs/baseline/.
# ==============================================================================

setwd("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation")
source("scripts/common.R"); source("scripts/ref_krebs.R")

CACHE <- "C:/Users/hp/AppData/Local/Temp/claude/baseline_fits"
OUT   <- file.path(EXP_HOME, "runs", "baseline")
p <- load_prep(); SYM <- gene_symbols()
sym <- function(g) ifelse(is.na(SYM[g]), g, SYM[g])

need <- c("cc_444_noli", "cc_444_li", "cc_474_noli",
          "spec_bp1_474", "spec_alldx3_444", "spec_allcase2_444")
miss <- need[!file.exists(file.path(CACHE, paste0(need, ".rds")))]
if (length(miss)) stop("missing cached fits: ", paste(miss, collapse = ", "),
                       "\nrun baseline.R and baseline_gap.R first")
G <- lapply(setNames(need, need), function(n) readRDS(file.path(CACHE, paste0(n, ".rds")))$tt)
tab <- function(x) { rownames(x) <- x$gene; x }
G <- lapply(G, tab)

# ---------------------------------------------------------------------------
# 1. the fraction of cases on lithium -- the predicted slope
# ---------------------------------------------------------------------------
qc <- p$meta$qc_pass; cases <- p$meta$dx != "Control"
f444 <- sum(qc & cases & p$meta$lithium == 1) / sum(qc & cases)
f474 <- sum(cases & p$meta$lithium == 1) / sum(cases)
cat(sprintf("== fraction of cases taking lithium ==\n  444 cohort: %d/%d = %.4f\n  474 cohort: %d/%d = %.4f\n",
  sum(qc & cases & p$meta$lithium == 1), sum(qc & cases), f444,
  sum(cases & p$meta$lithium == 1), sum(cases), f474))

# ---------------------------------------------------------------------------
# 2. does case logFC = f x lithium logFC?
# ---------------------------------------------------------------------------
cat("\n== is the case/control signature a lithium signature? ==\n")
pair <- function(ccname, liname, fexp, lab) {
  a <- G[[ccname]]; b <- G[[liname]]
  g <- intersect(rownames(a), rownames(b))
  x <- b[g, "logFC"]; y <- a[g, "logFC"]
  # through-origin slope: the prediction is y = f*x with no intercept
  s0 <- sum(x * y) / sum(x * x)
  m  <- lm(y ~ x)
  r2_0 <- 1 - sum((y - s0 * x)^2) / sum(y^2)
  data.frame(comparison = lab, n_genes = length(g),
    slope_through_origin = round(s0, 4), predicted_slope = round(fexp, 4),
    ratio_obs_pred = if (fexp > 0) round(s0 / fexp, 3) else NA_real_,
    slope_with_intercept = round(coef(m)[2], 4),
    intercept = signif(coef(m)[1], 3),
    r2 = round(summary(m)$r.squared, 4), r2_through_origin = round(r2_0, 4),
    pearson_r = round(cor(x, y), 4), spearman_rho = round(cor(x, y, method = "spearman"), 4))
}
pp <- rbind(
  pair("cc_444_noli", "spec_bp1_474",    f444, "case(444, unadj) ~ lithium(BP1 only)"),
  pair("cc_444_noli", "spec_alldx3_444", f444, "case(444, unadj) ~ lithium(full, 3-level dx)"),
  pair("cc_474_noli", "spec_bp1_474",    f474, "case(474, unadj) ~ lithium(BP1 only)"),
  pair("cc_444_li",   "spec_bp1_474",    0,    "case(444, lithium-ADJUSTED) ~ lithium(BP1)"))
print(pp, row.names = FALSE)
write.csv(pp, file.path(OUT, "confound_slope.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 3. what survives adjustment
# ---------------------------------------------------------------------------
cat("\n== what survives lithium adjustment ==\n")
a <- G$cc_444_noli; b <- G$cc_444_li
g <- intersect(rownames(a), rownames(b))
du <- rownames(a)[a$adj.P.Val < .05]; da <- rownames(b)[b$adj.P.Val < .05]
li <- rownames(G$spec_bp1_474)[G$spec_bp1_474$adj.P.Val < .05]
surv <- data.frame(
  metric = c("case DEG, lithium NOT adjusted", "case DEG, lithium adjusted",
             "of the unadjusted DEGs, still significant after adjustment",
             "median |logFC| unadjusted (all genes)", "median |logFC| adjusted (all genes)",
             "attenuation of median |logFC|",
             "unadjusted case DEGs that are also BP1 lithium DEGs",
             "as a percentage",
             "lithium DEGs (BP1 only)",
             "pi0 unadjusted", "pi0 adjusted"),
  value = c(length(du), length(da), length(intersect(du, da)),
            round(median(abs(a[g, "logFC"])), 4), round(median(abs(b[g, "logFC"])), 4),
            round(1 - median(abs(b[g, "logFC"])) / median(abs(a[g, "logFC"])), 4),
            length(intersect(du, li)),
            round(100 * length(intersect(du, li)) / length(du), 1),
            length(li),
            round(pi0_storey(a$P.Value), 4), round(pi0_storey(b$P.Value), 4)))
print(surv, row.names = FALSE)
write.csv(surv, file.path(OUT, "confound_survival.csv"), row.names = FALSE)

# hypergeometric: is the unadjusted-case / lithium overlap more than chance?
N <- length(intersect(rownames(a), rownames(G$spec_bp1_474)))
k <- length(intersect(du, li))
ph <- phyper(k - 1, length(li), N - length(li), length(du), lower.tail = FALSE)
cat(sprintf("\n  overlap %d of %d case DEGs with %d lithium DEGs in %d shared genes\n",
            k, length(du), length(li), N))
cat(sprintf("  hypergeometric p = %s   odds ratio = %.2f\n", format.pval(ph, digits = 3),
    (k / (length(du) - k)) / ((length(li) - k) / (N - length(du) - length(li) + k))))

# ---------------------------------------------------------------------------
# 4. the genes the case contrast keeps after adjustment
# ---------------------------------------------------------------------------
if (length(da)) {
  keep <- b[da, ]
  keep <- keep[order(keep$adj.P.Val), ]
  kt <- data.frame(gene = keep$gene, symbol = sym(keep$gene),
    logFC_adjusted = round(keep$logFC, 4), adjP_adjusted = signif(keep$adj.P.Val, 3),
    logFC_unadjusted = round(a[keep$gene, "logFC"], 4),
    adjP_unadjusted = signif(a[keep$gene, "adj.P.Val"], 3),
    lithium_logFC = round(G$spec_bp1_474[keep$gene, "logFC"], 4),
    lithium_adjP = signif(G$spec_bp1_474[keep$gene, "adj.P.Val"], 3))
  cat("\n== case/control DEGs surviving lithium adjustment ==\n"); print(kt, row.names = FALSE)
  write.csv(kt, file.path(OUT, "confound_surviving_genes.csv"), row.names = FALSE)
}
# ---------------------------------------------------------------------------
# 5. export the canonical topTables, so later experiments can compare against
#    the baseline without refitting it.
# ---------------------------------------------------------------------------
for (nm in c("lithium", "casecon", "lithium_krebs")) {
  f <- file.path(CACHE, paste0("canon_", nm, ".rds"))
  if (!file.exists(f)) next
  tt <- readRDS(f)$tt
  tt <- tt[order(tt$adj.P.Val, tt$P.Value), ]
  out <- data.frame(gene = tt$gene, symbol = sym(tt$gene),
    logFC = signif(tt$logFC, 6), AveExpr = signif(tt$AveExpr, 6),
    t = signif(tt$t, 6), P.Value = signif(tt$P.Value, 6),
    adj.P.Val = signif(tt$adj.P.Val, 6))
  gz <- gzfile(file.path(OUT, paste0("toptable_", nm, ".csv.gz")), "w")
  write.csv(out, gz, row.names = FALSE); close(gz)
  cat(sprintf("  exported toptable_%s.csv.gz  (%d genes)\n", nm, nrow(out)))
}

cat("\n[baseline_confound] done\n")
