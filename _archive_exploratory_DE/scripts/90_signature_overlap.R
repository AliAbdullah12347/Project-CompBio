#!/usr/bin/env Rscript
# ==============================================================================
# 90_signature_overlap.R -- is the "bipolar signature" separable from lithium?
#
# This is the project's question in its simplest possible form, and it does not
# need mediation machinery to ask. The cohort contains 74 bipolar I subjects who
# are NOT on lithium. If bipolar disorder has a blood transcriptome signature
# independent of medication, those 74 should differ from controls. If almost all
# of the reported signature lives in the medicated patients, the signature is
# substantially a lithium effect.
#
# Four contrasts, each against the same 234 controls where applicable:
#   A  BP1 off lithium (74) vs Control (234)   "unmedicated bipolar signature"
#   B  BP1 on  lithium (152) vs Control (234)  "medicated bipolar signature"
#   C  lithium within BP1 (74 vs 152)          "lithium signature"
#   D  all BD (240) vs Control (234)           "the published case-control contrast"
#
# Designs are built here rather than through build_design() because each needs
# its own exposure variable, which the shared CONTRASTS list does not define.
# Everything else -- filter, normalisation, estimator -- is identical to the
# canonical settings so the four are comparable to each other and to the
# baselines.
# ==============================================================================

source("scripts/common.R")
p <- load_prep()
set.seed(481)

OUT <- file.path(EXP_HOME, "runs", "signature_overlap")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

m <- p$meta
m$tob <- p$TOB[m$title, "tobacco_imp_01"]

# Each contrast: which samples, and a two-level exposure factor over them.
SPEC <- list(
  A_off_vs_ctrl = list(
    keep = (m$dx == "BP1" & m$lithium == 0) | m$dx == "Control",
    grp  = function(d) factor(ifelse(d$dx == "Control", "ctrl", "case"), c("ctrl", "case"))),
  B_on_vs_ctrl = list(
    keep = (m$dx == "BP1" & m$lithium == 1) | m$dx == "Control",
    grp  = function(d) factor(ifelse(d$dx == "Control", "ctrl", "case"), c("ctrl", "case"))),
  C_lithium = list(
    keep = m$dx == "BP1",
    grp  = function(d) factor(ifelse(d$lithium == 1, "case", "ctrl"), c("ctrl", "case"))),
  D_allBD_vs_ctrl = list(
    keep = rep(TRUE, nrow(m)),
    grp  = function(d) factor(ifelse(d$dx == "Control", "ctrl", "case"), c("ctrl", "case")))
)

fit_one <- function(keep, grpfun, label, adjust_ilr = FALSE) {
  d <- m[keep, ]
  d$grp <- grpfun(d)
  # Assessment group is constant among cases but varies once controls are in,
  # so it belongs in A, B and D and must be absent from C. Detected, not assumed.
  covars <- c("age", "sex", "tob", "rin", "plate", "seqpc1", "seqpc2", "seqpc3")
  if (nlevels(droplevels(d$group)) > 1) covars <- c(covars, "group")
  d$group <- droplevels(d$group); d$plate <- droplevels(d$plate); d$sex <- droplevels(d$sex)
  ok <- complete.cases(d[, c("grp", covars)])
  d <- d[ok, ]
  X <- model.matrix(as.formula(paste("~ grp +", paste(covars, collapse = " + "))), data = d)
  if (adjust_ilr) X <- cbind(X, p$ILR[d$title, , drop = FALSE])
  ne <- nonEstimable(X); if (!is.null(ne)) X <- X[, setdiff(colnames(X), ne), drop = FALSE]

  cnt <- p$counts[, d$title, drop = FALSE]
  cnt <- cnt[filter_genes(cnt, 10, 0.90), , drop = FALSE]
  dge <- calcNormFactors(DGEList(cnt), method = "TMM")
  v   <- voom(dge, X)
  fit <- eBayes(lmFit(v, X))
  tt  <- topTable(fit, coef = "grpcase", number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  cat(sprintf("  %-18s n=%3d (%d vs %d)  genes=%5d  DEG=%5d  pi0=%.3f\n",
              label, nrow(d), sum(d$grp == "ctrl"), sum(d$grp == "case"),
              nrow(tt), n_deg(tt), pi0_storey(tt$P.Value)))
  list(tt = tt, n = nrow(d), n_ctrl = sum(d$grp == "ctrl"), n_case = sum(d$grp == "case"))
}

cat("== unadjusted ==\n")
F1 <- lapply(names(SPEC), function(k) fit_one(SPEC[[k]]$keep, SPEC[[k]]$grp, k))
names(F1) <- names(SPEC)

cat("\n== additionally adjusted for the 4 ILR cell-composition balances ==\n")
F2 <- lapply(names(SPEC), function(k) fit_one(SPEC[[k]]$keep, SPEC[[k]]$grp, k, TRUE))
names(F2) <- names(SPEC)

## ---- summary table ---------------------------------------------------------
summ <- do.call(rbind, lapply(names(SPEC), function(k) data.frame(
  contrast = k, n = F1[[k]]$n, n_ctrl = F1[[k]]$n_ctrl, n_case = F1[[k]]$n_case,
  genes = nrow(F1[[k]]$tt),
  deg_raw  = n_deg(F1[[k]]$tt), deg_ilr = n_deg(F2[[k]]$tt),
  pct_lost = round(100 * (1 - n_deg(F2[[k]]$tt) / max(n_deg(F1[[k]]$tt), 1)), 1),
  pi0_raw  = round(pi0_storey(F1[[k]]$tt$P.Value), 3),
  pi0_ilr  = round(pi0_storey(F2[[k]]$tt$P.Value), 3),
  stringsAsFactors = FALSE)))
print(summ, row.names = FALSE)

## ---- pairwise overlap ------------------------------------------------------
# Jaccard is only meaningful between fits that tested the same genes, so every
# comparison is restricted to the intersection of the two gene universes and
# the universe size is reported alongside it.
pairs <- combn(names(SPEC), 2, simplify = FALSE)
ov <- do.call(rbind, lapply(pairs, function(pr) {
  a <- F1[[pr[1]]]$tt; b <- F1[[pr[2]]]$tt
  uni <- intersect(a$gene, b$gene)
  da <- intersect(deg_ids(a), uni); db <- intersect(deg_ids(b), uni)
  inter <- intersect(da, db)
  ja <- a[match(uni, a$gene), ]; jb <- b[match(uni, b$gene), ]
  sign_conc <- if (length(inter)) mean(sign(a$logFC[match(inter, a$gene)]) ==
                                       sign(b$logFC[match(inter, b$gene)])) else NA
  # Expected overlap if the two lists were independent draws from the universe.
  exp_ov <- length(da) * length(db) / length(uni)
  data.frame(pair = paste(pr, collapse = " vs "), universe = length(uni),
             deg_a = length(da), deg_b = length(db), overlap = length(inter),
             expected_if_independent = round(exp_ov, 1),
             fold_enrichment = round(length(inter) / max(exp_ov, 1e-9), 1),
             jaccard = round(length(inter) / length(union(da, db)), 3),
             sign_concordance = round(sign_conc, 3),
             logFC_spearman = round(cor(ja$logFC, jb$logFC, method = "spearman"), 3),
             stringsAsFactors = FALSE)
}))
cat("\n== pairwise DEG overlap (unadjusted) ==\n"); print(ov, row.names = FALSE)

write.csv(summ, file.path(OUT, "summary.csv"), row.names = FALSE)
write.csv(ov,   file.path(OUT, "pairwise_overlap.csv"), row.names = FALSE)
for (k in names(F1)) {
  write.csv(F1[[k]]$tt[order(F1[[k]]$tt$P.Value), ][1:200, ],
            file.path(OUT, paste0("top200_", k, ".csv")), row.names = FALSE)
}
saveRDS(list(F1 = F1, F2 = F2, summ = summ, ov = ov), file.path(OUT, "fits.rds"))
cat(sprintf("\nwrote runs/signature_overlap/ (%d files)\n",
            length(list.files(OUT))))
