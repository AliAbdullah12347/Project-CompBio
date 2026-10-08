#!/usr/bin/env Rscript
# ==============================================================================
# 01_prep.R -- assemble the mediation dataset, and refuse to continue if any
# assumption it rests on is violated.
#
# WHAT THIS ARM ASKS
#
#   How much of lithium's effect on whole-blood gene expression runs through the
#   cell mixture rather than through changed behaviour inside cells?
#
#   exposure  Lithium (0/1)
#   mediator  the 4 ILR balances, jointly, as ONE multivariate mediator
#   outcome   gene expression, one model per gene
#
# DECISIONS MADE HERE, AND WHY
#
# 1. 226 BIPOLAR I SUBJECTS ONLY. Lithium use varies only within bipolar I
#    (74 non-users vs 152 users); there is no lithium contrast among controls.
#    Nothing is gained and a great deal is risked by pooling.
#
# 2. ASSESSMENT GROUP IS EXCLUDED. All 226 bipolar I subjects are in group A,
#    so the variable is constant here and including it makes the design
#    singular. The script asserts this rather than trusting it. Every other
#    covariate from the DE arm is kept, so the two arms adjust identically.
#
# 3. THE 350 LM22 SIGNATURE GENES ARE DROPPED. This is the circularity control
#    required by CLAUDE.md 4.2: the cell fractions are estimated FROM the
#    expression matrix, and for the signature genes the mediator is nearly a
#    deterministic function of the outcome. 547 genes define LM22; 350 of them
#    survive our expression filter, leaving 12,018 for the primary analysis.
#    The 350 are kept aside, not deleted -- 05 runs them separately, and the
#    gap between the two is a direct measurement of how much circularity was
#    there to begin with.
#
# 4. voom LOG-CPM AS THE OUTCOME, with precision weights, exactly as in
#    de_final. Mediation on a different scale from the differential expression
#    arm would make the two incomparable. voom is recomputed here rather than
#    imported, because it must be fitted on these 226 subjects and this gene
#    set, not on all 474.
#
# 5. THE MEDIATOR IS TAKEN AS BUILT. The ILR balances come from prep.rds
#    unchanged; 00b_ilr_trace.R re-derives them from the CIBERSORTx output and
#    asserts they are identical, so there is no second construction to go
#    wrong here.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})
ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
HERE <- file.path(ROOT, "mediation")
setwd(HERE); set.seed(481)
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("data",    showWarnings = FALSE, recursive = TRUE)
options(width = 150)
say <- function(...) cat(sprintf(...))

p <- readRDS(file.path(ROOT, "de_analysis/data/prep.rds"))
say("== 01_prep ==\n")

## ---- subjects ---------------------------------------------------------------
s  <- p$sel$LI
mm <- p$meta[s$keep, , drop = FALSE]
mm$lithium01 <- as.integer(s$grp == "case")          # 1 = lithium user
stopifnot(nrow(mm) == 226)
stopifnot(sum(mm$lithium01) == 152, sum(1 - mm$lithium01) == 74)
say("subjects: %d bipolar I | %d lithium users vs %d non-users\n",
    nrow(mm), sum(mm$lithium01), sum(mm$lithium01 == 0))

# assessment group must be constant here -- assert, do not assume
ng <- length(unique(mm$group))
say("assessment group levels within these subjects: %d (%s) -> %s\n", ng,
    paste(sort(unique(as.character(mm$group))), collapse = ","),
    if (ng == 1) "EXCLUDED, as it must be" else "ERROR")
stopifnot(ng == 1)

## ---- covariates -------------------------------------------------------------
# Same set as the DE arm minus assessment group. tobacco is imputation 1,
# matching de_final; 05 repeats the analysis across all 20 and pools by Rubin.
COV <- c("age", "sex", "tobacco", "rin", "plate", "seqpc1", "seqpc2", "seqpc3")
stopifnot(all(COV %in% names(mm)))
mm$sex <- droplevels(factor(mm$sex)); mm$plate <- droplevels(factor(mm$plate))
say("covariates: %s\n", paste(COV, collapse = ", "))
say("  plate levels: %d | sex levels: %d | any NA in covariates: %s\n",
    nlevels(mm$plate), nlevels(mm$sex), any(is.na(mm[, COV])))
stopifnot(!any(is.na(mm[, COV])))

Ccols <- model.matrix(as.formula(paste("~", paste(COV, collapse = " + "))), data = mm)
Ccols <- Ccols[, -1, drop = FALSE]                   # drop intercept; added per model
say("  covariate design: %d columns after dummy coding\n", ncol(Ccols))

## ---- mediator ---------------------------------------------------------------
M <- p$ILR[mm$title, , drop = FALSE]
stopifnot(identical(rownames(M), mm$title), ncol(M) == 4, !any(is.na(M)))
say("\nmediator: %d x %d ILR balances (%s)\n", nrow(M), ncol(M), paste(colnames(M), collapse = ", "))
say("  means: %s\n", paste(sprintf("%s=%.3f", colnames(M), colMeans(M)), collapse = "  "))
say("  max |correlation| between balances: %.3f (non-redundant, not independent)\n",
    max(abs(cor(M)[upper.tri(cor(M))])))

## ---- outcome: split genes on LM22 membership --------------------------------
lm22 <- unique(toupper(trimws(read.delim(file.path(ROOT, "data/LM22.txt"),
                                         check.names = FALSE, stringsAsFactors = FALSE)[[1]])))
e2s <- read.delim(file.path(ROOT, "de_analysis/data/ens2sym.tsv"), header = FALSE,
                  col.names = c("ens", "sym"), stringsAsFactors = FALSE)
sym <- setNames(toupper(e2s$sym), e2s$ens)
allg <- rownames(p$counts)
stopifnot(!any(is.na(sym[allg])))
is_marker <- sym[allg] %in% lm22
gene_main <- allg[!is_marker]; gene_mark <- allg[is_marker]
say("\noutcome genes:\n  LM22 signature genes defined : %d\n", length(lm22))
say("  of our %d that are markers    : %d  -> held aside for the circularity check\n",
    length(allg), length(gene_mark))
say("  PRIMARY analysis set          : %d\n", length(gene_main))
stopifnot(length(gene_main) + length(gene_mark) == length(allg))

## ---- voom, fitted on these 226 subjects -------------------------------------
# The full design is needed so voom's weights reflect the model that will
# actually be fitted, including the interaction terms.
X_full <- cbind(`(Intercept)` = 1, lithium = mm$lithium01, M,
                M * mm$lithium01, Ccols)
colnames(X_full)[(2 + ncol(M) + 1):(2 + 2 * ncol(M))] <- paste0("li_x_", colnames(M))
stopifnot(qr(X_full)$rank == ncol(X_full))
say("\noutcome design: %d x %d | full rank: TRUE\n", nrow(X_full), ncol(X_full))
say("  terms: intercept, lithium, %d balances, %d interactions, %d covariate columns\n",
    ncol(M), ncol(M), ncol(Ccols))

mkvoom <- function(genes) {
  dge <- calcNormFactors(DGEList(p$counts[genes, mm$title, drop = FALSE]), method = "TMM")
  voom(dge, X_full)
}
v_main <- mkvoom(gene_main); v_mark <- mkvoom(gene_mark)
say("\nvoom: primary %s | markers %s\n",
    paste(dim(v_main$E), collapse = " x "), paste(dim(v_mark$E), collapse = " x "))
say("  weight range (primary): %.3f to %.3f | any non-finite: %s\n",
    min(v_main$weights), max(v_main$weights), any(!is.finite(v_main$weights)))
stopifnot(all(is.finite(v_main$E)), all(is.finite(v_main$weights)))

## ---- lithium is balanced against the covariates? ----------------------------
# Not an assumption of the method, but worth recording: a covariate strongly
# associated with lithium is where confounding would bite hardest.
say("\nassociation of each covariate with lithium (for the record):\n")
for (cv in COV) {
  v <- mm[[cv]]
  pv <- if (is.numeric(v)) t.test(v ~ mm$lithium01)$p.value else
        suppressWarnings(chisq.test(table(v, mm$lithium01))$p.value)
  say("  %-8s p = %.4f%s\n", cv, pv, if (pv < 0.05) "   <- imbalanced" else "")
}

saveRDS(list(meta = mm, lithium = mm$lithium01, M = M, C = Ccols, X_full = X_full,
             v_main = v_main, v_mark = v_mark,
             gene_main = gene_main, gene_mark = gene_mark,
             covariates = COV, lm22 = lm22),
        "data/med_input.rds")
say("\nwrote mediation/data/med_input.rds\n")
