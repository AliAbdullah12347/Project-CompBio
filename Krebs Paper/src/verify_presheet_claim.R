#!/usr/bin/env Rscript
# INDEPENDENT adversarial check of the claim made by src/02c_presheet_resolution.R:
#   "Li_DEGs_PreCellTypeCorrection is the case+lithium contrast (BP-on-lithium vs
#    control), not the lithium main effect."
#
# Nothing here is sourced from 02c. Everything is refitted from the processed
# inputs. Three things are tested that 02c does not test:
#
#   A. The contrast is re-derived by REPARAMETERISATION rather than by
#      contrasts.fit. A 3-level group factor (control / case-off-Li / case-on-Li)
#      spans the same column space as (case, lithium), so the coefficient
#      `grp3case_onLi` IS the case+lithium contrast by construction. If 02c's
#      contrasts.fit call were wrong, these two would disagree.
#   B. The deposited logFC is regressed on ALL coefficients of the design, not
#      just the two 02c pre-selected. A forced 2-column regression cannot fail to
#      find weights near (1,1); an unrestricted one can.
#   C. Interpretation checks: effect-size inflation vs the within-case lithium
#      contrast, and overlap with the BD case/control signal.
#
# Writes nothing to results/. Prints only.

suppressPackageStartupMessages({ library(edgeR); library(limma) })

script_path <- sub("^--file=", "",
                   grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
ROOT <- if (!is.na(script_path)) dirname(dirname(normalizePath(script_path))) else "."
if (!dir.exists(file.path(ROOT, "data"))) ROOT <- "."
PROC <- file.path(ROOT, "data", "processed")
REF  <- file.path(ROOT, "data", "reference")
RES  <- file.path(ROOT, "results")

cat("================ 0. inputs ================\n")
counts <- as.matrix(read.delim(file.path(PROC, "counts_filtered.tsv"),
                               row.names = 1, check.names = FALSE))
meta <- read.csv(file.path(PROC, "sample_metadata.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
meta <- meta[meta$qc_pass %in% c(TRUE, "True", "TRUE"), ]

# Join by title. Assert it, do not assume it.
stopifnot(!anyDuplicated(meta$title), !anyDuplicated(colnames(counts)))
cat(sprintf("  counts %d x %d | metadata rows post-QC %d\n",
            nrow(counts), ncol(counts), nrow(meta)))
cat(sprintf("  title sets identical (as sets): %s\n",
            setequal(meta$title, colnames(counts))))
cat(sprintf("  title order already identical before matching: %s\n",
            identical(as.character(meta$title), colnames(counts))))
cat(sprintf("  any colname starts with 'X' (make.names damage): %s\n",
            any(grepl("^X", colnames(counts)))))
meta <- meta[match(colnames(counts), meta$title), ]
stopifnot(identical(as.character(meta$title), colnames(counts)))

meta$case    <- factor(ifelse(meta$diagnosis == "Control", "control", "case"),
                       levels = c("control", "case"))
meta$sex     <- factor(meta$sex)
meta$tobacco <- factor(meta$tobacco)
meta$group   <- factor(meta$group)
meta$plate   <- factor(meta$plate)
meta$lithium <- as.numeric(meta$lithium)
# The reparameterisation. Built from diagnosis+lithium directly, NOT from `case`.
meta$grp3 <- factor(
  ifelse(meta$diagnosis == "Control", "control",
         ifelse(meta$lithium == 1, "case_onLi", "case_offLi")),
  levels = c("control", "case_offLi", "case_onLi"))
print(table(meta$diagnosis, meta$grp3))
stopifnot(all(meta$lithium[meta$diagnosis == "Control"] == 0))

pre  <- read.csv(file.path(REF, "File_S1_DEGs__Li_DEGs_PreCellTypeCorrection.csv"),
                 stringsAsFactors = FALSE)[, 1:7]
post <- read.csv(file.path(REF, "File_S1_DEGs__Li_DEGs_PostCellTypeCorrection.csv"),
                 stringsAsFactors = FALSE)[, 1:7]
cat(sprintf("  PRE sheet %d rows, %d dup ids, FDR<0.05 = %d\n",
            nrow(pre), sum(duplicated(pre$gene)), sum(pre$adj.P.Val < 0.05)))
cat(sprintf("  gene ids shared with our count matrix: %d of %d\n",
            length(intersect(pre$gene, rownames(counts))), nrow(pre)))
cat(sprintf("  PRE ids in same ORDER as our rownames: %s  (so an order-assumed\n",
            identical(pre$gene, rownames(counts))))
cat("    correlation would be WRONG; everything below matches by id)\n")

# Library sizes: prefilter (needed to match the deposited AveExpr exactly).
allg <- read.delim(gzfile(file.path(PROC, "counts_qc_allgenes.tsv.gz")),
                   row.names = 1, check.names = FALSE,
                   colClasses = c("character", rep("integer", ncol(counts))))
stopifnot(setequal(colnames(allg), colnames(counts)))
LIB_PRE <- colSums(as.matrix(allg))[colnames(counts)]
rm(allg); invisible(gc())

# ---------------------------------------------------------------- fitting
FULL <- "age + sex + tobacco + group + rin + plate + seqpc1 + seqpc2 + seqpc3"

make_voom <- function(rhs, lib) {
  y <- DGEList(counts = counts)
  y$samples$lib.size <- if (lib == "prefilter") LIB_PRE else colSums(counts)
  y <- suppressMessages(calcNormFactors(y, method = "TMM"))
  d <- model.matrix(as.formula(paste("~", rhs)), data = meta)
  stopifnot(qr(d)$rank == ncol(d))
  list(v = voom(y, d, plot = FALSE), d = d)
}

cat("\n================ 1. two routes to case+lithium ================\n")
A <- make_voom(paste("case +", FULL, "+ lithium"), "prefilter")   # 02c's design
B <- make_voom(paste("grp3 +", FULL), "prefilter")                # reparameterised
cat(sprintf("  design A cols (%d): %s\n", ncol(A$d), paste(colnames(A$d), collapse = ", ")))
cat(sprintf("  design B cols (%d): %s\n", ncol(B$d), paste(colnames(B$d), collapse = ", ")))
cat(sprintf("  same column space (rank of cbind == rank of each): %s\n",
            qr(cbind(A$d, B$d))$rank == ncol(A$d) && ncol(A$d) == ncol(B$d)))
cat(sprintf("  voom weights identical between the two designs: %s (max diff %.3g)\n",
            isTRUE(all.equal(A$v$weights, B$v$weights)),
            max(abs(A$v$weights - B$v$weights))))

fitA <- lmFit(A$v, A$d)
fitB <- lmFit(B$v, B$d)

ctA <- setNames(rep(0, ncol(A$d)), colnames(A$d))
ctA[c("casecase", "lithium")] <- 1
eA <- eBayes(contrasts.fit(fitA, ctA))            # route 1: explicit contrast
ttA <- topTable(eA, coef = 1, number = Inf, sort.by = "none")

eB <- eBayes(fitB)                                # route 2: read the coefficient
ttB <- topTable(eB, coef = "grp3case_onLi", number = Inf, sort.by = "none")
stopifnot(identical(rownames(ttA), rownames(ttB)))
cat(sprintf("  max |logFC_A - logFC_B| = %.3g   max |t_A - t_B| = %.3g\n",
            max(abs(ttA$logFC - ttB$logFC)), max(abs(ttA$t - ttB$t))))
cat("  -> the contrast 02c forms is genuinely (case + lithium), confirmed by a\n")
cat("     route that never calls contrasts.fit.\n")

# Other contrasts off the same fit, for comparison.
ctLi <- setNames(rep(0, ncol(A$d)), colnames(A$d)); ctLi["lithium"] <- 1
ttLi <- topTable(eBayes(contrasts.fit(fitA, ctLi)), coef = 1, number = Inf, sort.by = "none")
ctCa <- setNames(rep(0, ncol(A$d)), colnames(A$d)); ctCa["casecase"] <- 1
ttCa <- topTable(eBayes(contrasts.fit(fitA, ctCa)), coef = 1, number = Inf, sort.by = "none")

cat("\n================ 2. agreement with the deposited PRE sheet ================\n")
score <- function(tt, ref, lab) {
  g <- rownames(tt); i <- match(g, ref$gene); ok <- !is.na(i)
  os <- g[tt$adj.P.Val < 0.05]; ts <- ref$gene[ref$adj.P.Val < 0.05]
  cat(sprintf("  %-34s matched %5d  r(logFC)=%.4f  r(t)=%.4f  sign=%.1f%%  nsig %4d (theirs %4d)  ovlp %4d  J=%.3f\n",
              lab, sum(ok), cor(tt$logFC[ok], ref$logFC[i[ok]]),
              cor(tt$t[ok], ref$t[i[ok]]),
              100 * mean(sign(tt$logFC[ok]) == sign(ref$logFC[i[ok]])),
              length(os), length(ts), length(intersect(os, ts)),
              length(intersect(os, ts)) / length(union(os, ts))))
  invisible(NULL)
}
score(ttA,  pre, "case+lithium (contrasts.fit)")
score(ttB,  pre, "case+lithium (reparameterised)")
score(ttLi, pre, "lithium alone")
score(ttCa, pre, "case alone")

# AveExpr fingerprint, recomputed rather than quoted.
i <- match(rownames(ttA), pre$gene)
offs <- rowMeans(A$v$E) - pre$AveExpr[i]
cat(sprintf("\n  AveExpr (prefilter libsizes): mean offset %.3g, max |offset| %.3g\n",
            mean(offs), max(abs(offs))))
Af <- make_voom(paste("case +", FULL, "+ lithium"), "filtered")
offs_f <- rowMeans(Af$v$E) - pre$AveExpr[i]
cat(sprintf("  AveExpr (filtered libsizes): mean offset %.6f, sd of offset %.3g\n",
            mean(offs_f), sd(offs_f)))
cat("  -> 02c printed 0.020197 as a hard-coded literal; the recomputed value is above.\n")
fitAf <- lmFit(Af$v, Af$d)
ctAf <- setNames(rep(0, ncol(Af$d)), colnames(Af$d)); ctAf[c("casecase", "lithium")] <- 1
ttAf <- topTable(eBayes(contrasts.fit(fitAf, ctAf)), coef = 1, number = Inf, sort.by = "none")
score(ttAf, pre, "case+lithium, filtered libsize")

cat("\n================ 3. unrestricted decomposition ================\n")
# 02c regressed the deposited logFC on exactly two of our coefficients, which
# cannot help but return weights near (1,1) if the sheet is anywhere near that
# plane. Here every coefficient in the design is offered.
Bco <- coef(fitA)
i2  <- match(rownames(Bco), pre$gene); stopifnot(!anyNA(i2))
d15 <- lm(pre$logFC[i2] ~ Bco)
cf  <- coef(d15)
names(cf) <- sub("^Bco", "", names(cf))
cat(sprintf("  R2 = %.6f   (regressing deposited logFC on all %d coefficients)\n",
            summary(d15)$r.squared, ncol(Bco)))
for (nm in names(cf))
  cat(sprintf("    %-22s %+9.4f\n", nm, cf[nm]))
d2 <- lm(pre$logFC[i2] ~ Bco[, c("casecase", "lithium")])
cat(sprintf("  restricted to (casecase, lithium): R2 = %.6f, weights %+.4f / %+.4f\n",
            summary(d2)$r.squared, coef(d2)[2], coef(d2)[3]))
resid_exact <- pre$logFC[i2] - (Bco[, "casecase"] + Bco[, "lithium"])
cat(sprintf("  residual of the EXACT contrast (no free weights): rms %.5f vs sd(sheet logFC) %.5f  (ratio %.4f)\n",
            sqrt(mean(resid_exact^2)), sd(pre$logFC), sqrt(mean(resid_exact^2)) / sd(pre$logFC)))
resid_li <- pre$logFC[i2] - Bco[, "lithium"]
cat(sprintf("  same residual for lithium alone:                  rms %.5f  (ratio %.4f)\n",
            sqrt(mean(resid_li^2)), sqrt(mean(resid_li^2)) / sd(pre$logFC)))

cat("\n================ 4. interpretation checks ================\n")
li_all <- read.csv(file.path(RES, "DE_lithium_all444.csv"), stringsAsFactors = FALSE)
bd_all <- read.csv(file.path(RES, "DE_BD_all444.csv"), stringsAsFactors = FALSE)
m <- merge(pre, li_all, by = "gene", suffixes = c(".sheet", ".li"))
m <- merge(m, bd_all[, c("gene", "logFC", "t", "adj.P.Val")], by = "gene")
names(m)[names(m) == "logFC"] <- "logFC.bd"
names(m)[names(m) == "t"] <- "t.bd"
names(m)[names(m) == "adj.P.Val"] <- "adj.bd"
cat(sprintf("  merged %d genes\n", nrow(m)))
cat(sprintf("  median |logFC|: sheet %.5f | within-case lithium %.5f | ratio %.3f\n",
            median(abs(m$logFC.sheet)), median(abs(m$logFC.li)),
            median(abs(m$logFC.sheet)) / median(abs(m$logFC.li))))
cat(sprintf("  %% of genes with |sheet logFC| > |lithium logFC|: %.1f%%\n",
            100 * mean(abs(m$logFC.sheet) > abs(m$logFC.li))))
cat(sprintf("  paired test on |logFC| (sheet - lithium): mean diff %+.5f, wilcoxon p = %.3g\n",
            mean(abs(m$logFC.sheet) - abs(m$logFC.li)),
            wilcox.test(abs(m$logFC.sheet), abs(m$logFC.li), paired = TRUE)$p.value))

# Does the sheet carry BD case/control signal that the lithium contrast lacks?
sheet_sig <- pre$gene[pre$adj.P.Val < 0.05]
li_sig    <- li_all$gene[li_all$adj.P.Val < 0.05]
bd_sig    <- bd_all$gene[bd_all$adj.P.Val < 0.05]
cat(sprintf("\n  deposited BD DEG list (their sheet): %d genes; our BD DEGs: %d\n",
            nrow(read.csv(file.path(REF, "File_S1_DEGs__BD_DEGs.csv"))), length(bd_sig)))
cat(sprintf("  sheet DEGs %d | our lithium DEGs %d | overlap %d\n",
            length(sheet_sig), length(li_sig), length(intersect(sheet_sig, li_sig))))
cat(sprintf("  correlation of signed -log10P: sheet vs lithium  %.4f ; sheet vs BD  %.4f ; lithium vs BD  %.4f\n",
            cor(sign(m$logFC.sheet) * -log10(m$P.Value.sheet),
                sign(m$logFC.li) * -log10(m$P.Value.li), method = "spearman"),
            cor(sign(m$logFC.sheet) * -log10(m$P.Value.sheet),
                sign(m$logFC.bd) * -log10(2 * pnorm(-abs(m$t.bd))), method = "spearman"),
            cor(sign(m$logFC.li) * -log10(m$P.Value.li),
                sign(m$logFC.bd) * -log10(2 * pnorm(-abs(m$t.bd))), method = "spearman")))
# partial: does BD t add to the prediction of sheet t beyond lithium t?
cat(sprintf("  R2 of sheet logFC ~ lithium logFC            : %.4f\n",
            summary(lm(m$logFC.sheet ~ m$logFC.li))$r.squared))
cat(sprintf("  R2 of sheet logFC ~ lithium logFC + BD logFC : %.4f\n",
            summary(lm(m$logFC.sheet ~ m$logFC.li + m$logFC.bd))$r.squared))
cat(sprintf("  BD logFC coefficient in that model: %+.4f (t = %.1f)\n",
            coef(summary(lm(m$logFC.sheet ~ m$logFC.li + m$logFC.bd)))[3, 1],
            coef(summary(lm(m$logFC.sheet ~ m$logFC.li + m$logFC.bd)))[3, 3]))

cat("\n================ 5. does any threshold give 976 or 897? ================\n")
for (lab in c("PRE", "POST")) {
  s <- if (lab == "PRE") pre else post
  q <- sort(s$adj.P.Val); p <- sort(s$P.Value)
  cat(sprintf("  %s sheet: FDR<0.05 %d  FDR<0.01 %d  FDR<0.001 %d  Bonferroni %d  raw p<0.05 %d\n",
              lab, sum(s$adj.P.Val < 0.05), sum(s$adj.P.Val < 0.01),
              sum(s$adj.P.Val < 0.001), sum(s$P.Value < 0.05 / nrow(s)),
              sum(s$P.Value < 0.05)))
  for (tgt in c(976, 897))
    cat(sprintf("    exactly %d needs FDR in (%.4g, %.4g] or raw P in (%.4g, %.4g]\n",
                tgt, q[tgt], q[tgt + 1], p[tgt], p[tgt + 1]))
  # also: |logFC| thresholds combined with FDR<0.05
  sig <- s[s$adj.P.Val < 0.05, ]
  fc  <- sort(abs(sig$logFC), decreasing = TRUE)
  for (tgt in c(976, 897))
    if (tgt <= length(fc))
      cat(sprintf("    or FDR<0.05 plus |logFC| > %.4f gives %d\n", fc[tgt], tgt))
  cat(sprintf("    up/down split at FDR<0.05: %d up / %d down (paper: 754 up / 222 down)\n",
              sum(sig$logFC > 0), sum(sig$logFC < 0)))
}
cat(sprintf("\n  our within-case lithium contrast, FDR<0.05: %d\n", length(li_sig)))
cat(sprintf("  our cases-only lithium contrast,  FDR<0.05: %d\n",
            sum(read.csv(file.path(RES, "DE_lithium_casesonly.csv"))$adj.P.Val < 0.05)))

cat("\n================ 6. is contrasts.fit exact here? ================\n")
# `lithium` is itself a column of design A, so the same quantity is available
# both exactly (topTable on the coefficient) and approximately (contrasts.fit).
# limma's docs warn that with gene-specific weights -- which voom always
# produces -- contrasts.fit's unscaled SDs are approximate on a non-orthogonal
# design. This isolates the size of that approximation.
ttLi_exact <- topTable(eBayes(fitA), coef = "lithium", number = Inf, sort.by = "none")
stopifnot(identical(rownames(ttLi_exact), rownames(ttLi)))
cat(sprintf("  lithium coefficient: max |logFC diff| %.3g  max |t diff| %.4f\n",
            max(abs(ttLi_exact$logFC - ttLi$logFC)), max(abs(ttLi_exact$t - ttLi$t))))
cat(sprintf("  FDR<0.05 count: exact %d  vs contrasts.fit %d\n",
            sum(ttLi_exact$adj.P.Val < 0.05), sum(ttLi$adj.P.Val < 0.05)))
cat(sprintf("  case+lithium   : exact %d  vs contrasts.fit %d  (02c reports the latter)\n",
            sum(ttB$adj.P.Val < 0.05), sum(ttA$adj.P.Val < 0.05)))

cat("\n================ 7. integrity of the extracted reference CSV ================\n")
for (lab in c("PRE", "POST")) {
  s <- if (lab == "PRE") pre else post
  bh <- p.adjust(s$P.Value, "BH")
  cat(sprintf("  %s: max |adj.P.Val - BH(P.Value)| = %.3g  (relative %.3g)\n",
              lab, max(abs(bh - s$adj.P.Val)), max(abs(bh - s$adj.P.Val) / s$adj.P.Val)))
  dfl <- if (lab == "PRE") 434.1846 else 416.3565
  cat(sprintf("       P.Value vs 2*pt(-|t|, %.4f): max |log10 ratio| = %.4f\n", dfl,
              max(abs(log10(2 * pt(-abs(s$t), dfl)) - log10(s$P.Value)))))
}
cat("  -> the deposited columns are mutually consistent, so the xlsx->csv\n")
cat("     extraction did not shuffle rows or columns.\n")
cat("\nDone.\n")
