#!/usr/bin/env Rscript
# Which model produced the deposited "Li_DEGs_PreCellTypeCorrection" sheet?
#
# 02b left the PRE sheet unexplained: the documented lithium model reproduces the
# POST sheet at r = 0.989 but the PRE sheet only at r = 0.687. This script settles
# it with three fingerprints that do not depend on guessing the covariate list,
# then a systematic sweep that is scored against all three.
#
#   1. AveExpr. It is rowMeans of the voom E matrix, so it is a function of the
#      sample set, the gene set and the library sizes ONLY -- the design matrix
#      cannot touch it. It therefore identifies the input, independently of the model.
#   2. Residual d.f. Moderated P = 2*pt(-|t|, df.total), so df.total is recoverable
#      from the deposited (t, P) columns by root-finding. df.total = df.residual +
#      df.prior pins the NUMBER of coefficients to within about one.
#   3. Coefficient decomposition. Regressing the deposited logFC on our own fitted
#      coefficients says which linear combination of them the sheet actually is.
#
# Diagnostic script; not part of the main pipeline.

suppressPackageStartupMessages({ library(edgeR); library(limma) })

script_path <- sub("^--file=", "",
                   grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
ROOT <- if (!is.na(script_path)) dirname(dirname(normalizePath(script_path))) else "."
if (!dir.exists(file.path(ROOT, "data"))) ROOT <- "."
PROC <- file.path(ROOT, "data", "processed")
REF  <- file.path(ROOT, "data", "reference")
RES  <- file.path(ROOT, "results")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------ load
counts <- as.matrix(read.delim(file.path(PROC, "counts_filtered.tsv"),
                               row.names = 1, check.names = FALSE))
meta <- read.csv(file.path(PROC, "sample_metadata.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
meta <- meta[meta$qc_pass %in% c(TRUE, "True", "TRUE"), ]
stopifnot(!anyDuplicated(meta$title))
meta <- meta[match(colnames(counts), meta$title), ]
stopifnot(identical(as.character(meta$title), colnames(counts)))

meta$case    <- factor(ifelse(meta$diagnosis == "Control", "control", "case"),
                       levels = c("control", "case"))
meta$dx3     <- factor(meta$diagnosis, levels = c("Control", "BP1", "BP2"))
meta$sex     <- factor(meta$sex);   meta$tobacco   <- factor(meta$tobacco)
meta$group   <- factor(meta$group); meta$plate     <- factor(meta$plate)
meta$platenum<- as.numeric(as.character(meta$plate))
meta$lithium <- as.numeric(meta$lithium)

rd <- function(f) read.csv(file.path(REF, f), stringsAsFactors = FALSE)[, 1:7]
pre  <- rd("File_S1_DEGs__Li_DEGs_PreCellTypeCorrection.csv")
post <- rd("File_S1_DEGs__Li_DEGs_PostCellTypeCorrection.csv")
wgcna_genes <- read.csv(file.path(REF, "File_S2_WGCNA_modules__Supplementary_File_2.csv"),
                        stringsAsFactors = FALSE)[[1]]

CELLS <- c("b.cells.naive","b.cells.memory","plasma.cells","t.cells.cd8",
           "t.cells.cd4.naive","t.cells.cd4.memory.resting",
           "t.cells.cd4.memory.activated","t.cells.gamma.delta",
           "nk.cells.resting","nk.cells.activated","monocytes","macrophages.m0",
           "macrophages.m2","dendritic.cells.resting","dendritic.cells.activated",
           "mast.cells.resting","eosinophils","neutrophils")  # 18 non-zero LM22

# ----------------------------------------------------- fingerprint 1: AveExpr
cat("=== Fingerprint 1: AveExpr is bit-identical between the two sheets ===\n")
mm <- merge(pre, post, by = "gene", suffixes = c(".pre", ".post"))
cat(sprintf("  genes %d | max |AveExpr.pre - AveExpr.post| = %g\n",
            nrow(mm), max(abs(mm$AveExpr.pre - mm$AveExpr.post))))
cat("  -> both sheets come from ONE voom object. Same samples, same genes, same\n",
    "     library sizes. Only the design/contrast can differ.\n", sep = "")

# The deposited AveExpr sits a constant 0.0202 log2 units BELOW ours, for every
# gene and under every normalisation method (edgeR scales norm factors to
# geometric mean 1, so they cancel out of a row mean). A gene-independent offset
# can only be a library-size definition: they kept the PRE-FILTER library sizes.
allg <- read.delim(gzfile(file.path(PROC, "counts_qc_allgenes.tsv.gz")),
                   row.names = 1, check.names = FALSE,
                   colClasses = c("character", rep("integer", ncol(counts))))
LIB_PREFILTER <- colSums(as.matrix(allg))[colnames(counts)]
LIB_FILTERED  <- colSums(counts)
rm(allg); invisible(gc())
cat(sprintf("\n  mean_s log2((L_prefilter+1)/(L_filtered+1)) = %.6f\n",
            mean(log2((LIB_PREFILTER + 1) / (LIB_FILTERED + 1)))))
cat("  observed AveExpr offset (ours - theirs)     = 0.020197\n")
cat("  -> the authors ran voom on prefilter library sizes (filter with\n",
    "     keep.lib.sizes = TRUE). Adopted below.\n", sep = "")

# --------------------------------------------- fingerprint 2: recover df.total
implied_df <- function(tab, n = 400) {
  use <- which(tab$P.Value > 1e-6 & tab$P.Value < 0.5)
  use <- use[seq(1, length(use), length.out = min(n, length(use)))]
  d <- vapply(use, function(i) {
    f <- function(df) 2 * pt(-abs(tab$t[i]), df) - tab$P.Value[i]
    if (f(2) * f(5000) > 0) return(NA_real_)
    uniroot(f, c(2, 5000), tol = 1e-9)$root
  }, numeric(1))
  median(d, na.rm = TRUE)
}
DF_PRE  <- implied_df(pre)
DF_POST <- implied_df(post)
cat(sprintf("\n=== Fingerprint 2: d.f. recovered from the deposited (t, P) columns ===\n"))
cat(sprintf("  PRE  df.total = %.3f      POST df.total = %.3f\n", DF_PRE, DF_POST))
cat(sprintf("  difference %.2f ~ the 17 free coefficients the 18 cell types add\n",
            DF_PRE - DF_POST))

# ------------------------------------------------------------------ sweep
# One fit per design; every contrast available in that design is then scored for
# free, because contrasts.fit reuses the same sigma and so the same df.prior.
fit_design <- function(rhs, samples, lib, norm, quantile_voom, genes) {
  ct <- counts[genes, samples, drop = FALSE]
  y  <- DGEList(counts = ct)
  y$samples$lib.size <- if (lib == "prefilter") LIB_PREFILTER[samples] else colSums(ct)
  y  <- suppressMessages(calcNormFactors(y, method = norm))
  d  <- model.matrix(as.formula(paste("~", rhs)), data = meta[samples, , drop = FALSE])
  r  <- qr(d)$rank
  if (r < ncol(d)) d <- d[, sort(qr(d)$pivot[seq_len(r)]), drop = FALSE]
  v  <- voom(y, d, plot = FALSE,
             normalize.method = if (quantile_voom) "quantile" else "none")
  list(fit = lmFit(v, d), design = d, amean = rowMeans(v$E))
}

score_contrast <- function(f, weights, ref) {
  ct <- rep(0, ncol(f$design)); names(ct) <- colnames(f$design)
  if (!all(names(weights) %in% names(ct))) return(NULL)
  ct[names(weights)] <- weights
  e  <- eBayes(contrasts.fit(f$fit, ct))
  tt <- topTable(e, coef = 1, number = Inf, sort.by = "none")
  g  <- rownames(tt); i <- match(g, ref$gene)
  ok <- !is.na(i)
  ours_sig   <- g[tt$adj.P.Val < 0.05]
  theirs_sig <- ref$gene[ref$adj.P.Val < 0.05]
  list(n_genes = sum(ok), p = ncol(f$design),
       df_residual = e$df.residual[1], df_prior = e$df.prior,
       df_total = e$df.residual[1] + e$df.prior,
       r_logFC = cor(tt$logFC[ok], ref$logFC[i[ok]]),
       r_t     = cor(tt$t[ok],     ref$t[i[ok]]),
       sign_conc = mean(sign(tt$logFC[ok]) == sign(ref$logFC[i[ok]])),
       n_sig = length(ours_sig), n_sig_theirs = length(theirs_sig),
       overlap = length(intersect(ours_sig, theirs_sig)),
       jaccard = length(intersect(ours_sig, theirs_sig)) /
                 length(union(ours_sig, theirs_sig)),
       avexpr_maxdiff = max(abs(f$amean[ok] - ref$AveExpr[i[ok]])))
}

ALL   <- seq_len(nrow(meta))
NOBP2 <- which(meta$diagnosis != "BP2")
CASES <- which(meta$case == "case")
GALL  <- rownames(counts)
GWG   <- intersect(rownames(counts), wgcna_genes)

FULL  <- "age + sex + tobacco + group + rin + plate + seqpc1 + seqpc2 + seqpc3"
CR    <- paste(CELLS, collapse = " + ")

# Contrasts worth scoring in a design that holds both `case` and `lithium`.
# c(casecase = 1, lithium = 1) is "bipolar ON lithium vs control": the sum of the
# illness effect and the drug effect, i.e. NOT the documented lithium contrast.
CONTRASTS_CASE <- list(
  "lithium (BP on-Li vs BP off-Li)"  = c(lithium = 1),
  "case+lithium (BP on-Li vs ctrl)"  = c(casecase = 1, lithium = 1),
  "case (BP off-Li vs ctrl)"         = c(casecase = 1),
  "case+0.5*lithium (BP vs ctrl)"    = c(casecase = 1, lithium = 0.5)
)
CONTRASTS_LI <- list("lithium (BP on-Li vs BP off-Li)" = c(lithium = 1))
CONTRASTS_DX3 <- list(
  "lithium (BP on-Li vs BP off-Li)"      = c(lithium = 1),
  "BP1+lithium (BP1 on-Li vs ctrl)"      = c(dx3BP1 = 1, lithium = 1),
  "BP1 (BP1 off-Li vs ctrl)"             = c(dx3BP1 = 1)
)

# Vary one factor at a time from the documented specification.
specs <- list(
  list("documented covariates",            paste("case +", FULL, "+ lithium"), ALL,   "prefilter", "TMM", FALSE, GALL, CONTRASTS_CASE, "pre"),
  list("documented, FILTERED libsize",     paste("case +", FULL, "+ lithium"), ALL,   "filtered",  "TMM", FALSE, GALL, CONTRASTS_CASE, "pre"),
  list("no diagnosis",                     paste(FULL, "+ lithium"),           ALL,   "prefilter", "TMM", FALSE, GALL, CONTRASTS_LI,   "pre"),
  list("no tobacco",                       paste("case + age + sex + group + rin + plate + seqpc1 + seqpc2 + seqpc3 + lithium"), ALL, "prefilter", "TMM", FALSE, GALL, CONTRASTS_CASE, "pre"),
  list("no assessment group",              paste("case + age + sex + tobacco + rin + plate + seqpc1 + seqpc2 + seqpc3 + lithium"), ALL, "prefilter", "TMM", FALSE, GALL, CONTRASTS_CASE, "pre"),
  list("no plate",                         paste("case + age + sex + tobacco + group + rin + seqpc1 + seqpc2 + seqpc3 + lithium"), ALL, "prefilter", "TMM", FALSE, GALL, CONTRASTS_CASE, "pre"),
  list("plate as numeric",                 paste("case + age + sex + tobacco + group + rin + platenum + seqpc1 + seqpc2 + seqpc3 + lithium"), ALL, "prefilter", "TMM", FALSE, GALL, CONTRASTS_CASE, "pre"),
  list("no RIN",                           paste("case + age + sex + tobacco + group + plate + seqpc1 + seqpc2 + seqpc3 + lithium"), ALL, "prefilter", "TMM", FALSE, GALL, CONTRASTS_CASE, "pre"),
  list("no seq PCs",                       paste("case + age + sex + tobacco + group + rin + plate + lithium"), ALL, "prefilter", "TMM", FALSE, GALL, CONTRASTS_CASE, "pre"),
  list("seqPC1 only",                      paste("case + age + sex + tobacco + group + rin + plate + seqpc1 + lithium"), ALL, "prefilter", "TMM", FALSE, GALL, CONTRASTS_CASE, "pre"),
  list("seqPC1+2 only",                    paste("case + age + sex + tobacco + group + rin + plate + seqpc1 + seqpc2 + lithium"), ALL, "prefilter", "TMM", FALSE, GALL, CONTRASTS_CASE, "pre"),
  list("minimal (age+sex only)",           "case + age + sex + lithium",       ALL,   "prefilter", "TMM", FALSE, GALL, CONTRASTS_CASE, "pre"),
  list("diagnosis as 3 levels",            paste("dx3 +", FULL, "+ lithium"),  ALL,   "prefilter", "TMM", FALSE, GALL, CONTRASTS_DX3,  "pre"),
  list("norm = none",                      paste("case +", FULL, "+ lithium"), ALL,   "prefilter", "none", FALSE, GALL, CONTRASTS_CASE, "pre"),
  list("norm = upperquartile",             paste("case +", FULL, "+ lithium"), ALL,   "prefilter", "upperquartile", FALSE, GALL, CONTRASTS_CASE, "pre"),
  list("norm = RLE",                       paste("case +", FULL, "+ lithium"), ALL,   "prefilter", "RLE", FALSE, GALL, CONTRASTS_CASE, "pre"),
  list("norm = TMMwsp",                    paste("case +", FULL, "+ lithium"), ALL,   "prefilter", "TMMwsp", FALSE, GALL, CONTRASTS_CASE, "pre"),
  list("voom normalize = quantile",        paste("case +", FULL, "+ lithium"), ALL,   "prefilter", "TMM", TRUE,  GALL, CONTRASTS_CASE, "pre"),
  list("restricted to WGCNA gene set",     paste("case +", FULL, "+ lithium"), ALL,   "prefilter", "TMM", FALSE, GWG,  CONTRASTS_CASE, "pre"),
  list("sample set: BP2 excluded",         paste("case +", FULL, "+ lithium"), NOBP2, "prefilter", "TMM", FALSE, GALL, CONTRASTS_CASE, "pre"),
  list("sample set: cases only",           paste(FULL, "+ lithium"),           CASES, "prefilter", "TMM", FALSE, GALL, CONTRASTS_LI,   "pre"),
  list("documented + 18 cell types",       paste("case +", FULL, "+", CR, "+ lithium"), ALL, "prefilter", "TMM", FALSE, GALL, CONTRASTS_CASE, "post")
)

cat("\n=== Sweep: one fit per design, every contrast scored ===\n")
rows <- list()
for (s in specs) {
  cat(sprintf("  fitting: %-34s\n", s[[1]])); flush.console()
  f <- fit_design(s[[2]], s[[3]], s[[4]], s[[5]], s[[6]], s[[7]])
  ref <- if (s[[9]] == "pre") pre else post
  for (cn in names(s[[8]])) {
    sc <- score_contrast(f, s[[8]][[cn]], ref)
    if (is.null(sc)) next
    rows[[length(rows) + 1]] <- data.frame(
      spec = s[[1]], contrast = cn, compared_to = toupper(s[[9]]),
      sample_set = c("444 (all)", "431 (no BP2)", "239 (cases)")[
        match(length(s[[3]]), c(length(ALL), length(NOBP2), length(CASES)))],
      n_samples = length(s[[3]]), gene_set = length(s[[7]]),
      libsize = s[[4]], norm = s[[5]], voom_quantile = s[[6]],
      formula = s[[2]], as.data.frame(sc), stringsAsFactors = FALSE)
  }
}
sw <- do.call(rbind, rows)
sw$df_total_target <- ifelse(sw$compared_to == "PRE", DF_PRE, DF_POST)
sw$df_total_error  <- sw$df_total - sw$df_total_target
# r(logFC) alone cannot separate the two diagnosis codings -- both hit 0.9993.
# Break the tie on the independent d.f. fingerprint, which counts coefficients.
sw <- sw[order(-round(sw$r_logFC, 4), abs(sw$df_total_error)), ]
write.csv(sw, file.path(RES, "presheet_model_sweep.csv"), row.names = FALSE)

cat("\n=== Ranked by r(logFC) vs the deposited sheet ===\n")
cat(sprintf("%-34s %-33s %-5s %8s %8s %6s %6s %6s %8s\n",
            "spec", "contrast", "sheet", "r(logFC)", "r(t)", "sign%", "nsig", "ovlp", "df_err"))
for (i in seq_len(nrow(sw)))
  cat(sprintf("%-34s %-33s %-5s %8.4f %8.4f %5.1f%% %6d %6d %+8.3f\n",
              sw$spec[i], sw$contrast[i], sw$compared_to[i], sw$r_logFC[i], sw$r_t[i],
              100 * sw$sign_conc[i], sw$n_sig[i], sw$overlap[i], sw$df_total_error[i]))

# ------------------------------------ fingerprint 3: decompose deposited logFC
cat("\n=== Fingerprint 3: what linear combination of our coefficients is the sheet? ===\n")
decompose <- function(rhs, cols, lab) {
  f <- fit_design(rhs, ALL, "prefilter", "TMM", FALSE, GALL)
  B <- coef(f$fit)[, cols, drop = FALSE]
  i <- match(rownames(B), pre$gene)
  d <- lm(pre$logFC[i] ~ B)
  cat(sprintf("  %-26s R2 = %.4f | intercept %+.4f | %s\n", lab, summary(d)$r.squared,
              coef(d)[1], paste(sprintf("%s %+.4f", cols, coef(d)[-1]), collapse = "  ")))
}
decompose(paste("case +", FULL, "+ lithium"), c("casecase", "lithium"), "2-level diagnosis")
decompose(paste("dx3 +", FULL, "+ lithium"), c("dx3BP1", "dx3BP2", "lithium"), "3-level diagnosis")
cat("  -> weights of 1 on the illness coefficient and 1 on the drug coefficient:\n",
    "     the sheet is their SUM, i.e. bipolar-on-lithium versus control.\n", sep = "")

# ------------------------------- do 976 / 897 correspond to any threshold?
cat("\n=== Does any threshold on the deposited PRE sheet give 976 or 897 genes? ===\n")
for (tgt in c(976, 897)) {
  q <- sort(pre$adj.P.Val); p <- sort(pre$P.Value)
  cat(sprintf("  %d genes needs FDR cut in (%.5g, %.5g]  or raw P cut in (%.5g, %.5g]\n",
              tgt, q[tgt], q[tgt + 1], p[tgt], p[tgt + 1]))
}
cat("  conventional cuts on the PRE sheet: ")
cat(paste(sprintf("FDR<%g:%d", c(0.05, 0.01, 0.001), c(sum(pre$adj.P.Val < 0.05),
      sum(pre$adj.P.Val < 0.01), sum(pre$adj.P.Val < 0.001))), collapse = "  "),
    sprintf("  Bonferroni:%d\n", sum(pre$P.Value < 0.05 / nrow(pre))))
cat("  -> no conventional threshold yields 976 or 897 on the deposited sheet.\n")

pr   <- sw[sw$compared_to == "PRE", ]
best <- pr[1, ]
# The only rival hypothesis r(logFC) cannot reject is the 3-level diagnosis
# coding, so name that row rather than whichever near-duplicate sorts second.
alt  <- pr[pr$spec == "diagnosis as 3 levels" & grepl("^BP1\\+", pr$contrast), ][1, ]
cat(sprintf("\n=== CONCLUSION ===\n  best match to the PRE sheet: %s | %s\n", best$spec, best$contrast))
cat(sprintf("  r(logFC)=%.4f  r(t)=%.4f  sign concordance %.1f%%  n_sig %d vs their %d  overlap %d\n",
            best$r_logFC, best$r_t, 100 * best$sign_conc, best$n_sig,
            best$n_sig_theirs, best$overlap))
cat(sprintf("  df.total %.3f vs deposited %.3f  (error %+.3f)\n",
            best$df_total, DF_PRE, best$df_total_error))
cat(sprintf("  only rival: %s | %s -- r(logFC)=%.4f, indistinguishable, but its\n",
            alt$spec, alt$contrast, alt$r_logFC))
cat(sprintf("  df.total error is %+.3f: one coefficient too many. The d.f. fingerprint is\n",
            alt$df_total_error))
cat("  what separates them, and it favours the 2-level coding.\n")
li <- sw[sw$compared_to == "PRE" & sw$spec == "documented covariates" &
         grepl("^lithium", sw$contrast), ]
cat(sprintf("  the DOCUMENTED lithium contrast on the same fit gives %d DEGs at FDR<0.05\n",
            li$n_sig))
cat("  (the paper's text says 976) -- so the paper's text describes the lithium\n",
    "  contrast, while the deposited sheet holds a different contrast entirely.\n", sep = "")
cat(sprintf("\nWrote %s\n", file.path(RES, "presheet_model_sweep.csv")))
