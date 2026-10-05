#!/usr/bin/env Rscript
# ==============================================================================
# norm_sweep.R -- vary NORMALISATION ONLY, on the lithium contrast.
#
# Everything else is held fixed at the project baseline: bipolar-I-only subset
# (n=226), the same covariate set, the same gene filter (>10 counts in >=90% of
# samples), the same eBayes settings. The gene filter runs on RAW counts, so the
# gene set is identical across every normalisation method by construction --
# which is what makes the comparison clean. Any difference in the result is
# attributable to the scaling factors and to nothing else.
#
# Three questions:
#   Q1  How much does the normalisation method move the answer?
#   Q2  Do the normalisation factors themselves track lithium status? If they
#       do, normalisation is not a nuisance correction here -- it is a partial
#       erasure (or amplification) of the exposure effect, and the direction
#       matters.
#   Q3  If they track lithium, is that because lithium shifts the cell mixture?
#       TMM's reference quantile is set by the dominant transcripts. Neutrophils
#       are the dominant cell in whole blood and carry extremely high-abundance
#       transcripts, so a neutrophil shift is exactly the kind of shift a
#       scaling factor is built to notice.
#
# Writes only inside experimentation/.
# ==============================================================================

suppressPackageStartupMessages({library(edgeR); library(limma)})
source("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation/scripts/common.R")

OUTDIR <- file.path(EXP_HOME, "runs", "norm_sweep")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
RES <- file.path(EXP_HOME, "results")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

set.seed(20261218)
p <- load_prep()

cat("=========================================================\n")
cat("norm_sweep -- lithium contrast, normalisation varied only\n")
cat("=========================================================\n\n")

## ===========================================================================
## PART A -- the normalisation factors themselves
## ===========================================================================
# Rebuild exactly the DGEList de_fit() builds, so the factors examined here are
# the factors the models actually used. Doing this by hand rather than reaching
# into de_fit's internals keeps the engine untouched.
d   <- build_design(p, "lithium")
cnt <- p$counts[, d$samples, drop = FALSE]
gid <- filter_genes(cnt, 10, 0.90)
cnt <- cnt[gid, , drop = FALSE]
stopifnot(identical(colnames(cnt), d$meta$title))
cat(sprintf("subset: %d samples x %d genes after filter\n", ncol(cnt), nrow(cnt)))

METHODS <- c("TMM", "TMMwsp", "RLE", "upperquartile", "none")
NF <- matrix(NA_real_, ncol(cnt), length(METHODS),
             dimnames = list(colnames(cnt), METHODS))
for (mm in METHODS) {
  dge <- DGEList(counts = cnt)
  if (mm != "none") dge <- suppressMessages(calcNormFactors(dge, method = mm))
  NF[, mm] <- dge$samples$norm.factors
}
libsize <- colSums(cnt)

cat("\n-- normalisation factor spread (these multiply the library size) --\n")
nf_spread <- data.frame(
  method  = METHODS,
  min     = apply(NF, 2, min),
  q25     = apply(NF, 2, quantile, .25),
  median  = apply(NF, 2, median),
  q75     = apply(NF, 2, quantile, .75),
  max     = apply(NF, 2, max),
  sd      = apply(NF, 2, sd),
  iqr     = apply(NF, 2, IQR),
  range_ratio = apply(NF, 2, function(v) max(v) / min(v)),
  # the spread that actually matters is on the log scale, because the factor
  # enters the model as an offset
  sd_log2 = apply(NF, 2, function(v) sd(log2(v))),
  row.names = NULL)
print(format(nf_spread, digits = 4))

## ---- how much do the methods agree with each other? -----------------------
nfm <- METHODS[METHODS != "none"]
pair <- expand.grid(a = nfm, b = nfm, stringsAsFactors = FALSE)
pair <- pair[pair$a < pair$b, ]
pair$spearman <- mapply(function(a, b) cor(NF[, a], NF[, b], method = "spearman"),
                        pair$a, pair$b)
pair$pearson_log2 <- mapply(function(a, b) cor(log2(NF[, a]), log2(NF[, b])),
                            pair$a, pair$b)
cat("\n-- agreement between normalisation methods (factor vectors) --\n")
print(format(pair, digits = 4))

## ===========================================================================
## PART B -- do the factors track lithium? (the compositional-confound test)
## ===========================================================================
meta <- d$meta
lith <- meta$lithium
stopifnot(!any(is.na(lith)))

# Transcriptome dominance: what share of all counts sits in the few most
# abundant genes. This is the quantity a scaling factor is most sensitive to,
# so it is the natural mechanistic bridge between cell mixture and norm factor.
cpm_raw  <- sweep(cnt, 2, libsize, "/")
top_share <- function(k) apply(cpm_raw, 2, function(v) sum(sort(v, TRUE)[seq_len(k)]))
dom <- data.frame(top10 = top_share(10), top50 = top_share(50),
                  top100 = top_share(100))

# Globin is the classic whole-blood dominance problem: a handful of transcripts
# can take a double-digit share of the library and they are erythrocyte-derived,
# i.e. nothing to do with the leukocyte biology anyone is trying to measure.
GLOBIN <- c(HBB = "ENSG00000244734", HBA1 = "ENSG00000206172",
            HBA2 = "ENSG00000188536", HBD = "ENSG00000223609",
            HBG1 = "ENSG00000213934", HBG2 = "ENSG00000196565",
            HBM  = "ENSG00000206177", ALAS2 = "ENSG00000158578")
gl_present <- GLOBIN[GLOBIN %in% rownames(cnt)]
cat(sprintf("\nglobin/erythroid genes present in filtered set: %d of %d (%s)\n",
            length(gl_present), length(GLOBIN), paste(names(gl_present), collapse = ",")))
dom$globin <- if (length(gl_present)) colSums(cpm_raw[gl_present, , drop = FALSE]) else NA_real_

comp <- data.frame(p$lineage[meta$title, , drop = FALSE],
                   p$ILR[meta$title, , drop = FALSE], check.names = FALSE)

# Two tests per candidate. The unadjusted one answers "is it associated at all".
# The adjusted one uses the SAME covariate block as the DE model, so it answers
# the question that actually bears on the DE result: does the factor still move
# with lithium after everything the model already controls for?
Xcov <- d$X[, setdiff(colnames(d$X), "lithium"), drop = FALSE]
assoc <- function(y, label, group) {
  y  <- as.numeric(y)
  tt <- t.test(y ~ lith)
  m0 <- lm(y ~ Xcov - 1)
  m1 <- lm(y ~ lith + Xcov - 1)
  an <- anova(m0, m1)
  # partial correlation of y with lithium after the covariate block
  ry <- residuals(lm(y    ~ Xcov - 1))
  rl <- residuals(lm(lith ~ Xcov - 1))
  data.frame(group = group, variable = label,
             mean_nonuser = mean(y[lith == 0]), mean_user = mean(y[lith == 1]),
             diff = mean(y[lith == 1]) - mean(y[lith == 0]),
             t_unadj = unname(tt$statistic), p_unadj = tt$p.value,
             beta_adj = unname(coef(m1)["lith"]),
             p_adj = an$`Pr(>F)`[2],
             partial_r = cor(ry, rl),
             stringsAsFactors = FALSE)
}

rows <- list()
for (mm in METHODS[METHODS != "none"])
  rows[[length(rows) + 1]] <- assoc(log2(NF[, mm]), paste0("log2_normfactor_", mm), "normfactor")
rows[[length(rows) + 1]] <- assoc(log2(libsize), "log2_library_size", "depth")
for (mm in METHODS[METHODS != "none"])
  rows[[length(rows) + 1]] <- assoc(log2(libsize * NF[, mm]),
                                    paste0("log2_effective_lib_", mm), "effective_lib")
for (nm in names(dom))  rows[[length(rows) + 1]] <- assoc(dom[[nm]], paste0("dominance_", nm), "dominance")
for (nm in names(comp)) rows[[length(rows) + 1]] <- assoc(comp[[nm]], nm, "composition")
for (nm in c("Correlation", "RMSE"))
  rows[[length(rows) + 1]] <- assoc(p$cs_qc[meta$title, nm], paste0("cibersortx_", nm), "deconv_qc")
confound <- do.call(rbind, rows)
confound$q_adj <- p.adjust(confound$p_adj, "BH")

cat("\n-- association with lithium (unadjusted, and adjusted for the DE covariate block) --\n")
print(format(confound[, c("group", "variable", "diff", "p_unadj", "beta_adj",
                          "p_adj", "partial_r")], digits = 3), row.names = FALSE)

## ---- if the factor moves with lithium, is composition the reason? ---------
# Regress each log2 factor on the lineage balances. A high R^2 means the factor
# is largely a restatement of the cell mixture, which is the claim being tested.
nf_comp <- do.call(rbind, lapply(METHODS[METHODS != "none"], function(mm) {
  y <- log2(NF[, mm])
  r2 <- function(f) summary(lm(f))$r.squared
  data.frame(method = mm,
             r2_ILR4     = r2(y ~ as.matrix(p$ILR[meta$title, ])),
             r2_gran     = r2(y ~ comp$gran),
             r2_b1       = r2(y ~ comp$b1_myeloid_vs_lymphoid),
             r2_globin   = r2(y ~ dom$globin),
             r2_top50    = r2(y ~ dom$top50),
             r2_lib      = r2(y ~ log2(libsize)),
             r_gran      = cor(y, comp$gran),
             r_globin    = cor(y, dom$globin),
             r_top50     = cor(y, dom$top50),
             stringsAsFactors = FALSE)
}))
cat("\n-- how much of each log2 norm factor is explained by composition/dominance --\n")
print(format(nf_comp, digits = 3), row.names = FALSE)

## ===========================================================================
## PART C -- the differential expression sweep
## ===========================================================================
CFG <- list(
  list(id = "voom_TMM",            norm = "TMM",           method = "voom",  ref = TRUE),
  list(id = "voom_TMMwsp",         norm = "TMMwsp",        method = "voom"),
  list(id = "voom_RLE",            norm = "RLE",           method = "voom"),
  list(id = "voom_upperquartile",  norm = "upperquartile", method = "voom"),
  list(id = "voom_none_CPM",       norm = "none",          method = "voom"),
  list(id = "trend_TMM",           norm = "TMM",           method = "trend"),
  list(id = "trend_none_logCPM",   norm = "none",          method = "trend")
)

fits <- list()
for (cf in CFG) {
  t0 <- Sys.time()
  f <- de_fit(p, "lithium", norm = cf$norm, method = cf$method)
  fits[[cf$id]] <- f
  cat(sprintf("  %-20s  genes=%d  n=%d  DEG(0.05)=%4d  pi0=%.4f   [%.0fs]\n",
              cf$id, f$n_genes, f$n_samples, n_deg(f$tt),
              pi0_storey(f$tt$P.Value), as.numeric(Sys.time() - t0, units = "secs")))
}

# The gene filter is norm-independent, so this must hold. If it ever fails the
# comparison is no longer apples-to-apples and every metric below is void.
gsets <- lapply(fits, function(f) f$tt$gene)
stopifnot(all(vapply(gsets, function(g) identical(sort(g), sort(gsets[[1]])), logical(1))))
cat(sprintf("\nall %d configurations share the same %d genes -- comparison is clean\n",
            length(fits), fits[[1]]$n_genes))

ord <- order(fits[[1]]$tt$gene)
getv <- function(f, col) f$tt[[col]][order(f$tt$gene)]
REF <- "voom_TMM"
ref_deg <- deg_ids(fits[[REF]]$tt)
ref_t   <- getv(fits[[REF]], "t")
ref_lfc <- getv(fits[[REF]], "logFC")

jacc <- function(a, b) if (length(union(a, b)) == 0) NA else length(intersect(a, b)) / length(union(a, b))

summ <- do.call(rbind, lapply(names(fits), function(id) {
  f <- fits[[id]]; tt <- f$tt
  dg <- deg_ids(tt)
  tv <- getv(f, "t"); lf <- getv(f, "logFC")
  # a sign flip among genes significant in BOTH runs is the strongest possible
  # form of "the normalisation changed the conclusion"
  both <- intersect(dg, ref_deg)
  iboth <- match(both, sort(tt$gene))
  flips <- sum(sign(tv[iboth]) != sign(ref_t[iboth]))
  topref <- fits[[REF]]$tt$gene[order(fits[[REF]]$tt$P.Value)][1:100]
  topthis <- tt$gene[order(tt$P.Value)][1:100]
  data.frame(
    config = id, method = f$method, norm = f$norm,
    n_genes = f$n_genes, n_samples = f$n_samples,
    deg_fdr01 = n_deg(tt, 0.01), deg_fdr05 = n_deg(tt, 0.05), deg_fdr10 = n_deg(tt, 0.10),
    deg_fdr05_lfc0.1 = n_deg(tt, 0.05, 0.1),
    deg_fdr05_lfc0.2 = n_deg(tt, 0.05, 0.2),
    deg_up = sum(tt$adj.P.Val < 0.05 & tt$logFC > 0),
    deg_dn = sum(tt$adj.P.Val < 0.05 & tt$logFC < 0),
    pi0 = pi0_storey(tt$P.Value),
    p_lt_05_frac = mean(tt$P.Value < 0.05),
    jaccard_vs_TMM = jacc(dg, ref_deg),
    n_shared_deg = length(both),
    n_only_here  = length(setdiff(dg, ref_deg)),
    n_only_TMM   = length(setdiff(ref_deg, dg)),
    top100_overlap = length(intersect(topref, topthis)),
    spearman_t_vs_TMM = cor(tv, ref_t, method = "spearman"),
    pearson_t_vs_TMM  = cor(tv, ref_t),
    spearman_lfc_vs_TMM = cor(lf, ref_lfc, method = "spearman"),
    median_lfc_shift = median(lf - ref_lfc),
    sd_lfc_shift = sd(lf - ref_lfc),
    sign_flips_among_shared_deg = flips,
    nf_sd_log2 = if (f$norm == "none") 0 else nf_spread$sd_log2[match(f$norm, nf_spread$method)],
    stringsAsFactors = FALSE)
}))
cat("\n-- DE sweep summary --\n")
print(format(summ[, c("config", "deg_fdr05", "pi0", "jaccard_vs_TMM",
                      "spearman_t_vs_TMM", "median_lfc_shift",
                      "sign_flips_among_shared_deg", "nf_sd_log2")], digits = 4),
      row.names = FALSE)

## ---- full pairwise agreement, not just vs TMM -----------------------------
ids <- names(fits)
pw <- expand.grid(a = ids, b = ids, stringsAsFactors = FALSE)
pw <- pw[pw$a < pw$b, ]
pw$jaccard  <- mapply(function(a, b) jacc(deg_ids(fits[[a]]$tt), deg_ids(fits[[b]]$tt)), pw$a, pw$b)
pw$spearman_t <- mapply(function(a, b) cor(getv(fits[[a]], "t"), getv(fits[[b]], "t"),
                                           method = "spearman"), pw$a, pw$b)
cat("\n-- pairwise agreement across all configurations --\n")
print(format(pw, digits = 4), row.names = FALSE)

## ===========================================================================
## PART D -- can normalisation change a conclusion?
## ===========================================================================
# Three concrete operationalisations, because "changes the conclusion" is vague
# unless it is pinned to something a reader would actually assert.
cat("\n-- conclusion-stability probes --\n")

# (1) is the leading gene the same?
lead <- data.frame(config = ids,
                   top_gene = vapply(ids, function(i) fits[[i]]$tt$gene[which.min(fits[[i]]$tt$P.Value)], ""),
                   top_p    = vapply(ids, function(i) min(fits[[i]]$tt$P.Value), 0),
                   stringsAsFactors = FALSE, row.names = NULL)
print(format(lead, digits = 3), row.names = FALSE)

# (2) per-gene status churn: significant under some methods, not others
sigmat <- sapply(ids, function(i) {
  tt <- fits[[i]]$tt; setNames(tt$adj.P.Val < 0.05, tt$gene)[sort(fits[[1]]$tt$gene)] })
nsig <- rowSums(sigmat)
churn <- data.frame(
  n_sig_in_all   = sum(nsig == length(ids)),
  n_sig_in_none  = sum(nsig == 0),
  n_sig_unstable = sum(nsig > 0 & nsig < length(ids)),
  pct_of_ever_sig_that_are_unstable = 100 * sum(nsig > 0 & nsig < length(ids)) / sum(nsig > 0))
print(churn, row.names = FALSE)

# voom-only subset: the fair question is whether SWAPPING THE SCALING FACTOR
# (holding the estimator fixed) churns genes, since mixing in limma-trend
# confounds normalisation with estimator.
vids <- ids[grepl("^voom_", ids)]
sm_v <- sigmat[, vids]; nv <- rowSums(sm_v)
churn_v <- data.frame(
  scope = "voom only (normalisation varied alone)",
  n_sig_in_all = sum(nv == length(vids)), n_sig_in_none = sum(nv == 0),
  n_sig_unstable = sum(nv > 0 & nv < length(vids)),
  pct_of_ever_sig_that_are_unstable = 100 * sum(nv > 0 & nv < length(vids)) / sum(nv > 0))
print(churn_v, row.names = FALSE)

# (3) the direct test of the confound: put the log2 TMM factor into the design
# as an explicit covariate. If the factor carries lithium signal, adding it
# should absorb some of the exposure effect and the DEG count should fall.
adj_nf <- matrix(log2(NF[, "TMM"]), ncol = 1,
                 dimnames = list(rownames(NF), "log2_tmm_nf"))
f_nfadj <- de_fit(p, "lithium", adjust = adj_nf, norm = "TMM", method = "voom")
adj_ilr <- p$ILR
f_ilradj <- de_fit(p, "lithium", adjust = adj_ilr, norm = "TMM", method = "voom")
probe3 <- data.frame(
  model = c("baseline voom+TMM", "+ log2 TMM factor as covariate", "+ 4 ILR balances as covariates"),
  deg_fdr05 = c(n_deg(fits[[REF]]$tt), n_deg(f_nfadj$tt), n_deg(f_ilradj$tt)),
  pi0 = c(pi0_storey(fits[[REF]]$tt$P.Value), pi0_storey(f_nfadj$tt$P.Value),
          pi0_storey(f_ilradj$tt$P.Value)),
  jaccard_vs_baseline = c(1, jacc(deg_ids(f_nfadj$tt), ref_deg), jacc(deg_ids(f_ilradj$tt), ref_deg)),
  stringsAsFactors = FALSE)
cat("\n-- probe 3: is the scaling factor carrying exposure signal? --\n")
print(format(probe3, digits = 4), row.names = FALSE)

## ===========================================================================
## PART E -- write everything out
## ===========================================================================
write.csv(nf_spread, file.path(OUTDIR, "normfactor_spread.csv"), row.names = FALSE)
write.csv(pair,      file.path(OUTDIR, "normfactor_method_agreement.csv"), row.names = FALSE)
write.csv(confound,  file.path(OUTDIR, "lithium_confound_tests.csv"), row.names = FALSE)
write.csv(nf_comp,   file.path(OUTDIR, "normfactor_explained_by_composition.csv"), row.names = FALSE)
write.csv(pw,        file.path(OUTDIR, "pairwise_agreement.csv"), row.names = FALSE)
write.csv(probe3,    file.path(OUTDIR, "probe_normfactor_as_covariate.csv"), row.names = FALSE)
write.csv(rbind(cbind(scope = "all configs", churn), churn_v[, names(churn_v)[-1]] |>
                (\(x) cbind(scope = "voom only", x))()),
          file.path(OUTDIR, "gene_status_churn.csv"), row.names = FALSE)
nfout <- data.frame(title = rownames(NF), lithium = lith, lib = libsize,
                    NF, dom, comp, check.names = FALSE)
write.csv(nfout, file.path(OUTDIR, "per_sample_normfactors.csv"), row.names = FALSE)
saveRDS(list(summ = summ, NF = NF, confound = confound, fits_tt =
               lapply(fits, function(f) f$tt[, c("gene", "logFC", "t", "P.Value", "adj.P.Val")])),
        file.path(OUTDIR, "norm_sweep_objects.rds"))

write_run("norm_sweep", summ, notes = NULL)  # NOTES.md written separately
cat("\ndone.\n")
