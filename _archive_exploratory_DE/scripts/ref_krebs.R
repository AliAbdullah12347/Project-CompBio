#!/usr/bin/env Rscript
# ==============================================================================
# ref_krebs.R -- the published Krebs et al. (2020) results, loaded and made
# comparable to our fits. Read-only with respect to everything outside
# experimentation/.
#
# Two deposited DE tables exist (Supplementary File S1):
#   Li_DEGs_PreCellTypeCorrection   12,353 genes
#   Li_DEGs_PostCellTypeCorrection  12,353 genes
# Both carry versioned Ensembl IDs, so the "map symbols to Ensembl" step the
# task anticipated is not needed: after stripping the version suffix all
# 12,353/12,353 match our count-matrix row names exactly. The GENCODE v19 GTF
# is used only to attach symbols for human reading, never to join.
#
# source() this after common.R.
# ==============================================================================

REF_DIR <- file.path(dirname(EXP_HOME), "Krebs Paper", "data", "reference")
stopifnot(dir.exists(REF_DIR))

# Numbers the paper itself reports, kept in one place so that a comparison can
# never quietly drift from the published claim it is checking.
KREBS_PAPER <- list(
  n_genes_filtered = 12344,   # main text; the deposited tables hold 12,353
  li_ndeg          = 976,     # Results text
  li_ndeg_magma    = 897,     # Methods (MAGMA) text -- inconsistent with above
  li_nup           = 754,
  li_ndown         = 222,
  li_absFC_mean    = 0.20,
  li_absFC_max     = 0.82,
  li_absFC_sd      = 0.10,
  li_post_ndeg     = 233,     # after cell-type correction
  li_post_carry    = 194,     # of those 233, significant + concordant in the pre model
  li_post_carry_pct= 83.2,
  bd_ndeg          = 6,
  n_subjects       = 444, n_cases = 240, n_controls = 204,
  n_li_user = 152, n_li_nonuser = 88
)

read_ref <- function(f) {
  d <- read.csv(file.path(REF_DIR, f), stringsAsFactors = FALSE)
  d <- d[, intersect(c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"),
                     names(d))]
  d$gene <- sub("\\.\\d+$", "", d$gene)   # version suffix; ours are stripped
  stopifnot(!any(duplicated(d$gene)))
  d
}

krebs_tables <- function() {
  list(PRE  = read_ref("File_S1_DEGs__Li_DEGs_PreCellTypeCorrection.csv"),
       POST = read_ref("File_S1_DEGs__Li_DEGs_PostCellTypeCorrection.csv"),
       BD   = read_ref("File_S1_DEGs__BD_DEGs.csv"),
       ANAND = read.csv(file.path(REF_DIR, "File_S1_DEGs__Li_DEGs_AnandOverlap.csv"),
                        stringsAsFactors = FALSE),
       BREEN = read.csv(file.path(REF_DIR, "File_S1_DEGs__Li_DEGs_BreenOverlap.csv"),
                        stringsAsFactors = FALSE))
}

gene_symbols <- function() {
  f <- file.path(EXP_HOME, "data", "gene_symbols.csv")
  if (!file.exists(f)) return(NULL)
  s <- read.csv(f, stringsAsFactors = FALSE)
  setNames(s$symbol, s$gene)
}

# --------------------------------------------------------------- comparison
# Signed evidence: -log10(p) carrying the sign of the fold change. Rank
# correlation of this is the statistic Krebs themselves use (supplementary
# BMI analysis), so using it here keeps us on their scale.
signed_evidence <- function(logFC, p) sign(logFC) * -log10(p)

# Compare one of our topTables against a deposited sheet. Everything is
# computed on the genes the two have in common -- reported, not assumed.
compare_to_ref <- function(tt, ref, fdr = 0.05) {
  rownames(tt) <- tt$gene
  g <- intersect(ref$gene, tt$gene)
  a <- ref[match(g, ref$gene), ]; b <- tt[g, ]
  ours <- b$gene[b$adj.P.Val < fdr]
  theirs <- a$gene[a$adj.P.Val < fdr]
  sh <- intersect(ours, theirs)
  conc <- if (length(sh)) mean(sign(b[sh, "logFC"]) == sign(a[match(sh, a$gene), "logFC"])) else NA
  data.frame(
    n_shared_genes = length(g),
    n_ours = length(ours), n_theirs = length(theirs), n_overlap = length(sh),
    jaccard = length(sh) / length(union(ours, theirs)),
    # concordance among genes BOTH call significant
    sign_conc_overlap = conc,
    # concordance among all shared genes, and among their DEGs only
    sign_conc_all = mean(sign(b$logFC) == sign(a$logFC)),
    sign_conc_theirDEG = mean(sign(b[theirs, "logFC"]) ==
                              sign(a[match(theirs, a$gene), "logFC"])),
    r_logFC = cor(a$logFC, b$logFC),
    rho_logFC = cor(a$logFC, b$logFC, method = "spearman"),
    rho_signed_evidence = cor(signed_evidence(a$logFC, a$P.Value),
                              signed_evidence(b$logFC, b$P.Value),
                              method = "spearman"),
    # their DEGs recovered by us, and ours recovered by them
    recall_of_theirs = length(sh) / max(1, length(theirs)),
    precision_vs_theirs = length(sh) / max(1, length(ours))
  )
}

# The paper's own internal consistency check: of the 233 post-correction DEGs,
# how many are significant and sign-concordant in the pre-correction model.
# Krebs report 194 (83.2%). It is the sharpest available fingerprint of which
# model produced the 976, because it involves both tables at once.
carryover_from_post <- function(tt, POST, fdr = 0.05) {
  rownames(tt) <- tt$gene
  pd <- POST$gene[POST$adj.P.Val < fdr]
  sh <- intersect(pd, tt$gene)
  ps <- setNames(sign(POST$logFC), POST$gene)
  n <- sum(tt[sh, "adj.P.Val"] < fdr & sign(tt[sh, "logFC"]) == ps[sh])
  c(n_post_deg = length(pd), n_present = length(sh),
    n_carry = n, pct_carry = 100 * n / length(pd))
}

# Replication-list recovery. Krebs validated the 976 against two external
# studies; the deposited overlap files are the gene lists they matched.
external_recovery <- function(tt, lst, fdr = 0.05) {
  rownames(tt) <- tt$gene
  g <- intersect(lst$ensembl_gene_id, tt$gene)
  c(n_listed = nrow(lst), n_present = length(g),
    n_sig = sum(tt[g, "adj.P.Val"] < fdr))
}
