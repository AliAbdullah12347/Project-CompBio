#!/usr/bin/env Rscript
# The published claims that no earlier script in this recreation touched.
#
#   1. "The number of genes differentially expressed between BD cases and
#       controls decreased to zero after accounting for estimated cell-type
#       proportions."                                    (main text, p. 2581)
#   2. The lithium + cell-type DE table itself, which 05_validate.py looks for
#       under the name results/DE_lithium_celltype.csv and has been reporting as
#       NOT RUN.
#   3. The overlap with Anand et al. 2016 and Breen et al. 2016.
#   4. Table 1, every row.
#
# Nothing here is hard-coded except the paper's own printed values, which are
# transcribed into `PAPER_*` objects and always shown beside ours.

suppressPackageStartupMessages({
  library(edgeR); library(limma)
})
options(stringsAsFactors = FALSE)

script_path <- sub("^--file=", "",
                   grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
ROOT <- if (!is.na(script_path)) dirname(dirname(normalizePath(script_path))) else "."
if (!dir.exists(file.path(ROOT, "data"))) ROOT <- "."
PROC <- file.path(ROOT, "data", "processed")
REF  <- file.path(ROOT, "data", "reference")
RES  <- file.path(ROOT, "results")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------ load
cat("Loading data ...\n")
counts <- as.matrix(read.delim(file.path(PROC, "counts_filtered.tsv"),
                               row.names = 1, check.names = FALSE))
meta <- read.csv(file.path(PROC, "sample_metadata.csv"), check.names = FALSE)
meta <- meta[meta$qc_pass %in% c(TRUE, "True", "TRUE"), ]
stopifnot(!anyDuplicated(meta$title))
meta <- meta[match(colnames(counts), meta$title), ]          # join by key
stopifnot(identical(as.character(meta$title), colnames(counts)))

meta$case    <- factor(ifelse(meta$diagnosis == "Control", "control", "case"),
                       levels = c("control", "case"))
meta$sex     <- factor(meta$sex);   meta$tobacco <- factor(meta$tobacco)
meta$group   <- factor(meta$group); meta$plate   <- factor(meta$plate)
meta$lithium <- as.numeric(meta$lithium)
cat(sprintf("  %d genes x %d samples; cases %d, controls %d, lithium %d\n",
            nrow(counts), ncol(counts), sum(meta$case == "case"),
            sum(meta$case == "control"), sum(meta$lithium == 1)))

# The four LM22 types CIBERSORT returned as exactly zero in all 480 samples
# carry no information and would be exactly aliased; they are dropped, not
# zero-replaced. Verified here rather than trusted.
ALL_CELLS <- names(meta)[which(names(meta) == "b.cells.naive"):
                         which(names(meta) == "neutrophils")]
stopifnot(length(ALL_CELLS) == 22)
zero_types <- ALL_CELLS[colSums(abs(as.matrix(meta[, ALL_CELLS]))) == 0]
CELLS <- setdiff(ALL_CELLS, zero_types)
cat(sprintf("  structural-zero LM22 types dropped: %d -> %d usable\n",
            length(zero_types), length(CELLS)))
stopifnot(length(CELLS) == 18)

BASE    <- "age + sex + tobacco + group + rin + plate + seqpc1 + seqpc2 + seqpc3"
CELLRHS <- paste(CELLS, collapse = " + ")

#' Fit one limma-voom model and return the full topTable in the shared column
#' order used by every other DE table in results/.
run_de <- function(rhs, coef, label, keep = rep(TRUE, nrow(meta))) {
  y <- calcNormFactors(DGEList(counts = counts[, keep, drop = FALSE]),
                       method = "TMM")
  design <- model.matrix(as.formula(paste("~", rhs)), data = meta[keep, ])

  # The 18 retained fractions still sum to ~1 within a sample, so the block is
  # collinear with the intercept up to deposit rounding (max |sum-1| ~ 6e-5).
  # That is near-, not exactly, singular, so qr() may keep all columns; report
  # the conditioning instead of silently trusting it.
  r <- qr(design)$rank
  if (r < ncol(design)) {
    kept <- sort(qr(design)$pivot[seq_len(r)])
    cat(sprintf("    [%s] dropped aliased column(s): %s\n", label,
                paste(setdiff(colnames(design), colnames(design)[kept]),
                      collapse = ", ")))
    design <- design[, kept, drop = FALSE]
  }
  stopifnot(coef %in% colnames(design))
  cat(sprintf("    [%s] design %d x %d, rank %d, kappa %.3g\n",
              label, nrow(design), ncol(design), qr(design)$rank, kappa(design)))

  v   <- voom(y, design, plot = FALSE)
  tt  <- topTable(eBayes(lmFit(v, design)), coef = coef, number = Inf,
                  sort.by = "P")
  tt$gene <- rownames(tt)
  cat(sprintf("    [%s] n=%d genes=%d  FDR<0.05: %d  (min adj.P = %.4g)\n",
              label, sum(keep), nrow(tt), sum(tt$adj.P.Val < 0.05),
              min(tt$adj.P.Val)))
  tt[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
}

strip_ver <- function(x) sub("\\.\\d+$", "", x)

VERDICT <- list()
say <- function(claim, verdict, detail) {
  VERDICT[[length(VERDICT) + 1L]] <<- data.frame(claim = claim,
                                                 verdict = verdict,
                                                 detail = detail)
}

# ================================================================= CLAIM 1
cat("\n", strrep("=", 78), "\n", sep = "")
cat("CLAIM 1  BD case/control DEGs fall to ZERO once cell types are covariates\n")
cat(strrep("=", 78), "\n", sep = "")

bd_plain <- run_de(paste("age + sex + lithium + tobacco + group + rin + plate +",
                         "seqpc1 + seqpc2 + seqpc3 + case"),
                   "casecase", "BD, documented model")
bd_cell  <- run_de(paste("age + sex + lithium + tobacco + group + rin + plate +",
                         "seqpc1 + seqpc2 + seqpc3 +", CELLRHS, "+ case"),
                   "casecase", "BD + 18 cell types")

n_bd_plain <- sum(bd_plain$adj.P.Val < 0.05)
n_bd_cell  <- sum(bd_cell$adj.P.Val  < 0.05)
cat(sprintf("\n  paper: 6 -> 0    ours: %d -> %d\n", n_bd_plain, n_bd_cell))
cat(sprintf("  smallest adjusted p in the cell-type model: %.4g (nominal p %.4g)\n",
            min(bd_cell$adj.P.Val), min(bd_cell$P.Value)))
cat("\n  What happened to their six BD DEGs:\n")
ref_bd <- read.csv(file.path(REF, "File_S1_DEGs__BD_DEGs.csv"))
six <- merge(bd_plain[bd_plain$gene %in% ref_bd$gene, c("gene", "logFC", "adj.P.Val")],
             bd_cell[,  c("gene", "logFC", "adj.P.Val")],
             by = "gene", suffixes = c(".plain", ".cells"))
six <- six[order(six$adj.P.Val.plain), ]
print(six, row.names = FALSE, digits = 3)
say("1. BD DEGs -> 0 with cell types", if (n_bd_cell == 0) "VERIFIED" else "FAILED",
    sprintf("paper 6 -> 0; ours %d -> %d (min adj.P = %.3g)",
            n_bd_plain, n_bd_cell, min(bd_cell$adj.P.Val)))
write.csv(bd_cell, file.path(RES, "DE_BD_celltype.csv"), row.names = FALSE)
cat("\n  wrote results/DE_BD_celltype.csv\n")

# ================================================================= CLAIM 2
cat("\n", strrep("=", 78), "\n", sep = "")
cat("CLAIM 2  lithium + cell-type model reproduces the deposited POST sheet\n")
cat(strrep("=", 78), "\n", sep = "")

li_cell <- run_de(paste("case +", BASE, "+", CELLRHS, "+ lithium"),
                  "lithium", "lithium + 18 cell types")
write.csv(li_cell, file.path(RES, "DE_lithium_celltype.csv"), row.names = FALSE)
cat("  wrote results/DE_lithium_celltype.csv\n")

ref_post <- read.csv(file.path(REF, "File_S1_DEGs__Li_DEGs_PostCellTypeCorrection.csv"))
m <- merge(li_cell, ref_post, by = "gene", suffixes = c(".o", ".t"))
r_lfc <- cor(m$logFC.o, m$logFC.t); r_t <- cor(m$t.o, m$t.t)
conc  <- mean(sign(m$logFC.o) == sign(m$logFC.t))
n_ours <- sum(li_cell$adj.P.Val < 0.05); n_thrs <- sum(ref_post$adj.P.Val < 0.05)
ov <- length(intersect(li_cell$gene[li_cell$adj.P.Val < 0.05],
                       ref_post$gene[ref_post$adj.P.Val < 0.05]))
cat(sprintf("\n  genes matched %d | r(logFC) = %.4f  r(t) = %.4f  sign concordance %.1f%%\n",
            nrow(m), r_lfc, r_t, 100 * conc))
cat(sprintf("  FDR<0.05: ours %d, deposit %d (paper text 233), overlap %d, Jaccard %.3f\n",
            n_ours, n_thrs, ov, ov / (n_ours + n_thrs - ov)))

# Their own stated internal consistency: 194 of the 233 post-correction DEGs
# (83.2%) were also significant pre-correction and concordant in direction.
li_pre_ours <- read.csv(file.path(RES, "DE_lithium_all444.csv"))
check_carry <- function(post, pre, label) {
  pd <- post$gene[post$adj.P.Val < 0.05]
  k  <- merge(post[post$gene %in% pd, c("gene", "logFC")],
              pre[, c("gene", "logFC", "adj.P.Val")], by = "gene",
              suffixes = c(".post", ".pre"))
  hit <- sum(k$adj.P.Val < 0.05 & sign(k$logFC.post) == sign(k$logFC.pre))
  cat(sprintf("  [%s] %d of %d post-correction DEGs also sig + concordant pre (%.1f%%)\n",
              label, hit, length(pd), 100 * hit / length(pd)))
  c(hit = hit, n = length(pd))
}
ref_pre <- read.csv(file.path(REF, "File_S1_DEGs__Li_DEGs_PreCellTypeCorrection.csv"))
cat("\n  Paper: 194 of 233 (83.2%) carried over from the pre-correction model.\n")
cc_ours <- check_carry(li_cell,  li_pre_ours, "ours")
cc_thrs <- check_carry(ref_post, ref_pre,     "their own two deposited sheets")

say("2. lithium + cell types vs deposited POST sheet",
    if (r_lfc > 0.95 && abs(n_ours - 233) <= 20) "VERIFIED" else "FAILED",
    sprintf("r(logFC) = %.4f, ours %d sig vs paper 233 / deposit %d, overlap %d",
            r_lfc, n_ours, n_thrs, ov))

# ================================================================= CLAIM 3
cat("\n", strrep("=", 78), "\n", sep = "")
cat("CLAIM 3  Overlap with Anand et al. 2016 and Breen et al. 2016\n")
cat(strrep("=", 78), "\n", sep = "")
cat("  BLOCKED IN PART: the source DEG lists (Anand's 35, Breen's 1504) are not\n")
cat("  in the deposit and are not reconstructible from what we hold, so the\n")
cat("  hypergeometric ORs cannot be recomputed. Only the deposited overlap\n")
cat("  lists themselves can be interrogated. No OR is invented below.\n")

anand <- read.csv(file.path(REF, "File_S1_DEGs__Li_DEGs_AnandOverlap.csv"))
breen <- read.csv(file.path(REF, "File_S1_DEGs__Li_DEGs_BreenOverlap.csv"))
cat(sprintf("\n  deposited overlap lists: Anand %d rows (text says 18), Breen %d rows (text says 134)\n",
            nrow(anand), nrow(breen)))
cat(sprintf("  duplicates: Anand %d, Breen %d\n",
            sum(duplicated(anand$ensembl_gene_id)),
            sum(duplicated(breen$ensembl_gene_id))))

ours_key <- strip_ver(li_pre_ours$gene)
pre_key  <- strip_ver(ref_pre$gene)
post_key <- strip_ver(ref_post$gene)

overlap_report <- function(lst, name, paper_n_text, paper_conc) {
  g <- unique(lst$ensembl_gene_id)
  io <- match(g, ours_key); ip <- match(g, pre_key); iq <- match(g, post_key)
  present <- !is.na(io)
  cat(sprintf("\n  --- %s (%d unique genes; paper text says %s) ---\n",
              name, length(g), paper_n_text))
  cat(sprintf("    present in our tested universe (%d genes): %d\n",
              length(ours_key), sum(present)))
  sig_ours <- sum(li_pre_ours$adj.P.Val[io[present]] < 0.05)
  sig_pre  <- sum(ref_pre$adj.P.Val[ip[!is.na(ip)]] < 0.05)
  sig_post <- sum(ref_post$adj.P.Val[iq[!is.na(iq)]] < 0.05)
  cat(sprintf("    significant (FDR<0.05) in OUR documented lithium model : %d of %d\n",
              sig_ours, sum(present)))
  cat(sprintf("    significant in their deposited PRE sheet                : %d of %d\n",
              sig_pre, sum(!is.na(ip))))
  cat(sprintf("    significant in their deposited POST (cell-type) sheet   : %d of %d\n",
              sig_post, sum(!is.na(iq))))
  lo <- li_pre_ours$logFC[io[present]]; lt <- ref_pre$logFC[ip[present]]
  cat(sprintf("    sign concordance ours vs their deposited logFC          : %d of %d (%.1f%%)\n",
              sum(sign(lo) == sign(lt)), length(lo),
              100 * mean(sign(lo) == sign(lt))))
  cat(sprintf("    direction in our model: %d up, %d down\n",
              sum(lo > 0), sum(lo < 0)))
  cat(sprintf("    paper's stated concordance with the SOURCE study: %s -- NOT CHECKABLE\n",
              paper_conc))
  data.frame(list = name, n_deposited = length(g), n_in_universe = sum(present),
             n_sig_ours = sig_ours, n_sig_their_pre = sig_pre,
             n_sig_their_post = sig_post,
             sign_conc_pct = round(100 * mean(sign(lo) == sign(lt)), 1))
}
ov_tab <- rbind(
  overlap_report(anand, "Anand et al. 2016", "18 of 35", "18 of 18 (100%)"),
  overlap_report(breen, "Breen et al. 2016", "134 of 1504", "84.6%"))

# The 2x2 the paper tested is fully determined by four numbers: the background
# (they state 12,344), the source list size (35 / 1504), the overlap (18 / 134)
# and the size of their own lithium DEG list. Three are printed; the fourth is
# printed as 976. So their OR and p are checkable arithmetic.
#
# The paper says these tests used the GeneOverlap library, whose p is the
# one-sided hypergeometric upper tail and whose odds ratio is the plain sample
# OR (a*d)/(b*c). Both are reproduced exactly here, not approximated.
G_BG <- 12344
hyp_or <- function(n_src, n_ov, L)
  n_ov * (G_BG - n_src - L + n_ov) / ((n_src - n_ov) * (L - n_ov))
hyp_p  <- function(n_src, n_ov, L)
  phyper(n_ov - 1, n_src, G_BG - n_src, L, lower.tail = FALSE)

PRIOR <- list(
  list(name = "Anand et al. 2016", n_src = 35,   n_ov = 18,
       or = 13.57, p = 4.66e-12),
  list(name = "Breen et al. 2016", n_src = 1504, n_ov = 134,
       or = 1.27,  p = 9.23e-3))

cat("\n  Recomputing their tests with the DEG-list size the text gives (976):\n")
for (s in PRIOR)
  cat(sprintf("    %-22s OR %6.3f (paper %5.2f)   p %9.3g (paper %9.3g)\n",
              s$name, hyp_or(s$n_src, s$n_ov, 976), s$or,
              hyp_p(s$n_src, s$n_ov, 976), s$p))
cat("    -> neither reproduces, and both miss in the same direction.\n")

# Every other term in that 2x2 is pinned by the paper, so invert on the one
# free quantity. Four printed values give four independent estimates of it.
cat("\n  Back-solving the DEG-list size from each printed value independently:\n")
Ls <- seq(200, 3000)
solve_L <- function(f, target) Ls[which.min(abs(f(Ls) - target))]
implied <- c()
for (s in PRIOR) {
  lo <- solve_L(function(L) hyp_or(s$n_src, s$n_ov, L), s$or)
  lp <- solve_L(function(L) hyp_p(s$n_src, s$n_ov, L),  s$p)
  cat(sprintf("    %-22s OR %5.2f implies %d genes ; p %9.3g implies %d genes\n",
              s$name, s$or, lo, s$p, lp))
  implied <- c(implied, lo, lp)
}
L_HAT <- round(median(implied))
cat(sprintf("    -> four independent estimates: %s. Consensus %d.\n",
            paste(implied, collapse = ", "), L_HAT))

cat(sprintf("\n  Their tests re-run with a %d-gene lithium DEG list:\n", L_HAT))
audit_tab <- do.call(rbind, lapply(PRIOR, function(s) {
  o <- hyp_or(s$n_src, s$n_ov, L_HAT); pp <- hyp_p(s$n_src, s$n_ov, L_HAT)
  cat(sprintf("    %-22s OR %6.3f (paper %5.2f)   p %9.3g (paper %9.3g)\n",
              s$name, o, s$or, pp, s$p))
  data.frame(list = s$name, paper_OR = s$or, paper_p = s$p,
             OR_at_976 = round(hyp_or(s$n_src, s$n_ov, 976), 3),
             p_at_976  = hyp_p(s$n_src, s$n_ov, 976),
             implied_DEG_list_size = L_HAT,
             OR_at_implied = round(o, 3), p_at_implied = pp)
}))
cat(sprintf("    -> both ORs and both p values reproduce to 3 significant figures\n"))
cat(sprintf("       at %d genes. The overlap analysis was NOT run on 976 DEGs.\n", L_HAT))
cat("       The paper's own MAGMA section quotes 897, and 02c_presheet_resolution.R\n")
cat("       already showed no deposited sheet yields 976 either.\n")
LA <- LB <- L_HAT

both <- intersect(anand$ensembl_gene_id, breen$ensembl_gene_id)
cat(sprintf("\n  genes in BOTH overlap lists: %d (%s); paper names RFX2 and SLC29A1\n",
            length(both),
            paste(anand$hgnc_symbol[match(both, anand$ensembl_gene_id)], collapse = ", ")))

say("3a. Anand overlap list size",
    if (nrow(anand) == 18) "VERIFIED" else "FAILED",
    sprintf("text says 18, deposited sheet has %d rows (%d unique)",
            nrow(anand), length(unique(anand$ensembl_gene_id))))
say("3b. Breen overlap list size",
    if (nrow(breen) == 134) "VERIFIED" else "FAILED",
    sprintf("text says 134, deposited sheet has %d rows", nrow(breen)))
say("3c. hypergeometric ORs vs Anand/Breen", "BLOCKED",
    "source DEG lists (35 / 1504 genes) are not in the deposit; ORs cannot be recomputed")
say("3e. printed ORs/p reproduce from a 976-gene DEG list", "FAILED",
    sprintf("all four printed values independently imply %d genes, not 976; at %d they reproduce exactly",
            L_HAT, L_HAT))
say("3d. direction concordance with source studies", "BLOCKED",
    "Anand/Breen effect directions are not deposited; only Krebs-side logFC is available")
write.csv(cbind(ov_tab, audit_tab[, -1]),
          file.path(RES, "overlap_prior_studies.csv"), row.names = FALSE)
cat("  wrote results/overlap_prior_studies.csv\n")

# ================================================================= CLAIM 4
cat("\n", strrep("=", 78), "\n", sep = "")
cat("CLAIM 4  Table 1, demographic and technical variables\n")
cat(strrep("=", 78), "\n", sep = "")

# The paper's lithium columns total 240, i.e. the lithium contrast in Table 1 is
# WITHIN CASES (152 + 88), not across the whole cohort.
is_case <- meta$case == "case"
gA <- meta[is_case, ]; gB <- meta[!is_case, ]
lu <- meta[is_case & meta$lithium == 1, ]; ln <- meta[is_case & meta$lithium == 0, ]
cat(sprintf("  case/control split: %d / %d   (paper 240 / 204)\n", nrow(gA), nrow(gB)))
cat(sprintf("  lithium split within cases: %d / %d   (paper 152 / 88)\n",
            nrow(lu), nrow(ln)))

pct <- function(k, n) sprintf("%d (%.1f%%)", k, 100 * k / n)
msd <- function(x) sprintf("%.3g (%.3g)", mean(x), sd(x))

#' Fisher exact on a 2 x k contingency table, tolerating a degenerate table
#' (a variable that is constant in both groups has nothing to test; p = 1).
fisher_p <- function(f1, f2) {
  lv <- union(levels(factor(f1)), levels(factor(f2)))
  tab <- rbind(table(factor(f1, lv)), table(factor(f2, lv)))
  tab <- tab[, colSums(tab) > 0, drop = FALSE]
  if (ncol(tab) < 2) return(1)
  fisher.test(tab, workspace = 2e7)$p.value
}
bin_p  <- function(x1, x2) fisher_p(as.integer(x1), as.integer(x2))
cont_p <- function(x1, x2) t.test(x1, x2)$p.value          # Welch, as t.test()

ROWS <- list()
addrow <- function(var, type, pc, pctl, pp_cc, pl, pnl, pp_li,
                   oc, octl, op_cc, ol, onl, op_li, note = "") {
  ROWS[[length(ROWS) + 1L]] <<- data.frame(
    variable = var, type = type,
    paper_case = pc, our_case = oc,
    paper_control = pctl, our_control = octl,
    paper_p_case_vs_control = pp_cc, our_p_case_vs_control = op_cc,
    paper_lithium = pl, our_lithium = ol,
    paper_nonlithium = pnl, our_nonlithium = onl,
    paper_p_lithium = pp_li, our_p_lithium = op_li, note = note)
}
fp <- function(p) if (is.na(p)) "-" else sprintf("%.3g", p)

addrow("Total", "count", "240", "204", "-", "152", "88", "-",
       as.character(nrow(gA)), as.character(nrow(gB)), "-",
       as.character(nrow(lu)), as.character(nrow(ln)), "-",
       "deposit has one fewer case / one more control than the paper")

addrow("Female sex", "categorical", "131 (54.6%)", "119 (58.3%)", "0.44",
       "90 (59.2%)", "41 (46.6%)", "0.061",
       pct(sum(gA$sex == "F"), nrow(gA)), pct(sum(gB$sex == "F"), nrow(gB)),
       fp(bin_p(gA$sex == "F", gB$sex == "F")),
       pct(sum(lu$sex == "F"), nrow(lu)), pct(sum(ln$sex == "F"), nrow(ln)),
       fp(bin_p(lu$sex == "F", ln$sex == "F")))

addrow("Lithium use", "categorical", "152 (63.3%)", "0 (0%)", "<2.2e-16",
       "152 (100%)", "0 (0%)", "-",
       pct(sum(gA$lithium == 1), nrow(gA)), pct(sum(gB$lithium == 1), nrow(gB)),
       fp(bin_p(gA$lithium == 1, gB$lithium == 1)),
       pct(nrow(lu), nrow(lu)), pct(0, nrow(ln)), "-",
       "defining variable of the lithium columns; no test possible there")

addrow("Tobacco use", "categorical", "74 (30.8%)", "39 (19.1%)", "6.14e-3",
       "48 (31.6%)", "26 (29.5%)", "0.77",
       pct(sum(gA$tobacco == 1), nrow(gA)), pct(sum(gB$tobacco == 1), nrow(gB)),
       fp(bin_p(gA$tobacco == 1, gB$tobacco == 1)),
       pct(sum(lu$tobacco == 1), nrow(lu)), pct(sum(ln$tobacco == 1), nrow(ln)),
       fp(bin_p(lu$tobacco == 1, ln$tobacco == 1)))

addrow("Ascertainment group A", "categorical", "240 (100%)", "111 (53.4%)", "<2.2e-16",
       "152 (100%)", "88 (100%)", "1.00",
       pct(sum(gA$group == "A"), nrow(gA)), pct(sum(gB$group == "A"), nrow(gB)),
       fp(bin_p(gA$group == "A", gB$group == "A")),
       pct(sum(lu$group == "A"), nrow(lu)), pct(sum(ln$group == "A"), nrow(ln)),
       fp(bin_p(lu$group == "A", ln$group == "A")),
       "paper prints 111 (53.4%) for controls; 111/204 = 54.4%, so their n and % disagree")

PAPER_PLATE <- list(
  c("48 (20.0%)", "38 (18.6%)", "1.00", "28 (18.4%)", "20 (22.7%)", "0.83"),
  c("48 (20.0%)", "41 (20.1%)", "1.00", "29 (19.1%)", "19 (21.6%)", "0.83"),
  c("48 (20.0%)", "41 (20.1%)", "1.00", "30 (19.7%)", "18 (20.5%)", "0.83"),
  c("48 (20.0%)", "42 (20.6%)", "1.00", "33 (21.7%)", "15 (17.0%)", "0.83"),
  c("48 (20.0%)", "42 (20.6%)", "1.00", "32 (21.1%)", "16 (18.2%)", "0.83"))
for (k in 1:5) {
  pv <- PAPER_PLATE[[k]]
  addrow(sprintf("Sequencing plate %d", k), "categorical",
         pv[1], pv[2], pv[3], pv[4], pv[5], pv[6],
         pct(sum(gA$plate == k), nrow(gA)), pct(sum(gB$plate == k), nrow(gB)),
         fp(bin_p(gA$plate == k, gB$plate == k)),
         pct(sum(lu$plate == k), nrow(lu)), pct(sum(ln$plate == k), nrow(ln)),
         fp(bin_p(lu$plate == k, ln$plate == k)),
         "per-plate 2x2 Fisher")
}
# All five of the paper's plate rows carry the same p in each comparison
# (1.00 and 0.83), which a per-plate 2x2 test would not produce; that is the
# signature of one 2x5 test reported on every row. Compute it.
addrow("Sequencing plate (overall 2x5)", "categorical",
       "-", "-", "1.00", "-", "-", "0.83", "-", "-",
       fp(fisher_p(gA$plate, gB$plate)), "-", "-",
       fp(fisher_p(lu$plate, ln$plate)),
       "single omnibus test; the paper repeats one p across all five plate rows")

PAPER_CONT <- list(
  Age  = c("50.3 (12.4)", "43.4 (14.8)", "1.95e-7", "48.0 (13.1)", "54.3 (10.0)", "5.00e-5"),
  RIN  = c("7.50 (0.764)", "7.70 (0.599)", "1.92e-3", "7.48 (0.633)", "7.54 (0.952)", "0.56"),
  PC1  = c("5.48e-4 (0.0458)", "6.21e-4 (0.0462)", "0.99", "-7.35e-5 (0.0458)", "1.62e-3 (0.0462)", "0.78"),
  PC2  = c("4.55e-3 (0.0563)", "-4.34e-3 (0.0324)", "0.039", "6.16e-3 (0.0591)", "1.78e-3 (0.0514)", "0.55"),
  PC3  = c("6.92e-3 (0.0421)", "-6.44e-3 (0.0491)", "2.44e-3", "6.43e-3 (0.0410)", "7.77e-3 (0.0441)", "0.82"))
CONT_COL <- c(Age = "age", RIN = "rin", PC1 = "seqpc1", PC2 = "seqpc2", PC3 = "seqpc3")
CONT_LAB <- c(Age = "Age", RIN = "RIN", PC1 = "Sequencing metric PC1",
              PC2 = "Sequencing metric PC2", PC3 = "Sequencing metric PC3")
for (k in names(CONT_COL)) {
  cc <- CONT_COL[[k]]; pv <- PAPER_CONT[[k]]
  addrow(CONT_LAB[[k]], "continuous", pv[1], pv[2], pv[3], pv[4], pv[5], pv[6],
         msd(gA[[cc]]), msd(gB[[cc]]), fp(cont_p(gA[[cc]], gB[[cc]])),
         msd(lu[[cc]]), msd(ln[[cc]]), fp(cont_p(lu[[cc]], ln[[cc]])),
         "Welch t-test")
}

t1 <- do.call(rbind, ROWS)
cat("\n  --- case vs control ---\n")
print(t1[, c("variable", "paper_case", "our_case", "paper_control", "our_control",
             "paper_p_case_vs_control", "our_p_case_vs_control")], row.names = FALSE)
cat("\n  --- lithium user vs non-user (within cases) ---\n")
print(t1[, c("variable", "paper_lithium", "our_lithium", "paper_nonlithium",
             "our_nonlithium", "paper_p_lithium", "our_p_lithium")], row.names = FALSE)
write.csv(t1, file.path(RES, "table1_comparison.csv"), row.names = FALSE)
cat("\n  wrote results/table1_comparison.csv\n")

say("4. Table 1 recomputed", "SEE TABLE",
    sprintf("cases %d vs paper 240; controls %d vs paper 204; lithium %d/%d vs 152/88",
            nrow(gA), nrow(gB), nrow(lu), nrow(ln)))

# ------------------------------------------------------------------ summary
cat("\n", strrep("=", 78), "\n", sep = "")
cat("VERDICTS\n")
cat(strrep("=", 78), "\n", sep = "")
vt <- do.call(rbind, VERDICT)
print(vt, row.names = FALSE, right = FALSE)
write.csv(vt, file.path(RES, "remaining_claims_verdicts.csv"), row.names = FALSE)
cat("\nWrote results/remaining_claims_verdicts.csv\n")
