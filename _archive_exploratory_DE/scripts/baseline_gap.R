#!/usr/bin/env Rscript
# ==============================================================================
# baseline_gap.R -- stage two of the baseline: account for what is left of the
# distance between our fits and the published numbers.
#
# Three questions, each answered by fitting rather than by argument:
#
#   A. The case/control arm. The engine's `casecon` returns 1,816 DEGs; the
#      published BD model returns 6. Two things differ -- the sample set and
#      whether lithium use is a covariate. Which one does the work?
#
#   B. The residual 1,023 vs 976 gap on the lithium arm. Krebs' Table 1 counts
#      240 cases / 204 controls where the deposit has 239 / 205, so one subject
#      is a case in the paper and a control in the deposit. Is a single subject
#      enough to move the DEG count by 47? Answered by relabelling one control
#      as a case, repeatedly, and reading off the spread.
#
#   C. Estimator choice. Could the published spec return 976 under a different
#      limma setting? Checked narrowly, on the published spec only.
#
# Writes only inside experimentation/runs/baseline/.
# ==============================================================================

setwd("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation")
source("scripts/common.R"); source("scripts/ref_krebs.R")

CACHE <- "C:/Users/hp/AppData/Local/Temp/claude/baseline_fits"
OUT   <- file.path(EXP_HOME, "runs", "baseline")
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT,   showWarnings = FALSE, recursive = TRUE)
SEED <- 20261218   # fixed here, and nowhere else in this script

p <- load_prep(); K <- krebs_tables(); SYM <- gene_symbols()
sym <- function(g) ifelse(is.na(SYM[g]), g, SYM[g])
qc <- p$meta$qc_pass
CV <- c("age", "sex", "tob", "group", "rin", "plate", "seqpc1", "seqpc2", "seqpc3")
allof <- rep(TRUE, nrow(p$meta))
p$CONTRASTS$li_case     <- list(keep = allof, exposure = "lithium", covars = c("case", CV), note = "")
p$CONTRASTS$casecon_li  <- list(keep = allof, exposure = "case", covars = c("lithium", CV), note = "")
p$CONTRASTS$casecon_nol <- list(keep = allof, exposure = "case", covars = CV, note = "")

fit <- function(key, contrast, ...) {
  f <- file.path(CACHE, paste0(key, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  r <- de_fit(p, contrast, ...); saveRDS(r, f); r
}

# ===========================================================================
# 0. is tobacco even in play on the 444 subset?
# ===========================================================================
cat("== 0. tobacco on the Krebs 444 ==\n")
cat(sprintf("  tobacco observed for all 444 QC-pass subjects: %s\n",
            all(!is.na(p$meta$tob_obs[qc]))))
cat(sprintf("  the %d subjects with missing tobacco are all outside the 444: %s\n",
            sum(is.na(p$meta$tob_obs)), all(!qc[is.na(p$meta$tob_obs)])))
cat("  -> on the 444 the imputation is inert; it cannot explain any 444-based gap.\n")

# ===========================================================================
# A. CASE/CONTROL: what turns 1,816 into 6?
# ===========================================================================
cat("\n== A. case/control decomposition ==\n")
ccA <- list(
  cc_474_noli = fit("cc_474_noli", "casecon_nol"),
  cc_444_noli = fit("cc_444_noli", "casecon_nol", samples = qc),
  cc_474_li   = fit("cc_474_li",   "casecon_li"),
  cc_444_li   = fit("cc_444_li",   "casecon_li", samples = qc))
lab <- c(cc_474_noli = "474, lithium NOT a covariate (engine `casecon`)",
         cc_444_noli = "444, lithium NOT a covariate",
         cc_474_li   = "474, lithium a covariate",
         cc_444_li   = "444, lithium a covariate  <- PUBLISHED BD SPEC")
ccrows <- lapply(names(ccA), function(nm) {
  f <- ccA[[nm]]; tt <- f$tt
  ours <- deg_ids(tt, 0.05); rownames(tt) <- tt$gene
  sh <- intersect(K$BD$gene, tt$gene)
  data.frame(spec = nm, label = lab[[nm]], n = f$n_samples, n_genes = f$n_genes,
    DEG10 = n_deg(tt, .10), DEG05 = length(ours), DEG01 = n_deg(tt, .01),
    pi0 = round(pi0_storey(tt$P.Value), 4),
    overlap_their6 = length(intersect(ours, K$BD$gene)),
    their6_medFDR = signif(median(tt[sh, "adj.P.Val"]), 3))
})
cc <- do.call(rbind, ccrows); print(cc[, -2], row.names = FALSE)
write.csv(cc, file.path(OUT, "casecon_decomposition.csv"), row.names = FALSE)

# which six genes do we call, and are they theirs?
tt6 <- ccA$cc_444_li$tt
o6  <- tt6[order(tt6$adj.P.Val), ][1:8, ]
six <- data.frame(rank = 1:8, gene = o6$gene, symbol = sym(o6$gene),
  logFC = round(o6$logFC, 4), adj.P.Val = signif(o6$adj.P.Val, 3),
  in_krebs6 = o6$gene %in% K$BD$gene,
  krebs_logFC = round(K$BD$logFC[match(o6$gene, K$BD$gene)], 4))
cat("\n  our top 8 in the published BD spec (6 clear FDR .05):\n"); print(six, row.names = FALSE)
kb <- data.frame(gene = K$BD$gene, symbol = sym(K$BD$gene),
  krebs_logFC = round(K$BD$logFC, 4), krebs_adjP = signif(K$BD$adj.P.Val, 3),
  our_logFC = round(tt6$logFC[match(K$BD$gene, tt6$gene)], 4),
  our_adjP  = signif(tt6$adj.P.Val[match(K$BD$gene, tt6$gene)], 3))
cat("\n  their six, as we fit them:\n"); print(kb, row.names = FALSE)
write.csv(rbind(cbind(set = "ours_top8", six[, c("gene","symbol","logFC","adj.P.Val")]),
                cbind(set = "krebs_six", data.frame(gene = kb$gene, symbol = kb$symbol,
                      logFC = kb$our_logFC, adj.P.Val = kb$our_adjP))),
          file.path(OUT, "bd_six_genes.csv"), row.names = FALSE)

# ===========================================================================
# B. ONE-SUBJECT SENSITIVITY on the published lithium spec
# ===========================================================================
# The paper counts one more case and one fewer control than the deposit, and
# its extra case is a lithium non-user (152 users both ways, 88 vs 87
# non-users). Emulate exactly that: promote one control to a non-using case.
cat("\n== B. one-subject sensitivity, published lithium spec ==\n")
BFILE <- file.path(CACHE, "relabel_sens.rds")
if (file.exists(BFILE)) {
  rel <- readRDS(BFILE)
} else {
  set.seed(SEED)
  ctrl <- which(qc & p$meta$dx == "Control")
  picks <- sample(ctrl, 25)
  rel <- data.frame()
  for (i in seq_along(picks)) {
    q <- p; q$meta$case[picks[i]] <- "BD"; q$meta$dx[picks[i]] <- "BP1"
    f <- de_fit(q, "li_case", samples = qc)
    rel <- rbind(rel, data.frame(iter = i, subject = p$meta$title[picks[i]],
      n = f$n_samples, n_genes = f$n_genes, DEG05 = n_deg(f$tt, .05)))
    cat(sprintf("   relabel %2d/%d  %s  DEG=%d\n", i, length(picks),
                p$meta$title[picks[i]], tail(rel$DEG05, 1)))
  }
  saveRDS(rel, BFILE)
}
base1023 <- n_deg(readRDS(file.path(CACHE, "spec_allcase2_444.rds"))$tt, .05)
cat(sprintf("\n  unmodified published spec        : %d DEG\n", base1023))
cat(sprintf("  one control promoted to a case   : median %d, range %d-%d, sd %.1f (n=%d)\n",
    median(rel$DEG05), min(rel$DEG05), max(rel$DEG05), sd(rel$DEG05), nrow(rel)))
cat(sprintf("  paper's figure                   : %d\n", KREBS_PAPER$li_ndeg))
cat(sprintf("  is 976 inside the observed range : %s\n",
    KREBS_PAPER$li_ndeg >= min(rel$DEG05) && KREBS_PAPER$li_ndeg <= max(rel$DEG05)))
write.csv(rel, file.path(OUT, "one_subject_sensitivity.csv"), row.names = FALSE)

# ===========================================================================
# C. ESTIMATOR CHOICE on the published spec
# ===========================================================================
cat("\n== C. estimator variants, published spec only ==\n")
VAR <- list(
  voom            = list(method = "voom",      robust = FALSE),
  voom_robust     = list(method = "voom",      robust = TRUE),
  voomWQW         = list(method = "voomWQW",   robust = FALSE),
  voomWQW_robust  = list(method = "voomWQW",   robust = TRUE),
  trend           = list(method = "trend",     robust = FALSE),
  QLF             = list(method = "QLF",       robust = FALSE))
vrows <- lapply(names(VAR), function(nm) {
  f <- fit(paste0("var_", nm), "li_case", samples = qc,
           method = VAR[[nm]]$method, robust = VAR[[nm]]$robust)
  tt <- f$tt; d <- tt[which(tt$adj.P.Val < .05), ]
  co <- carryover_from_post(tt, K$POST)
  data.frame(variant = nm, DEG10 = n_deg(tt, .10), DEG05 = nrow(d),
    DEG01 = n_deg(tt, .01), pi0 = round(pi0_storey(tt$P.Value), 4),
    absFC_mean = round(mean(abs(d$logFC)), 3), absFC_max = round(max(abs(d$logFC)), 3),
    pct_up = round(100 * mean(d$logFC > 0), 1),
    carry_n = co[["n_carry"]], carry_pct = round(co[["pct_carry"]], 1))
})
vt <- do.call(rbind, vrows)
cat("  paper: DEG05=976  absFC_mean=.20  absFC_max=.82  up=77.3  carry=194 (83.2%)\n")
print(vt, row.names = FALSE)
write.csv(vt, file.path(OUT, "estimator_variants_published_spec.csv"), row.names = FALSE)

cat("\n[baseline_gap] done\n")
