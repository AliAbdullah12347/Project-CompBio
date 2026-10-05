#!/usr/bin/env Rscript
# ==============================================================================
# baseline.R -- the reference differential-expression fits, and a replication
# against Krebs et al. (2020), the paper this dataset comes from.
#
# Four things happen here, in order:
#   1. Canonical fits for the three engine contrasts, with every summary the
#      rest of the experiments will be compared against.
#   2. Replication against the two deposited Krebs DE tables.
#   3. A factorial decomposition of the gap between our count and the
#      published 976, varying one modelling choice at a time.
#   4. Design diagnostics -- in particular whether `assessment group` really is
#      constant among cases, which decides whether the published model could
#      have contained it.
#
# Writes only inside experimentation/runs/baseline/.
# Fits are cached in the scratchpad so a rerun is cheap; delete the cache to
# force a refit.
# ==============================================================================

setwd("C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation/experimentation")
source("scripts/common.R")
source("scripts/ref_krebs.R")

CACHE <- "C:/Users/hp/AppData/Local/Temp/claude/baseline_fits"
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
OUT <- file.path(EXP_HOME, "runs", "baseline")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

p   <- load_prep()
K   <- krebs_tables()
SYM <- gene_symbols()
sym <- function(g) ifelse(is.na(SYM[g]), g, SYM[g])

# ---------------------------------------------------------------------------
# Extra contrast definitions. These live in the in-memory copy of `prep` only;
# data/prep.rds is never rewritten, so the other experiments are unaffected.
#
# The engine ships `lithium_krebs` as "all cases, diagnosis retained". The
# supplementary methods say something different: "Although there were no
# controls being treated with lithium, diagnosis was included in the
# lithium-use comparison to account for BD-effects within the non-lithium
# using group." A sentence remarking that no control took lithium is only
# meaningful if controls were in the model. So the published lithium fit was
# run on the WHOLE cohort, not on cases. `li_case` and `li_dx` below are that
# model under the two ways "BD diagnosis" can be coded.
# ---------------------------------------------------------------------------
CV <- c("age", "sex", "tob", "group", "rin", "plate", "seqpc1", "seqpc2", "seqpc3")
allof <- rep(TRUE, nrow(p$meta))
p$CONTRASTS$li_case <- list(keep = allof, exposure = "lithium",
  covars = c("case", CV), note = "published spec, diagnosis as 2-level case/control")
p$CONTRASTS$li_dx <- list(keep = allof, exposure = "lithium",
  covars = c("dx", CV), note = "published spec, diagnosis as 3-level Control/BP1/BP2")
# Krebs' case/control model lists lithium use among its covariates; the engine's
# `casecon` omits it. Both are fitted so the 6 published BD DEGs can be checked.
p$CONTRASTS$casecon_li <- list(keep = allof, exposure = "case",
  covars = c("lithium", CV), note = "published BD spec: lithium retained as covariate")

qc    <- p$meta$qc_pass          # Krebs' own 444-subject flag
noBP2 <- p$meta$dx != "BP2"

fit <- function(key, contrast, ...) {
  f <- file.path(CACHE, paste0(key, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  r <- de_fit(p, contrast, ...)
  saveRDS(r, f); r
}

# ===========================================================================
# 1. CANONICAL FITS
# ===========================================================================
cat("== 1. canonical fits ==\n")
CANON <- list(
  lithium       = fit("canon_lithium",       "lithium"),
  casecon       = fit("canon_casecon",       "casecon"),
  lithium_krebs = fit("canon_lithium_krebs", "lithium_krebs"))

fc_dist <- function(x) {
  q <- quantile(x, c(0, .01, .25, .5, .75, .99, 1))
  setNames(round(q, 4), c("min", "p01", "q25", "median", "q75", "p99", "max"))
}

canon_rows <- lapply(names(CANON), function(nm) {
  f <- CANON[[nm]]; tt <- f$tt; d <- tt[which(tt$adj.P.Val < 0.05), ]
  q  <- fc_dist(tt$logFC); qa <- fc_dist(abs(tt$logFC))
  data.frame(contrast = nm, n_samples = f$n_samples, n_genes = f$n_genes,
    n_params = ncol(f$design$X), resid_df = f$n_samples - ncol(f$design$X),
    DEG_fdr10 = n_deg(tt, 0.10), DEG_fdr05 = n_deg(tt, 0.05),
    DEG_fdr01 = n_deg(tt, 0.01),
    DEG_fdr05_lfc0.1 = n_deg(tt, 0.05, 0.1),
    DEG_fdr05_lfc0.2 = n_deg(tt, 0.05, 0.2),
    pi0 = round(pi0_storey(tt$P.Value), 4),
    est_true_pos = round((1 - pi0_storey(tt$P.Value)) * f$n_genes),
    pct_up_DEG = round(100 * mean(d$logFC > 0), 1),
    absFC_mean_DEG = round(mean(abs(d$logFC)), 4),
    absFC_max_DEG  = round(max(abs(d$logFC)), 4),
    absFC_sd_DEG   = round(sd(abs(d$logFC)), 4),
    logFC_min = q[["min"]], logFC_q25 = q[["q25"]], logFC_median = q[["median"]],
    logFC_q75 = q[["q75"]], logFC_max = q[["max"]],
    absFC_median_all = qa[["median"]], absFC_p99_all = qa[["p99"]],
    dropped_constant = paste(f$design$dropped_constant, collapse = ";"),
    dropped_aliased  = paste(f$design$dropped_aliased,  collapse = ";"))
})
canon <- do.call(rbind, canon_rows)
print(canon[, c("contrast","n_samples","n_genes","resid_df","DEG_fdr10","DEG_fdr05",
                "DEG_fdr01","pi0","pct_up_DEG","absFC_max_DEG")], row.names = FALSE)
write.csv(canon, file.path(OUT, "canonical_summary.csv"), row.names = FALSE)

# top 30 by adjusted p, per contrast
top30 <- do.call(rbind, lapply(names(CANON), function(nm) {
  tt <- CANON[[nm]]$tt
  o  <- order(tt$adj.P.Val, tt$P.Value)[1:30]
  data.frame(contrast = nm, rank = 1:30, gene = tt$gene[o], symbol = sym(tt$gene[o]),
    logFC = round(tt$logFC[o], 4), AveExpr = round(tt$AveExpr[o], 3),
    t = round(tt$t[o], 3), P.Value = signif(tt$P.Value[o], 3),
    adj.P.Val = signif(tt$adj.P.Val[o], 3))
}))
write.csv(top30, file.path(OUT, "top30_per_contrast.csv"), row.names = FALSE)

# full logFC distribution, binned, so the shape is on record without a figure
bins <- seq(-1.1, 1.1, by = 0.05)
hist_rows <- do.call(rbind, lapply(names(CANON), function(nm) {
  tt <- CANON[[nm]]$tt
  h  <- table(cut(tt$logFC, bins))
  data.frame(contrast = nm, bin_lo = head(bins, -1), bin_hi = tail(bins, -1),
             n = as.integer(h))
}))
write.csv(hist_rows, file.path(OUT, "logFC_histogram.csv"), row.names = FALSE)

# ===========================================================================
# 2. REPLICATION AGAINST THE DEPOSITED KREBS TABLES
# ===========================================================================
cat("\n== 2. replication against deposited tables ==\n")
cat(sprintf("deposited tables: %d genes each, identical gene sets: %s\n",
            nrow(K$PRE), identical(sort(K$PRE$gene), sort(K$POST$gene))))
cat(sprintf("  PRE  DEG @ FDR .05/.01/.10 : %d / %d / %d\n",
  sum(K$PRE$adj.P.Val < .05), sum(K$PRE$adj.P.Val < .01), sum(K$PRE$adj.P.Val < .10)))
cat(sprintf("  POST DEG @ FDR .05/.01/.10 : %d / %d / %d   (paper reports 233)\n",
  sum(K$POST$adj.P.Val < .05), sum(K$POST$adj.P.Val < .01), sum(K$POST$adj.P.Val < .10)))

# All candidate specifications, one modelling choice varied at a time.
SPECS <- list(
  bp1_474        = list(c = "lithium",       s = NULL,          lab = "BP1 only (engine `lithium`)"),
  cases_474      = list(c = "lithium_krebs", s = NULL,          lab = "cases only (engine `lithium_krebs`)"),
  cases_444      = list(c = "lithium_krebs", s = qc,            lab = "cases only, Krebs 444"),
  allcase2_474   = list(c = "li_case",       s = NULL,          lab = "full cohort, diagnosis 2-level"),
  allcase2_444   = list(c = "li_case",       s = qc,            lab = "full cohort, diagnosis 2-level, Krebs 444  <- PUBLISHED SPEC"),
  alldx3_474     = list(c = "li_dx",         s = NULL,          lab = "full cohort, diagnosis 3-level"),
  alldx3_444     = list(c = "li_dx",         s = qc,            lab = "full cohort, diagnosis 3-level, Krebs 444"),
  allcase2_444_noBP2 = list(c = "li_case",   s = qc & noBP2,    lab = "full cohort minus BP2, 2-level"),
  alldx3_444_noBP2   = list(c = "li_dx",     s = qc & noBP2,    lab = "full cohort minus BP2, 3-level (must equal above)")
)
FITS <- lapply(names(SPECS), function(nm)
  fit(paste0("spec_", nm), SPECS[[nm]]$c, samples = SPECS[[nm]]$s))
names(FITS) <- names(SPECS)

rep_rows <- lapply(names(FITS), function(nm) {
  f <- FITS[[nm]]; tt <- f$tt; d <- tt[which(tt$adj.P.Val < 0.05), ]
  cp <- compare_to_ref(tt, K$PRE); cq <- compare_to_ref(tt, K$POST)
  co <- carryover_from_post(tt, K$POST)
  an <- external_recovery(tt, K$ANAND); br <- external_recovery(tt, K$BREEN)
  data.frame(spec = nm, label = SPECS[[nm]]$lab, n = f$n_samples,
    n_genes = f$n_genes, resid_df = f$n_samples - ncol(f$design$X),
    DEG05 = nrow(d),
    absFC_mean = round(mean(abs(d$logFC)), 3), absFC_max = round(max(abs(d$logFC)), 3),
    absFC_sd = round(sd(abs(d$logFC)), 3), pct_up = round(100 * mean(d$logFC > 0), 1),
    # vs PRE sheet
    PRE_overlap = cp$n_overlap, PRE_jaccard = round(cp$jaccard, 3),
    PRE_r_logFC = round(cp$r_logFC, 3),
    PRE_rho_signed = round(cp$rho_signed_evidence, 3),
    PRE_sign_conc_overlap = round(100 * cp$sign_conc_overlap, 1),
    PRE_recall = round(100 * cp$recall_of_theirs, 1),
    # vs POST sheet
    POST_overlap = cq$n_overlap, POST_r_logFC = round(cq$r_logFC, 3),
    POST_rho_signed = round(cq$rho_signed_evidence, 3),
    POST_sign_conc_overlap = round(100 * cq$sign_conc_overlap, 1),
    # the paper's own internal check: 194/233 (83.2%)
    carry_n = co[["n_carry"]], carry_pct = round(co[["pct_carry"]], 1),
    anand_sig = an[["n_sig"]], anand_n = an[["n_present"]],
    breen_sig = br[["n_sig"]], breen_n = br[["n_present"]])
})
repl <- do.call(rbind, rep_rows)
print(repl[, c("spec","n","n_genes","DEG05","absFC_mean","absFC_max","pct_up",
               "carry_n","carry_pct","PRE_r_logFC","POST_r_logFC")], row.names = FALSE)
write.csv(repl, file.path(OUT, "replication_vs_krebs.csv"), row.names = FALSE)

# ---- two checks that are stronger than the summary statistics ---------------
# (a) Matching their gene COUNT is weak; matching the gene SET is the real test
#     of whether our filter is their filter. (b) Sign concordance is quoted as
#     100% above, which is the sort of number that is usually a bug, so it is
#     recomputed here from scratch and the disagreement count is printed.
cat("\n-- filter reproduces their gene set? --\n")
setrows <- lapply(names(FITS), function(nm) {
  g <- FITS[[nm]]$tt$gene
  data.frame(spec = nm, n = length(g),
    identical_to_deposited = setequal(g, K$PRE$gene),
    ours_not_theirs = length(setdiff(g, K$PRE$gene)),
    theirs_not_ours = length(setdiff(K$PRE$gene, g)))
})
st <- do.call(rbind, setrows); print(st, row.names = FALSE)
write.csv(st, file.path(OUT, "filter_gene_set_identity.csv"), row.names = FALSE)

cat("\n-- sign concordance, recomputed with disagreements counted --\n")
scrows <- list()
for (nm in names(FITS)) for (refnm in c("PRE", "POST")) {
  tt <- FITS[[nm]]$tt; rownames(tt) <- tt$gene
  R <- K[[refnm]]; g <- intersect(R$gene, tt$gene)
  a <- R[match(g, R$gene), ]; b <- tt[g, ]
  sh <- intersect(b$gene[b$adj.P.Val < .05], a$gene[a$adj.P.Val < .05])
  same <- sign(b[sh, "logFC"]) == sign(a[match(sh, a$gene), "logFC"])
  scrows[[paste(nm, refnm)]] <- data.frame(spec = nm, sheet = refnm,
    n_co_significant = length(sh), n_agree = sum(same), n_disagree = sum(!same),
    pct_agree = round(100 * mean(same), 2),
    n_all_shared = length(g),
    pct_agree_all_shared = round(100 * mean(sign(b$logFC) == sign(a$logFC)), 1))
}
sc <- do.call(rbind, scrows)
print(sc[sc$sheet == "POST", ], row.names = FALSE)
cat(sprintf("  total sign disagreements across all specs and both sheets: %d\n",
            sum(sc$n_disagree)))
write.csv(sc, file.path(OUT, "sign_concordance.csv"), row.names = FALSE)

# ---- is the PRE sheet the lithium main effect at all? ----------------------
# REPRODUCIBILITY.md concludes it is not: it is case+lithium summed. Checked
# here independently by asking which of OUR coefficients the sheet tracks.
cat("\n-- what does the deposited PRE sheet actually contain? --\n")
fD <- FITS$alldx3_444
rn <- fD$tt; rownames(rn) <- rn$gene
g  <- intersect(K$PRE$gene, rn$gene)
cat(sprintf("  PRE vs our lithium coef      r = %.3f\n",
            cor(K$PRE$logFC[match(g, K$PRE$gene)], rn[g, "logFC"])))
cat(sprintf("  PRE max|logFC| = %.3f   paper says its 976 DEGs max at 0.82\n",
            max(abs(K$PRE$logFC))))
cat(sprintf("  POST max|logFC| = %.3f  POST n at FDR.05 = %d (paper: 233)\n",
            max(abs(K$POST$logFC)), sum(K$POST$adj.P.Val < .05)))

# ---- the published BD case/control model ----------------------------------
cat("\n-- BD case/control replication (paper: 6 DEGs) --\n")
bd <- list(casecon_engine = fit("bd_engine", "casecon"),
           casecon_publspec = fit("bd_publspec", "casecon_li", samples = qc))
bd_rows <- lapply(names(bd), function(nm) {
  f <- bd[[nm]]; tt <- f$tt
  ours <- deg_ids(tt, 0.05); theirs <- K$BD$gene
  rownames(tt) <- tt$gene; sh <- intersect(theirs, tt$gene)
  data.frame(spec = nm, n = f$n_samples, n_genes = f$n_genes,
    DEG05 = length(ours), overlap_with_their6 = length(intersect(ours, theirs)),
    their6_present = length(sh),
    their6_median_ourFDR = signif(median(tt[sh, "adj.P.Val"]), 3),
    their6_sign_conc = round(100 * mean(sign(tt[sh, "logFC"]) ==
                        sign(K$BD$logFC[match(sh, K$BD$gene)])), 1))
})
bdt <- do.call(rbind, bd_rows); print(bdt, row.names = FALSE)
write.csv(bdt, file.path(OUT, "bd_replication.csv"), row.names = FALSE)

# ===========================================================================
# 3. DESIGN DIAGNOSTICS -- the assessment-group question
# ===========================================================================
cat("\n== 3. design diagnostics ==\n")
m <- p$meta
tab <- table(m$dx, m$group)
cat("diagnosis x assessment group (474):\n"); print(tab)
cases <- m$dx != "Control"
cat(sprintf("\ngroup levels among the %d cases: %s  -> CONSTANT: %s\n",
  sum(cases), paste(names(which(table(droplevels(m$group[cases])) > 0)), collapse = ","),
  length(unique(m$group[cases])) == 1))
cat(sprintf("group levels among cases, Krebs 444 subset (%d): %s\n",
  sum(cases & qc), paste(unique(as.character(m$group[cases & qc])), collapse = ",")))
cat(sprintf("group B occurs only among controls: %s  (B: %d control, %d case)\n",
  all(m$dx[m$group == "B"] == "Control"), sum(m$group == "B" & !cases),
  sum(m$group == "B" & cases)))
cat(sprintf("\nlithium by diagnosis:\n")); print(table(m$dx, m$lithium))

diag_tab <- data.frame(
  check = c("group constant among all cases (474)",
            "group constant among cases in Krebs 444",
            "group B exists only among controls",
            "lithium=1 exists only among BP1",
            "n cases 474", "n cases 444", "n BP2 474", "n BP2 444",
            "BP1 lithium users", "BP1 lithium non-users"),
  value = c(length(unique(m$group[cases])) == 1,
            length(unique(m$group[cases & qc])) == 1,
            all(m$dx[m$group == "B"] == "Control"),
            all(m$dx[m$lithium == 1] == "BP1"),
            sum(cases), sum(cases & qc), sum(m$dx == "BP2"),
            sum(m$dx == "BP2" & qc),
            sum(m$dx == "BP1" & m$lithium == 1),
            sum(m$dx == "BP1" & m$lithium == 0)))
write.csv(diag_tab, file.path(OUT, "design_diagnostics.csv"), row.names = FALSE)

# ===========================================================================
# 4. GAP DECOMPOSITION -- write the ladder as a table
# ===========================================================================
ladder <- data.frame(
  step = c("Krebs published (Results text)",
           "Krebs published (MAGMA methods text)",
           "deposited PRE sheet @ FDR .05",
           "deposited POST sheet @ FDR .05",
           "our reproduction of the published spec",
           "  + 3-level diagnosis instead of 2-level",
           "  + drop BP2 entirely",
           "  + our 474 cohort instead of Krebs 444",
           "engine `lithium_krebs` (cases only, 474)",
           "engine `lithium` (BP1 only, 474)"),
  n_DEG = c(KREBS_PAPER$li_ndeg, KREBS_PAPER$li_ndeg_magma,
            sum(K$PRE$adj.P.Val < .05), sum(K$POST$adj.P.Val < .05),
            n_deg(FITS$allcase2_444$tt), n_deg(FITS$alldx3_444$tt),
            n_deg(FITS$allcase2_444_noBP2$tt), n_deg(FITS$alldx3_474$tt),
            n_deg(FITS$cases_474$tt), n_deg(FITS$bp1_474$tt)))
cat("\n== 4. gap ladder ==\n"); print(ladder, row.names = FALSE)
write.csv(ladder, file.path(OUT, "gap_ladder.csv"), row.names = FALSE)

saveRDS(list(canon = canon, repl = repl, ladder = ladder, diag = diag_tab),
        file.path(CACHE, "baseline_objects.rds"))
cat("\n[baseline] wrote runs/baseline/\n")
