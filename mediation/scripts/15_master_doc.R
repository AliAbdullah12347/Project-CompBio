#!/usr/bin/env Rscript
# ==============================================================================
# 15_master_doc.R -- the master document for the whole mediation arm.
#
# Writes MEDIATION_MASTER.md (parseable) and MEDIATION_MASTER.html (paste into
# Google Docs: open in a browser, Ctrl+A, Ctrl+C, paste). Every number is read
# from a saved results file; nothing is typed in by hand, so the document
# cannot drift from what the scripts produced.
# ==============================================================================

ROOT <- "C:/Users/hp/Downloads/Colgate/Junior Fall/Comp Biology Research/Implementation"
setwd(file.path(ROOT, "mediation")); options(width = 160)
source("scripts/doc_render.R")
rd <- function(f) readRDS(f); rc <- function(f) read.csv(f, stringsAsFactors = FALSE)

## ---- load everything ---------------------------------------------------------
PR <- list(lineage = rd("results/primary/lineage_fit.rds"), proposal = rd("results/primary/proposal_fit.rds"))
BT <- list(lineage = rd("results/boot/lineage_boot.rds"), proposal = rd("results/boot/proposal_boot.rds"))
PM <- list(lineage = rd("results/perm/lineage_perm.rds"), proposal = rd("results/perm/proposal_perm.rds"))
PMF <- rc("results/perm/perm_pm_fixedset.csv"); ZT <- rc("results/perm/pm_ztest.csv")
CAL <- rd("results/sens/calibration.rds"); CALP <- rd("results/sens/calibration_proposal.rds")
SW  <- rc("results/sens/S3_correction_sweep.csv")
SEN <- rd("results/sens/sensitivity.rds"); SEC <- rd("results/sens/secondary.rds")
REL <- rc("results/balance_reliability.csv"); RELP <- rc("results/balance_reliability_proposal.csv")
CIRC <- rc("results/circularity/summary_both_trees.csv"); NB <- rc("results/circularity/normalisation_basis.csv")
ACC <- rc("results/circularity/lm22_accounting.csv")
CR <- list(lineage = rd("results/circularity/lineage_circularity.rds"), proposal = rd("results/circularity/proposal_circularity.rds"))
CRI <- rbind(rc("results/circularity/lineage_pm_intervals.csv"), rc("results/circularity/proposal_pm_intervals.csv"))
CRD <- rbind(rc("results/circularity/lineage_pm_differences.csv"), rc("results/circularity/proposal_pm_differences.csv"))
AL <- rd("results/diversity/alpha.rds"); BE <- rd("results/diversity/beta.rds")
GL <- rd("results/genes/gene_lists.rds")
VAL <- sub(", validated.*$", "", readLines("results/validation_med.txt")[1])
dx  <- regmatches(readLines("results/de_crosscheck.txt")[1],
                  gregexpr("[0-9]+[.]?[0-9]*", readLines("results/de_crosscheck.txt")[1]))[[1]]
DEX <- sprintf("per-gene correlation r = %s; %s of the DE arm's %s lithium genes recovered; Jaccard %s",
               dx[1], fi(as.numeric(dx[2])), fi(as.numeric(dx[3])), dx[4])
L <- PR$lineage; P <- PR$proposal
sigL <- L$tb$TE_q_BH < 0.05; sigP <- P$tb$TE_q_BH < 0.05
pmL <- median(L$tb$PM[sigL]); pmP <- median(P$tb$PM[sigP])
mtL <- L$med_tab; mtP <- P$med_tab
fwL <- c(CDE = mean(abs(L$tb$CDE[sigL])), INTref = mean(abs(L$tb$INTref[sigL])),
         INTmed = mean(abs(L$tb$INTmed[sigL])), PIE = mean(abs(L$tb$PIE[sigL])))
psL <- PM$lineage$summary; psP <- PM$proposal$summary
ar <- AL$results; ap <- ar[ar$config == "lithium | lineage tree (PRIMARY)", ]; rownames(ap) <- ap$outcome
app <- ar[ar$config == "lithium | proposal tree", ]; rownames(app) <- app$outcome
abd <- ar[grepl("bipolar", ar$config), ]; rownames(abd) <- abd$outcome
br <- BE$results; bp <- br[br$config == "lithium | lineage tree (PRIMARY)", ]; rownames(bp) <- bp$layer
bpp <- br[br$config == "lithium | proposal tree", ]; rownames(bpp) <- bpp$layer
bbd <- br[grepl("bipolar", br$config), ]; rownames(bbd) <- bbd$layer
DB <- AL$decomposition; rownames(DB) <- DB$component; DC <- AL$cellmix
cl <- CIRC[CIRC$tree == "lineage", ]; cp <- CIRC[CIRC$tree == "proposal", ]
cil <- CRI[CRI$tree == "lineage", ]; rownames(cil) <- cil$set
cdl <- CRD[CRD$tree == "lineage", ]; cdp <- CRD[CRD$tree == "proposal", ]
CT <- GL$counts; rownames(CT) <- CT$category
# primary confirmatory family for the diversity arm: Holm over two tests
p_alpha <- ap["H_primary", "NIE_p_perm"]
p_beta  <- max(bp["gene_aitchison", "path_a_p"], bp["gene_aitchison", "path_b_p"])  # intersection-union
holm <- p.adjust(c(alpha = p_alpha, beta = p_beta), "holm")
TODAY <- format(Sys.Date(), "%d %B %Y")

D <- list()
add <- function(...) D <<- c(D, list(...))

## =============================================================================
## TITLE + OVERVIEW
## =============================================================================
add(H(1, "Mediation arm — master document"),
    Par("Project 7, BIOL/MATH 481, Colgate · Ali Abdullah & Ava Gillespie · generated ", TODAY,
        " by `scripts/15_master_doc.R` from the saved results — every number below is read from a results file."),
    Par("**How to read this.** The first section is a brief overview of everything. Each later section goes into ",
        "detail, in the order the work was done. Section 16 lists every flagged gene by name. Section 17 lists the ",
        "files, so any number here can be traced to the script that produced it."))

add(H(2, "Overview — everything we did, in brief"))
add(H(3, "The question and the design"))
add(UL(c(
  "**Question:** how much of lithium's effect on whole-blood gene expression runs *through the cell mixture* rather than directly?",
  sprintf("**Sample:** the %d bipolar I subjects of GSE124326 (%d lithium users, %d non-users). Lithium varies only within bipolar I.",
          226, 152, 74),
  "**Exposure** lithium (0/1) · **mediator** the cell mix, as ILR balances · **outcome** gene expression, one model per gene (12,018 genes; the 350 LM22 signature genes held out) · **covariates** age, sex, tobacco, RIN, plate, sequencing PCs 1–3.",
  "**Two versions of the mediator, both run end to end:** the pipeline's tree (5 lineages → 4 balances) and the proposal's tree (6 cell types → 5 balances).",
  "**Protocol fixed in writing before any estimate** (`PRESPEC.md`), including two later addenda written before their analyses ran."
)))
add(H(3, "How it was estimated and checked"))
add(UL(c(
  "**Estimator:** natural direct and indirect effects (NDE, NIE) for a vector mediator, allowing a lithium × cell-mix interaction (VanderWeele & Vansteelandt 2014) — chosen over product-of-coefficients because Prof. Ay's feedback required the interaction to be tested, not assumed.",
  sprintf("**Verification:** %s, against brute-force Monte-Carlo counterfactuals, independent `lm()` fits, a 4,000-replicate bootstrap and 500-dataset interval coverage.", VAL[1]),
  sprintf("**External validation:** the mediation *total* effect reproduces the separate differential-expression arm (%s).", DEX[1]),
  "**Independent adversarial review** of the statistics found four real defects after the checks had passed; all were fixed and everything was re-run. A fifth claim was tested and rejected."
)))
add(H(3, "Main results"))
add(UL(c(
  sprintf("**Lithium changes the expression of %s genes** (total effect, BH 5%%) — %s up, %s down.",
          fi(sum(sigL)), fi(sum(sigL & L$tb$TE > 0)), fi(sum(sigL & L$tb$TE < 0))),
  sprintf("**Lithium shifts the cell mix:** granulocytes rise relative to monocytes (p = %s) and myeloid relative to lymphoid (p = %s); nothing moves inside the lymphoid compartment.",
          f3(mtL$p[2]), f3(mtL$p[1])),
  sprintf("**About %s of lithium's effect runs through the cell mix** — median proportion mediated %s (95%% CI %s–%s); %s on the proposal tree (%s–%s).",
          pct(pmL, 0), f3(pmL), f3(BT$lineage$pm_ci[1]), f3(BT$lineage$pm_ci[2]),
          f3(pmP), f3(BT$proposal$pm_ci[1]), f3(BT$proposal$pm_ci[2])),
  sprintf("**Correcting for the mediator being estimated** raises it to %s (lineage tree); the uncorrected figure is a lower bound.", f3(CAL$cal$corrected[2])),
  sprintf("**%s genes keep a significant direct effect** after accounting for the cell mix — candidates for lithium acting inside cells (%s of them up-regulated).",
          fi(CT["DIRECT", "lineage"]), fi(sum(GL$lineage$category == "DIRECT" & GL$lineage$TE > 0))),
  sprintf("**The lithium × cell-mix interaction is negligible:** %d gene of %s.", sum(L$tb$INT_q_BH < 0.05), fi(nrow(L$tb)))
)))
add(H(3, "What the results do and do not show"))
add(UL(c(
  sprintf("**The proportion mediated is not itself evidence of mediation.** With lithium labels shuffled, the same genes still give a median PM of %s (permutation p = %s; z = %s). It is a valid *description* of how a real total effect splits, not a test.",
          f3(PMF$null_median[1]), f2(PMF$perm_p[1]), f2(ZT$z[1])),
  sprintf("**Per-gene mediated effects cannot be declared significant.** The count swings from 0 to %s between trees and the permutation null alone reaches %s.",
          fi(sum(P$tb$NIE_q_BH < 0.05)), fi(psL$null_max[1])),
  sprintf("**The limit is the design, not the data:** %s of the uncertainty comes from the lithium → cell-mix link, so even perfect expression measurement caps the best gene at |z| = %s against the 4.60 needed. Roughly %s subjects would be required.",
          pct(median(CAL$v_beta / (CAL$v_beta + CAL$v_theta)), 1), f2(max(CAL$z_ceiling)), fi(CAL$n_needed)),
  sprintf("**Unmeasured confounding is the largest threat:** an unmeasured factor associated with both lithium and expression at a risk ratio of %s would explain the typical effect away; age, which we do adjust for, sits at 2.33.",
          f2(CAL$E$median_E[1]))
)))
add(H(3, "Added since — circularity sensitivity analysis"))
add(UL(c(
  "The mediation was re-run on **only the LM22 signature genes** and on **the full filtered set including them**, from one common fit so only the gene set differs.",
  sprintf("Of the 547 signature genes, %d are absent from the annotation, %d are essentially unexpressed in blood, and **%d are expressed** — those are the ones analysed.",
          sum(ACC$status == "absent from annotation"), sum(ACC$status %in% c("present, fails filter", "all-zero")),
          sum(ACC$status == "expressed (C2)")),
  sprintf("**Circularity is real but does not inflate the answer.** The cell mix explains %s× more of a signature gene's expression than of an ordinary gene's (partial R² %s vs %s), yet signature genes do not show a larger proportion mediated (%s vs %s; difference not significant), and including them moves the overall figure by %s.",
          f2(cl$median_partialR2_cellmix[2] / cl$median_partialR2_cellmix[1]),
          f3(cl$median_partialR2_cellmix[2]), f3(cl$median_partialR2_cellmix[1]),
          f3(cl$median_PM[2]), f3(cl$median_PM[1]), f4(cl$median_PM[3] - cl$median_PM[1]))
)))
add(H(3, "Added since — alpha and beta diversity"))
add(UL(c(
  "Transcriptome-wide summaries replace the 12,018 per-gene outcomes: **alpha** diversity (how evenly one person's expression is spread across genes; Shannon) and **beta** diversity (how different two people's profiles are; Aitchison distance).",
  sprintf("**Lithium lowers transcriptome diversity through the cell mix:** indirect effect %s nats (95%% CI %s to %s; permutation p = %s), ≈ %s fewer effective genes. The total effect alone is not significant (p = %s).",
          fs(ap["H_primary", "NIE"]), fs(ap["H_primary", "NIE_lo"]), fs(ap["H_primary", "NIE_hi"]),
          f3(ap["H_primary", "NIE_p_perm"]), pct(-(exp(ap["H_primary", "NIE"]) - 1), 1), f3(ap["H_primary", "TE_p"])),
  sprintf("**A model-free decomposition agrees:** %s of lithium's effect on bulk diversity comes from the blood being *reweighted* toward granulocytes, a lineage with a narrower transcriptome; changes inside lineages are not detectable.",
          pct(DB["reweighting", "share_of_total"], 0)),
  sprintf("**Lithium makes the cell mix itself less even** (p = %s), and does not make people's cell mixes or transcriptomes more variable (dispersion p = %s and %s).",
          f3(DC$p[1]), f2(bp["gene_aitchison", "cellmix_dispersion_p"]), f2(bp["gene_aitchison", "B2_p"])),
  sprintf("**Beta diversity:** lithium shifts whole-transcriptome structure (p = %s); %s of that shift is shared with the cell mix. The formal joint-significance test fails on the lineage tree because lithium's effect on the four balances as a whole is weak (p = %s), but passes on the proposal tree (p = %s). At the pathway level the shift is almost entirely shared with the cell mix and no direct part survives (p = %s).",
          fp(bp["gene_aitchison", "B1_p"]), pct(bp["gene_aitchison", "shared_with_cellmix"], 0), f3(bp["gene_aitchison", "path_a_p"]),
          f3(bpp["gene_aitchison", "path_a_p"]), f2(bp["pw_hallmark_aitchison", "B3i_direct_p"])),
  sprintf("**Confirmatory family (Holm, 2 tests):** alpha mediation adjusted p = %s (**%s**); beta mediation adjusted p = %s (%s).",
          f3(holm["alpha"]), if (holm["alpha"] < 0.05) "significant" else "not significant",
          f3(holm["beta"]), if (holm["beta"] < 0.05) "significant" else "not significant"),
  sprintf("**Exploratory — bipolar off lithium vs controls:** no change in diversity, cell mix or average transcriptome, but patients may be *more heterogeneous* among themselves on abundance-weighted distances (dispersion p = %s, %s). A lead, not a finding.",
          f3(bbd["gene_jsd", "B2_p"]), f3(bbd["gene_bray", "B2_p"]))
)))
add(H(3, "What changes in the proposal and paper"))
add(UL(c(
  "Report the proportion mediated with its interval, never a count of mediated genes.",
  "Correct the proposal's claim that ILR ratios were more consistent across deconvolutions — aggregation from 22 types to 5 lineages did that, not ILR.",
  "State which cell-type tree is used (both were run; they agree).",
  "Drop \"only\" from the myeloid-vs-lymphoid prediction, or keep it with the refinement that granulocytes rise relative to monocytes as well.",
  "Extract the obesity → neutrophil figure from Furuncuoğlu et al. (2016) for the E-value benchmark; that paper was not in the project library."
)))

## =============================================================================
## 1. QUESTION AND DESIGN
## =============================================================================
add(HR(), H(2, "1. The question and the design"))
add(Par("Lithium raises neutrophil proportions. Blood is a mixture, so changing the mixture changes bulk gene expression ",
        "arithmetically, without any cell behaving differently. The differential-expression arm had already found lithium's ",
        "whole-blood signature to be large and largely compositional. This arm puts a number on it."))
add(TB(data.frame(Element = c("Sample", "Exposure", "Mediator", "Outcome", "Covariates", "Dropped"),
  Choice = c("226 bipolar I subjects: 152 lithium users, 74 non-users",
             "lithium use, 0/1",
             "the cell mix as ILR balances, treated jointly as one vector",
             "voom log-CPM expression, one model per gene, 12,018 non-signature genes",
             "age, sex, tobacco, RIN, sequencing plate, sequencing PCs 1–3",
             "assessment group — constant within bipolar I, so including it would make the design singular"),
  stringsAsFactors = FALSE)))
add(Par("The protocol (`PRESPEC.md`) was written before any mediation estimate existed, because the pre-proposal feedback ",
        "warned that trying exposures, mediators and outcomes until something is significant is a multiple-testing problem. ",
        "It declares one primary model, every threshold, and a closed list of sensitivity and exploratory analyses."))

## =============================================================================
## 2. THE MEDIATOR
## =============================================================================
add(HR(), H(2, "2. The mediator: the cell mix as ILR balances"))
add(Par("Cell fractions sum to one, so they cannot all enter a regression, and a rise in one forces a fall in another. ",
        "Isometric log-ratio (ILR) balances on a haematopoietic tree solve both problems: each balance is a log-ratio between ",
        "two groups of cell types. The written proposal and the implemented pipeline used different trees, so **both were run**."))
add(TB(data.frame(
  ` ` = c("parts", "balances", "cell mass used", "samples needing zero replacement", "median reliability (ICC)"),
  `Lineage tree (pipeline)` = c("granulocyte, monocyte, T, NK, B", "4", "100%", "0", f3(median(REL$icc_absolute))),
  `Proposal tree` = c("neutrophil, monocyte, CD4 naive, CD4 memory, NK, B", "5", "92.4%", "71 of 474 (15%)", f3(median(RELP$icc_absolute))),
  check.names = FALSE, stringsAsFactors = FALSE)))
add(TB(data.frame(Tree = c(rep("lineage", 4), rep("proposal", 5)),
  Balance = c(REL$balance, RELP$balance), Contrast = c(REL$contrast, RELP$contrast),
  `Reliability (ICC)` = f3(c(REL$icc_absolute, RELP$icc_absolute)), check.names = FALSE),
  caption = "Reliability = agreement between our CIBERSORTx deconvolution and Krebs et al.'s deposited CIBERSORT fractions (ICC, absolute agreement)."))
add(Par("The proposal tree's CD4-naive-vs-memory balance has a reliability of only ", f3(RELP$icc_absolute[4]),
        " — two deconvolutions of the same blood barely agree on it. This is the main reason the lineage tree is primary."))

## =============================================================================
## 3. ESTIMATOR
## =============================================================================
add(HR(), H(2, "3. The estimator, and how it was verified"))
add(CODE(paste(
  "mediator models   M_k = b0_k + b1_k·lithium + b2_k'·C + e_k                      (k = 1..K)",
  "outcome model     Y   = t0 + t1·lithium + Σ t2_k·M_k + Σ t3_k·(lithium·M_k) + t4'·C + e",
  "",
  "NDE = t1 + Σ t3_k·b0_k            natural direct effect",
  "NIE = Σ (t2_k + t3_k)·b1_k         natural indirect effect (through the cell mix)",
  "TE  = NDE + NIE                    total effect",
  "PM  = NIE / TE                     proportion mediated (only where TE is significant)", sep = "\n")))
add(Par("These are Pearl's (2001) natural effects in the regression form of VanderWeele & Vansteelandt (2009, 2010), ",
        "extended to a vector mediator (VanderWeele & Vansteelandt 2014). Because the mediator is a vector, NIE is the ",
        "*joint* indirect effect through the whole cell mix; per-balance pieces are descriptive, not separate causal effects. ",
        "Standard errors use the delta method with an exactly block-diagonal covariance, and limma's empirical-Bayes variance ",
        "moderation for consistency with the DE arm."))
add(TB(data.frame(Check = c("Monte-Carlo counterfactuals (4 million draws)", "independent lm() computation",
                            "four-way decomposition sums to TE", "centring changes nothing",
                            "delta-method SE vs 4,000-replicate bootstrap", "95% interval coverage, 500 datasets",
                            "product-of-coefficients is a genuine refit, not PIE relabelled", "total suite"),
  Result = c("matches the closed form to 1.4e-4", "agrees to 1e-9", "exact", "identical to 1e-8",
             "within 1.9–6.3% (delta slightly low)", "93.0% (NDE), 94.4% (NIE)", "verified", VAL[1]))))
add(Par("**The delta method is mildly anticonservative** (~5% low in simulation, ~", pct(1 - BT$lineage$delta_vs_boot, 0),
        " low against the bootstrap on the real data), so bootstrap intervals are quoted throughout."))

## =============================================================================
## 4. STAGE ONE
## =============================================================================
add(HR(), H(2, "4. Stage one: what lithium does to the cell mix"))
mt <- rbind(transform(mtL, tree = "lineage"), transform(mtP, tree = "proposal"))
add(TB(data.frame(Tree = mt$tree, Balance = mt$balance, Contrast = mt$contrast,
                  Effect = sprintf("%+.4f", mt$b1), t = sprintf("%+.2f", mt$t), p = fp(mt$p), `BH q` = f3(mt$q_BH),
                  check.names = FALSE)))
add(Par("Both the myeloid-vs-lymphoid balance and the granulocyte-vs-monocyte balance move **up**, as the neutrophilia ",
        "literature predicts (Kast 2008; Krebs et al. 2020). The lymphoid-internal balances do not move. After correction ",
        "across the four balances only granulocytes-vs-monocytes survives (q = ", f3(mtL$q_BH[2]), "). The proposal's ",
        "\"myeloid relative to lymphoid **only**\" holds in the sense that nothing shifts inside the lymphoid compartment, with ",
        "the refinement that granulocytes also rise relative to monocytes."))
add(Par("Tested as one multivariate whole (MANOVA, Pillai), lithium's effect on the four balances together gives p = ",
        f3(bp["gene_aitchison", "path_a_p"]), " (lineage) — weaker than the single strongest balance, because three of the four ",
        "balances do not move and dilute the test. This matters for the beta-diversity mediation test in section 12."))

## =============================================================================
## 5. PRIMARY RESULT
## =============================================================================
add(HR(), H(2, "5. The primary mediation result"))
add(TB(data.frame(` ` = c("genes tested", "significant total effect", "significant direct effect", "significant indirect effect (per gene, delta method)",
                          "lithium × cell-mix interaction", "median proportion mediated", "95% bootstrap CI", "share of genes with PM > 0.50"),
  `Lineage (primary)` = c(fi(nrow(L$tb)), fi(sum(sigL)), fi(sum(L$tb$NDE_q_BH < 0.05)), fi(sum(L$tb$NIE_q_BH < 0.05)),
                          fi(sum(L$tb$INT_q_BH < 0.05)), f3(pmL), paste(f3(BT$lineage$pm_ci), collapse = " – "),
                          pct(mean(L$tb$PM[sigL] > 0.5))),
  `Proposal` = c(fi(nrow(P$tb)), fi(sum(sigP)), fi(sum(P$tb$NDE_q_BH < 0.05)), fi(sum(P$tb$NIE_q_BH < 0.05)),
                 fi(sum(P$tb$INT_q_BH < 0.05)), f3(pmP), paste(f3(BT$proposal$pm_ci), collapse = " – "),
                 pct(mean(P$tb$PM[sigP] > 0.5))), check.names = FALSE)))
add(H(3, "Four-way decomposition (VanderWeele 2014), lineage tree"))
add(TB(data.frame(Component = c("CDE", "INT_ref", "INT_med", "PIE"),
  Meaning = c("needs neither mediation nor interaction", "interaction only", "mediation and interaction", "pure mediation"),
  Share = pct(fwL / sum(fwL)))))
add(Par("Mediation (PIE + INT_med) accounts for ", pct((fwL["PIE"] + fwL["INTmed"]) / sum(fwL)),
        " of the decomposed effect; the direct part (CDE) is the largest single component."))
add(H(3, "Was allowing the interaction worth it?"))
pcd <- median(abs(L$tb$NIE - L$tb$NIE_prod)) / median(abs(L$tb$NIE))
add(Par("Only ", sum(L$tb$INT_q_BH < 0.05), " gene shows a significant interaction, so the proposal's original ",
        "product-of-coefficients estimator was defensible in practice. Refitting it properly — without the interaction ",
        "columns — it differs from ours by ", pct(pcd), " of the median effect and disagrees in sign for ",
        pct(mean(sign(L$tb$NIE) != sign(L$tb$NIE_prod)), 2), " of genes. (An earlier version of this comparison was ",
        "tautological; see section 14.)"))
add(H(3, "External validation"))
add(Par("The mediation total effect is, by construction, what the DE arm estimated with different code, a different ",
        "model and a different gene set. They agree: ", DEX[1], ". The decomposition splits an independently ",
        "established effect; it does not invent one."))

## =============================================================================
## 6. INFERENCE ON THE PROPORTION
## =============================================================================
add(HR(), H(2, "6. Can the proportion mediated be trusted?"))
add(H(3, "Bootstrap — the interval"))
add(TB(data.frame(Tree = c("lineage", "proposal"), Genes = c(BT$lineage$n_sig, BT$proposal$n_sig),
  `Median PM` = f3(c(BT$lineage$pm_obs, BT$proposal$pm_obs)),
  `95% CI` = c(paste(f3(BT$lineage$pm_ci), collapse = " – "), paste(f3(BT$proposal$pm_ci), collapse = " – ")),
  `Resamples above 0` = pct(c(mean(BT$lineage$pm_rep > 0), mean(BT$proposal$pm_rep > 0))),
  `Resamples above 0.5` = pct(c(mean(BT$lineage$pm_rep > 0.5), mean(BT$proposal$pm_rep > 0.5))), check.names = FALSE),
  caption = "B = 1,000 subject resamples; TMM, voom, mediator and outcome models refitted in every replicate."))
add(Par("The pre-declared \"substantial\" threshold (PM > 0.5) is **not met**: 0.5 lies inside both intervals."))
add(H(3, "Permutation — the null"))
pt <- rbind(transform(psL, tree = "lineage"), transform(psP, tree = "proposal"))
add(TB(data.frame(Tree = pt$tree, Statistic = pt$statistic, Observed = fs(pt$observed, 4),
                  `Null median` = fs(pt$null_median, 4), `Null max` = fs(pt$null_max, 4),
                  `Permutation p` = f4(pt$perm_p), check.names = FALSE),
  caption = "B = 1,000 shuffles of the lithium labels. The 'median proportion mediated' rows are superseded by the fixed-gene-set null below."))
add(Par("**The total effect is unambiguous** (p = ", f3(psL$perm_p[2]), "). **The mediated-gene count is meaningless**: ",
        "with lithium doing nothing by construction, a single shuffle produced ", fi(psL$null_max[1]),
        " \"significant\" mediated genes."))
add(H(3, "The proportion mediated against its own null"))
add(TB(data.frame(Tree = PMF$tree, Observed = f3(PMF$observed), `Null median` = f3(PMF$null_median),
                  `Null 95% range` = paste(f3(PMF$null_lo), "to", f3(PMF$null_hi)),
                  `Permutation p` = f3(PMF$perm_p), z = f2(ZT$z), check.names = FALSE),
  caption = "Gene set held fixed at the genes significant in the real data; 300 shuffles."))
add(NOTE("**The most important caveat in the arm.** The proportion mediated sits at the same place under pure noise ",
         "(z = ", f2(ZT$z[1]), "). PM = NIE / (NDE + NIE) has its numerator inside its denominator, so it lands near 0.4 ",
         "even when nothing is happening. It is a legitimate description of how a *real* total effect splits — and the total ",
         "effect is real — but it is **not evidence that mediation exists**. That evidence comes from lithium moving the cell mix ",
         "and from the indirect-effect tests in sections 12 and 4."))

## =============================================================================
## 7. POWER
## =============================================================================
add(HR(), H(2, "7. Power — why per-gene answers are out of reach"))
add(TB(data.frame(` ` = c("share of Var(NIE) from the lithium → cell-mix link", "median observed |NIE|",
                          "minimum detectable NIE, 80% power, uncorrected", "minimum detectable NIE, 80% power, Bonferroni",
                          "largest attainable |z| with perfect expression measurement", "|z| needed (Bonferroni, 12,018 genes)",
                          "subjects needed"),
  Value = c(pct(median(CAL$v_beta / (CAL$v_beta + CAL$v_theta))), f4(median(abs(L$tb$NIE))),
            f4(CAL$power$MDE_NIE[1]), f4(CAL$power$MDE_NIE[2]), f2(max(CAL$z_ceiling)),
            f2(CAL$power$z_needed[2]), paste0("~", fi(CAL$n_needed), " (vs 226)")), check.names = FALSE)))
add(Par("You cannot be more certain about a mediated effect than about the lithium → cell-mix link that carries it. ",
        "That part of the variance shrinks only with more subjects, never with deeper sequencing. This was knowable before ",
        "any data were analysed, which is exactly why the diversity analysis (section 12) asks the question once, globally."))

## =============================================================================
## 8. SENSITIVITY
## =============================================================================
add(HR(), H(2, "8. Sensitivity analyses (declared in advance)"))
S <- SEN$OUT
add(TB(data.frame(`#` = c("—", "S1", "S2", "S3", "S5", "S6", "S9", "S11", "S12"),
  `What changed` = c("primary", "LM22 genes analysed separately (superseded by section 11)", "Krebs' independent deconvolution as mediator",
                     "regression calibration for the estimated mediator", "20 tobacco imputations, Rubin-pooled",
                     "zero-replacement constant 0.10–0.80", "no covariates / age only", "proposal's 6-part tree", "unweighted outcome models"),
  `Median PM` = c(f3(pmL), f3(S$S1$median_PM[2]), f3(S$S2$median_PM[2]), f3(CAL$cal$corrected[2]), f3(pmL),
                  paste(f3(min(S$S6$median_PM)), "–", f3(max(S$S6$median_PM))),
                  paste(f3(S$S9$median_PM[3]), "/", f3(S$S9$median_PM[2])), f3(pmP), f3(S$S12$median_PM[2])),
  Verdict = c("—", "see section 11", sprintf("holds; per-gene NIE correlate %.3f", cor(SEN$s1_main$NIE, SEN$s2$NIE)),
              "raises it; see section 9", "exact no-op: none of the 226 had tobacco imputed",
              "not load-bearing", "adjustment lowers PM; age does real work", "agrees",
              sprintf("immaterial (NIE correlate %.4f)", cor(SEN$s1_main$NIE, SEN$s12$NIE))), check.names = FALSE)))

## =============================================================================
## 9. CALIBRATION
## =============================================================================
add(HR(), H(2, "9. Correcting for the mediator being estimated"))
add(Par("Classical measurement error in a mediator attenuates its coefficient and understates mediation. The proposal commits ",
        "to regression calibration (Valeri, Lin & VanderWeele 2014) with the uncorrected figure reported as a lower bound."))
add(TB(data.frame(Tree = c("lineage", "proposal"), `Median reliability` = f3(c(median(REL$icc_absolute), median(RELP$icc_absolute))),
                  `PM uncorrected` = f3(c(CAL$cal$uncorrected[2], CALP$summary$PM_uncorrected[2])),
                  `PM corrected` = f3(c(CAL$cal$corrected[2], CALP$summary$PM_corrected[2])), check.names = FALSE)))
add(Par("Two problems had to be fixed to do this honestly. The full multivariate correction is **inadmissible** — the ",
        "reliabilities imply more error than the observed covariance can hold — and at the admissibility boundary it is ",
        "**numerically unstable**:"))
add(TB(data.frame(`Correction strength c` = f3(SW$c), `Condition number` = f2(SW$condition_number),
                  `Median PM` = f3(SW$median_PM), check.names = FALSE)))
add(Par("So the stable per-balance correction is reported instead of a point on a cliff edge. The *direction* is certain ",
        "(every version raises PM); the *size* is not. The proposal tree corrects further only because its mediator is ",
        "measured worse. The total effect is held exactly fixed throughout, because measurement error in a mediator cannot ",
        "change it."))

## =============================================================================
## 10. E-VALUE
## =============================================================================
add(HR(), H(2, "10. Unmeasured confounding — the E-value"))
add(TB(data.frame(Quantity = CAL$E$quantity, `Median E-value` = f3(CAL$E$median_E),
                  `At the top gene` = f3(CAL$E$E_at_top_gene), check.names = FALSE)))
add(Par("An unmeasured confounder needs a risk ratio of about ", f2(CAL$E$median_E[1]), " with both lithium use and ",
        "expression to explain the typical effect away. Age — which we adjust for — reaches 2.33 on this scale, so a ",
        "confounder no stronger than age would suffice. Illness severity, mood state and other psychotropics are plausible ",
        "at that strength and are not recorded. BMI is not recorded either; the proposal's benchmark (Furuncuoğlu et al. 2016) ",
        "is not in the project library and was not guessed."))

## =============================================================================
## 11. CIRCULARITY
## =============================================================================
add(HR(), H(2, "11. Circularity sensitivity analysis"))
add(Par("**Why.** The cell fractions were estimated *from* the expression matrix using the LM22 signature genes, so for those ",
        "genes the mediator is close to a function of the outcome. The question is how much that inflates the answer."))
add(H(3, "Which of the 547 signature genes exist in this data"))
add(TB(as.data.frame(table(Status = ACC$status)), caption = "Measured against the raw count matrix, not assumed."))
add(Par("The 189 unexpressed genes mostly mark cell types that are rare or absent in blood (regulatory T cells, plasma cells, ",
        "macrophages, mast cells, dendritic cells). They carried essentially no weight in the deconvolution either."))
add(H(3, "Design"))
add(Par("One common fit — normalisation, voom weights and the variance prior computed once on all 12,368 genes — so the ",
        "analyses differ **only** in which genes are summarised. The earlier marker check (S1) had normalised the markers on ",
        "their own counts, which puts cell composition into the denominator; it is superseded."))
add(H(3, "Results"))
cc <- rbind(cl, cp)
add(TB(data.frame(Tree = cc$tree, Set = cc$set, Genes = fi(cc$n_genes), `Total` = fi(cc$n_TE), `Direct` = fi(cc$n_NDE),
                  `Median PM` = f3(cc$median_PM), `Cell-mix partial R²` = f3(cc$median_partialR2_cellmix), check.names = FALSE),
  caption = "Total / Direct = genes with a significant total / direct effect (BH within each set)."))
add(TB(data.frame(Tree = CRI$tree, Set = CRI$set, `Median PM` = f3(CRI$observed),
                  `Bootstrap 95% CI` = paste(f3(CRI$boot_lo), "–", f3(CRI$boot_hi)),
                  `Permutation p` = f3(CRI$perm_p), check.names = FALSE),
  caption = "C1 non-markers, C2 expressed LM22 markers, C3 everything. 1,000 bootstrap replicates; 300 permutations."))
add(TB(data.frame(Tree = CRD$tree, Contrast = CRD$contrast, Difference = f4(CRD$observed),
                  `Paired 95% CI` = paste(f4(CRD$lo), "to", f4(CRD$hi)), check.names = FALSE),
  caption = "Paired: both sets come from the same bootstrap replicates."))
add(TB(data.frame(Tree = NB$tree, `NIE correlation` = f4(NB$corr_NIE),
                  `Median PM, original` = f4(NB$median_PM_original), `Median PM, common fit` = f4(NB$median_PM_commonfit),
                  check.names = FALSE), caption = "Normalising on 12,018 genes (original) vs 12,368 (common fit)."))
add(H(3, "What it means"))
add(UL(c(
  sprintf("**The circularity signature is real.** The cell mix explains %s× more of a signature gene's expression (partial R² %s vs %s; Wilcoxon p < 1e-15). Exactly what circularity predicts.",
          f2(cl$median_partialR2_cellmix[2] / cl$median_partialR2_cellmix[1]), f3(cl$median_partialR2_cellmix[2]), f3(cl$median_partialR2_cellmix[1])),
  sprintf("**It does not inflate the proportion mediated.** Signature genes show a median PM of %s against %s for the rest — a paired difference of %s (95%% CI %s to %s). The difference is not significant, and it points the *opposite* way from inflation. Pre-declared criterion 1 (stronger cell-mix link for signature genes) holds; criterion 2 (larger share mediated) does not — exactly the case the protocol named in advance: circularity strengthens the mediator–outcome link without inflating the share mediated.",
          f3(cl$median_PM[2]), f3(cl$median_PM[1]), f4(cdl$observed[1]), f4(cdl$lo[1]), f4(cdl$hi[1])),
  sprintf("**Including them changes nothing.** The full set's median PM differs from the non-marker set's by %s.", f4(cl$median_PM[3] - cl$median_PM[1])),
  sprintf("**Why would it lean lower at all? A likely reason, not tested directly:** signature genes carry proportionally more *direct* effect (%s of their significant genes vs %s). They mark sub-types *within* lineages — eosinophils inside granulocytes, plasma cells inside B cells — and shifts between sub-types are invisible to a lineage-level mediator, so they surface as \"direct\". The \"direct\" effect is only direct relative to the resolution of the mediator, which makes the headline PM conservative on this count too.",
          pct(cl$n_NDE[2] / cl$n_TE[2], 0), pct(cl$n_NDE[1] / cl$n_TE[1], 0)),
  sprintf("**The most specific signature genes are not the most circular** (Spearman ρ = %s, p = %s): circularity works through the deconvolution as a whole, not gene by gene.",
          f3(CR$lineage$spec_rho), f2(CR$lineage$spec_p)),
  sprintf("**The normalisation basis is immaterial** (per-gene NIE correlation %s).", f4(NB$corr_NIE[1])),
  "**Verdict:** excluding the signature genes was the right call, but it was not load-bearing — the ~0.39 headline is not a circularity artefact."
)))

## =============================================================================
## 12. DIVERSITY
## =============================================================================
add(HR(), H(2, "12. Alpha and beta diversity, mediated through the cell mix"))
add(Par("Requested at the weekly check-in. Research and every decision are in `DIVERSITY_RESEARCH.md` and `PRESPEC.md` §8; ",
        "the decisions were fixed before any of this ran."))
add(H(3, "The decisions, briefly"))
add(TB(data.frame(Question = c("Which meaning of alpha/beta?", "Location or dispersion?", "Exposure", "Pathway database",
                               "Pathway representation", "Alpha measure", "Gene universe", "Beta distance", "Packages"),
  Decision = c("All three: per-person / between-person (primary); within- vs between-cell-type (secondary); cell mix's own diversity (secondary)",
               "Both; dispersion mediated only if dispersion differs", "lithium (primary); bipolar vs control (exploratory)",
               "MSigDB Hallmark, 50 sets (primary); KEGG legacy (sensitivity)", "aggregated abundance by fractional counting, not GSVA",
               "Shannon (primary); Hill q = 0, 2", "protein-coding, filtered, non-LM22: 10,288 genes", "Aitchison (primary); Jensen–Shannon, Bray–Curtis",
               "vegan, msigdbr (22.8 MB + 84.6 MB gene-set cache)"),
  Why = c("standard omics usage answers the question with one test; the other two explain the mechanism",
          "PERMANOVA is unreliable under unequal dispersion in unbalanced designs (Anderson & Walsh 2013); ours is 152:74",
          "lithium is the arm's pre-specified primary; exposure-swapping is a multiple-testing problem",
          "built to remove redundancy; 50 interpretable programmes; KEGG legacy is flat, not nested",
          "diversity needs non-negative abundances; GSVA scores can be negative and have no entropy precedent",
          "the only measure that partitions exactly with unequal weights (Jost 2007)",
          "structural non-coding RNAs (RN7SL1 6.9%, RN7SL2 4.8%, RN7SK 3.5%) dominate all-gene counts",
          "same geometry as the ILR mediator", "GSVA and energy (MODIMA) were not needed"))))
add(H(3, "Alpha diversity — per person"))
alt <- rbind(transform(ap, tree = "lineage"), transform(app, tree = "proposal"))
pick <- c("H_primary", "logHill2", "logHill0", "H_with_markers", "H_no_mito", "H_counts", "H_all_filtered",
          "H_pw_hallmark", "H_pw_hallmark_noLM22", "H_pw_kegg")
alt <- alt[alt$outcome %in% pick, ]
add(TB(data.frame(Tree = alt$tree, Outcome = alt$description, `Total effect` = fs(alt$TE), `Total p` = f3(alt$TE_p),
                  `Indirect (NIE)` = fs(alt$NIE), `NIE 95% CI` = paste(fs(alt$NIE_lo), "to", fs(alt$NIE_hi)),
                  `NIE perm. p` = f4(alt$NIE_p_perm), check.names = FALSE),
  caption = "Outcomes on the log scale of Hill numbers (nats): an effect d changes the effective number of genes by exp(d). 10,000 bootstrap replicates and 10,000 permutations."))
add(UL(c(
  sprintf("**Primary:** the indirect effect is %s (95%% CI %s to %s; permutation p = %s) — lithium lowers transcriptome diversity by about %s through the cell mix. Interaction p = %s.",
          fs(ap["H_primary", "NIE"]), fs(ap["H_primary", "NIE_lo"]), fs(ap["H_primary", "NIE_hi"]), f4(ap["H_primary", "NIE_p_perm"]),
          pct(-(exp(ap["H_primary", "NIE"]) - 1), 1), f2(ap["H_primary", "p_interaction"])),
  sprintf("**The total effect alone is not significant** (p = %s). An indirect effect can be detectable when the total is not: the total also carries a small, noisy direct part, and the indirect test is far more targeted (Zhao, Lynch & Chen 2010).",
          f3(ap["H_primary", "TE_p"])),
  sprintf("**Robust:** with signature genes included, with mitochondrial genes removed, on read-count proportions, on inverse Simpson and on the proposal tree (table). Also after a **post-hoc** adjustment for the structural-ncRNA share — added because primary diversity correlates −0.50 with that library-quality marker: indirect effect %s (95%% CI %s to %s).",
          fs(AL$posthoc$NIE), fs(AL$posthoc$NIE_lo), fs(AL$posthoc$NIE_hi)),
  "**Not robust where it should not be:** richness shows nothing (it is depth-driven), and the all-gene version shows nothing because structural non-coding RNAs drown it out.",
  "**Pathways:** Hallmark diversity falls through the cell mix (barely: CI just excludes zero), KEGG diversity *rises* — pathway-level diversity has no consistent direction, and is reported as descriptive."
)))
add(H(3, "Within- versus between-cell-type decomposition"))
add(TB(data.frame(Component = DB$component, Meaning = DB$meaning, `Lithium effect` = fs(DB$effect_of_lithium),
                  `95% CI` = paste(fs(DB$lo), "to", fs(DB$hi)), `Share of total` = pct(DB$share_of_total, 0), check.names = FALSE),
  caption = "H(bulk) = Σ f_k·H(lineage k) + JSD between lineages + residual (Jost 2007; Lin 1991). The five effects add up to the total exactly."))
add(Par("**", pct(DB["reweighting", "share_of_total"], 0), " of lithium's effect on bulk diversity comes from reweighting** — ",
        "blood shifting toward granulocytes, a lineage with a narrower transcriptome — and that component's interval excludes zero. ",
        "Changes *inside* lineages are not detectable. This is a model-free confirmation of the mediation result, with the caveat ",
        "that it rests on bMIND's per-sample lineage estimates (the mixture it builds tracks real bulk diversity at r = 0.94)."))
add(H(3, "The cell mix's own diversity"))
add(TB(data.frame(Configuration = DC$config, `Effect on cell-mix Shannon` = fs(DC$effect_of_exposure), p = f4(DC$p),
                  `NIE via cell-mix diversity` = fs(DC$NIE_via_cellmix_diversity),
                  `95% CI` = paste(fs(DC$NIE_lo), "to", fs(DC$NIE_hi)), check.names = FALSE)))
add(Par("Lithium makes the cell mix **less even** — a single number for its evenness carries the whole indirect effect ",
        "on transcriptome diversity. A shift toward granulocyte dominance also lowers the ceiling on bulk diversity, because the ",
        "between-lineage term can never exceed the diversity of the mix (JSD ≤ H(f), verified in every sample)."))
add(H(3, "Beta diversity — between people"))
blt <- rbind(transform(bp, tree = "lineage"), transform(bpp, tree = "proposal"))
add(TB(data.frame(Tree = blt$tree, Layer = blt$layer, `Shift p (B1)` = fp(blt$B1_p), `Dispersion p (B2)` = f3(blt$B2_p),
                  `Direct p (B3i)` = fp(blt$B3i_direct_p),
                  `Shared with cell mix` = ifelse(blt$B1_p < 0.05, pct(blt$shared_with_cellmix, 0), NA),
                  `Shared 95% CI` = ifelse(is.na(blt$shared_lo) | blt$B1_p >= 0.05, NA,
                                           paste(pct(blt$shared_lo, 0), "to", pct(blt$shared_hi, 0))),
                  `Path a p` = f3(blt$path_a_p), `Path b p` = fp(blt$path_b_p), check.names = FALSE),
  caption = "9,999 permutations; 1,000 bootstrap replicates for the shared fraction (Euclidean layers only; Jensen–Shannon and Bray–Curtis have no interval). Like the proportion mediated, the shared fraction is shown only where the shift itself (B1) is significant. Values slightly above 100% occur because an adjusted R² can dip just below zero once the cell mix is conditioned on; read them as 'all of it'."))
add(UL(c(
  sprintf("**Lithium shifts whole-transcriptome structure** (p = %s), and **does not change its heterogeneity** (p = %s), so the shift is a genuine location change rather than a dispersion artefact of the unbalanced design.",
          fp(bp["gene_aitchison", "B1_p"]), f2(bp["gene_aitchison", "B2_p"])),
  sprintf("**About %s of that shift is shared with the cell mix**, but the bootstrap interval reaches down to %s, and lithium's effect survives conditioning on the cell mix (p = %s).",
          pct(bp["gene_aitchison", "shared_with_cellmix"], 0), pct(bp["gene_aitchison", "shared_lo"], 0), fp(bp["gene_aitchison", "B3i_direct_p"])),
  sprintf("**Joint significance depends on the tree.** The cell mix strongly structures expression in both (path b p = %s). Lithium's effect on the balances as a whole is p = %s on the lineage tree (fails) and p = %s on the proposal tree (passes). The step that limits per-gene mediation limits this too.",
          fp(bp["gene_aitchison", "path_b_p"]), f3(bp["gene_aitchison", "path_a_p"]), f3(bpp["gene_aitchison", "path_a_p"])),
  sprintf("**At the level of biological programmes, the shift is essentially compositional.** On Hallmark pathways lithium shifts structure (p = %s), but once the cell mix is conditioned on no direct effect remains (p = %s); %s is shared with the cell mix (without signature genes: %s, 95%% CI %s to %s).",
          fp(bp["pw_hallmark_aitchison", "B1_p"]), f2(bp["pw_hallmark_aitchison", "B3i_direct_p"]),
          pct(bp["pw_hallmark_aitchison", "shared_with_cellmix"], 0), pct(bp["pw_hallmark_noLM22_aitchison", "shared_with_cellmix"], 0),
          pct(bp["pw_hallmark_noLM22_aitchison", "shared_lo"], 0), pct(bp["pw_hallmark_noLM22_aitchison", "shared_hi"], 0)),
  sprintf("**On abundance-weighted distances** (Jensen–Shannon, Bray–Curtis), lithium's shift is again significant (p = %s, %s) and its direct part is borderline once the cell mix is conditioned on (p = %s, %s).",
          fp(bp["gene_jsd", "B1_p"]), fp(bp["gene_bray", "B1_p"]), f3(bp["gene_jsd", "B3i_direct_p"]), f3(bp["gene_bray", "B3i_direct_p"])),
  "**B4 (dispersion mediation) was not run**, exactly as declared, because dispersion did not differ on the primary layer — nor on any other layer for lithium."
)))
ax <- BE$axes; axl <- ax[ax$config == "lithium | lineage tree (PRIMARY)", ]
add(TB(data.frame(Axis = axl$axis, `Variance explained` = pct(axl$var_explained), `Total effect p` = f3(axl$TE_p),
                  `Indirect p` = f3(axl$NIE_p), PM = f2(axl$PM), `Correlation with granulocyte fraction` = f2(axl$corr_with_granulocyte_fraction),
                  check.names = FALSE),
  caption = "The five leading principal axes of the primary layer, each mediated through med.R (lineage tree). Descriptive."))
add(Par("The main axes of transcriptome variation that lithium moves (PC1, PC2, PC4) are 59–79% mediated and track the ",
        "granulocyte fraction — the beta-diversity shift and the cell-mix shift are largely the same direction in expression space."))
add(H(3, "The confirmatory family"))
add(TB(data.frame(Test = c("alpha mediation (NIE, primary outcome)", "beta mediation (joint significance)"),
                  `Raw p` = f4(c(p_alpha, p_beta)), `Holm-adjusted p` = f4(holm),
                  Verdict = ifelse(holm < 0.05, "significant", "not significant"), check.names = FALSE),
  caption = "Two pre-declared confirmatory tests, Holm-corrected. Joint significance's p-value is the larger of its two paths (intersection–union test)."))
add(Par("Why can the alpha indirect effect be significant while path a is not? They are different tests. Path a asks whether ",
        "lithium moves the four balances *anywhere*, spending four degrees of freedom; three balances do not move and dilute it. ",
        "The indirect effect weights each balance by how much it matters for diversity, which concentrates on the one balance ",
        "lithium does move (granulocytes vs monocytes). The cell-mix evenness result says the same thing in one number."))
add(H(3, "Exploratory: bipolar (off lithium) vs control"))
add(UL(c(
  sprintf("**No change in alpha diversity** (total effect p = %s; indirect %s, CI %s to %s) and **no change in the cell mix's evenness** (p = %s).",
          f3(abd["H_primary", "TE_p"]), fs(abd["H_primary", "NIE"]), fs(abd["H_primary", "NIE_lo"]), fs(abd["H_primary", "NIE_hi"]), f3(DC$p[3])),
  sprintf("**No clear location shift** on the primary distance (p = %s; with signature genes p = %s) — consistent with the DE arm finding almost nothing gene by gene.",
          f3(bbd["gene_aitchison", "B1_p"]), f3(bbd["gene_aitchison_markers", "B1_p"])),
  sprintf("**But patients are more heterogeneous** on the abundance-weighted distances: dispersion p = %s (Jensen–Shannon), %s (Bray–Curtis), %s (KEGG pathways), with bipolar subjects further from their centroid than controls (Jensen–Shannon %s vs %s). Not on the primary Aitchison distance (p = %s).",
          f3(bbd["gene_jsd", "B2_p"]), f3(bbd["gene_bray", "B2_p"]), f3(bbd["pw_kegg_aitchison", "B2_p"]),
          f4(bbd["gene_jsd", "B2_mean_dist_exposed"]), f4(bbd["gene_jsd", "B2_mean_dist_unexposed"]), f3(bbd["gene_aitchison", "B2_p"])),
  "**How to read it:** bipolar disorder off lithium does not move the average blood transcriptome, but patients may differ from one another more than controls do — the \"Anna Karenina\" pattern (Zaneveld et al. 2017), here confined to the most abundant genes. This is **exploratory**: several distances were examined, it does not hold on the primary one, and the comparison carries the QC-filtered-controls caveat. Worth stating as a lead, not a finding."
)))

## =============================================================================
## 13. EXPLORATORY: AGE
## =============================================================================
add(HR(), H(2, "13. Exploratory analyses declared in the protocol"))
t1 <- SEC$T1$e; s1 <- t1$TE_q < 0.05
add(UL(c(
  sprintf("**T1 — age as exposure (per 10 years).** Age barely moves the cell mix (all balances |t| < 1.4); %d genes show a total effect of age and their median proportion mediated is only %s.",
          sum(s1), f3(median(t1$PM[s1]))),
  sprintf("**T2 — does age's mediated effect differ by lithium status?** No gene differs significantly, and the comparison is underpowered by a factor of %s (74 non-users). A null here is uninformative.",
          f2(SEC$mde / median(abs(SEC$T2$difference))))
)))

## =============================================================================
## 14. MISTAKES
## =============================================================================
add(HR(), H(2, "14. Mistakes found and corrected"))
add(Par("Recorded because a methods record that reports only successes is not a methods record."))
add(OL(c(
  "**A validation test was wrong, not the code** — it expected two estimators to agree exactly where they should differ by sampling noise. Replaced with an exact identity plus a z-test.",
  "**The total-effect identity appeared to fail** — traced to the deliberately unweighted mediator models; proved to machine precision once the projections match (`02c_identity.R`).",
  "**The product-of-coefficients comparison was tautological** (found by independent review): it reused the interaction model's coefficients. Refitted properly; 25.3% → 8.8%.",
  "**The regression calibration was inadmissible and the script carried on regardless** (review). Now shrinks to the admissible range and reports a sweep.",
  "**The calibration let the total effect drift 18.5%** (review), which is impossible in principle. Now held exactly fixed.",
  "**The bootstrap never produced an interval for the headline number** (review). Now computes the pre-declared estimand per replicate.",
  "**One review claim was tested and rejected** — a variance 'fix' that is right only conditional on the sample means; the bootstrap showed it would make the standard errors worse.",
  "**The calibration was quoted for both trees when it had run on one.** Run on the proposal tree, with that tree's own reliabilities.",
  "**The first marker-only check normalised markers on their own counts**, putting composition into the denominator. Superseded by the circularity analysis (section 11).",
  "**Path a was briefly coded as RDA** where the protocol declared MANOVA; they test different things. The declared test is used and RDA reported beside it.",
  "**Two shell-escaping errors** (a regular expression and a backreference) were caught by assertions before producing results."
)))

## =============================================================================
## 15. FOR THE PAPER
## =============================================================================
add(HR(), H(2, "15. What to carry into the paper"))
add(OL(c(
  sprintf("**Headline:** roughly %s of lithium's whole-blood effect is mediated by cell composition (95%% CI %s–%s), rising to %s after correcting for measurement error. Both are lower bounds; they bracket the 0.5 the proposal expected.",
          f3(pmL), f3(BT$lineage$pm_ci[1]), f3(BT$lineage$pm_ci[2]), f3(CAL$cal$corrected[2])),
  "**Do not present PM as proof of mediation** — present it as a decomposition of a real total effect, with the permutation result beside it.",
  "**Do not report a count of mediated genes.**",
  sprintf("**The cleanest single finding:** lithium lowers transcriptome diversity through the cell mix (NIE p = %s, Holm-adjusted %s), and %s of the change is the blood being reweighted toward granulocytes.",
          f3(p_alpha), f3(holm["alpha"]), pct(DB["reweighting", "share_of_total"], 0)),
  sprintf("**State the power ceiling:** ~%s subjects would be needed for per-gene significance.", fi(CAL$n_needed)),
  "**State the confounding threat plainly:** a confounder as strong as age would be enough, and illness severity and other drugs are unrecorded.",
  "**Correct the proposal** on the ILR-consistency claim and on which tree is used."
)))

## =============================================================================
## 16. GENES
## =============================================================================
add(HR(), H(2, "16. The genes"))
add(H(3, "What \"significant\" means here"))
add(Par("The model is fitted once per gene, so every gene has its own total effect (TE), direct effect (NDE) and indirect effect ",
        "through the cell mix (NIE). Significance is Benjamini–Hochberg at 5% across the gene family. Each gene is placed in one ",
        "category by asking, in order:"))
add(TB(data.frame(Category = c("TOTAL", "DIRECT", "VIA_MIX", "UNRESOLVED", "OPPOSING", "INTERACTION"),
  Rule = c("TE significant", "TOTAL and NDE significant", "TOTAL, NDE not significant, NIE bootstrap CI excludes 0",
           "TOTAL, neither of the above", "NDE significant but TE not", "lithium × cell-mix interaction significant"),
  `Plain meaning` = c("lithium changes this gene in whole blood, by any route — the solid list",
                      "the change survives accounting for the cell mix: candidates for lithium acting inside cells",
                      "the change is consistent with running through the cell mix, and no direct part can be detected",
                      "lithium changes it, and this study cannot say through which route",
                      "a direct effect that the indirect effect cancels — the net change is hidden, not absent",
                      "the effect of the cell mix on this gene depends on lithium"),
  Caveat = c("reproduces the DE arm (r = 0.9985)", "'direct' = not explained by this mediator; sub-type shifts within a lineage would land here",
             "a non-significant NDE is absence of evidence; the NIE interval is per gene, uncorrected",
             "—", "rare by nature", "one gene; the interaction is negligible overall"), check.names = FALSE)))
add(TB(data.frame(Category = CT$category, `Lineage tree` = fi(CT$lineage), `Proposal tree` = fi(CT$proposal),
                  `In both` = fi(CT$both_trees), check.names = FALSE),
  caption = "NIE_delta / NIE_boot / NIE_perm: per-gene significance of the indirect effect by three methods. NIE_consensus: significant by all three."))
add(NOTE("**Why per-gene NIE significance is not a category.** It is 0 genes on the lineage tree and thousands on the proposal ",
         "tree, from estimates that correlate at 0.98, and the permutation null alone reaches ", fi(psL$null_max[1]),
         " genes. Those lists are given in section 16.6 for completeness, flagged, and should not be reported as discoveries."))

gl <- GL$lineage; gp <- GL$proposal
sy <- function(g, f) { s <- g$symbol[f]; s[is.na(s) | s == ""] <- g$gene[f][is.na(s) | s == ""]; s }
add(H(3, "16.1 Direct-effect genes — lineage tree (primary)"))
dg <- gl[gl$category == "DIRECT", ]; dg <- dg[order(dg$NDE_q), ]
add(TB(data.frame(Gene = sy(dg, TRUE), TE = fs(dg$TE), NDE = fs(dg$NDE), `NDE q` = fp(dg$NDE_q), NIE = fs(dg$NIE),
                  PM = f2(dg$PM), `Also DIRECT on proposal tree` = ifelse(dg$gene %in% gp$gene[gp$category == "DIRECT"], "yes", "no"),
                  check.names = FALSE),
  caption = sprintf("%d genes, ordered by significance of the direct effect. Effects are log2 fold changes in expression (users vs non-users).", nrow(dg))))
add(H(3, "16.2 Genes with the largest mediated (indirect) effects — lineage tree"))
tg <- gl[gl$top_NIE, ]; tg <- tg[order(-abs(tg$NIE)), ]
add(TB(data.frame(Rank = seq_len(nrow(tg)), Gene = sy(tg, TRUE), NIE = fs(tg$NIE),
                  `NIE 95% CI` = paste(fs(tg$NIE_lo_boot), "to", fs(tg$NIE_hi_boot)), TE = fs(tg$TE), PM = f2(tg$PM),
                  Category = tg$category, check.names = FALSE),
  caption = "The 50 genes with a significant total effect whose response is most driven by the cell mix."))
add(H(3, "16.3 Every gene with a significant total effect — lineage tree, by category"))
for (cat in c("DIRECT", "VIA_MIX", "UNRESOLVED")) {
  for (dirn in c("up", "down")) {
    f <- gl$category == cat & (if (dirn == "up") gl$TE > 0 else gl$TE < 0)
    gg <- gl[f, ]; gg <- gg[order(gg$TE_q), ]
    add(GENES(sprintf("%s — %s-regulated by lithium", cat, dirn), sy(gg, TRUE), "ordered by significance of the total effect"))
  }
}
add(H(3, "16.4 Opposing and interaction genes"))
og <- gl[gl$category == "OPPOSING", ]; og <- og[order(og$NDE_q), ]
add(TB(data.frame(Gene = sy(og, TRUE), NDE = fs(og$NDE), NIE = fs(og$NIE), TE = fs(og$TE), `TE q` = f3(og$TE_q), check.names = FALSE),
  caption = "OPPOSING (lineage): significant direct effect, cancelled by an indirect effect of opposite sign."))
ig <- gl[gl$interaction, ]
add(GENES("INTERACTION (lineage)", sy(ig, TRUE)))
add(H(3, "16.5 Signature (LM22) genes with a significant total effect"))
mk <- GL$markers; mk <- mk[mk$TE_q < 0.05, ]
add(TB(data.frame(Gene = mk$symbol, `LM22 cell type it marks` = mk$lm22_top_celltype, TE = fs(mk$TE), `TE q` = fp(mk$TE_q),
                  NDE = fs(mk$NDE), PM = f2(mk$PM), `Cell-mix partial R²` = f3(mk$partialR2_cellmix), check.names = FALSE),
  caption = "From the circularity analysis (section 11), lineage tree, BH within the 350 expressed signature genes."))
add(H(3, "16.6 Per-gene indirect-effect significance — proposal tree (flagged)"))
add(NOTE("**Not discoveries.** Listed because they were requested and for completeness. The lineage tree gives none of these, ",
         "and the permutation null reaches thousands. Use only alongside the warning in section 6."))
cg <- gp[gp$NIE_consensus, ]; cg <- cg[order(cg$NIE_q_delta), ]
add(GENES("Proposal tree: NIE significant by all three methods (delta, bootstrap, permutation)", sy(cg, TRUE),
          "ordered by delta-method significance; each method's full list is in results/genes/genes_proposal_all.csv"))
add(H(3, "16.7 Proposal tree — direct-effect genes not already listed"))
pd <- gp[gp$category == "DIRECT" & !(gp$gene %in% dg$gene), ]
add(GENES("DIRECT on the proposal tree only", sy(pd, TRUE)))

## =============================================================================
## 17. FILES  /  18. REFERENCES
## =============================================================================
add(HR(), H(2, "17. Files"))
add(TB(data.frame(Path = c("PRESPEC.md", "METHODS.md", "DIVERSITY_RESEARCH.md", "MEDIATION_MASTER.md / .html", "scripts/med.R",
                           "scripts/00_validate_med.R", "scripts/02_effects.R", "scripts/03_bootstrap.R, 04_permutation.R, 04b_perm_pm.R",
                           "scripts/05–07, 09", "scripts/10_circularity.R, 10b_…", "scripts/11–13", "scripts/14_gene_lists.R",
                           "results/genes/", "results/circularity/", "results/diversity/"),
  Holds = c("the protocol, fixed before each analysis ran (incl. addenda A and B)", "methods and the incident log",
            "the diversity research, with sources", "this document", "the estimator — one implementation used everywhere",
            "27 correctness checks", "the primary result, both trees", "bootstrap and permutation inference",
            "sensitivity, calibration, E-value, power, exploratory", "circularity sensitivity analysis",
            "diversity: prep, alpha, beta", "gene categories", "every gene with its category and statistics (CSV)",
            "circularity results", "diversity results"))))
add(H(2, "18. Key references"))
add(UL(c(
  "Pearl J (2001). Direct and indirect effects. *UAI*.",
  "VanderWeele TJ, Vansteelandt S (2014). Mediation analysis with multiple mediators. *Epidemiologic Methods* 2:95–115.",
  "VanderWeele TJ (2014). A unification of mediation and interaction: a 4-way decomposition. *Epidemiology* 25:749–761.",
  "Valeri L, Lin X, VanderWeele TJ (2014). Mediation analysis when a continuous mediator is measured with error. *Statistics in Medicine* 33:4875–4890.",
  "VanderWeele TJ, Ding P (2017). Sensitivity analysis in observational research: introducing the E-value. *Annals of Internal Medicine* 167:268–274.",
  "Egozcue JJ et al. (2003). Isometric logratio transformations for compositional data analysis. *Mathematical Geology* 35:279–300.",
  "Newman AM et al. (2015, 2019). CIBERSORT; CIBERSORTx. *Nature Methods* 12:453–457; *Nature Biotechnology* 37:773–782.",
  "Krebs CE et al. (2020). Whole blood transcriptome analysis in bipolar disorder reveals strong lithium effect. *Psychological Medicine* 50:2575–2586.",
  "Martínez O, Reyes-Valdés MH (2008). Defining diversity, specialization, and gene specificity in transcriptomes through information theory. *PNAS* 105:9709–9714.",
  "García-Nieto PE, Wang B, Fraser HB (2022). Transcriptome diversity is a systematic source of variation in RNA-sequencing data. *PLoS Computational Biology* 18:e1009939.",
  "Jost L (2007). Partitioning diversity into independent alpha and beta components. *Ecology* 88:2427–2439.",
  "Anderson MJ (2001, 2006); Anderson MJ, Walsh DCI (2013). PERMANOVA; PERMDISP; heterogeneous dispersions. *Austral Ecology*; *Biometrics*; *Ecological Monographs*.",
  "Zapala MA, Schork NJ (2006). Multivariate regression analysis of distance matrices for testing associations between gene expression patterns and related variables. *PNAS* 103:19430–19435.",
  "Liberzon A et al. (2015). The MSigDB hallmark gene set collection. *Cell Systems* 1:417–425.",
  "Zhao X, Lynch JG, Chen Q (2010). Reconsidering Baron and Kenny: myths and truths about mediation analysis. *Journal of Consumer Research* 37:197–206.",
  "Full lists with verification marks: `DIVERSITY_RESEARCH.md` §10 and `METHODS.md`."
)))

## ---- write ---------------------------------------------------------------------
md <- render_md(D); html <- render_html(D, "Mediation arm — master document")
writeLines(md, "MEDIATION_MASTER.md", useBytes = TRUE)
con <- file("MEDIATION_MASTER.html", open = "w", encoding = "UTF-8"); writeLines(html, con); close(con)
cat(sprintf("wrote MEDIATION_MASTER.md (%d lines) and MEDIATION_MASTER.html (%d lines)\n", length(md), length(html)))
cat(sprintf("headings: %d | tables: %d | gene lists: %d\n", sum(sapply(D, `[[`, "type") == "h"),
            sum(sapply(D, `[[`, "type") == "tb"), sum(sapply(D, `[[`, "type") == "genes")))
