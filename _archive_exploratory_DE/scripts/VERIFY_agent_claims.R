#!/usr/bin/env Rscript
# ==============================================================================
# VERIFY_agent_claims.R -- independent re-computation of every agent-sourced
# number that FINDINGS.md leans on.
#
# The 14 parallel agents died at their reporting step, so no adversarial pass
# ran. This is that pass. Each claim is recomputed by a DIFFERENT route than the
# agent used -- own code, own parameters, own random seeds where relevant -- and
# the verdict is recorded against a stated tolerance rather than by eye.
#
# Claims under test:
#   1. case-control with lithium as a covariate collapses to ~6 DEGs
#   2. permutation null for the lithium contrast is ~0 DEGs
#   3. skipping normalisation inflates the count to ~4,013 with Jaccard ~0.21
#   4. TMM normalisation factors are explained by cell composition, r2 ~0.24
#   5. the spike-in calibration: composition adjustment preserves transcriptional
#      effects and removes compositional ones
#   6. split-half DEG reproducibility is ~0.32
# ==============================================================================

source("scripts/common.R")
p <- load_prep()
set.seed(20260929)               # deliberately NOT the project seed
OUT <- file.path(EXP_HOME, "runs", "VERIFY")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
V <- list()
vd <- function(claim, agent_val, mine, tol, note = "") {
  ok <- if (is.na(agent_val)) NA else abs(mine - agent_val) <= tol
  V[[length(V) + 1]] <<- data.frame(claim = claim, agent = agent_val, verified = mine,
                                    tolerance = tol, holds = ok, note = note,
                                    stringsAsFactors = FALSE)
  cat(sprintf("  [%s] %-46s agent %-9s mine %-9s\n",
              if (isTRUE(ok)) "HOLDS" else if (is.na(ok)) " -- " else "FAILS",
              claim, format(agent_val), format(round(mine, 4))))
  if (nzchar(note)) cat(sprintf("         %s\n", note))
}

## ---- 1. case-control adjusted for lithium ---------------------------------
cat("\n=== 1. case-control, lithium as a covariate ===\n")
# Different route: build the design by hand rather than through the agent's
# helper, and run it on both cohorts.
m <- p$meta; m$tob <- p$TOB[m$title, "tobacco_imp_01"]
cc_fit <- function(keep, with_li) {
  d <- m[keep, ]
  d$case <- factor(ifelse(d$dx == "Control", "ctrl", "case"), c("ctrl", "case"))
  cv <- c("age", "sex", "tob", "group", "rin", "plate", "seqpc1", "seqpc2", "seqpc3")
  if (with_li) cv <- c(cv, "lithium")
  d$plate <- droplevels(d$plate); d$sex <- droplevels(d$sex); d$group <- droplevels(d$group)
  d <- d[complete.cases(d[, c("case", cv)]), ]
  X <- model.matrix(as.formula(paste("~ case +", paste(cv, collapse = " + "))), data = d)
  ne <- nonEstimable(X); if (!is.null(ne)) X <- X[, setdiff(colnames(X), ne), drop = FALSE]
  cn <- p$counts[, d$title, drop = FALSE]; cn <- cn[filter_genes(cn, 10, 0.90), , drop = FALSE]
  tt <- topTable(eBayes(lmFit(voom(calcNormFactors(DGEList(cn), "TMM"), X), X)),
                 coef = "casecase", number = Inf, sort.by = "none")
  tt$gene <- rownames(tt); list(tt = tt, n = nrow(d))
}
a1 <- cc_fit(rep(TRUE, nrow(m)), FALSE)
a2 <- cc_fit(rep(TRUE, nrow(m)), TRUE)
a3 <- cc_fit(m$qc_pass, TRUE)
cat(sprintf("  474 without lithium: %d DEG | 474 with lithium: %d | 444 with lithium: %d\n",
            n_deg(a1$tt), n_deg(a2$tt), n_deg(a3$tt)))
vd("casecon no-lithium DEG", 1816, n_deg(a1$tt), 60)
vd("casecon +lithium, n=474", 1, n_deg(a2$tt), 5)
vd("casecon +lithium, n=444 (published spec)", 6, n_deg(a3$tt), 5)

## ---- 2. permutation null ---------------------------------------------------
cat("\n=== 2. permutation null, lithium contrast (20 permutations, new seed) ===\n")
B <- 20
nullv <- vapply(seq_len(B), function(i) {
  perm <- sample(226)
  r <- de_fit(p, "lithium", permute_exposure = perm)
  c(deg = n_deg(r$tt), pi0 = pi0_storey(r$tt$P.Value))
}, numeric(2))
cat(sprintf("  null DEG: mean %.2f  max %d  nonzero in %d/%d | null pi0 mean %.3f\n",
            mean(nullv["deg", ]), max(nullv["deg", ]), sum(nullv["deg", ] > 0), B,
            mean(nullv["pi0", ])))
vd("permutation null mean DEG", 0, mean(nullv["deg", ]), 2)
vd("permutation null mean pi0", 0.971, mean(nullv["pi0", ]), 0.05)

## ---- 3 & 4. normalisation --------------------------------------------------
cat("\n=== 3. normalisation ===\n")
r_tmm  <- de_fit(p, "lithium", norm = "TMM")
r_none <- de_fit(p, "lithium", norm = "none")
u <- intersect(r_tmm$tt$gene, r_none$tt$gene)
da <- intersect(deg_ids(r_tmm$tt), u); db <- intersect(deg_ids(r_none$tt), u)
jac <- length(intersect(da, db)) / length(union(da, db))
vd("DEG with TMM", 1382, n_deg(r_tmm$tt), 30)
vd("DEG with no normalisation", 4013, n_deg(r_none$tt), 120)
vd("Jaccard none-vs-TMM", 0.212, jac, 0.03)

cat("\n=== 4. are TMM factors explained by cell composition? ===\n")
mb <- p$meta[p$meta$dx == "BP1", ]
cntb <- p$counts[, mb$title, drop = FALSE]
cntb <- cntb[filter_genes(cntb, 10, 0.90), , drop = FALSE]
nf <- calcNormFactors(DGEList(cntb), method = "TMM")$samples$norm.factors
ILRb <- p$ILR[mb$title, , drop = FALSE]
r2_ilr  <- summary(lm(nf ~ ILRb))$r.squared
r2_gran <- summary(lm(nf ~ p$lineage[mb$title, "gran"]))$r.squared
r_gran  <- cor(nf, p$lineage[mb$title, "gran"])
cat(sprintf("  r2(normfactors ~ 4 ILR balances) = %.3f | r2(~granulocyte%%) = %.3f | r = %+.3f\n",
            r2_ilr, r2_gran, r_gran))
vd("r2 normfactors ~ ILR balances", 0.237, r2_ilr, 0.05)
vd("r normfactors vs granulocyte %", -0.476, r_gran, 0.06)

## ---- 5. spike-in calibration, own implementation --------------------------
cat("\n=== 5. spike-in: does the adjustment tell composition from transcription? ===\n")
# Independent construction. Take the real 226-sample matrix, add synthetic
# signal to 1,200 randomly chosen genes, and ask what survives adjustment.
#   TRANS  : add c * lithium_i          -- a true exposure effect, b1-independent
#   COMPOS : add c * (b1_i - mean(b1))  -- acts ONLY through the mixture, so the
#            gene's lithium association arises entirely via b1
# c is chosen per arm so both arms produce a comparable MARGINAL lithium effect
# (the COMPOS coefficient is scaled by 1/delta, delta = 0.156).
mb$tob <- p$TOB[mb$title, "tobacco_imp_01"]
covars <- c("age", "sex", "tob", "rin", "plate", "seqpc1", "seqpc2", "seqpc3")
mb <- mb[complete.cases(mb[, c("lithium", covars)]), ]
mb$plate <- droplevels(mb$plate); mb$sex <- droplevels(mb$sex)
b1 <- p$ILR[mb$title, "b1_myeloid_vs_lymphoid"]; b1c <- b1 - mean(b1)
delta <- coef(lm(as.formula(paste("b1 ~ lithium +", paste(covars, collapse = " + "))),
                 data = cbind(mb, b1 = b1)))["lithium"]
cnt2 <- p$counts[, mb$title, drop = FALSE]
cnt2 <- cnt2[filter_genes(cnt2, 10, 0.90), , drop = FALSE]
E <- cpm(calcNormFactors(DGEList(cnt2), "TMM"), log = TRUE, prior.count = 3)
Xm <- model.matrix(as.formula(paste("~ lithium +", paste(covars, collapse = " + "))), data = mb)
Xa <- cbind(Xm, b1 = b1)
spk <- sample(rownames(E), 1200)

spike_run <- function(arm, target = 0.20) {
  Y <- E
  add <- if (arm == "TRANS") target * mb$lithium else (target / delta) * b1c
  Y[spk, ] <- Y[spk, ] + matrix(rep(add, each = length(spk)), nrow = length(spk))
  t1 <- topTable(eBayes(lmFit(Y, Xm), trend = TRUE), coef = "lithium", number = Inf, sort.by = "none")
  t2 <- topTable(eBayes(lmFit(Y, Xa), trend = TRUE), coef = "lithium", number = Inf, sort.by = "none")
  det <- intersect(spk, rownames(t1)[t1$adj.P.Val < 0.05])
  sur <- intersect(det, rownames(t2)[t2$adj.P.Val < 0.05])
  shr <- mean(abs(t2$logFC[match(spk, rownames(t2))])) /
         mean(abs(t1$logFC[match(spk, rownames(t1))]))
  cat(sprintf("  %-7s detected %4d/1200 (%.0f%%) | survive adjustment %4d (%.0f%% of detected) | |logFC| ratio %.3f\n",
              arm, length(det), 100 * length(det) / 1200, length(sur),
              100 * length(sur) / max(length(det), 1), shr))
  c(det = length(det), sur = length(sur),
    rate = length(sur) / max(length(det), 1), shrink = shr)
}
tr <- spike_run("TRANS"); co <- spike_run("COMPOS")
vd("TRANS survival rate", 0.907, tr["rate"], 0.10,
   "transcriptional spike should SURVIVE composition adjustment")
vd("TRANS logFC shrinkage", 1.013, tr["shrink"], 0.10, "should be ~1, i.e. none")
vd("COMPOS survival rate", 0.000, co["rate"], 0.10,
   "compositional spike should be REMOVED by adjustment")
vd("COMPOS logFC shrinkage", 0.174, co["shrink"], 0.15, "should be far below 1")

## ---- 6. split-half reproducibility ----------------------------------------
cat("\n=== 6. split-half DEG reproducibility (10 splits, new seed) ===\n")
idx <- which(p$meta$dx == "BP1")
rep_v <- vapply(1:10, function(i) {
  s <- sample(idx); h1 <- rep(FALSE, nrow(p$meta)); h2 <- h1
  h1[s[1:113]] <- TRUE; h2[s[114:226]] <- TRUE
  t1 <- de_fit(p, "lithium", samples = h1)$tt; t2 <- de_fit(p, "lithium", samples = h2)$tt
  uu <- intersect(t1$gene, t2$gene)
  d1 <- intersect(deg_ids(t1), uu); d2 <- intersect(deg_ids(t2), uu)
  if (!length(union(d1, d2))) return(c(repro = NA, sign = NA))
  sh <- intersect(d1, d2)
  c(repro = length(sh) / min(length(d1), length(d2)),
    sign = if (length(sh)) mean(sign(t1$logFC[match(sh, t1$gene)]) ==
                                sign(t2$logFC[match(sh, t2$gene)])) else NA)
}, numeric(2))
cat(sprintf("  reproducibility mean %.3f (SD %.3f) | sign agreement %.3f\n",
            mean(rep_v["repro", ], na.rm = TRUE), sd(rep_v["repro", ], na.rm = TRUE),
            mean(rep_v["sign", ], na.rm = TRUE)))
vd("split-half reproducibility", 0.320, mean(rep_v["repro", ], na.rm = TRUE), 0.12)
vd("split-half sign agreement", 0.967, mean(rep_v["sign", ], na.rm = TRUE), 0.06)

## ---- verdict ---------------------------------------------------------------
R <- do.call(rbind, V)
cat("\n=============================== VERDICT ===============================\n")
print(R[, c("claim", "agent", "verified", "holds")], row.names = FALSE, digits = 4)
cat(sprintf("\n  %d of %d claims hold within tolerance\n", sum(R$holds, na.rm = TRUE), nrow(R)))
if (any(!R$holds, na.rm = TRUE))
  cat(sprintf("  FAILED: %s\n", paste(R$claim[!R$holds], collapse = "; ")))
write.csv(R, file.path(OUT, "verdict.csv"), row.names = FALSE)
cat("\nwrote runs/VERIFY/verdict.csv\n")
