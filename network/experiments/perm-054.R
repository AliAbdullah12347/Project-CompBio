#!/usr/bin/env Rscript
# ===========================================================================
# perm-054  A legitimate null for the BPD module enrichment (NEXT_STEPS item 2)
#
# THE CLAIM UNDER TEST
#
# analyses/019 found 27 of 69 modules rank-enriched for BPD-associated genes,
# against 21 for lithium -- even though whole blood yields only 4 BPD DEGs at
# BH<0.05 against 1,426 for lithium. That is the surprising part: a contrast
# with almost no significant genes lights up as many modules as one with
# fourteen hundred.
#
# WHY THE EXISTING NULL DOES NOT ESTABLISH IT
#
# 019's null draws random size-matched gene sets. Those are UNCORRELATED,
# while a module's genes are co-expressed by construction, so their DE
# statistics move together. A Wilcoxon that assumes independence is
# anti-conservative on correlated genes, and a null built from uncorrelated
# sets cannot see that. It rules out module size and nothing else.
#
# THE HONEST NULL
#
# Permute the diagnosis label across the 308 WB_BPD samples, recompute the
# whole DE fit, and re-score every module. The expression matrix is never
# touched, so the gene-gene correlation structure that makes the Wilcoxon
# anti-conservative is preserved exactly; only the group signal is destroyed.
# The null distribution of "how many modules reach BH<0.05" then absorbs the
# correlation inflation, and the observed 27 is read against it.
#
#   p = (1 + #{null >= observed}) / (1 + K)
#
# Add-one (Phipson & Smyth 2010), not mean(). Three scripts in this arm --
# perm-031, perm-032 and hub-037 -- published counts off mean(null >= obs),
# which returns exactly 0 when the observed beats every draw. A p of 0 is not
# a probability.
#
# POWER. This is ONE statistic, not 69, so no multiple-testing correction
# applies and the floor 1/(K+1) is the only constraint. At K=500 that is
# 0.002, far below any threshold we would use.
#
# COVARIATES. The label moves; the covariates stay with their samples,
# because they are properties of the subject. Permuting therefore also
# destroys the real association between diagnosis and tobacco (19% of
# controls, 32.7% of BP1), which makes each permuted design slightly
# better conditioned than the real one. That makes this null mildly
# OPTIMISTIC -- it is easier to beat than the truth -- so a negative here is
# conservative and a positive is the weaker reading. Said plainly because it
# cannot be fixed without a restricted permutation scheme that the NEXT_STEPS
# specification does not ask for.
#
# SLICING. PERM_FROM/PERM_TO run a contiguous slice of draws so the work
# spreads over pods; PERM_AGGREGATE=1 skips the draws and only aggregates.
# Every draw checkpoints, so a pod lost mid-slice still contributes.
# ===========================================================================

# One script, two rows. The lithium contrast needs exactly the same null --
# 019 scored both against the same invalid random-gene-set null, so the
# BPD-vs-lithium comparison is not symmetric until both have been redone the
# same way, on the same software. Separate run directories keep their
# checkpoints apart; everything between them is identical by construction.
CONTRASTS <- list("perm-054" = c(sel = "BPD", published = "WB_BPD"),
                  "perm-059" = c(sel = "LI",  published = "WB_LI"))
id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "perm-054"
stopifnot(id %in% names(CONTRASTS))

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

run_dir <- file.path("network/runs", id)
cfg_dir <- "network/config"
for (d_ in c(run_dir, cfg_dir)) dir.create(d_, showWarnings = FALSE, recursive = TRUE)

K_TOTAL   <- 500L
# One pod per slice. SLICE_ID/N_SLICES are what network/cluster/mkjobs.py
# exports to each pod, so the same two variables drive every sliced row in
# this arm rather than each inventing its own pair.
SLICE_ID  <- as.integer(Sys.getenv("SLICE_ID",  "0"))
N_SLICES  <- as.integer(Sys.getenv("N_SLICES",  "1"))
if (SLICE_ID >= 1L && N_SLICES >= 1L) {
  per    <- ceiling(K_TOTAL / N_SLICES)
  K_FROM <- (SLICE_ID - 1L) * per + 1L
  K_TO   <- min(SLICE_ID * per, K_TOTAL)
} else {
  K_FROM <- as.integer(Sys.getenv("PERM_FROM", "1"))
  K_TO   <- as.integer(Sys.getenv("PERM_TO",   as.character(K_TOTAL)))
}
AGG_ONLY  <- Sys.getenv("PERM_AGGREGATE", "0") == "1"
SEED_BASE <- 20261009L + 54000L
# de_final keys its selections "LI"/"BPD" and names its output files
# "WB_LI"/"WB_BPD" -- run_de.R's JOBS table maps between them. Both are needed
# here: the first to pull the sample set out of prep.rds, the second to find
# the published table this fit is checked against.
CONTRAST  <- unname(CONTRASTS[[id]]["sel"])        # key into p$sel
PUBLISHED <- unname(CONTRASTS[[id]]["published"])  # de_final/results/DE_<...>.csv.gz
MIN_MOD_N <- 20L          # 019's threshold; it is what makes the count 69 and not 71

suppressPackageStartupMessages({ library(edgeR); library(limma) })
options(stringsAsFactors = FALSE)

## ---- inputs -----------------------------------------------------------------
p <- readRDS("de_analysis/data/prep.rds")

g <- read.csv("network/runs/grp-051/genes_control_A.csv")
ref <- list(genes = g$gene, colors = g$module)
modules_all <- setdiff(unique(ref$colors), "grey")
msg("partition: ", length(ref$genes), " genes, ", length(modules_all), " non-grey modules")

## ---- the DE fit, copied from de_final/run_de.R ------------------------------
# Copied rather than sourced: run_de.R executes all four comparisons top to
# bottom and setwd()s to an absolute path on the iMac. A copy can drift from
# the original, so the unpermuted fit below is checked against the published
# DE_WB_BPD.csv.gz before any permutation runs. If the copy has diverged, that
# assertion fails and this row stops rather than producing a null for a model
# nobody reported.
design_for <- function(contrast, grp_override = NULL) {
  s  <- p$sel[[contrast]]
  mm <- p$meta[s$keep, , drop = FALSE]
  mm$grp <- if (is.null(grp_override)) s$grp else grp_override
  mm$plate <- droplevels(factor(mm$plate)); mm$sex <- droplevels(factor(mm$sex))
  if ("group" %in% s$cov) mm$group <- droplevels(factor(mm$group))
  cv <- s$cov[vapply(s$cov, function(v) length(unique(mm[[v]])) > 1, logical(1))]
  X  <- model.matrix(as.formula(paste("~ grp +", paste(cv, collapse = " + "))), data = mm)
  ne <- nonEstimable(X); if (!is.null(ne)) X <- X[, setdiff(colnames(X), ne), drop = FALSE]
  stopifnot(qr(X)$rank == ncol(X), "grpcase" %in% colnames(X))
  list(X = X, samples = mm$title)
}

fit_wb_p <- function(contrast, grp_override = NULL) {
  d   <- design_for(contrast, grp_override)
  dge <- calcNormFactors(DGEList(p$counts[, d$samples, drop = FALSE]), method = "TMM")
  v   <- voom(dge, d$X)
  tt  <- topTable(eBayes(lmFit(v, d$X)), coef = "grpcase", number = Inf, sort.by = "none")
  setNames(tt$P.Value, rownames(tt))
}

## ---- 019's scoring, applied to any p-value vector ---------------------------
# Identical to network/analyses/019-bpd-module-connection/run.R: Wilcoxon
# rank-sum on -log10(p), module genes against every other matched gene,
# one-sided greater, BH across modules, modules with fewer than 20 matched
# genes skipped.
count_enriched <- function(pvals) {
  pv    <- pvals[names(pvals) %in% ref$genes]
  score <- setNames(-log10(pmax(pv, 1e-300)), names(pv))
  res <- lapply(modules_all, function(m) {
    gin <- intersect(ref$genes[ref$colors == m], names(score))
    if (length(gin) < MIN_MOD_N) return(NULL)
    w <- suppressWarnings(wilcox.test(score[gin], score[setdiff(names(score), gin)],
                                      alternative = "greater"))
    data.frame(module = m, n = length(gin), p = w$p.value)
  })
  R <- do.call(rbind, res)
  R$padj <- p.adjust(R$p, "BH")
  list(n_tested = nrow(R), n_sig = sum(R$padj < 0.05), tab = R)
}

## ---- config -----------------------------------------------------------------
cfg_path <- file.path(cfg_dir, paste0(id, ".json"))
if (!file.exists(cfg_path)) {
  sha <- tryCatch(trimws(system("git rev-parse --short HEAD", intern = TRUE)[1]),
                  error = function(e) "unknown", warning = function(w) "unknown")
  writeLines(c("{",
    paste0('  "id": "', id, '",'),
    '  "family": "permutation",',
    '  "spec": "sample-label permutation null for the 019 BPD module rank-enrichment count",',
    paste0('  "contrast": "', PUBLISHED, '",'),
    paste0('  "k_total": ', K_TOTAL, ','),
    paste0('  "min_module_n": ', MIN_MOD_N, ','),
    '  "estimator": "add-one (1 + sum(null >= obs)) / (1 + K)",',
    '  "partition": "network/runs/grp-051/genes_control_A.csv",',
    paste0('  "seed_base": ', SEED_BASE, ','),
    paste0('  "git_sha": "', if (length(sha) && !is.na(sha)) sha else "unknown", '",'),
    paste0('  "r_version": "', paste(R.version$major, R.version$minor, sep = "."), '"'),
    "}"), cfg_path)
  msg("config written: ", cfg_path)
}

## ---- observed, and the check that the copied fit is the published one -------
obs_file <- file.path(run_dir, "observed.rds")
if (!file.exists(obs_file)) {
  msg("refitting ", PUBLISHED, " unpermuted, to check this copy against de_final ...")
  pv_obs <- fit_wb_p(CONTRAST)

  pub <- read.csv(gzfile(sprintf("de_final/results/DE_%s.csv.gz", PUBLISHED)))
  shared <- intersect(names(pv_obs), pub$gene)
  stopifnot(length(shared) > 10000)
  dif <- max(abs(pv_obs[shared] - pub$P[match(shared, pub$gene)]))
  msg("  ", length(shared), " genes shared with the published table; max |dp| = ",
      signif(dif, 3))
  if (dif > 1e-8) stop("the copied fit does not reproduce de_final/results/DE_",
                       PUBLISHED, ".csv.gz (max |dp| = ", signif(dif, 3),
                       "). Fix the copy before running a null for it.")
  msg("  copy reproduces the published fit exactly")

  o <- count_enriched(pv_obs)
  saveRDS(list(n_tested = o$n_tested, n_sig = o$n_sig, tab = o$tab), obs_file)
  write.csv(o$tab, file.path(run_dir, "observed_module_enrichment.csv"), row.names = FALSE)
  msg("observed: ", o$n_sig, " of ", o$n_tested, " modules at BH<0.05")
}
obs <- readRDS(obs_file)

## ---- the permutation draws --------------------------------------------------
s_bpd  <- p$sel[[CONTRAST]]
grp_in <- s_bpd$grp
msg("sample set: n=", length(grp_in), " (",
    sum(grp_in == "case"), " case, ", sum(grp_in == "ctrl"), " ctrl)")

for (k in if (AGG_ONLY) integer(0) else seq(K_FROM, K_TO)) {
  ck <- file.path(run_dir, sprintf("perm_%03d.rds", k))
  if (file.exists(ck)) next
  set.seed(SEED_BASE + k)
  gperm <- sample(grp_in)              # labels move, covariates stay with samples
  pv <- tryCatch(fit_wb_p(CONTRAST, gperm),
                 error = function(e) { msg("draw ", k, " FAILED: ", conditionMessage(e)); NULL })
  if (is.null(pv)) next
  r <- count_enriched(pv)
  saveRDS(list(draw = k, n_tested = r$n_tested, n_sig = r$n_sig), ck)
  msg("  draw ", k, "/", K_TO, ": ", r$n_sig, " of ", r$n_tested, " modules")
}

## ---- aggregate --------------------------------------------------------------
fs <- sort(list.files(run_dir, "^perm_[0-9]+[.]rds$", full.names = TRUE))
if (!length(fs)) { msg("no draws yet; nothing to aggregate"); quit(save = "no", status = 0L) }
draws <- do.call(rbind, lapply(fs, function(f) {
  x <- readRDS(f); data.frame(draw = x$draw, n_sig = x$n_sig, n_tested = x$n_tested)
}))
draws <- draws[order(draws$draw), ]
write.csv(draws, file.path(run_dir, "null_draws.csv"), row.names = FALSE)

K <- nrow(draws)
pval <- (1 + sum(draws$n_sig >= obs$n_sig)) / (1 + K)
res <- data.frame(contrast = PUBLISHED, observed_n_sig = obs$n_sig,
                  n_modules_tested = obs$n_tested, K = K,
                  null_mean = mean(draws$n_sig), null_sd = sd(draws$n_sig),
                  null_median = median(draws$n_sig),
                  null_q95 = unname(quantile(draws$n_sig, 0.95)),
                  null_max = max(draws$n_sig),
                  n_null_ge_obs = sum(draws$n_sig >= obs$n_sig),
                  empirical_p_addone = pval, p_floor = 1 / (1 + K))
write.csv(res, file.path(run_dir, "results_enrichment_null.csv"), row.names = FALSE)

msg("")
msg("observed  : ", obs$n_sig, " of ", obs$n_tested, " modules at BH<0.05")
msg("null (K=", K, "): mean ", sprintf("%.1f", mean(draws$n_sig)),
    ", median ", median(draws$n_sig),
    ", 95th pct ", sprintf("%.0f", quantile(draws$n_sig, 0.95)),
    ", max ", max(draws$n_sig))
msg("draws >= observed: ", sum(draws$n_sig >= obs$n_sig), " of ", K)
msg("empirical p (add-one): ", signif(pval, 3), "   floor ", signif(1 / (1 + K), 3))
# The verdict line used to hardcode "27 of 69", which is the BPD observed
# count. On the lithium row that printed a sentence about a number the row
# had never computed. Say what this row actually found.
msg(if (pval < 0.05)
      sprintf("  -> %s: %d of %d exceeds what co-expression alone produces (p = %.3f)",
              PUBLISHED, obs$n_sig, obs$n_tested, pval)
    else
      sprintf("  -> %s: %d of %d is WITHIN what co-expression alone produces (p = %.3f); the claim must be withdrawn",
              PUBLISHED, obs$n_sig, obs$n_tested, pval))

if (K < K_TOTAL) {
  msg("NOTE: ", K, " of ", K_TOTAL, " draws present; this is a partial but valid null.")
  quit(save = "no", status = 1L)   # requeue until the slices land
}
msg("=== ", id, " COMPLETE ===")
