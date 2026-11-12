#!/usr/bin/env Rscript
# ===========================================================================
# eos-057  Is the eosinophil module's lithium association driven by a handful
#          of subjects?
#
# THE WORRY
#
# mediumpurple2 has 35 genes. An eigengene built from 35 genes is the first
# principal component of a small matrix, and a small number of unusual
# subjects can move it a long way. Every other caveat on this result has been
# checked -- it replicates across three decompositions, survives the gene
# filter, is not diagnosis, and tracks the measured eosinophil fraction -- but
# nobody has asked whether a few people are carrying it.
#
# THREE TESTS, WEAKEST TO STRONGEST
#
# 1. Leave-one-out. Drop each subject in turn, refit. If one subject's removal
#    pushes p above 0.05, the result rests on that subject.
#
# 2. Bootstrap over subjects, 2000 resamples, REBUILDING THE EIGENGENE inside
#    each one. Rebuilding matters: an eigengene computed once on all samples
#    and then resampled tests only the regression's stability, not the
#    module summary's. The eigengene's sign is arbitrary (it is a PC), so
#    each resample is sign-aligned to the full-sample eigengene before use.
#
# 3. Fragility. Remove the single most influential subject, refit, repeat.
#    How many have to go before p crosses 0.05? A result that survives one
#    removal but not three is not the same as one that survives twenty.
#    This is the number that actually answers the question.
#
# THE CONTRAST. trait-053's M_051b on BP1 only: lithium users against
# non-users inside bipolar I, adjusted for age, sex, RIN, plate and the four
# ILR balances. Diagnosis is constant there and drops out. Observed
# p = 4.23e-05. BP1-only is the clean comparison -- the all-sample version
# mixes in the case-control contrast, which runs against a QC-filtered
# control group.
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "eos-057"
stopifnot(id == "eos-057")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

run_dir <- file.path("network/runs", id)
cfg_dir <- "network/config"
for (d_ in c(run_dir, cfg_dir)) dir.create(d_, showWarnings = FALSE, recursive = TRUE)

MODULE  <- "mediumpurple2"
N_BOOT  <- 2000L
SEED    <- 20261009L + 57000L
ALPHA   <- 0.05

suppressPackageStartupMessages({ library(WGCNA) })
enableWGCNAThreads()
options(stringsAsFactors = FALSE)
source("network/00_load.R")
d <- load_project(cache = TRUE, verbose = FALSE)

## ---- the module's genes -----------------------------------------------------
gm <- read.csv("network/runs/ann-044/reference_gene_modules.csv")
genes_mod <- gm$gene[gm$module == MODULE]
msg(MODULE, ": ", length(genes_mod), " genes")
stopifnot(length(genes_mod) >= 10)

E <- expr_for(d, group = NULL)                       # samples x genes, offset-removed
stopifnot(all(genes_mod %in% colnames(E)))
Emod <- E[, genes_mod, drop = FALSE]

## ---- covariates, exactly as trait-053 ---------------------------------------
md <- d$meta
md$li <- as.integer(d$groups == "bp_lith")
ilr <- d$ilr[match(md$title, d$ilr$sample), c("b1", "b2", "b3", "b4"), drop = FALSE]
keep <- !is.na(d$groups) & d$groups %in% c("bp_lith", "bp_nolith")   # BP1 only
msg("BP1 only: n=", sum(keep), " (", sum(md$li[keep]), " on lithium, ",
    sum(1 - md$li[keep]), " off)")

covars <- data.frame(li = md$li, age = md$age, sex = md$sex, rin = md$rin,
                     plate = factor(md$plate), b1 = ilr$b1, b2 = ilr$b2,
                     b3 = ilr$b3, b4 = ilr$b4)[keep, ]
Ekeep <- Emod[keep, , drop = FALSE]
n <- nrow(Ekeep)
stopifnot(n == sum(keep))

# The eigengene is the first PC of the module's genes. Recomputing it inside a
# resample is the point of the exercise, so it is a function, not a constant.
eigengene <- function(M) {
  me <- moduleEigengenes(M, colors = rep(MODULE, ncol(M)), excludeGrey = FALSE)
  as.numeric(me$eigengenes[[1]])
}
align <- function(y, ref) if (sum(y * ref) < 0) -y else y   # PCs have arbitrary sign

fit_p <- function(y, cv) {
  f <- lm(y ~ li + age + sex + rin + plate + b1 + b2 + b3 + b4,
          data = cbind(cv, y = y))
  s <- summary(f)$coefficients
  if ("li" %in% rownames(s)) c(beta = s["li", 1], p = s["li", 4]) else c(beta = NA, p = NA)
}

## ---- observed ---------------------------------------------------------------
ME_full <- eigengene(Ekeep)
obs <- fit_p(ME_full, covars)
msg(sprintf("observed: beta = %+.5f   p = %.3g", obs["beta"], obs["p"]))

cfg_path <- file.path(cfg_dir, paste0(id, ".json"))
if (!file.exists(cfg_path)) {
  sha <- tryCatch(trimws(system("git rev-parse --short HEAD", intern = TRUE)[1]),
                  error = function(e) "unknown", warning = function(w) "unknown")
  writeLines(c("{",
    paste0('  "id": "', id, '",'), '  "family": "sensitivity",',
    '  "spec": "leave-one-out, subject bootstrap and fragility for the eosinophil module lithium association",',
    paste0('  "module": "', MODULE, '", "n_genes": ', length(genes_mod), ','),
    '  "model": "ME ~ li + age + sex + rin + plate + b1..b4, BP1 only",',
    paste0('  "n_boot": ', N_BOOT, ', "alpha": ', ALPHA, ','),
    '  "eigengene": "rebuilt inside every resample, sign-aligned to the full-sample PC",',
    paste0('  "seed": ', SEED, ','),
    paste0('  "git_sha": "', if (length(sha) && !is.na(sha)) sha else "unknown", '",'),
    paste0('  "r_version": "', paste(R.version$major, R.version$minor, sep = "."), '"'),
    "}"), cfg_path)
  msg("config written: ", cfg_path)
}

## ---- 1. leave-one-out -------------------------------------------------------
msg("leave-one-out over ", n, " subjects ...")
loo <- do.call(rbind, lapply(seq_len(n), function(i) {
  y  <- align(eigengene(Ekeep[-i, , drop = FALSE]), ME_full[-i])
  r  <- fit_p(y, covars[-i, , drop = FALSE])
  data.frame(dropped = i, beta = unname(r["beta"]), p = unname(r["p"]))
}))
write.csv(loo, file.path(run_dir, "leave_one_out.csv"), row.names = FALSE)
msg(sprintf("  p range %.3g to %.3g; subjects whose removal lifts p above %.2f: %d",
            min(loo$p), max(loo$p), ALPHA, sum(loo$p > ALPHA)))

## ---- 2. bootstrap over subjects --------------------------------------------
msg("bootstrap, ", N_BOOT, " resamples, eigengene rebuilt each time ...")
set.seed(SEED)
boot <- do.call(rbind, lapply(seq_len(N_BOOT), function(b) {
  idx <- sample.int(n, n, replace = TRUE)
  cv  <- covars[idx, , drop = FALSE]
  # A resample can lose a plate level or leave lithium constant; both make the
  # design rank-deficient. Those are resampling artifacts, not evidence, so
  # they are recorded as NA and counted rather than silently dropped.
  if (length(unique(cv$li)) < 2) return(data.frame(b = b, beta = NA, p = NA))
  cv$plate <- droplevels(factor(cv$plate))
  y <- tryCatch(align(eigengene(Ekeep[idx, , drop = FALSE]), ME_full[idx]),
                error = function(e) NULL)
  if (is.null(y)) return(data.frame(b = b, beta = NA, p = NA))
  r <- tryCatch(fit_p(y, cv), error = function(e) c(beta = NA, p = NA))
  data.frame(b = b, beta = unname(r["beta"]), p = unname(r["p"]))
}))
write.csv(boot, file.path(run_dir, "bootstrap.csv"), row.names = FALSE)
ok <- boot[!is.na(boot$beta), ]
ci <- quantile(ok$beta, c(0.025, 0.975))
msg(sprintf("  %d of %d usable; beta 95%% CI [%+.5f, %+.5f]; crosses zero: %s",
            nrow(ok), N_BOOT, ci[1], ci[2],
            if (ci[1] <= 0 && ci[2] >= 0) "YES" else "no"))
msg(sprintf("  resamples with p < %.2f: %.1f%%; same sign as observed: %.1f%%",
            ALPHA, 100 * mean(ok$p < ALPHA),
            100 * mean(sign(ok$beta) == sign(obs["beta"]))))

## ---- 3. fragility -----------------------------------------------------------
# Greedy: at each step drop the subject whose removal lowers the evidence most,
# and stop when p crosses alpha. This is deliberately the most hostile way to
# attack the result -- a real adversary choosing who to delete.
msg("fragility: removing the most damaging subject at each step ...")
alive <- seq_len(n); steps <- list(); cur <- obs
for (s in seq_len(40)) {
  cand <- vapply(alive, function(i) {
    k <- setdiff(alive, i)
    y <- align(eigengene(Ekeep[k, , drop = FALSE]), ME_full[k])
    cv <- covars[k, , drop = FALSE]; cv$plate <- droplevels(factor(cv$plate))
    if (length(unique(cv$li)) < 2) return(NA_real_)
    unname(fit_p(y, cv)["p"])
  }, numeric(1))
  if (all(is.na(cand))) break
  worst <- alive[which.max(replace(cand, is.na(cand), -Inf))]
  cur_p <- max(cand, na.rm = TRUE)
  alive <- setdiff(alive, worst)
  steps[[s]] <- data.frame(step = s, removed_subject = worst, p_after = cur_p,
                           n_remaining = length(alive))
  if (cur_p > ALPHA) break
}
FR <- do.call(rbind, steps)
write.csv(FR, file.path(run_dir, "fragility.csv"), row.names = FALSE)
n_needed <- if (any(FR$p_after > ALPHA)) min(FR$step[FR$p_after > ALPHA]) else NA_integer_
msg(sprintf("  subjects that must be removed to push p above %.2f: %s (of %d)",
            ALPHA, if (is.na(n_needed)) paste0(">", nrow(FR)) else n_needed, n))

## ---- summary ----------------------------------------------------------------
summ <- data.frame(
  module = MODULE, n_genes = length(genes_mod), n_subjects = n,
  observed_beta = unname(obs["beta"]), observed_p = unname(obs["p"]),
  loo_p_max = max(loo$p), loo_n_breaking = sum(loo$p > ALPHA),
  boot_usable = nrow(ok), boot_ci_lo = unname(ci[1]), boot_ci_hi = unname(ci[2]),
  boot_frac_sig = mean(ok$p < ALPHA), boot_frac_same_sign = mean(sign(ok$beta) == sign(obs["beta"])),
  fragility_n = n_needed,
  fragility_pct = if (is.na(n_needed)) NA_real_ else 100 * n_needed / n)
write.csv(summ, file.path(run_dir, "results_robustness.csv"), row.names = FALSE)

msg("")
msg("=== verdict ===")
msg(sprintf("  leave-one-out : worst p = %.3g; %d of %d subjects break it",
            max(loo$p), sum(loo$p > ALPHA), n))
msg(sprintf("  bootstrap     : %.1f%% of resamples significant, CI %s zero",
            100 * mean(ok$p < ALPHA), if (ci[1] <= 0 && ci[2] >= 0) "INCLUDES" else "excludes"))
msg(sprintf("  fragility     : %s subject(s) must be deleted to break it, %s of the sample",
            if (is.na(n_needed)) paste0(">", nrow(FR)) else as.character(n_needed),
            if (is.na(n_needed)) "more than tested" else sprintf("%.1f%%", 100 * n_needed / n)))
msg("=== ", id, " COMPLETE ===")
