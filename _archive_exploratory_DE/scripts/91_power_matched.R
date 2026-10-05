#!/usr/bin/env Rscript
# ==============================================================================
# 91_power_matched.R -- is "4 DEGs off lithium vs 2,603 on lithium" just power?
#
# 90_signature_overlap.R found that bipolar I subjects NOT on lithium are almost
# indistinguishable from controls (4 DEGs) while those ON lithium differ
# strongly (2,603). The obvious objection is sample size: 74 cases against 152.
#
# This settles it by brute force. Subsample the on-lithium group down to 74 --
# exactly matching the off-lithium group, against the identical 234 controls --
# and refit. If the power-matched on-lithium contrast still yields far more than
# 4 DEGs, sample size is not the explanation.
#
# Two calibrations run alongside:
#   * control-vs-control: 74 random controls against the remaining 160. Should
#     give ~0. If it does not, the pipeline is anticonservative and every count
#     in this folder is suspect.
#   * off-lithium subsampled to 74 (i.e. itself) is the observed value, 4.
# ==============================================================================

source("scripts/common.R")
p <- load_prep()
set.seed(481)

OUT <- file.path(EXP_HOME, "runs", "power_matched")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

N_DRAWS <- 40          # reported exactly; nothing is silently truncated
m <- p$meta
m$tob <- p$TOB[m$title, "tobacco_imp_01"]

ctrl <- m$title[m$dx == "Control"]
off  <- m$title[m$dx == "BP1" & m$lithium == 0]
on   <- m$title[m$dx == "BP1" & m$lithium == 1]
cat(sprintf("controls %d | BP1 off lithium %d | BP1 on lithium %d\n",
            length(ctrl), length(off), length(on)))

fit_two_group <- function(case_ids, ctrl_ids) {
  d <- m[match(c(ctrl_ids, case_ids), m$title), ]
  d$grp <- factor(c(rep("ctrl", length(ctrl_ids)), rep("case", length(case_ids))),
                  c("ctrl", "case"))
  covars <- c("age", "sex", "tob", "rin", "plate", "seqpc1", "seqpc2", "seqpc3")
  if (nlevels(droplevels(d$group)) > 1) covars <- c(covars, "group")
  d$group <- droplevels(d$group); d$plate <- droplevels(d$plate); d$sex <- droplevels(d$sex)
  d <- d[complete.cases(d[, c("grp", covars)]), ]
  X <- model.matrix(as.formula(paste("~ grp +", paste(covars, collapse = " + "))), data = d)
  ne <- nonEstimable(X); if (!is.null(ne)) X <- X[, setdiff(colnames(X), ne), drop = FALSE]
  cnt <- p$counts[, d$title, drop = FALSE]
  cnt <- cnt[filter_genes(cnt, 10, 0.90), , drop = FALSE]
  v <- voom(calcNormFactors(DGEList(cnt), method = "TMM"), X)
  tt <- topTable(eBayes(lmFit(v, X)), coef = "grpcase", number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  list(deg = n_deg(tt), genes = nrow(tt), pi0 = pi0_storey(tt$P.Value), tt = tt)
}

cat(sprintf("\n== A. ON-LITHIUM subsampled to n=74, vs all %d controls (%d draws) ==\n",
            length(ctrl), N_DRAWS))
on74 <- vapply(seq_len(N_DRAWS), function(i) {
  r <- fit_two_group(sample(on, length(off)), ctrl)
  cat(sprintf("  draw %2d  DEG=%5d  genes=%5d  pi0=%.3f\n", i, r$deg, r$genes, r$pi0))
  c(deg = r$deg, pi0 = r$pi0)
}, numeric(2))

cat(sprintf("\n== B. CONTROL vs CONTROL null: 74 controls vs the other %d (%d draws) ==\n",
            length(ctrl) - 74, N_DRAWS))
cc <- vapply(seq_len(N_DRAWS), function(i) {
  s <- sample(ctrl); r <- fit_two_group(s[1:74], s[75:length(s)])
  cat(sprintf("  draw %2d  DEG=%5d  pi0=%.3f\n", i, r$deg, r$pi0))
  c(deg = r$deg, pi0 = r$pi0)
}, numeric(2))

obs_off <- fit_two_group(off, ctrl)

res <- data.frame(
  comparison = c("BP1 off lithium (74) vs control  [OBSERVED]",
                 "BP1 on lithium subsampled to 74 vs control",
                 "control (74) vs control (160)  [NULL]"),
  n_case = 74, n_draws = c(1, N_DRAWS, N_DRAWS),
  deg_median = c(obs_off$deg, median(on74["deg", ]), median(cc["deg", ])),
  deg_min = c(obs_off$deg, min(on74["deg", ]), min(cc["deg", ])),
  deg_max = c(obs_off$deg, max(on74["deg", ]), max(cc["deg", ])),
  pi0_median = round(c(obs_off$pi0, median(on74["pi0", ]), median(cc["pi0", ])), 3),
  stringsAsFactors = FALSE)
cat("\n== POWER-MATCHED COMPARISON ==\n"); print(res, row.names = FALSE)

p_emp <- mean(on74["deg", ] <= obs_off$deg)
cat(sprintf("\nDraws of the power-matched ON-lithium contrast giving <= %d DEGs: %d/%d (p = %.3f)\n",
            obs_off$deg, sum(on74["deg", ] <= obs_off$deg), N_DRAWS, p_emp))
cat(sprintf("Control-vs-control null median %g DEG -- pipeline is %s\n",
            median(cc["deg", ]),
            if (median(cc["deg", ]) <= 5) "calibrated" else "ANTICONSERVATIVE, investigate"))

write.csv(res, file.path(OUT, "summary.csv"), row.names = FALSE)
write.csv(data.frame(draw = seq_len(N_DRAWS), on_lithium_74 = on74["deg", ],
                     on_lithium_pi0 = on74["pi0", ],
                     control_null = cc["deg", ], control_null_pi0 = cc["pi0", ]),
          file.path(OUT, "draws.csv"), row.names = FALSE)
saveRDS(list(on74 = on74, cc = cc, obs_off = obs_off$deg, res = res),
        file.path(OUT, "fits.rds"))
cat("\nwrote runs/power_matched/\n")
