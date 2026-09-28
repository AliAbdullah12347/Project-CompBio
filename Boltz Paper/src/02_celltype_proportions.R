#!/usr/bin/env Rscript
# Cell-type proportions: Boltz et al. (2024) Table 1, the ">0.02" cell-type selection, and the three
# logistic-regression contrasts of Figure S6 / p. 330.
#
# Paper, Subjects and methods: "We built logistic regression models to evaluate the effect of
# cell-type proportion on case or control status, as well as lithium-use status within only the BP
# individuals. These models included the proportion of one cell type at a time, along with
# covariates including age, sex, RNA concentration, and RNA integrity number (RIN) ... glm() ...
# family = binomial."
#
# Deviations (see REPRODUCIBILITY.md): fractions are the CIBERSORT LM22 values deposited in
# GSE124326, not CIBERSORTx with batch correction; RNA concentration is not in the GEO deposit and
# is omitted; cases are BP-I + BP-II (the cohort has no SCZ).

suppressPackageStartupMessages(library(stats))
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
dir.create(file.path(root, "results"), showWarnings = FALSE)

meta <- read.csv(file.path(root, "data/processed/sample_metadata.csv"), check.names = FALSE,
                 colClasses = c(title = "character"))
stopifnot(nrow(meta) == 444)

LM22 <- c("b.cells.naive", "b.cells.memory", "plasma.cells", "t.cells.cd8", "t.cells.cd4.naive",
          "t.cells.cd4.memory.resting", "t.cells.cd4.memory.activated", "t.cells.follicular.helper",
          "t.cells.regulatory", "t.cells.gamma.delta", "nk.cells.resting", "nk.cells.activated",
          "monocytes", "macrophages.m0", "macrophages.m1", "macrophages.m2",
          "dendritic.cells.resting", "dendritic.cells.activated", "mast.cells.resting",
          "mast.cells.activated", "eosinophils", "neutrophils")
# Boltz Table 1: the eight types they kept, in their order.
BOLTZ8 <- c("b.cells.naive", "b.cells.memory", "t.cells.cd8", "t.cells.cd4.naive",
            "t.cells.cd4.memory.resting", "nk.cells.resting", "monocytes", "neutrophils")
ref <- read.csv(file.path(root, "data/reference/paper_targets.csv"))
tgt <- setNames(ref$value, ref$quantity)

# ---- Table 1 --------------------------------------------------------------------------------
cat("=== Table 1: mean (s.d.) cell-type proportion, all 444 samples ===\n")
tab1 <- data.frame(cell_type = LM22,
                   mean = sapply(LM22, function(k) mean(meta[[k]])),
                   sd = sapply(LM22, function(k) sd(meta[[k]])),
                   pct_nonzero = sapply(LM22, function(k) 100 * mean(meta[[k]] > 0)))
tab1$passes_0.02 <- tab1$mean > 0.02
tab1$in_boltz_table1 <- tab1$cell_type %in% BOLTZ8
tab1$boltz_mean <- unname(tgt[paste0("prop_mean_", tab1$cell_type)])
tab1$boltz_sd <- unname(tgt[paste0("prop_sd_", tab1$cell_type)])
tab1 <- tab1[order(-tab1$mean), ]
print(format(tab1, digits = 3), row.names = FALSE)

passing <- tab1$cell_type[tab1$passes_0.02]
cat(sprintf("\nTypes with mean > 0.02 here: %d  [%s]\n", length(passing), paste(passing, collapse = ", ")))
cat(sprintf("In Boltz's eight but below 0.02 here: [%s]\n", paste(setdiff(BOLTZ8, passing), collapse = ", ")))
cat(sprintf("Above 0.02 here but not in Boltz's eight: [%s]\n", paste(setdiff(passing, BOLTZ8), collapse = ", ")))
# bMIND and the adjusted DE use the types this rule selects (the method as written); Boltz's own
# eight are kept in the regressions too so each can be matched to a published p-value.
CELLS <- union(BOLTZ8, passing)
cat("Regressions below cover both sets; bMIND (script 04) uses the rule-selected set.\n")
write.csv(tab1, file.path(root, "results/table1_celltype_proportions.csv"), row.names = FALSE)
writeLines(passing, file.path(root, "data/processed/celltypes_selected.txt"))

# ---- Logistic regressions -------------------------------------------------------------------
meta$sex <- factor(meta$sex)
contrasts <- list(
  case_vs_control     = list(rows = rep(TRUE, nrow(meta)), y = meta$case,
                             note = "BP-I + BP-II (1) vs controls (0), all 444"),
  lithium_user_vs_non = list(rows = meta$case == 1, y = meta$lithium,
                             note = "lithium users (1) vs non-users (0), cases only"),
  nonuser_vs_control  = list(rows = meta$case == 0 | meta$lithium == 0, y = meta$case,
                             note = "lithium non-user cases (1) vs controls (0)"))

fit_one <- function(d, y, k) {
  d$y <- y
  d$p <- d[[k]]
  m <- glm(y ~ p + age + sex + rin, data = d, family = binomial)
  s <- summary(m)$coefficients["p", ]
  # beta is per unit proportion; report per 0.1 (10 percentage points) as the readable scale
  c(beta = unname(s["Estimate"]), se = unname(s["Std. Error"]), z = unname(s["z value"]),
    p = unname(s["Pr(>|z|)"]), or_per_0.1 = exp(0.1 * unname(s["Estimate"])),
    mean_y1 = mean(d$p[y == 1]), mean_y0 = mean(d$p[y == 0]))
}

res <- do.call(rbind, lapply(names(contrasts), function(cn) {
  cc <- contrasts[[cn]]
  d <- meta[cc$rows, ]
  y <- cc$y[cc$rows]
  cat(sprintf("\n=== %s: %s  (n1 = %d, n0 = %d) ===\n", cn, cc$note, sum(y == 1), sum(y == 0)))
  out <- do.call(rbind, lapply(CELLS, function(k) {
    r <- fit_one(d, y, k)
    data.frame(contrast = cn, cell_type = k, n1 = sum(y == 1), n0 = sum(y == 0), t(r))
  }))
  out$bonferroni_p <- pmin(1, out$p * length(CELLS))
  out$higher_in <- ifelse(out$beta > 0, "group 1", "group 0")
  print(format(out[, c("cell_type", "beta", "p", "bonferroni_p", "mean_y1", "mean_y0")], digits = 3),
        row.names = FALSE)
  out
}))

# Attach the paper's p-values. Boltz's "CD4 T cells" in the case/control sentence is unspecified,
# so it is matched to both CD4 types.
res$boltz_p <- NA
key <- paste(res$contrast, res$cell_type)
put <- function(cn, k, q) res$boltz_p[key == paste(cn, k)] <<- tgt[q]
put("case_vs_control", "t.cells.cd4.naive", "p_casecontrol_cd4_t")
put("case_vs_control", "t.cells.cd4.memory.resting", "p_casecontrol_cd4_t")
put("case_vs_control", "nk.cells.resting", "p_casecontrol_nk.cells.resting")
put("case_vs_control", "neutrophils", "p_casecontrol_neutrophils")
for (k in c("t.cells.cd4.naive", "t.cells.cd4.memory.resting", "nk.cells.resting", "neutrophils"))
  put("lithium_user_vs_non", k, paste0("p_lithium_", k))
write.csv(res, file.path(root, "results/celltype_logistic_regressions.csv"), row.names = FALSE)

# ---- Figure S6 analogue: neutrophil proportion adjusted for covariates, three contrasts --------
png(file.path(root, "results/figS6_neutrophils.png"), width = 1500, height = 520, res = 130)
par(mfrow = c(1, 3), mar = c(4, 4.5, 3, 1))
for (cn in names(contrasts)) {
  cc <- contrasts[[cn]]
  d <- meta[cc$rows, ]
  y <- cc$y[cc$rows]
  adj <- resid(lm(neutrophils ~ age + sex + rin, data = d)) + mean(d$neutrophils)
  lab <- switch(cn, case_vs_control = c("Control", "Case"),
                lithium_user_vs_non = c("Non-user", "Li user"),
                nonuser_vs_control = c("Control", "Non-user"))
  p <- res$p[res$contrast == cn & res$cell_type == "neutrophils"]
  boxplot(adj ~ factor(y, levels = 0:1, labels = lab), col = c("grey85", "grey60"),
          ylab = "Neutrophil proportion (covariate-adjusted)", xlab = "",
          main = sprintf("%s\nlogistic p = %.2g", gsub("_", " ", cn), p))
}
dev.off()
cat("\nWrote results/table1_celltype_proportions.csv, celltype_logistic_regressions.csv, figS6_neutrophils.png\n")
