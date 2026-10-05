#!/usr/bin/env Rscript
# ==============================================================================
# 97_composition_by_group.R -- does blood composition differ by GROUP, or only
# by lithium?
#
# The expression results say the bipolar signature lives almost entirely in the
# medicated patients. If that signature is compositional, the composition itself
# must show the same pattern:
#
#     control  ~=  bipolar I OFF lithium  !=  bipolar I ON lithium
#
# If instead bipolar-off-lithium already differs compositionally from controls,
# the compositional account of the expression result is in trouble.
#
# Tested on the ILR balances rather than the raw fractions because fractions are
# constrained to sum to 1 and cannot all enter a linear model; the balances are
# unconstrained real coordinates and are what the project analyses.
# ==============================================================================

source("scripts/common.R")
p <- load_prep()
OUT <- file.path(EXP_HOME, "runs", "composition_by_group")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

m <- p$meta
m$tob <- p$TOB[m$title, "tobacco_imp_01"]
m$grp <- factor(ifelse(m$dx == "Control", "control",
                ifelse(m$dx == "BP1" & m$lithium == 0, "BP1_off",
                ifelse(m$dx == "BP1" & m$lithium == 1, "BP1_on", NA))),
                levels = c("control", "BP1_off", "BP1_on"))
m <- m[!is.na(m$grp), ]
cat(sprintf("groups: %s\n", paste(sprintf("%s=%d", levels(m$grp), table(m$grp)), collapse = "  ")))

B <- cbind(p$ILR[m$title, , drop = FALSE], p$lineage[m$title, , drop = FALSE])
covars <- c("age", "sex", "tob", "rin", "plate", "seqpc1", "seqpc2", "seqpc3", "group")
ok <- complete.cases(m[, covars]); m <- m[ok, ]; B <- B[ok, , drop = FALSE]
m$plate <- droplevels(m$plate); m$sex <- droplevels(m$sex); m$group <- droplevels(m$group)

cat(sprintf("n = %d after covariate completeness\n\n", nrow(m)))

res <- do.call(rbind, lapply(colnames(B), function(v) {
  y <- B[, v]
  d <- cbind(m, y = y)
  f  <- lm(as.formula(paste("y ~ grp +", paste(covars, collapse = " + "))), data = d)
  f0 <- update(f, . ~ . - grp)
  an <- anova(f0, f)
  cf <- summary(f)$coefficients
  # Adjusted group means, holding covariates at their observed distribution.
  adj <- tapply(residuals(f0) , d$grp, mean)
  sdy <- sd(residuals(f0))
  data.frame(
    variable = v,
    p_overall_group = an$`Pr(>F)`[2],
    off_vs_control_beta = cf["grpBP1_off", 1], off_vs_control_p = cf["grpBP1_off", 4],
    on_vs_control_beta  = cf["grpBP1_on", 1],  on_vs_control_p  = cf["grpBP1_on", 4],
    off_d = unname((adj["BP1_off"] - adj["control"]) / sdy),
    on_d  = unname((adj["BP1_on"]  - adj["control"]) / sdy),
    stringsAsFactors = FALSE)
}))
res$fdr_overall <- p.adjust(res$p_overall_group, "BH")

cat("=== composition by group, covariate-adjusted ===\n")
cat("   off = bipolar I NOT on lithium vs control | on = bipolar I on lithium vs control\n\n")
print(res[, c("variable", "off_vs_control_p", "off_d", "on_vs_control_p", "on_d",
              "p_overall_group", "fdr_overall")], row.names = FALSE, digits = 3)

## ---- the direct contrast: on vs off ---------------------------------------
cat("\n=== bipolar I ON vs OFF lithium (the within-diagnosis contrast) ===\n")
bp <- m[m$grp != "control", ]; bp$grp <- droplevels(bp$grp)
Bb <- B[m$grp != "control", , drop = FALSE]
bp$plate <- droplevels(bp$plate)
cv2 <- setdiff(covars, "group")     # constant among cases
onoff <- do.call(rbind, lapply(colnames(Bb), function(v) {
  d <- cbind(bp, y = Bb[, v])
  f <- lm(as.formula(paste("y ~ grp +", paste(cv2, collapse = " + "))), data = d)
  cf <- summary(f)$coefficients["grpBP1_on", ]
  f0 <- update(f, . ~ . - grp)
  data.frame(variable = v, beta = cf[1], p = cf[4],
             d = unname((mean(residuals(f0)[d$grp == "BP1_on"]) -
                         mean(residuals(f0)[d$grp == "BP1_off"])) / sd(residuals(f0))),
             stringsAsFactors = FALSE)
}))
onoff$fdr <- p.adjust(onoff$p, "BH")
print(onoff, row.names = FALSE, digits = 3)

cat("\n=== verdict ===\n")
nb <- sum(res$off_vs_control_p < 0.05); na_ <- sum(res$on_vs_control_p < 0.05)
cat(sprintf("  balances/lineages differing from control at p<.05:  OFF lithium %d/%d | ON lithium %d/%d\n",
            nb, nrow(res), na_, nrow(res)))
cat(sprintf("  mean |Cohen d| vs control:  OFF %.3f | ON %.3f\n",
            mean(abs(res$off_d)), mean(abs(res$on_d))))
cat("  If OFF-lithium patients resemble controls compositionally while ON-lithium\n")
cat("  patients do not, the compositional account of the expression result holds.\n")

write.csv(res, file.path(OUT, "group_vs_control.csv"), row.names = FALSE)
write.csv(onoff, file.path(OUT, "on_vs_off.csv"), row.names = FALSE)
cat("\nwrote runs/composition_by_group/\n")
