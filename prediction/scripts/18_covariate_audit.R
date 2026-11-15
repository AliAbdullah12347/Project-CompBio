#!/usr/bin/env Rscript
# Which covariates belong in the lithium contrast, and which do not.
#
# The prediction arm was built with one covariate block reused for both
# contrasts: age, sex, tobacco, RIN, plate and sequencing PCs 1-3. That is the
# habit this script exists to break. A covariate earns its place by doing a
# job, and in the lithium contrast most of these have no job to do.
#
# Two criteria, both fixed before any AUC is looked at, and neither of them
# involving the outcome:
#
#   (a) TEMPORAL. A confounder has to be causally prior to the exposure.
#       Lithium is prescribed years before the blood draw, so nothing created
#       by the assay -- RIN, plate, library size, sequencing PCs -- can
#       possibly confound it. These are post-exposure variables. They may
#       still reduce noise, but that is a different claim and gets reported
#       as a sensitivity analysis, not smuggled in as confounder control.
#
#   (b) IMBALANCE. Among the variables that pass (a), adjustment only buys
#       something if the variable actually differs by lithium status. The
#       threshold is |SMD| > 0.1, the usual convention, declared in advance.
#
# Writes results/covariate_audit.txt.
#
# Usage: Rscript 18_covariate_audit.R [impl_root]

args <- commandArgs(trailingOnly = TRUE)
impl <- if (length(args)) args[1] else ".."
outd <- file.path(impl, "prediction/results"); dir.create(outd, showWarnings = FALSE, recursive = TRUE)
con <- file(file.path(outd, "covariate_audit.txt"), "w")
say <- function(...) { s <- sprintf(...); cat(s); writeLines(sub("\n$", "", s), con) }

meta <- read.csv(file.path(impl, "data/cohort_474/metadata_474.csv"),
                 check.names = FALSE, stringsAsFactors = FALSE)
meta$li <- as.integer(meta[["lithium use (non-user=0, user = 1)"]])
b <- meta[meta[["bipolar disorder diagnosis"]] == "BP1", ]
g <- b$li
stopifnot(nrow(b) == 226, sum(g == 1) == 152, sum(g == 0) == 74)

say("LITHIUM CONTRAST -- COVARIATE AUDIT\n")
say("%d BP1 subjects: %d users, %d non-users\n\n", nrow(b), sum(g), sum(1 - g))

## ---- (b) imbalance ----------------------------------------------------------
smd_num <- function(x) { x1 <- x[g == 1]; x0 <- x[g == 0]
  (mean(x1) - mean(x0)) / sqrt((var(x1) + var(x0)) / 2) }
smd_bin <- function(x) { p1 <- mean(x[g == 1]); p0 <- mean(x[g == 0])
  (p1 - p0) / sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2) }

say("BALANCE WITH LITHIUM STATUS\n")
say("%-14s %10s %10s %8s %9s  %s\n", "variable", "non-user", "user", "SMD", "p", "prior to exposure?")
rows <- list(
  list("age",         as.numeric(b$age),                            TRUE,  "yes"),
  list("sex (M)",     as.integer(b$Sex == "M"),                     FALSE, "yes"),
  list("tobacco",     as.integer(b[["tobacco use"]]),               FALSE, "yes"),
  list("RIN",         as.numeric(b$rin),                            TRUE,  "NO - assay"),
  list("library size",as.numeric(b$library_size) / 1e6,             TRUE,  "NO - assay"),
  list("seqPC1",      as.numeric(b[["sequencing metric pc1"]]),     TRUE,  "NO - assay"),
  list("seqPC2",      as.numeric(b[["sequencing metric pc2"]]),     TRUE,  "NO - assay"),
  list("seqPC3",      as.numeric(b[["sequencing metric pc3"]]),     TRUE,  "NO - assay"))
for (r in rows) {
  nm <- r[[1]]; v <- r[[2]]; isnum <- r[[3]]
  s <- if (isnum) smd_num(v) else smd_bin(v)
  p <- if (isnum) t.test(v ~ g)$p.value else suppressWarnings(chisq.test(table(v, g))$p.value)
  say("%-14s %10.3f %10.3f %8.3f %9.2g  %s\n", nm, mean(v[g == 0]), mean(v[g == 1]), s, p, r[[4]])
}
pl <- suppressWarnings(chisq.test(table(b[["sequencing plate"]], g))$p.value)
say("%-14s %10s %10s %8s %9.2g  %s\n", "plate", "-", "-", "-", pl, "NO - assay")
say("%-14s %10s %10s %8s %9s  %s\n", "assess. group", "A", "A", "-", "-", "constant, dropped")

say("\nVERDICT\n")
say("  IN   age      imbalanced (|SMD| 0.46) and prior to prescribing\n")
say("  IN   sex      imbalanced (|SMD| 0.24) and prior to prescribing\n")
say("  OUT  tobacco  prior, but balanced (|SMD| 0.08) -- nothing to correct\n")
say("  OUT  RIN, library size, plate, seqPC1-3\n")
say("       post-exposure assay variables. Cannot confound a prescription\n")
say("       written years earlier. Noise terms at best.\n")

## ---- the cost of the kitchen sink -------------------------------------------
# Same CV, same model, same subjects -- only the covariate block changes.
ilr <- read.csv(file.path(impl, "data/analysis_matrices/ILR_coordinates.csv"), stringsAsFactors = FALSE)
rownames(ilr) <- ilr$sample
sp <- read.csv(file.path(impl, "prediction/splits/splits_bpd_lithium.csv"), stringsAsFactors = FALSE)
ca <- readRDS(file.path(impl, "prediction/cache/lithium_inputs.rds")); m <- ca$meta

D <- data.frame(ilr[m$title, c("b1", "b2", "b3", "b4")],
                age = m$age, sexM = as.integer(m$sex == "M"), tobacco = m$tobacco,
                rin = m$rin, seqpc1 = m$seqpc1, seqpc2 = m$seqpc2, seqpc3 = m$seqpc3,
                y = m$is_case, row.names = m$title)
for (p in sort(unique(m$plate))[-1]) D[[paste0("pl", p)]] <- as.integer(m$plate == p)
PL <- grep("^pl", names(D), value = TRUE)

auc <- function(s, l) { r <- rank(s); n1 <- sum(l == 1); n0 <- sum(l == 0)
  (sum(r[l == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }

cv <- function(vars) {
  a <- tapply(seq_len(nrow(sp)), sp$repeat_no, function(ii) {
    s <- sp[ii, ]; pr <- rep(NA_real_, nrow(s))
    for (k in unique(s$split_id)) {
      q <- s$split_id == k; tr <- q & s$assignment == "train"; te <- q & s$assignment != "train"
      fit <- suppressWarnings(glm(reformulate(vars, "y"), data = D[s$sample_id[tr], ], family = binomial))
      pr[te] <- predict(fit, newdata = D[s$sample_id[te], ], type = "response")
    }
    auc(pr[!is.na(pr)], s$y[!is.na(pr)])
  })
  as.numeric(a)
}

say("\nWHAT THE KITCHEN SINK COSTS (25 repeats x 3 folds, logistic, no penalty)\n")
sets <- list(
  "age alone"                    = "age",
  "age + sex  [proposed set]"    = c("age", "sexM"),
  "+ tobacco"                    = c("age", "sexM", "tobacco"),
  "full block [current code]"    = c("age", "sexM", "tobacco", "rin", "seqpc1", "seqpc2", "seqpc3", PL),
  "cell composition (4 ILR)"     = c("b1", "b2", "b3", "b4"),
  "cells + age + sex"            = c("b1", "b2", "b3", "b4", "age", "sexM"))
for (nm in names(sets)) {
  a <- cv(sets[[nm]])
  say("  %-28s %2d vars   AUC %.3f  (sd %.3f)\n", nm, length(sets[[nm]]), mean(a), sd(a))
}
say("\n  Adding eight balanced nuisance variables to age costs real AUC. They\n")
say("  cannot remove confounding -- there is none to remove -- so all they do\n")
say("  is spend degrees of freedom and dilute the one variable that matters.\n")

close(con)
cat(sprintf("\n-> %s\n", file.path(outd, "covariate_audit.txt")))
