#!/usr/bin/env Rscript
# trait-053  Why does finding 001 say 0/33 lithium associations survive composition
#            adjustment while grp-051b says 28/69?
#
# The two differ in THREE ways at once, so neither result refutes the other:
#
#   |          | finding 001                   | grp-051b                        |
#   | modules  | mod-000 labels_p12_ds4.rds    | grp-051 control_A, 69 modules   |
#   |          | (the RAW-matrix 33-module set,|                                 |
#   |          |  built before the 10-18 fix)  |                                 |
#   | samples  | BP1 only, n~226               | all labelled, n~460             |
#   | model    | ME ~ li + b1..b4              | ME ~ li + dx + technical + b1..b4|
#
# This holds the module set FIXED at the verified 71-module offset-removed
# reference and crosses the other two factors, so the discrepancy is attributed
# rather than argued about.
#
# Cheap: eigengenes plus lm, no network building. Seconds, not hours.
#
#   Rscript network/experiments/trait-053-composition-reconcile.R trait-053
suppressPackageStartupMessages(library(WGCNA))
id <- commandArgs(TRUE)[1]; if (is.na(id) || id=="") id <- "trait-053"
run_dir <- file.path("network/runs", id); dir.create(run_dir, recursive=TRUE, showWarnings=FALSE)
msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), paste0(...)))

source("network/00_load.R")
d <- load_project(cache=TRUE, verbose=FALSE)

REF <- "network/cache/ref_control_all_bicor_signed_p12_ds4_mch00_offsetrm.rds"
if (!ref_cache_ok(REF)) stop("verified reference missing or unstamped: ", REF)
ref <- readRDS(REF)
msg("reference: ", ref$n_modules, " modules, expr_mode=", ref$expr_mode)

E  <- expr_for(d, group=NULL)                      # offset-removed, ALL samples
stopifnot(ncol(E) == length(ref$moduleColors))
ME <- moduleEigengenes(E, colors=ref$moduleColors, excludeGrey=TRUE)$eigengenes
msg("eigengenes: ", nrow(ME), " x ", ncol(ME))

md  <- d$meta
md$li <- as.integer(d$groups == "bp_lith")
md$dx <- as.integer(d$groups %in% c("bp_nolith","bp_lith"))
ilr <- d$ilr[match(md$title, d$ilr$sample), c("b1","b2","b3","b4"), drop=FALSE]

sample_sets <- list(
  BP1_only = !is.na(d$groups) & d$groups %in% c("bp_lith","bp_nolith"),
  all      = !is.na(d$groups))

# M_001 reproduces finding 001's adjustment; M_051b reproduces grp-051b's.
fit_and_count <- function(keep, model) {
  p <- vapply(colnames(ME), function(m) {
    y <- ME[[m]][keep]
    df <- data.frame(y=y, li=md$li[keep], dx=md$dx[keep], age=md$age[keep],
                     sex=md$sex[keep], rin=md$rin[keep], plate=factor(md$plate[keep]),
                     ilr[keep,])
    f <- switch(model,
      marginal = lm(y ~ li, data=df),
      M_001    = lm(y ~ li + b1 + b2 + b3 + b4, data=df),
      M_051b   = if (all(df$dx == df$dx[1]))   # dx constant inside BP1: drop it
                   lm(y ~ li + age + sex + rin + plate + b1 + b2 + b3 + b4, data=df)
                 else
                   lm(y ~ li + dx + age + sex + rin + plate + b1 + b2 + b3 + b4, data=df))
    s <- summary(f)$coefficients
    if ("li" %in% rownames(s)) s["li", 4] else NA_real_
  }, numeric(1))
  list(p=p, n_bh=sum(p.adjust(p,"BH") < 0.05, na.rm=TRUE))
}

out <- list(); per_mod <- list()
for (ss in names(sample_sets)) {
  keep <- sample_sets[[ss]]
  for (mm in c("marginal","M_001","M_051b")) {
    r <- fit_and_count(keep, mm)
    # Counts alone cannot be followed up. "1 of 71 survives" is the headline,
    # but which one decides whether it is the eosinophil module replicating
    # across decompositions or a different module each time, and the first
    # version of this row did not record it.
    per_mod[[length(per_mod)+1]] <- data.frame(
      sample_set=ss, n=sum(keep), model=mm,
      module=sub("^ME", "", names(r$p)), p=unname(r$p),
      p_bh=unname(p.adjust(r$p, "BH")))
    out[[length(out)+1]] <- data.frame(sample_set=ss, n=sum(keep), model=mm,
                                       n_modules=ncol(ME), n_BH05=r$n_bh)
    msg(sprintf("  %-9s n=%-3d  %-9s -> %d of %d modules at BH<0.05",
                ss, sum(keep), mm, r$n_bh, ncol(ME)))
  }
}
R <- do.call(rbind, out)
write.csv(R, file.path(run_dir, "reconciliation.csv"), row.names=FALSE)

PM <- do.call(rbind, per_mod)
write.csv(PM, file.path(run_dir, "per_module_p.csv"), row.names=FALSE)
surv <- PM[PM$p_bh < 0.05 & PM$model != "marginal", ]
for (ss in unique(surv$sample_set)) for (mm in unique(surv$model)) {
  g <- surv[surv$sample_set == ss & surv$model == mm, ]
  if (nrow(g) && nrow(g) <= 6)
    msg(sprintf("  survivors %-9s %-7s: %s", ss, mm, paste(g$module, collapse=", ")))
}

msg("")
msg("=== READ THIS TABLE AS ===")
msg("  marginal vs M_001  : how much the ILR balances alone remove")
msg("  M_001 vs M_051b    : how much dx + technical change it on the SAME modules")
msg("  BP1_only vs all    : whether the sample set drives it")
msg("Any residual gap against grp-051b's 28/69 is then attributable to the module set")
msg("alone (71 verified reference modules here vs 69 group-specific there).")
msg("=== trait-053 COMPLETE ===")
