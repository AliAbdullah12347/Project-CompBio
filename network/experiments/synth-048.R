#!/usr/bin/env Rscript
# ===========================================================================
# synth-048  Cross-run synthesis: which conclusions are stable across all variants?
#
# Spec: cross-run synthesis: which conclusions are stable across all
#       variants? write to RESULTS.md
#
# Reads all completed results CSVs from network/runs/*/
# Synthesizes: for each experiment, what was the key Zsummary finding?
# Stability criterion: a conclusion is "stable" if it appears in ≥80% of
# sensitivity/robustness variants that have completed.
#
# Writes a markdown summary section to network/RESULTS.md
# (appends a "## Cross-run synthesis" section with timestamp).
# ===========================================================================

id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(id) || id == "") id <- "synth-048"
stopifnot(id == "synth-048")

now_ <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
msg  <- function(...) message(sprintf("[%s] %s", now_(), paste0(..., collapse = "")))

msg("=== ", id, " starting ===")

run_dir <- file.path("network/runs", id)
cfg_dir <- "network/config"
dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cfg_dir, showWarnings = FALSE, recursive = TRUE)

RESULTS_MD    <- "network/RESULTS.md"
BASE001_LITH  <- "network/runs/base-001/results_lith_alldraws.csv"
BASE001_NOLITH <- "network/runs/base-001/results_nolith.csv"
STABILITY_THRESHOLD <- 0.80   # 80% of completed variants must agree

cfg_path <- file.path(cfg_dir, paste0(id, ".json"))
if (!file.exists(cfg_path)) {
  scl <- function(v) {
    if (is.logical(v)) return(if (isTRUE(v)) "true" else "false")
    if (is.integer(v) || is.numeric(v)) return(as.character(v[1]))
    paste0('"', gsub('"', '\\"', as.character(v[1]), fixed = TRUE), '"')
  }
  jsn <- function(nm, v) paste0('  "', nm, '": ', scl(v))
  sha <- tryCatch(trimws(system("git rev-parse --short HEAD 2>/dev/null", intern=TRUE)[1]),
                  error = function(e) "unknown")
  cfg <- list(id=id, family="synthesis",
    spec="cross-run synthesis: which conclusions are stable across all variants? write to RESULTS.md",
    stability_threshold=STABILITY_THRESHOLD,
    output=RESULTS_MD,
    git_sha=if(length(sha)&&!is.na(sha))sha else "unknown",
    r_version=paste(R.version$major,R.version$minor,sep="."))
  lines <- mapply(jsn, names(cfg), cfg, SIMPLIFY=TRUE)
  writeLines(c("{", paste(lines, collapse=",\n"), "}"), cfg_path)
  msg("config written: ", cfg_path)
}

options(stringsAsFactors = FALSE)
`%||%` <- function(x, y) if (!is.null(x)) x else y

# ---- helper functions --------------------------------------------------------
z_interp <- function(z) {
  if (is.na(z)) return("NA")
  if (z < 2)  return("not preserved (Z<2)")
  if (z < 10) return("weakly preserved (2≤Z<10)")
  return("strongly preserved (Z≥10)")
}
z_bracket <- function(z) {
  if (is.na(z)) return("NA")
  if (z < 2)  return("Z<2")
  if (z < 10) return("2≤Z<10")
  return("Z≥10")
}
safe_read <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(read.csv(path, stringsAsFactors=FALSE), error=function(e) NULL)
}
summarise_lith <- function(df) {
  if (is.null(df)||nrow(df)==0) return(NULL)
  modules <- unique(df$module)
  do.call(rbind, lapply(modules, function(m) {
    s <- df[df$module==m,]
    data.frame(module=m, n_draws=nrow(s),
               Zsummary_mean=mean(s$Zsummary,na.rm=T),
               Zsummary_sd=sd(s$Zsummary,na.rm=T),
               Zsummary_min=min(s$Zsummary,na.rm=T),
               Zsummary_max=max(s$Zsummary,na.rm=T),
               stringsAsFactors=FALSE)
  }))
}

msg("=== Reading experiment results ===")

# ---- base-001 (primary) ------------------------------------------------------
base001_lith   <- safe_read(BASE001_LITH)
base001_nolith <- safe_read(BASE001_NOLITH)
base001_lith_sum   <- summarise_lith(base001_lith)
msg("base-001 lith: ", if(!is.null(base001_lith_sum)) nrow(base001_lith_sum) else "MISSING", " modules")

# ---- robustness experiments: lith Zsummary per experiment --------------------
robustness_exps <- list(
  "pres-030" = "network/runs/pres-030/results_lith_summary.csv",
  "inp-021"  = "network/runs/inp-021/results_lith_alldraws.csv",
  "inp-022"  = "network/runs/inp-022/results_lith_alldraws.csv",
  "inp-023"  = "network/runs/inp-023/results_lith_alldraws.csv",
  "inp-024"  = "network/runs/inp-024/results_lith_alldraws.csv",
  "gene-025" = "network/runs/gene-025/results_lith_alldraws.csv",
  "gene-026" = "network/runs/gene-026/results_lith_alldraws.csv",
  "gene-027" = "network/runs/gene-027/results_lith_alldraws.csv",
  "gene-028" = "network/runs/gene-028/results_lith_alldraws.csv",
  "corr-005" = "network/runs/corr-005/results_lith_alldraws.csv",
  "corr-006" = "network/runs/corr-006/results_lith_alldraws.csv",
  "net-007"  = "network/runs/net-007/results_lith_alldraws.csv",
  "net-008"  = "network/runs/net-008/results_lith_alldraws.csv",
  "pow-010"  = "network/runs/pow-010/results_lith_alldraws.csv",
  "pow-011"  = "network/runs/pow-011/results_lith_alldraws.csv",
  "ref-002"  = "network/runs/ref-002/results_lith_alldraws.csv",
  "ref-003"  = "network/runs/ref-003/results_lith_alldraws.csv"
)
mod_experiments <- list(
  "mod-013"  = "network/runs/mod-013/results_lith_alldraws.csv",
  "mod-014"  = "network/runs/mod-014/results_lith_alldraws.csv",
  "mod-015"  = "network/runs/mod-015/results_lith_alldraws.csv",
  "mod-016"  = "network/runs/mod-016/results_lith_alldraws.csv",
  "mod-017"  = "network/runs/mod-017/results_lith_alldraws.csv",
  "mod-018"  = "network/runs/mod-018/results_lith_alldraws.csv",
  "mod-019"  = "network/runs/mod-019/results_lith_alldraws.csv",
  "mod-020"  = "network/runs/mod-020/results_lith_alldraws.csv"
)

# Read all robustness experiments
rob_results <- lapply(c(robustness_exps, mod_experiments), function(path) {
  df <- safe_read(path)
  if (is.null(df)) return(NULL)
  summarise_lith(df)
})
completed_rob <- names(rob_results)[!vapply(rob_results, is.null, logical(1))]
msg("completed robustness experiments: ", length(completed_rob), "/", length(rob_results))

# ---- nolith results -----------------------------------------------------------
nolith_exps <- list(
  "base-001"=BASE001_NOLITH,
  "pres-030"="network/runs/pres-030/results_nolith.csv"
)
nolith_results <- lapply(nolith_exps, safe_read)
nolith_completed <- names(nolith_results)[!vapply(nolith_results, is.null, logical(1))]

# ---- other results -----------------------------------------------------------
perm031 <- safe_read("network/runs/perm-031/results_permutation_test.csv")
perm032 <- safe_read("network/runs/perm-032/results_permutation_test.csv")
ann041  <- safe_read("network/runs/ann-041/results_deg_overlap.csv")
ann042  <- safe_read("network/runs/ann-042/results_geneset_overlap.csv")
trait033 <- safe_read("network/runs/trait-033/results_lm_traits.csv")
trait034 <- safe_read("network/runs/trait-034/results_lm_traits.csv")
trait035 <- safe_read("network/runs/trait-035/results_lm_traits.csv")
hub036  <- safe_read("network/runs/hub-036/results_kme_stability.csv")
hub037  <- safe_read("network/runs/hub-037/results_delta_kme.csv")
samp043 <- safe_read("network/runs/samp-043/results_sample_clustering.csv")
filt044 <- safe_read("network/runs/filt-044/results_module_survival_summary.csv")
tob045  <- safe_read("network/runs/tob-045/results_sensitivity_summary.csv")
bmind047_lith <- safe_read("network/runs/bmind-047/results_lith_alldraws.csv")

# ---- stability analysis ------------------------------------------------------
# For each module in base-001, collect Zsummary_mean across all completed variants
# that share that module name. A module conclusion is "stable" if ≥80% of variants
# report the same Z-bracket (Z<2, 2≤Z<10, Z≥10).
msg("computing stability ...")
if (!is.null(base001_lith_sum)) {
  modules <- base001_lith_sum$module
  stability_rows <- do.call(rbind, lapply(modules, function(mod) {
    base_z <- base001_lith_sum$Zsummary_mean[base001_lith_sum$module==mod]
    # collect from all completed robustness experiments
    variant_z <- vapply(completed_rob, function(exp_id) {
      s <- rob_results[[exp_id]]
      if (is.null(s)) return(NA_real_)
      row <- s[s$module==mod,]
      if (nrow(row)==0) return(NA_real_)
      row$Zsummary_mean[1]
    }, numeric(1))
    variant_z <- variant_z[!is.na(variant_z)]
    base_bracket <- z_bracket(base_z)
    n_agree <- sum(vapply(variant_z, function(z) z_bracket(z)==base_bracket, logical(1)))
    n_total <- length(variant_z)
    frac_agree <- if (n_total>0) n_agree/n_total else NA_real_
    data.frame(module=mod, n_variants=n_total,
               base001_Zsummary_mean=base_z, base001_bracket=base_bracket,
               n_agree=n_agree, frac_agree=frac_agree,
               stable=!is.na(frac_agree) && frac_agree >= STABILITY_THRESHOLD,
               stringsAsFactors=FALSE)
  }))
  write.csv(stability_rows, file.path(run_dir,"stability_analysis.csv"), row.names=FALSE)
  n_stable <- sum(stability_rows$stable, na.rm=TRUE)
  msg(n_stable, "/", nrow(stability_rows), " modules stable across ≥",
      round(STABILITY_THRESHOLD*100), "% of variants")
} else {
  stability_rows <- NULL
}

# ---- write RESULTS.md section ------------------------------------------------
msg("writing cross-run synthesis to RESULTS.md ...")

mk_line <- function(...) paste0(..., collapse="")
fmt_z <- function(x) if(!is.null(x) && !is.na(x)) sprintf("%.2f", x) else "pending"

build_section <- function() {
  ts   <- now_()
  lines <- character(0)
  lines <- c(lines, "",
             "---",
             "",
             paste0("## Cross-run synthesis (", ts, ")"),
             "",
             paste0("Generated by synth-048. Stability threshold: ",
                    round(STABILITY_THRESHOLD*100), "% of completed variants."),
             "")

  # Completed experiments
  all_exp_ids <- c("base-001", names(robustness_exps), names(mod_experiments),
                   "perm-031","perm-032","ann-041","ann-042",
                   "trait-033","trait-034","trait-035","hub-036","hub-037",
                   "cent-040","comm-038","comm-039","samp-043","filt-044","tob-045",
                   "bmind-046","bmind-047","synth-048")
  completed_flags <- vapply(all_exp_ids, function(exp) {
    if (exp=="base-001") return(!is.null(base001_lith))
    if (exp=="perm-031") return(!is.null(perm031))
    if (exp=="perm-032") return(!is.null(perm032))
    if (exp=="ann-041")  return(!is.null(ann041))
    if (exp=="ann-042")  return(!is.null(ann042))
    if (exp=="filt-044") return(!is.null(filt044))
    if (exp=="tob-045")  return(!is.null(tob045))
    if (exp=="bmind-047") return(!is.null(bmind047_lith))
    if (exp %in% names(rob_results)) return(!is.null(rob_results[[exp]]))
    file.exists(file.path("network/runs", exp))
  }, logical(1))
  n_done  <- sum(completed_flags)
  n_total <- length(all_exp_ids)
  lines <- c(lines,
             paste0("**Progress**: ", n_done, "/", n_total, " experiments with results"),
             "",
             "| Experiment | Status |",
             "|---|---|",
             vapply(all_exp_ids, function(exp)
               sprintf("| %s | %s |", exp, if(completed_flags[exp]) "complete" else "pending"),
               character(1)),
             "")

  # Primary result: base-001
  lines <- c(lines, "### Primary result (base-001)", "")
  if (!is.null(base001_lith_sum)) {
    lines <- c(lines, "bp_lith preservation (mean across 20 draws):", "",
               "| Module | Zsummary mean (±SD) | medianRank mean | Interpretation |",
               "|---|---|---|---|")
    for (i in seq_len(nrow(base001_lith_sum))) {
      r <- base001_lith_sum[i,]
      lines <- c(lines, sprintf("| %s | %.2f ± %.2f | — | %s |",
                                r$module, r$Zsummary_mean, r$Zsummary_sd,
                                z_interp(r$Zsummary_mean)))
    }
    lines <- c(lines, "")
  } else {
    lines <- c(lines, "*base-001 lith results not yet available*", "")
  }
  if (!is.null(base001_nolith)) {
    lines <- c(lines, "bp_nolith preservation (Zsummary):", "")
    bp_no <- base001_nolith[base001_nolith$module!="grey",]
    for (i in seq_len(nrow(bp_no))) {
      lines <- c(lines, sprintf("- %s: Z=%.2f (%s)",
                                bp_no$module[i], bp_no$Zsummary[i],
                                z_interp(bp_no$Zsummary[i])))
    }
    lines <- c(lines, "")
  }

  # Stability
  lines <- c(lines, "### Stability across sensitivity variants", "")
  if (!is.null(stability_rows) && nrow(stability_rows)>0) {
    lines <- c(lines, paste0("Stability threshold: ≥", round(STABILITY_THRESHOLD*100),
                             "% of ", length(completed_rob), " completed variants agree on Z-bracket."),
               "",
               "| Module | Base-001 Z (mean) | Z-bracket | Variants | Agree | Stable |",
               "|---|---|---|---|---|---|")
    for (i in seq_len(nrow(stability_rows))) {
      r <- stability_rows[i,]
      lines <- c(lines, sprintf("| %s | %.2f | %s | %d | %d | %s |",
                                r$module, r$base001_Zsummary_mean, r$base001_bracket,
                                r$n_variants, r$n_agree, if(r$stable) "YES" else "NO"))
    }
    lines <- c(lines, "")
  } else {
    lines <- c(lines, "*Stability analysis unavailable (base-001 pending)*", "")
  }

  # Annotation results
  lines <- c(lines, "### Annotation (ann-041, ann-042)", "")
  if (!is.null(ann041)) {
    sig <- ann041[!is.na(ann041$padj) & ann041$padj<0.05,]
    lines <- c(lines, paste0("ann-041 DEG overlap: ", nrow(sig),
                             " significant module×DEG-list pairs (BH<5%)"))
    for (dl in unique(sig$deg_list)) {
      ms <- sig$module[sig$deg_list==dl]
      if (length(ms)>0) lines <- c(lines, paste0("  - ", dl, ": ", paste(ms,collapse=",")))
    }
    lines <- c(lines, "")
  } else { lines <- c(lines, "*ann-041 pending*", "") }

  if (!is.null(ann042)) {
    sig <- ann042[!is.na(ann042$padj) & ann042$padj<0.05,]
    lines <- c(lines, paste0("ann-042 LM22/CT-DE overlap: ", nrow(sig), " significant pairs (BH<5%)"), "")
  } else { lines <- c(lines, "*ann-042 pending*", "") }

  # Filter robustness (filt-044)
  lines <- c(lines, "### Filter robustness (filt-044)", "")
  if (!is.null(filt044)) {
    for (i in seq_len(nrow(filt044))) {
      r <- filt044[i,]
      lines <- c(lines, sprintf("- %s: %.0f%% survive all 31 filters (worst filter: %.0f%%)",
                                r$module,
                                r$frac_survive_all_filters*100,
                                r$min_frac_any_filter*100))
    }
    lines <- c(lines, "")
  } else { lines <- c(lines, "*filt-044 pending*", "") }

  # Tobacco sensitivity (tob-045)
  lines <- c(lines, "### Tobacco sensitivity (tob-045)", "")
  if (!is.null(tob045)) {
    lines <- c(lines, "Zsummary range across 5 tobacco imputations (sensitivity, 1 draw each):", "")
    for (i in seq_len(nrow(tob045))) {
      r <- tob045[i,]
      lines <- c(lines, sprintf("- %s: Z mean=%.2f, range=%.2f (min=%.2f, max=%.2f)",
                                r$module, r$Zsummary_mean, r$Zsummary_range,
                                r$Zsummary_min, r$Zsummary_max))
    }
    lines <- c(lines, "")
  } else { lines <- c(lines, "*tob-045 pending*", "") }

  # bMIND per-lineage (bmind-047)
  lines <- c(lines, "### Per-lineage networks (bmind-047)", "")
  if (!is.null(bmind047_lith)) {
    by_lin <- split(bmind047_lith, bmind047_lith$lineage)
    for (lin in names(by_lin)) {
      s <- by_lin[[lin]]
      mods <- unique(s$module)
      z_means <- vapply(mods, function(m) mean(s$Zsummary[s$module==m],na.rm=T), numeric(1))
      lines <- c(lines, sprintf("- %s: %d modules; Z range [%.2f, %.2f]",
                                lin, length(mods), min(z_means), max(z_means)))
    }
    lines <- c(lines, "")
  } else { lines <- c(lines, "*bmind-047 pending (requires bmind-046)*", "") }

  # Main conclusions
  lines <- c(lines, "### Conclusions", "")
  lines <- c(lines,
    "1. **Reference network**: 4 modules detected in whole-blood controls (power=12, bicor, signed)",
    "   after mergeCloseModules (cutHeight=0.25). Biologically plausible: cell-type composition",
    "   dominates whole-blood co-expression.",
    "",
    "2. **Scale-free topology**: NOT achieved at any tested power (max signed-R²=0.955 at power=2,",
    "   positive slope). Power=12 used per WGCNA FAQ for signed networks with n>40.",
    "",
    "3. **Zsummary interpretation**: DESCRIPTIVE ONLY. No p-values computed from Zsummary",
    "   (no CDF transformation, no BH across draws). Draws are not independent.",
    "",
    "4. **Stability**: see table above (pending completion of robustness experiments).",
    "")

  lines
}

synthesis_lines <- build_section()

# Write to run_dir CSV for reference
write.csv(stability_rows %||% data.frame(note="base-001 pending"),
          file.path(run_dir, "stability_analysis.csv"), row.names=FALSE)
writeLines(synthesis_lines, file.path(run_dir, "synthesis_section.md"))
msg("synthesis section written to ", file.path(run_dir, "synthesis_section.md"), " (", length(synthesis_lines), " lines)")

# Append to RESULTS.md
if (file.exists(RESULTS_MD)) {
  existing <- readLines(RESULTS_MD)
  # Remove any prior cross-run synthesis section
  synth_start <- grep("^## Cross-run synthesis", existing)
  if (length(synth_start) > 0) {
    cutoff <- synth_start[1] - 2   # remove the --- separator too
    if (cutoff > 0) existing <- existing[seq_len(cutoff)]
    else existing <- character(0)
  }
  writeLines(c(existing, synthesis_lines), RESULTS_MD)
  msg("appended synthesis to RESULTS.md (total ", length(existing)+length(synthesis_lines), " lines)")
} else {
  writeLines(synthesis_lines, RESULTS_MD)
  msg("wrote RESULTS.md (", length(synthesis_lines), " lines)")
}

msg("=== ", id, " COMPLETE ===")
