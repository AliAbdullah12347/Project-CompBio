#!/usr/bin/env Rscript
# ===========================================================================
# 011-darkred-annotation
#
# Analysis 009 showed that 31/33 modules differ between controls and
# untreated BP patients after composition adjustment. Darkred (164 genes)
# and one other module are the exceptions: ctrl-nolith bonf=1.0 and
# ctrl-lith bonf=1.0. Why?
#
# This analysis annotates darkred the same way 006 annotated darkturquoise:
# GO over-representation (Fisher exact, one-sided, FDR<0.05, min 10 genes,
# universe = all 12,368 network genes), hub gene kME in control samples,
# and comparison of darkred ME distribution vs darkturquoise.
#
# Pre-specified: same protocol as 006. Control samples used for kME
# (n=234) since darkred is not significantly different between groups,
# controls provide the most stable reference.
# ===========================================================================
setwd("/Users/Aaylab/Project-CompBio")
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(WGCNA)
  library(AnnotationDbi)
  library(GO.db)
})
enableWGCNAThreads()
OUT <- "network/analyses/011-darkred-annotation/results"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), paste0(...,collapse="")))
set.seed(20261014L)

msg("Loading project data ...")
source("network/00_load.R"); d <- load_project(cache = TRUE, verbose = FALSE)

all_genes <- rownames(d$logtpm)   # unversioned Ensembl IDs, n=12368
msg(length(all_genes), " genes in network")

msg("Loading module labels ...")
lab    <- readRDS("network/runs/mod-000-resolution/labels_p12_ds4.rds")
colors <- WGCNA::labels2colors(lab)
names(colors) <- all_genes

dr_genes <- all_genes[colors == "darkred"]
dt_genes <- all_genes[colors == "darkturquoise"]
msg("Darkred: ", length(dr_genes), " genes; darkturquoise: ", length(dt_genes), " genes")

write.csv(data.frame(ensembl_id = dr_genes, module = "darkred"),
          file.path(OUT, "darkred_genes.csv"), row.names = FALSE)

# ---- GO enrichment (same run_go() as 006) -----------------------------------
org_ok <- suppressPackageStartupMessages(
  tryCatch({ library(org.Hs.eg.db); TRUE }, error = function(e) FALSE)
)
msg("org.Hs.eg.db available: ", org_ok)

if (org_ok) {
  eg2go <- suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db, keys = all_genes,
                          columns = c("GOALL", "ONTOLOGYALL"),
                          keytype = "ENSEMBL")
  )
  eg2go <- eg2go[!is.na(eg2go$GOALL), ]
  msg(length(unique(eg2go$ENSEMBL)), " genes with GO annotation")

  run_go <- function(module_genes, ns) {
    sub_df  <- eg2go[eg2go$ONTOLOGYALL == ns, ]
    go_terms <- split(sub_df$ENSEMBL, sub_df$GOALL)
    go_terms <- go_terms[sapply(go_terms, function(g) length(intersect(g, all_genes)) >= 10)]
    if (length(go_terms) == 0) return(NULL)
    n_all <- length(all_genes)
    n_mod <- length(module_genes)
    res <- do.call(rbind, lapply(names(go_terms), function(term) {
      in_term  <- go_terms[[term]]
      m <- matrix(c(
        length(intersect(module_genes, in_term)),
        length(module_genes) - length(intersect(module_genes, in_term)),
        length(intersect(all_genes, in_term)) - length(intersect(module_genes, in_term)),
        n_all - n_mod - length(intersect(all_genes, in_term)) + length(intersect(module_genes, in_term))
      ), 2, 2)
      pval <- fisher.test(m, alternative = "greater")$p.value
      data.frame(go_id = term,
                 n_mod = length(intersect(module_genes, in_term)),
                 n_all = length(intersect(all_genes, in_term)),
                 p_fisher = pval, stringsAsFactors = FALSE)
    }))
    res$fdr <- p.adjust(res$p_fisher, "BH")
    res$term_name <- tryCatch(Term(res$go_id), error = function(e) NA_character_)
    res <- res[order(res$p_fisher), ]
    res
  }

  go_results <- list()
  for (ns in c("BP","CC","MF")) {
    msg("GO ", ns, " enrichment ...")
    r <- run_go(dr_genes, ns)
    if (is.null(r) || nrow(r) == 0) {
      msg("  no terms tested")
      next
    }
    sig <- r[!is.na(r$fdr) & r$fdr < 0.05, ]
    msg("  ", nrow(sig), " terms FDR < 0.05 / ", nrow(r), " tested")
    write.csv(r, file.path(OUT, paste0("go_", tolower(ns), "_all.csv")), row.names = FALSE)
    go_results[[ns]] <- r
    if (nrow(sig) > 0) {
      cat(sprintf("\n======= GO %s (FDR<0.05, top 15) =======\n", ns))
      print(head(sig[, c("go_id","term_name","n_mod","n_all","p_fisher","fdr")], 15),
            row.names = FALSE)
    } else {
      cat(sprintf("\n======= GO %s: 0 significant terms (top 5 by p) =======\n", ns))
      print(head(r[, c("go_id","term_name","n_mod","n_all","p_fisher","fdr")], 5),
            row.names = FALSE)
    }
  }
} else {
  msg("org.Hs.eg.db not available — skipping GO enrichment.")
}

# ---- kME in control samples (n=234) -----------------------------------------
msg("Computing kME for darkred genes in control samples ...")
ctrl_idx <- !is.na(d$groups) & d$groups == "control"
stopifnot(sum(ctrl_idx) == 234L)
E_ctrl <- t(d$logtpm[, ctrl_idx])   # 234 x 12368

ME_dr_ctrl <- WGCNA::moduleEigengenes(E_ctrl, colors = colors,
                                       excludeGrey = TRUE)$eigengenes[, "MEdarkred"]
E_dr_ctrl  <- E_ctrl[, dr_genes, drop = FALSE]
kME_dr <- apply(E_dr_ctrl, 2, function(g) cor(g, ME_dr_ctrl, use = "pairwise.complete.obs"))
kME_df <- data.frame(ensembl_id = dr_genes,
                     kME_darkred = round(kME_dr, 4),
                     stringsAsFactors = FALSE)
kME_df <- kME_df[order(-kME_df$kME_darkred), ]
write.csv(kME_df, file.path(OUT, "darkred_kME.csv"), row.names = FALSE)
cat("\n======= Top 20 hub genes (kME) in darkred (control samples) =======\n")
print(head(kME_df, 20), row.names = FALSE)

# ---- Contrast: darkred vs darkturquoise composition correlations -------------
msg("Checking darkred ME vs composition (all 460 samples) ...")
samp <- read.csv("network/analyses/007-darkturquoise-across-groups/results/sample_data.csv",
                 stringsAsFactors = FALSE)
me_all <- read.csv("network/analyses/009-all-modules-across-groups/results/module_eigengenes_all.csv",
                   stringsAsFactors = FALSE)
stopifnot(nrow(me_all) == 460L)

ME_dr_all <- me_all[["MEdarkred"]]
ME_dt_all <- me_all[["MEdarkturquoise"]]

cat("\n======= Module ME correlations with ILR balances (460 samples) =======\n")
cat("Module\t\tr_b1\t\tr_b2\t\tr_b3\t\tr_b4\n")
for (mod_name in c("darkred","darkturquoise")) {
  me_vec <- if (mod_name == "darkred") ME_dr_all else ME_dt_all
  cat(sprintf("%-15s\t%.3f\t\t%.3f\t\t%.3f\t\t%.3f\n",
              mod_name,
              cor(me_vec, samp$b1), cor(me_vec, samp$b2),
              cor(me_vec, samp$b3), cor(me_vec, samp$b4)))
}

# Save comparison
comp_df <- data.frame(
  module      = c("darkred","darkturquoise"),
  n_genes     = c(length(dr_genes), length(dt_genes)),
  r_b1        = c(cor(ME_dr_all, samp$b1), cor(ME_dt_all, samp$b1)),
  r_b2        = c(cor(ME_dr_all, samp$b2), cor(ME_dt_all, samp$b2)),
  r_b3        = c(cor(ME_dr_all, samp$b3), cor(ME_dt_all, samp$b3)),
  r_b4        = c(cor(ME_dr_all, samp$b4), cor(ME_dt_all, samp$b4))
)
write.csv(comp_df, file.path(OUT, "module_composition_comparison.csv"), row.names = FALSE)

# ---- Group means of darkred ME ----------------------------------------------
msg("Darkred ME by group ...")
grp_means <- do.call(rbind, lapply(unique(samp$group), function(g) {
  s <- ME_dr_all[samp$group == g]
  data.frame(group=g, n=length(s), mean=round(mean(s),6), sd=round(sd(s),4))
}))
cat("\n======= Darkred ME by group =======\n")
print(grp_means[order(grp_means$group), ], row.names = FALSE)
write.csv(grp_means, file.path(OUT, "darkred_me_by_group.csv"), row.names = FALSE)

cat("\nANALYSIS_011_DONE\n")
