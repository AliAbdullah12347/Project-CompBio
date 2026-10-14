#!/usr/bin/env Rscript
# ===========================================================================
# 006-darkturquoise-annotation
#
# Analysis 005 showed darkturquoise (140 genes) dominates PC2, the
# lithium-composition dimension. What biological processes do these genes
# represent?
#
# Pre-specified: GO over-representation test (Fisher exact, one-sided),
# FDR < 0.05, min 10 annotated genes per term. All three GO namespaces
# (BP, CC, MF). Universe = all 12,368 network genes.
#
# Uses GO.db and AnnotationDbi from the approved package list.
# Gene IDs in the network are Ensembl (unversioned; 00_load.R strips versions).
# ===========================================================================
setwd("/Users/Aaylab/Project-CompBio")
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(WGCNA)
  library(AnnotationDbi)
  library(GO.db)
})
enableWGCNAThreads()
OUT <- "network/analyses/006-darkturquoise-annotation/results"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), paste0(...,collapse="")))
set.seed(20261014L)

msg("Loading project data ...")
source("network/00_load.R"); d <- load_project(cache = TRUE, verbose = FALSE)

# Gene IDs in the network (rows of expression matrix)
all_genes <- rownames(d$logtpm)   # unversioned Ensembl IDs, n=12368
msg(length(all_genes), " genes in network")

# Darkturquoise module membership
msg("Loading module labels ...")
lab    <- readRDS("network/runs/mod-000-resolution/labels_p12_ds4.rds")
colors <- WGCNA::labels2colors(lab)
names(colors) <- all_genes
stopifnot(length(colors) == length(all_genes))

dt_genes <- all_genes[colors == "darkturquoise"]
msg("Darkturquoise: ", length(dt_genes), " genes")

# Save gene list
write.csv(data.frame(ensembl_id = dt_genes, module = "darkturquoise"),
          file.path(OUT, "darkturquoise_genes.csv"), row.names = FALSE)

# GO annotation via org.Hs.eg.db if available, otherwise GO.db only
# Attempt to load org.Hs.eg.db (Bioconductor human annotation package)
org_ok <- suppressPackageStartupMessages(
  tryCatch({ library(org.Hs.eg.db); TRUE }, error = function(e) FALSE)
)
msg("org.Hs.eg.db available: ", org_ok)

if (org_ok) {
  # Map Ensembl -> GO terms
  eg2go <- suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db, keys = all_genes,
                          columns = c("GOALL", "ONTOLOGYALL"),
                          keytype = "ENSEMBL")
  )
  eg2go <- eg2go[!is.na(eg2go$GOALL), ]
  msg(length(unique(eg2go$ENSEMBL)), " genes with GO annotation")

  run_go <- function(ns) {
    sub_df  <- eg2go[eg2go$ONTOLOGYALL == ns, ]
    go_terms <- split(sub_df$ENSEMBL, sub_df$GOALL)
    go_terms <- go_terms[sapply(go_terms, function(g) length(intersect(g, all_genes)) >= 10)]
    if (length(go_terms) == 0) return(NULL)
    n_all <- length(all_genes)
    n_dt  <- length(dt_genes)
    res <- do.call(rbind, lapply(names(go_terms), function(term) {
      in_term  <- go_terms[[term]]
      m <- matrix(c(
        length(intersect(dt_genes, in_term)),
        length(dt_genes) - length(intersect(dt_genes, in_term)),
        length(intersect(all_genes, in_term)) - length(intersect(dt_genes, in_term)),
        n_all - n_dt - length(intersect(all_genes, in_term)) + length(intersect(dt_genes, in_term))
      ), 2, 2)
      pval <- fisher.test(m, alternative = "greater")$p.value
      n_ann_dt <- length(intersect(dt_genes, in_term))
      n_ann_all <- length(intersect(all_genes, in_term))
      data.frame(go_id = term, n_dt = n_ann_dt, n_all = n_ann_all,
                 p_fisher = pval, stringsAsFactors = FALSE)
    }))
    res$fdr <- p.adjust(res$p_fisher, "BH")
    res$term_name <- tryCatch(Term(res$go_id), error = function(e) NA_character_)
    res <- res[order(res$p_fisher), ]
    res
  }

  for (ns in c("BP","CC","MF")) {
    msg("GO ", ns, " enrichment ...")
    r <- run_go(ns)
    if (is.null(r) || nrow(r) == 0) {
      msg("  no terms tested")
      next
    }
    sig <- r[!is.na(r$fdr) & r$fdr < 0.05, ]
    msg("  ", nrow(sig), " terms FDR < 0.05 / ", nrow(r), " tested")
    write.csv(r, file.path(OUT, paste0("go_", tolower(ns), "_all.csv")), row.names = FALSE)
    if (nrow(sig) > 0) {
      cat(sprintf("\n======= GO %s (FDR<0.05, top 15) =======\n", ns))
      print(head(sig[, c("go_id","term_name","n_dt","n_all","p_fisher","fdr")], 15),
            row.names = FALSE)
    }
  }
} else {
  msg("org.Hs.eg.db not available — skipping GO enrichment. Gene list saved.")
  msg("Install: BiocManager::install('org.Hs.eg.db') (requires user approval per arm rules)")
}

# kME (module membership) within darkturquoise
# kME = Pearson correlation of gene expression with module eigengene
msg("Computing kME for darkturquoise genes ...")
bp <- !is.na(d$groups) & d$groups %in% c("bp_nolith","bp_lith")
stopifnot(sum(bp) == 226L)
E   <- t(d$logtpm[, bp])
lab_bp <- readRDS("network/runs/mod-000-resolution/labels_p12_ds4.rds")
col_bp <- WGCNA::labels2colors(lab_bp)
names(col_bp) <- all_genes

ME_dt <- WGCNA::moduleEigengenes(E, colors = col_bp, excludeGrey = TRUE)$eigengenes[, "MEdarkturquoise"]
E_dt  <- E[, dt_genes, drop = FALSE]   # 226 x 140
kME_dt <- apply(E_dt, 2, function(g) cor(g, ME_dt, use = "pairwise.complete.obs"))
kME_df <- data.frame(ensembl_id = dt_genes,
                     kME_darkturquoise = round(kME_dt, 4),
                     stringsAsFactors = FALSE)
kME_df <- kME_df[order(-kME_df$kME_darkturquoise), ]
write.csv(kME_df, file.path(OUT, "darkturquoise_kME.csv"), row.names = FALSE)
cat("\n======= Top 20 hub genes (kME) in darkturquoise =======\n")
print(head(kME_df, 20), row.names = FALSE)

cat("\nANALYSIS_006_DONE\n")
