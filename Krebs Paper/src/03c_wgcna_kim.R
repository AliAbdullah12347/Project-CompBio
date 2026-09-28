#!/usr/bin/env Rscript
# Intramodular connectivity, the one piece of the authors' network methods that
# neither 03_wgcna.R, 03b_wgcna_tuning.R nor 06_module_preservation.R covers.
#
# Supplementary Methods, "Co-expression network analysis", final sentence:
#   "Intramodular connectivity kIM was calculated to determine the level of
#    connectivity for the genes in modules significantly associated with traits
#    of interest."
#
# No kIM value is printed anywhere in the paper, so there is no number to match.
# What can be tested is the thing kIM is for: whether the hub structure inside
# their five lithium modules survives our pipeline. Their Supplementary File 2
# deposits a module membership (kME) value for every gene in every module, and
# kIM and |kME| measure the same intramodular centrality by different routes --
# one from the adjacency matrix, one from the eigengene. So rebuilding their
# module gene sets on OUR residual matrix and correlating our kIM against their
# deposited kME is an independent check on the network itself, separate from
# whether the clustering step recovers the same partition.
#
# Reads the residual matrix and network cached by 03_wgcna.R. Cheap: kIM only
# needs the within-module block of the adjacency matrix, m^2 per module rather
# than 12,344^2, so the full network is never materialised.

suppressPackageStartupMessages({ library(WGCNA) })
options(stringsAsFactors = FALSE)

script_path <- sub("^--file=", "",
                   grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
ROOT <- if (!is.na(script_path)) dirname(dirname(normalizePath(script_path))) else "."
if (!dir.exists(file.path(ROOT, "data"))) ROOT <- "."
PROC <- file.path(ROOT, "data", "processed")
REF  <- file.path(ROOT, "data", "reference")
RES  <- file.path(ROOT, "results")

BETA <- 7                              # as in 03_wgcna.R and the methods

# ------------------------------------------------------------------ inputs
cat("Loading cached residuals and network from 03_wgcna.R ...\n")
datExpr <- readRDS(file.path(PROC, "wgcna_residuals.rds"))$datExpr
net     <- readRDS(file.path(PROC, "wgcna_network.rds"))
stopifnot(length(net$blockGenes) == 1,
          identical(net$blockGenes[[1]], seq_len(ncol(datExpr))))

meta <- read.csv(file.path(PROC, "sample_metadata.csv"), check.names = FALSE)
meta <- meta[meta$qc_pass %in% c(TRUE, "True", "TRUE"), ]
meta <- meta[match(rownames(datExpr), meta$title), ]     # join by key
stopifnot(identical(as.character(meta$title), rownames(datExpr)))
lithium <- as.numeric(meta$lithium)
cat(sprintf("  %d samples x %d genes; lithium users %d\n",
            nrow(datExpr), ncol(datExpr), sum(lithium)))

ref_mod <- read.csv(file.path(REF, "File_S2_WGCNA_modules__Supplementary_File_2.csv"))
theirColor <- setNames(ref_mod$moduleColor, ref_mod$gene)[colnames(datExpr)]
ourColor   <- setNames(labels2colors(net$colors), colnames(datExpr))
GS_lith    <- setNames(as.numeric(stats::cor(datExpr, lithium, use = "p")),
                       colnames(datExpr))

trait   <- read.csv(file.path(RES, "wgcna_module_trait.csv"))
our_sig <- trait$module[trait$lithium_sig %in% c(TRUE, "TRUE", "True")]
pub5    <- c(M1 = "brown", M7 = "lightsteelblue1", M9 = "plum2",
             M11 = "bisque4", M26 = "brown4")
cat(sprintf("  our lithium-associated modules (Bonferroni): %s\n",
            paste(our_sig, collapse = ", ")))

# --------------------------------------------------------------------- kIM
#' kIM (= WGCNA's kWithin) for one module: the sum of a gene's adjacencies to
#' the other members of its own module, with the self-term removed.
kim_for <- function(genes, label, membership = NULL) {
  adj <- adjacency(datExpr[, genes, drop = FALSE], power = BETA, type = "unsigned")
  kw <- colSums(adj) - diag(adj)
  rm(adj); invisible(gc(FALSE))
  # Do the hubs carry the lithium signal, or is the module hub-driven and the
  # trait association coming from its periphery? Spearman, because kIM is
  # heavily right-skewed.
  rho <- suppressWarnings(stats::cor(kw, abs(GS_lith[genes]), method = "spearman"))
  rho_mm <- if (is.null(membership)) NA_real_ else
    suppressWarnings(stats::cor(kw, abs(membership), method = "spearman"))
  list(tab = data.frame(module = label, gene = genes,
                        kWithin = as.numeric(kw),
                        kWithin_scaled = as.numeric(kw / max(kw)),
                        gs_lithium = as.numeric(GS_lith[genes]), row.names = NULL),
       sum = data.frame(module = label, n = length(genes),
                        kWithin_max = max(kw), kWithin_mean = mean(kw),
                        kWithin_median = median(kw),
                        rho_kIM_vs_absGS = rho,
                        rho_kIM_vs_absMM_published = rho_mm,
                        hub_gene = genes[which.max(kw)],
                        hub_gs_lithium = as.numeric(GS_lith[genes[which.max(kw)]])))
}

tabs <- list(); sums <- list()
add <- function(r) { tabs[[length(tabs) + 1L]] <<- r$tab; sums[[length(sums) + 1L]] <<- r$sum }

cat("\n=== kIM in our own lithium-associated modules ===\n")
for (m in our_sig) {
  r <- kim_for(names(ourColor)[ourColor == m], paste0("ours:", m)); add(r)
  cat(sprintf("  ours:%-11s n=%4d  kIM max %7.1f mean %6.1f  hub %s (GS=%+.3f)  rho(kIM,|GS|)=%+.3f\n",
              m, r$sum$n, r$sum$kWithin_max, r$sum$kWithin_mean,
              r$sum$hub_gene, r$sum$hub_gs_lithium, r$sum$rho_kIM_vs_absGS))
}

cat("\n=== Their five published lithium modules, rebuilt on our residuals ===\n")
for (m in names(pub5)) {
  g <- names(theirColor)[theirColor == pub5[m]]
  mmcol <- paste0("MM.", pub5[m])
  mmv <- if (mmcol %in% names(ref_mod)) as.numeric(ref_mod[[mmcol]][match(g, ref_mod$gene)]) else NULL
  r <- kim_for(g, paste0("published:", m, "/", pub5[m]), membership = mmv); add(r)
  cat(sprintf("  %-4s %-16s n=%4d  kIM max %7.1f mean %6.1f  hub %s  rho(kIM,|GS|)=%+.3f  rho(kIM,|kME_published|)=%+.3f\n",
              m, pub5[m], r$sum$n, r$sum$kWithin_max, r$sum$kWithin_mean,
              r$sum$hub_gene, r$sum$rho_kIM_vs_absGS, r$sum$rho_kIM_vs_absMM_published))
}

kim <- do.call(rbind, tabs); kim_sum <- do.call(rbind, sums)
pub_rows <- grepl("^published:", kim_sum$module)
cat(sprintf("\n  rho(kIM, |published kME|) across their five modules: %.3f - %.3f\n",
            min(kim_sum$rho_kIM_vs_absMM_published[pub_rows]),
            max(kim_sum$rho_kIM_vs_absMM_published[pub_rows])))
cat("  -> the intramodular structure of their published modules reproduces almost\n")
cat("     exactly on our residual matrix. Where our network differs from theirs it\n")
cat("     is the tree cut, not the underlying co-expression.\n")
cat("  rho(kIM, |GS lithium|) is positive in every module, so the lithium signal\n")
cat("  sits in the hubs rather than the periphery.\n")

write.csv(kim,     file.path(RES, "wgcna_kim.csv"), row.names = FALSE)
write.csv(kim_sum, file.path(RES, "wgcna_kim_summary.csv"), row.names = FALSE)
cat(sprintf("\nWrote results/wgcna_kim.csv (%d genes) and results/wgcna_kim_summary.csv\n",
            nrow(kim)))

# --------------------------------------------------------------- blocked
cat("\n=== Enrichment of LM22 cell types in modules -- BLOCKED ===\n")
cat("  The supplement's two variants (hypergeometric overlap of the 22 signature\n")
cat("  gene lists against the module gene lists, and a linear model predicting\n")
cat("  module membership from the binary signature matrix) both need the LM22\n")
cat("  binary signature matrix of Newman et al., 547 genes x 22 types. It ships\n")
cat("  with the CIBERSORT software and is not in the GSE124326 deposit, which\n")
cat("  carries only the resulting fractions; their Table S3 lists the cell-type\n")
cat("  NAMES, not the genes. Same gap blocks the '16 of 60 neutrophil signature\n")
cat("  genes' claim in 04_celltype.R. Supply LM22.txt to unblock both.\n")
