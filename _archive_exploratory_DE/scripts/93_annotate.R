#!/usr/bin/env Rscript
# ==============================================================================
# 93_annotate.R -- map Ensembl IDs to HGNC symbols and interpret the genes.
#
# A result stated in Ensembl IDs is not interpretable by anyone, including us.
# The mapping comes from the same GENCODE v19 GTF the counts were generated
# against, so it is the annotation the data was actually built with rather than
# a current-release lookup that would silently rename or merge genes.
# ==============================================================================

source("scripts/common.R")
GTF <- file.path(dirname(EXP_HOME), "cibersortx", "data", "gencode.v19.annotation.gtf.gz")

cat("parsing GENCODE v19 gene records ...\n")
con <- gzfile(GTF, "rt")
ids <- character(0); nms <- character(0); tys <- character(0)
repeat {
  l <- readLines(con, n = 100000)
  if (!length(l)) break
  l <- l[!startsWith(l, "#")]
  if (!length(l)) next
  fld <- vapply(strsplit(l, "\t", fixed = TRUE), function(x) x[3], character(1))
  l <- l[fld == "gene"]
  if (!length(l)) next
  ids <- c(ids, sub("\\..*$", "", sub('.*gene_id "([^"]+)".*', "\\1", l)))
  nms <- c(nms, sub('.*gene_name "([^"]+)".*', "\\1", l))
  tys <- c(tys, sub('.*gene_type "([^"]+)".*', "\\1", l))
}
close(con)
ann <- data.frame(gene = ids, symbol = nms, type = tys, stringsAsFactors = FALSE)
ann <- ann[!duplicated(ann$gene), ]
cat(sprintf("  %d gene records\n", nrow(ann)))

d <- read.csv(file.path(EXP_HOME, "runs", "just_composition", "per_gene.csv"),
              stringsAsFactors = FALSE)
i <- match(d$gene, ann$gene)
d$symbol <- ann$symbol[i]; d$type <- ann$type[i]
write.csv(d, file.path(EXP_HOME, "runs", "just_composition", "per_gene.csv"), row.names = FALSE)
cat(sprintf("  annotated %d of %d tested genes\n", sum(!is.na(d$symbol)), nrow(d)))

show <- function(df, cols, n = 20) print(head(df[, cols], n), row.names = FALSE, digits = 3)

cat("\n=== 20 genes LEAST explained by cell composition ===\n")
cat("    large lithium effect, small predicted-from-composition -> direct-effect candidates\n\n")
show(d[order(-abs(d$residual)), ],
     c("symbol", "type", "beta_lithium", "predicted_from_composition", "residual", "fdr_direct"))

cat("\n=== 15 lithium DEGs MOST explained by composition ===\n")
e <- d[d$fdr_marginal < 0.05, ]
show(e[order(-abs(e$predicted_from_composition)), ],
     c("symbol", "type", "beta_lithium", "predicted_from_composition", "residual"), 15)

cat("\n=== top 20 lithium DEGs by marginal significance ===\n")
show(d[order(d$fdr_marginal), ],
     c("symbol", "type", "beta_lithium", "fdr_marginal", "fdr_direct"))

cat("\n=== gene-type make-up ===\n")
tt <- table(d$type); dd <- table(d$type[d$fdr_direct < 0.05])
cmp <- data.frame(type = names(tt), tested = as.integer(tt),
                  direct = as.integer(dd[names(tt)]), stringsAsFactors = FALSE)
cmp$direct[is.na(cmp$direct)] <- 0
cmp$pct_tested <- round(100 * cmp$tested / sum(cmp$tested), 1)
cmp$pct_direct <- round(100 * cmp$direct / sum(cmp$direct), 1)
print(cmp[order(-cmp$tested), ][1:8, ], row.names = FALSE)

# Neutrophil-identity genes. If the lithium signature is compositional, these
# should be up in users and their effect should be well predicted by gamma*delta.
NEUT <- c("FCGR3B", "CSF3R", "CEACAM3", "CXCR2", "S100A8", "S100A9", "S100A12",
          "MMP8", "MMP9", "DEFA4", "LTF", "CAMP", "ELANE", "MPO", "LCN2",
          "ARG1", "OLFM4", "CEACAM8", "BPI", "PGLYRP1")
LYMPH <- c("CD3E", "CD3D", "CD2", "IL7R", "CCR7", "LEF1", "TCF7", "CD27",
           "MS4A1", "CD79A", "CD79B", "KLRD1", "GZMK", "CCL5", "GNLY")
for (set in list(list("neutrophil identity", NEUT), list("lymphocyte identity", LYMPH))) {
  s <- d[!is.na(d$symbol) & d$symbol %in% set[[2]], ]
  if (!nrow(s)) next
  cat(sprintf("\n=== %s genes present (%d of %d) ===\n", set[[1]], nrow(s), length(set[[2]])))
  s <- s[order(-s$beta_lithium), ]
  print(s[, c("symbol", "beta_lithium", "predicted_from_composition", "residual",
              "fdr_marginal", "fdr_direct")], row.names = FALSE, digits = 3)
  cat(sprintf("  mean beta %+.3f | DE at FDR .05: %d | still direct after balance: %d\n",
              mean(s$beta_lithium), sum(s$fdr_marginal < 0.05), sum(s$fdr_direct < 0.05)))
}
