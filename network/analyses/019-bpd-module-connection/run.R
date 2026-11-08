# Do any modules connect to BPD diagnosis at all?
#
# grp-051 tested module x DEG overlap for WB_LI and WB_BPD, but WB_BPD was
# silently skipped by `if (length(degs) < 5) next`: whole blood yields only 4
# BPD DEGs at BH<0.05 against 1426 for lithium. A Fisher test on 4 genes is
# meaningless, so this asks the question a way that does not require
# significance: is each module ENRICHED AT THE TOP of the BPD ranking?
#
# Wilcoxon rank-sum on -log10(p) for module genes vs all others. No threshold,
# uses the whole ranking, and is the standard fallback when a DE contrast is
# underpowered. Run identically for lithium so the two are directly comparable.
OUT <- "network/analyses/019-bpd-module-connection/results"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
say <- function(...) cat(sprintf(...), "\n")

g <- read.csv("network/runs/grp-051/genes_control_A.csv")
gcol <- names(g)[grepl("gene|ensembl|id", names(g), ignore.case=TRUE)][1]
ccol <- names(g)[grepl("color|module", names(g), ignore.case=TRUE)][1]
if (is.na(gcol)) gcol <- names(g)[1]; if (is.na(ccol)) ccol <- names(g)[2]
ref <- list(genes = g[[gcol]], colors = g[[ccol]])
say("reference: %d genes, %d modules", length(ref$genes), length(setdiff(unique(ref$colors),"grey")))

res <- list()
for (lab in c("WB_BPD","WB_LI")) {
  t <- read.csv(gzfile(sprintf("de_final/results/DE_%s.csv.gz", lab)))
  gc_ <- intersect(c("gene","Gene","ensembl","ID"), names(t))[1]
  # "P" is the raw p-value here. An earlier version listed only lowercase "p";
  # intersect() is case-sensitive, so it silently fell back to BH, which has
  # 34% unique values and tie groups up to 71 genes - a degraded ranking.
  pc  <- intersect(c("P","P.Value","pvalue","PValue"), names(t))[1]
  if (is.na(pc)) stop("no raw p-value column in ", lab)
  t <- t[t[[gc_]] %in% ref$genes, ]
  score <- setNames(-log10(pmax(t[[pc]], 1e-300)), t[[gc_]])
  say("")
  say("=== %s (ranking on %s, %d genes matched) ===", lab, pc, length(score))
  for (m in setdiff(unique(ref$colors), "grey")) {
    gin <- intersect(ref$genes[ref$colors==m], names(score))
    if (length(gin) < 20) next
    w <- suppressWarnings(wilcox.test(score[gin], score[setdiff(names(score), gin)], alternative="greater"))
    res[[length(res)+1]] <- data.frame(deg_list=lab, module=m, n=length(gin),
      median_in=median(score[gin]), median_out=median(score[setdiff(names(score),gin)]),
      p=w$p.value)
  }
}
R <- do.call(rbind, res)
R$padj <- ave(R$p, R$deg_list, FUN=function(x) p.adjust(x,"BH"))
write.csv(R[order(R$deg_list, R$p),], file.path(OUT,"module_rank_enrichment.csv"), row.names=FALSE)

for (lab in unique(R$deg_list)) {
  s <- R[R$deg_list==lab,]
  say("")
  say("%s: %d of %d modules enriched at BH<0.05", lab, sum(s$padj<0.05), nrow(s))
  s <- s[order(s$p),]
  print(head(s[,c("module","n","median_in","median_out","padj")], 6), row.names=FALSE)
}
say("")
say("ANALYSIS_019_DONE")
