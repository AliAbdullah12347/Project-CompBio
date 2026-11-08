set.seed(11)
g <- read.csv("network/runs/grp-051/genes_control_A.csv")
t <- read.csv(gzfile("de_final/results/DE_WB_BPD.csv.gz")); t <- t[t$gene %in% g$gene,]
score <- setNames(-log10(pmax(t$P,1e-300)), t$gene)
mods <- setdiff(unique(g$module),"grey")
sizes <- sapply(mods, function(m) length(intersect(g$gene[g$module==m], names(score))))
sizes <- sizes[sizes>=20]
nullA <- replicate(60, {
  p <- sapply(sizes, function(n){ gin <- sample(names(score), n)
    suppressWarnings(wilcox.test(score[gin], score[setdiff(names(score),gin)], alternative="greater")$p.value) })
  sum(p.adjust(p,"BH") < 0.05) })
cat(sprintf("random size-matched sets (60 reps): median %.1f significant, 95th pct %.1f, max %.1f\n",
            median(nullA), quantile(nullA,.95), max(nullA)))
cat(sprintf("observed for WB_BPD: 27 of %d modules\n", length(sizes)))
cat(sprintf("verdict: %s\n", if (27 > quantile(nullA,.95)) "27 EXCEEDS the size-only null" else "27 is WITHIN the size-only null"))
cat("CAVEAT: random sets are uncorrelated; real module genes are co-expressed, so this\n")
cat("null is too permissive. The honest null permutes SAMPLE labels and recomputes DE.\n")
