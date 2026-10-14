#!/usr/bin/env Rscript
# ===========================================================================
# 011-module-annotation  — covers catalogue items 15, 16 and 11 at once
#
#  15. module overlap with the lithium / bipolar DEG lists and LM22 markers
#      (Fisher exact, BH across modules)
#  16. sample-level structure: do samples separate by group in eigengene space?
#  11. hub genes per module by kME
#
# All at the usable decomposition (deepSplit=4, mergeCutHeight=0, 33 modules).
# Uses cached labels + cached project data; no network rebuild. Writes only here.
# ===========================================================================
setwd("/Users/Aaylab/Project-CompBio")
suppressPackageStartupMessages(library(WGCNA)); options(stringsAsFactors=FALSE); enableWGCNAThreads()
OUT <- "network/analyses/011-module-annotation/results"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), paste0(...,collapse="")))
set.seed(20261014L)
source("network/00_load.R"); d <- load_project(cache=TRUE, verbose=FALSE)
lab <- readRDS("network/runs/mod-000-resolution/labels_p12_ds4.rds"); colors <- labels2colors(lab)
genes <- rownames(d$logtpm); sym <- d$sym[genes]
mods <- setdiff(unique(colors), "grey")

# ---------- 15. module vs DEG lists (Fisher exact, BH across modules) ------
de <- read.csv("de_final/results/significant_genes.csv")
de <- de[de$significant %in% c(TRUE,"TRUE",1), ]
cmps <- unique(de$comparison)
msg("DE comparisons: ", paste(cmps, collapse=" | "))
sets <- lapply(cmps, function(cp) unique(de$symbol[de$comparison==cp & !is.na(de$symbol)]))
names(sets) <- cmps
sets[["LM22_markers"]] <- d$lm22

ann <- do.call(rbind, lapply(names(sets), function(sn) {
  tgt <- !is.na(sym) & sym %in% sets[[sn]]
  if (sum(tgt) < 5) return(NULL)
  r <- do.call(rbind, lapply(mods, function(m) {
    inm <- colors == m
    tb <- matrix(c(sum(inm & tgt), sum(inm & !tgt), sum(!inm & tgt), sum(!inm & !tgt)), 2)
    ft <- fisher.test(tb)
    data.frame(set=sn, set_size=sum(tgt), module=m, n_genes=sum(inm),
               overlap=tb[1,1], expected=round(sum(inm)*sum(tgt)/length(tgt),1),
               odds_ratio=round(unname(ft$estimate),2), p=ft$p.value)
  }))
  r$p_BH <- p.adjust(r$p, "BH"); r[order(r$p_BH), ]
}))
write.csv(ann, file.path(OUT,"module_vs_degs_fisher.csv"), row.names=FALSE)
cat("\n=============== 15. MODULE / DEG OVERLAP (Fisher, BH across modules) ===============\n")
for (sn in unique(ann$set)) {
  s <- ann[ann$set==sn, ]; sig <- s[s$p_BH < 0.05, ]
  cat(sprintf("\n  %s (%d genes mapped): %d of %d modules enriched at BH<0.05\n", sn, s$set_size[1], nrow(sig), nrow(s)))
  if (nrow(sig)) print(head(sig[,c("module","n_genes","overlap","expected","odds_ratio","p_BH")], 6), row.names=FALSE)
}

# ---------- 16. sample structure in eigengene space -------------------------
keep <- !is.na(d$groups)
E <- t(d$logtpm[, keep]); grp <- droplevels(d$groups[keep])
me <- moduleEigengenes(E, colors=colors, excludeGrey=TRUE)$eigengenes
pc <- prcomp(me, scale.=TRUE)
cat("\n=============== 16. SAMPLE STRUCTURE IN EIGENGENE SPACE ===============\n")
sep <- do.call(rbind, lapply(1:4, function(i) {
  a <- summary(aov(pc$x[,i] ~ grp))[[1]]
  data.frame(PC=paste0("PC",i), var_pct=round(100*pc$sdev[i]^2/sum(pc$sdev^2),1),
             F=round(a[1,"F value"],2), p=signif(a[1,"Pr(>F)"],3))
}))
sep$p_BH <- signif(p.adjust(sep$p,"BH"),3)
print(sep, row.names=FALSE)
cat(sprintf("  PCs separating the three groups at BH<0.05: %d of 4\n", sum(sep$p_BH<0.05)))
# how well can group be told apart at all?
cen <- aggregate(pc$x[,1:3], by=list(group=grp), FUN=mean)
cat("\n  group centroids (PC1-3):\n"); print(cen, row.names=FALSE, digits=3)
write.csv(sep, file.path(OUT,"sample_structure_anova.csv"), row.names=FALSE)

# ---------- 11. hub genes by kME -------------------------------------------
ctrl <- !is.na(d$groups) & d$groups=="control"
Ec <- t(d$logtpm[, ctrl])
mec <- moduleEigengenes(Ec, colors=colors, excludeGrey=TRUE)$eigengenes
kME <- signedKME(Ec, mec)
hubs <- do.call(rbind, lapply(mods, function(m) {
  idx <- which(colors == m); col <- paste0("kME", m)
  if (!col %in% colnames(kME)) return(NULL)
  k <- kME[idx, col]; o <- order(-k)[1:min(5,length(idx))]
  data.frame(module=m, n_genes=length(idx),
             hub_rank=seq_along(o), gene=genes[idx][o],
             symbol=ifelse(is.na(sym[idx][o]), "", sym[idx][o]), kME=round(k[o],3))
}))
write.csv(hubs, file.path(OUT,"hub_genes_kme.csv"), row.names=FALSE)
cat("\n=============== 11. TOP HUB GENE PER MODULE (kME, control reference) ===============\n")
h1 <- hubs[hubs$hub_rank==1, c("module","n_genes","symbol","kME")]
print(head(h1[order(-h1$n_genes),], 12), row.names=FALSE)
saveRDS(list(ann=ann, sep=sep, hubs=hubs), file.path(OUT,"summary.rds"))
cat("\nANALYSIS_005_DONE\n")
