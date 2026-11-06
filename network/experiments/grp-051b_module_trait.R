suppressPackageStartupMessages(library(WGCNA))
source("network/00_load.R"); d <- load_project(cache=TRUE, verbose=FALSE)
g <- read.csv("network/runs/grp-051/genes_control_A.csv")
stopifnot(identical(g$gene, rownames(d$logtpm)))
E  <- expr_for(d, group=NULL)
ME <- moduleEigengenes(E, colors=g$module, excludeGrey=TRUE)$eigengenes
md <- d$meta
md$li <- as.integer(d$groups=="bp_lith")          # on lithium
md$dx <- as.integer(d$groups %in% c("bp_nolith","bp_lith"))  # bipolar
ilr <- d$ilr[,c("b1","b2","b3","b4")]
keep <- !is.na(d$groups)
cat(sprintf("n = %d  (li=%d, dx=%d)\n", sum(keep), sum(md$li[keep]), sum(md$dx[keep])))

# JOINT model: lithium and diagnosis mutually adjusted, plus the technical set
# trait-033 used. Without dx in the model, "lithium" is partly just "bipolar",
# because every lithium user is a patient.
base <- data.frame(li=md$li, dx=md$dx, age=md$age, sex=md$sex,
                   rin=md$rin, plate=factor(md$plate))[keep,]
full <- cbind(base, ilr[keep,])
get <- function(fit, t) { s <- summary(fit)$coefficients
  if (t %in% rownames(s)) s[t,c(1,4)] else c(NA,NA) }
rows <- lapply(colnames(ME), function(m) {
  y <- ME[[m]][keep]
  f1 <- lm(y ~ ., data=base); f2 <- lm(y ~ ., data=full)
  f3 <- lm(y ~ ., data=cbind(base[,c("age","sex","rin","plate")], ilr[keep,]))
  data.frame(module=sub("^ME","",m),
             li_b=get(f1,"li")[1],  li_p=get(f1,"li")[2],
             dx_b=get(f1,"dx")[1],  dx_p=get(f1,"dx")[2],
             li_adj_b=get(f2,"li")[1], li_adj_p=get(f2,"li")[2],
             dx_adj_b=get(f2,"dx")[1], dx_adj_p=get(f2,"dx")[2],
             ilr_F_p=anova(f3)["b1","Pr(>F)"])
})
R <- do.call(rbind, rows)
for (p in c("li_p","dx_p","li_adj_p","dx_adj_p","ilr_F_p"))
  R[[sub("_p$","_padj",p)]] <- p.adjust(R[[p]], "BH")
write.csv(R[order(R$li_p),], "network/runs/grp-051/module_trait_joint.csv", row.names=FALSE)
s <- function(k) sum(R[[k]] < 0.05, na.rm=TRUE)
cat(sprintf("\nmodules: %d\n", nrow(R)))
cat(sprintf("  lithium (adj for diagnosis + technical)        : %d\n", s("li_padj")))
cat(sprintf("  diagnosis (adj for lithium + technical)        : %d\n", s("dx_padj")))
cat(sprintf("  composition (ILR b1)                           : %d\n", s("ilr_padj")))
cat(sprintf("  lithium, ALSO adjusted for composition          : %d\n", s("li_adj_padj")))
cat(sprintf("  diagnosis, ALSO adjusted for composition        : %d\n", s("dx_adj_padj")))
li <- R$module[R$li_padj<0.05]; cat(sprintf("\n  lithium-significant modules: %s\n", paste(head(li,15), collapse=", ")))
