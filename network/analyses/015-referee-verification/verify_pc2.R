# Does the HEADLINE (finding 002: PC2 = lithium dimension, mostly composition)
# survive the library-offset attack that destroys PC1?
OUT <- "network/analyses/015-referee-verification/results"
say <- function(...) cat(sprintf(...), "\n")

me  <- read.csv("network/analyses/009-all-modules-across-groups/results/module_eigengenes_all.csv",
                stringsAsFactors=FALSE, check.names=FALSE)
tmm <- read.csv("data/analysis_matrices/TMM_normalisation_factors.csv", stringsAsFactors=FALSE)
cov <- read.csv("data/analysis_matrices/covariates_474.csv", stringsAsFactors=FALSE)
dn  <- readRDS("network/analyses/015-referee-verification/results/density_check.rds")
offs <- dn$offs

m <- merge(me, tmm, by="sample"); m <- merge(m, cov, by.x="sample", by.y="title", suffixes=c("",".c"))
mc <- grep("^ME", names(m), value=TRUE)
X  <- scale(as.matrix(m[, mc])); pc <- prcomp(X, center=FALSE, scale.=FALSE)
ve <- pc$sdev^2/sum(pc$sdev^2)

for (k in 1:3) {
  s <- pc$x[,k]
  r <- cor(s, m$tmm_norm_factor)
  say("PC%d (%.1f%% var):  r(PC%d, TMM) = %+.3f   R2 on TMM = %.3f", k, 100*ve[k], k, r, r^2)
}

# the headline test, restricted to BP1 as 002 did, with and without the offset
bp <- m[m$group %in% c("bp_lith","bp_nolith"), ]
bp$lith <- as.numeric(bp$group == "bp_lith")
Xb <- scale(as.matrix(bp[, mc])); pcb <- prcomp(Xb, center=FALSE, scale.=FALSE)
veb <- pcb$sdev^2/sum(pcb$sdev^2)
say("")
say("--- within BP1 (n=%d) ---", nrow(bp))
bal <- grep("^b[1-4]$", names(bp), value=TRUE)
for (k in 1:2) {
  s <- pcb$x[,k]
  m0 <- summary(lm(s ~ bp$lith))
  f1 <- as.formula(paste("s ~ bp$lith +", paste(paste0("bp$",bal), collapse=" + ")))
  m1 <- summary(lm(f1))
  m2 <- summary(lm(s ~ bp$lith + bp$tmm_norm_factor))
  f3 <- as.formula(paste("s ~ bp$lith + bp$tmm_norm_factor +", paste(paste0("bp$",bal), collapse=" + ")))
  m3 <- summary(lm(f3))
  gp <- function(z) { cf <- coef(z); i <- grep("lith", rownames(cf))[1]; c(cf[i,1], cf[i,4]) }
  say("PC%d (%.1f%% var)  r(PC%d,TMM)=%+.3f", k, 100*veb[k], k, cor(s, bp$tmm_norm_factor))
  say("   marginal            beta=%+.4f p=%.3g", gp(m0)[1], gp(m0)[2])
  say("   + composition       beta=%+.4f p=%.3g", gp(m1)[1], gp(m1)[2])
  say("   + library offset    beta=%+.4f p=%.3g", gp(m2)[1], gp(m2)[2])
  say("   + offset AND comp   beta=%+.4f p=%.3g", gp(m3)[1], gp(m3)[2])
}
