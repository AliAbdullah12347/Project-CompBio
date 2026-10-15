# The structural claim: is the reference network a near-complete graph created
# by a per-sample global offset in log2(TPM+1)?
source("network/00_load.R")
OUT <- "network/analyses/015-referee-verification/results"
say <- function(...) cat(sprintf(...), "\n")
set.seed(1)

d <- load_project(cache = TRUE)
E <- expr_for(d, group = "control", transpose = FALSE)   # genes x samples
say("control expression matrix: %d genes x %d samples", nrow(E), ncol(E))

idx <- sample(nrow(E), min(1500, nrow(E)))
S   <- E[idx, , drop = FALSE]

offs <- colMeans(E)                       # per-sample global offset
say("per-sample mean log2TPM: mean %.3f  sd %.3f  range [%.3f, %.3f]",
    mean(offs), sd(offs), min(offs), max(offs))

C1 <- cor(t(S)); diag(C1) <- NA
say("RAW          mean gene-gene r = %+.3f   %% pairs r>0.5 = %.1f%%",
    mean(C1, na.rm=TRUE), 100*mean(C1 > 0.5, na.rm=TRUE))

Sc <- sweep(S, 2, offs, "-")              # remove the offset
C2 <- cor(t(Sc)); diag(C2) <- NA
say("OFFSET-REMOVED mean gene-gene r = %+.3f   %% pairs r>0.5 = %.1f%%",
    mean(C2, na.rm=TRUE), 100*mean(C2 > 0.5, na.rm=TRUE))

tf <- read.csv("data/analysis_matrices/TMM_normalisation_factors.csv", stringsAsFactors=FALSE)
sn <- colnames(E); mt <- match(sn, tf$sample)
ok <- !is.na(mt)
say("r(per-sample offset, tmm_norm_factor) = %+.3f  (n=%d)",
    cor(offs[ok], tf$tmm_norm_factor[mt[ok]]), sum(ok))

grp <- d$groups[match(colnames(d$logtpm), colnames(d$logtpm))]
alloff <- colMeans(d$logtpm)
g <- d$groups
say("per-sample offset by group:")
print(round(tapply(alloff, g, mean), 4))
aw <- oneway.test(alloff ~ factor(g))
say("Welch one-way p = %.3g", aw$p.value)

saveRDS(list(raw_mean=mean(C1,na.rm=TRUE), adj_mean=mean(C2,na.rm=TRUE),
             raw_pct=mean(C1>0.5,na.rm=TRUE), adj_pct=mean(C2>0.5,na.rm=TRUE),
             offs=offs), file.path(OUT,"density_check.rds"))
