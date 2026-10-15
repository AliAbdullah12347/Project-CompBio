# Independent verification of the referee workflow's load-bearing claims.
# Uses only saved CSVs + the cached project object. Writes nothing outside
# network/analyses/015-referee-verification/.
suppressMessages({library(stats)})
OUT <- "network/analyses/015-referee-verification/results"
say <- function(...) cat(sprintf(...), "\n")

me <- read.csv("network/analyses/009-all-modules-across-groups/results/module_eigengenes_all.csv",
               stringsAsFactors = FALSE, check.names = FALSE)
tmm <- read.csv("data/analysis_matrices/TMM_normalisation_factors.csv", stringsAsFactors = FALSE)
say("eigengene table: %d x %d", nrow(me), ncol(me))
say("tmm table cols: %s", paste(names(tmm), collapse=", "))

idc  <- intersect(c("sample","sample_id","Sample","id"), names(me))[1]
idt  <- intersect(c("sample","sample_id","Sample","id"), names(tmm))[1]
say("join keys: me$%s  tmm$%s", idc, idt)
mcols <- grep("^ME", names(me), value = TRUE)
say("ME columns: %d", length(mcols))

m <- merge(me, tmm, by.x = idc, by.y = idt)
say("merged n = %d", nrow(m))

# --- CLAIM 1: PC1 of the ME matrix is a global scaling axis -------------------
X  <- as.matrix(m[, mcols]); X <- scale(X, center = TRUE, scale = TRUE)
pc <- prcomp(X, center = FALSE, scale. = FALSE)
ve <- pc$sdev^2 / sum(pc$sdev^2)
say("PC1 var explained = %.1f%%   PC2 = %.1f%%   PC3 = %.1f%%", 100*ve[1], 100*ve[2], 100*ve[3])
ld <- pc$rotation[,1]
say("PC1 loadings: range [%.3f, %.3f]  all same sign = %s",
    min(ld), max(ld), all(sign(ld)==sign(ld[1])))

tcol <- grep("tmm|factor|norm", names(m), value = TRUE, ignore.case = TRUE)
say("candidate TMM columns: %s", paste(tcol, collapse=", "))
pc1 <- pc$x[,1]
for (tc in tcol) {
  if (!is.numeric(m[[tc]])) next
  say("  r(PC1, %s) = %+.3f", tc, cor(pc1, m[[tc]], use="complete.obs"))
}

gcol <- intersect(c("group","grp","dx","diagnosis"), names(m))[1]
if (!is.na(gcol)) {
  say("group column: %s   levels: %s", gcol, paste(unique(m[[gcol]]), collapse=", "))
  dis <- as.numeric(m[[gcol]] != "control")
  say("  r(PC1, disease[0/1]) = %+.3f", cor(pc1, dis))
  for (tc in tcol) if (is.numeric(m[[tc]])) {
    say("  group means of %s:", tc)
    print(round(tapply(m[[tc]], m[[gcol]], mean), 4))
  }
  # partial out the TMM factor
  tcn <- tcol[sapply(tcol, function(z) is.numeric(m[[z]]))][1]
  if (!is.na(tcn)) {
    rp <- residuals(lm(pc1 ~ m[[tcn]]))
    rd <- residuals(lm(dis ~ m[[tcn]]))
    say("  r(PC1, disease | %s) = %+.3f   (was %+.3f)", tcn, cor(rp, rd), cor(pc1, dis))
  }
}
saveRDS(list(ve=ve, loadings=ld), file.path(OUT,"pc1_check.rds"))
write.csv(data.frame(sample=m[[idc]], PC1=pc1), file.path(OUT,"pc1_scores.csv"), row.names=FALSE)
say("--- claim 1 done ---")
