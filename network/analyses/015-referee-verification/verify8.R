# Does finding #8 ("31/33 modules differ control vs bp_nolith") survive the
# project's own standard covariate set, and the per-sample library offset?
OUT <- "network/analyses/015-referee-verification/results"
say <- function(...) cat(sprintf(...), "\n")

me  <- read.csv("network/analyses/009-all-modules-across-groups/results/module_eigengenes_all.csv",
                stringsAsFactors = FALSE, check.names = FALSE)
cov <- read.csv("data/analysis_matrices/covariates_474.csv", stringsAsFactors = FALSE)
tmm <- read.csv("data/analysis_matrices/TMM_normalisation_factors.csv", stringsAsFactors = FALSE)

idc <- "sample"
key <- intersect(c("title","sample","sample_id"), names(cov))[1]
say("covariate join key: %s", key)
m <- merge(me, cov, by.x = idc, by.y = key, suffixes = c("", ".cov"))
m <- merge(m, tmm, by.x = idc, by.y = "sample")
say("merged n = %d", nrow(m))

gcol <- "group"
mcols <- grep("^ME", names(m), value = TRUE)
say("modules: %d", length(mcols))

bal <- grep("^b[1-4]$", names(m), value = TRUE)
say("ILR balances present: %s", paste(bal, collapse=", "))

sub <- m[m[[gcol]] %in% c("control","bp_nolith"), ]
sub[[gcol]] <- factor(sub[[gcol]], levels = c("control","bp_nolith"))
say("control vs bp_nolith: n = %d (%s)", nrow(sub),
    paste(names(table(sub[[gcol]])), table(sub[[gcol]]), collapse=" / "))

fitp <- function(df, mod, extra) {
  f <- as.formula(paste(mod, "~", gcol, if (nzchar(extra)) paste("+", extra) else ""))
  s <- summary(lm(f, data = df))
  cf <- coef(s)
  r <- grep(paste0("^", gcol), rownames(cf))
  if (!length(r)) return(NA_real_)
  cf[r[1], 4]
}

models <- list(
  "009 as published (composition only)"      = paste(bal, collapse=" + "),
  "+ technical (age sex rin plate seqpc1-3)" = paste(c(bal,"age","sex","rin","factor(plate)","seqpc1","seqpc2","seqpc3"), collapse=" + "),
  "+ library offset (tmm_norm_factor)"       = paste(c(bal,"tmm_norm_factor"), collapse=" + "),
  "+ technical AND library offset"           = paste(c(bal,"age","sex","rin","factor(plate)","seqpc1","seqpc2","seqpc3","tmm_norm_factor"), collapse=" + ")
)

res <- list()
for (nm in names(models)) {
  p <- sapply(mcols, function(mm) tryCatch(fitp(sub, mm, models[[nm]]), error=function(e) NA_real_))
  b3 <- pmin(p * 3, 1)   # 009 rule: each module p x 3 pairwise contrasts
  bh <- p.adjust(p, method="BH")                  # project convention: BH across modules
  res[[nm]] <- data.frame(module=mcols, p=p, bonf3=b3, bh=bh)
  say("%-42s  raw<.05: %2d/33   bonf-n3<.05: %2d/33   BH<.05: %2d/33   min BH = %.3f",
      nm, sum(p<0.05,na.rm=TRUE), sum(b3<0.05,na.rm=TRUE), sum(bh<0.05,na.rm=TRUE), min(bh,na.rm=TRUE))
}
saveRDS(res, file.path(OUT,"finding8_recheck.rds"))
write.csv(do.call(rbind, lapply(names(res), function(n) cbind(model=n, res[[n]]))),
          file.path(OUT,"finding8_recheck.csv"), row.names=FALSE)

say("")
say("=== same four models, control vs bp_lith (the other half of #8) ===")
sub2 <- m[m[[gcol]] %in% c("control","bp_lith"), ]
sub2[[gcol]] <- factor(sub2[[gcol]], levels=c("control","bp_lith"))
for (nm in names(models)) {
  p <- sapply(mcols, function(mm) tryCatch(fitp(sub2, mm, models[[nm]]), error=function(e) NA_real_))
  say("%-42s  bonf-n3<.05: %2d/33   BH<.05: %2d/33", nm,
      sum(pmin(p*3,1)<0.05,na.rm=TRUE), sum(p.adjust(p,"BH")<0.05,na.rm=TRUE))
}
