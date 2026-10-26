# THE matched comparison: one reference, two test groups, one instrument.
#
# base-002-ceil scores the ctrl160 network against:
#   ctrl_ceiling  74 held-out CONTROLS   500 perms   <- the null ceiling
#   nolith        74 bp_nolith           500 perms   <- the observed difference
#   lith_draw_*   74 bp_lith x 20         50 perms   <- different instrument
#
# The first two are the arm's central question. Everything this session learned
# is enforced as a precondition rather than assumed:
#   - same reference? (module sets must be identical, not just same-named)
#   - same instrument? (nperms and n_test must match)
#   - is Zsummary behaving, or is it just module size again?
# It refuses the comparison rather than producing a confident wrong number.
OUT <- "network/analyses/018-matched-ceiling-comparison/results"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf(...), "\n"); sd_ <- stats::sd
RUN <- "network/runs/base-002-ceil"

need <- c(ceiling = file.path(RUN, "ceiling_ctrl_preservation.rds"),
          nolith  = file.path(RUN, "nolith_preservation.rds"))
miss <- need[!file.exists(need)]
if (length(miss)) { say("not ready: %s", paste(basename(miss), collapse=", ")); quit(save="no", status=1) }

C <- readRDS(need[["ceiling"]]); N <- readRDS(need[["nolith"]])

say("=== PRECONDITION 1: same reference? ===")
say("  ceiling modules %d   nolith modules %d", nrow(C), nrow(N))
if (!setequal(C$module, N$module)) {
  say("  REFUSED: module sets differ. WGCNA names modules by size rank within a run,")
  say("  so these are different gene sets and Zsummary is not on one scale.")
  quit(save="no", status=1)
}
say("  identical module sets: %d shared. OK", length(intersect(C$module, N$module)))

say("")
say("=== PRECONDITION 2: same instrument? ===")
say("  ceiling nperms=%s n_test=%s | nolith nperms=%s n_test=%s",
    paste(unique(C$nperms),collapse=","), paste(unique(C$n_test),collapse=","),
    paste(unique(N$nperms),collapse=","), paste(unique(N$n_test),collapse=","))
if (!identical(unique(C$nperms), unique(N$nperms)) || !identical(unique(C$n_test), unique(N$n_test))) {
  say("  REFUSED: instruments differ; Zsummary scales with both."); quit(save="no", status=1)
}
say("  matched. OK")

m <- merge(C[,c("module","Zsummary","medianRank")], N[,c("module","Zsummary","medianRank")],
           by="module", suffixes=c("_ceil","_nolith"))
gC <- m$Zsummary_ceil[m$module=="gold"]; gN <- m$Zsummary_nolith[m$module=="gold"]
r  <- m[!m$module %in% c("gold","grey"), ]

say("")
say("=== SANITY: is Zsummary just module size again? ===")
szf <- file.path(RUN, "reference_module_sizes.csv")
if (file.exists(szf)) {
  sz <- read.csv(szf); nmc <- names(sz)[1]; nsc <- names(sz)[sapply(sz, is.numeric)][1]
  x <- merge(r, sz[,c(nmc,nsc)], by.x="module", by.y=nmc)
  cc <- cor(log10(x[[nsc]]), x$Zsummary_ceil); cn <- cor(log10(x[[nsc]]), x$Zsummary_nolith)
  say("  cor(log10 size, Z): ceiling %+.3f   nolith %+.3f", cc, cn)
  if (abs(cc - cn) > 0.4)
    say("  WARNING: size dependence differs sharply between arms. On 2026-10-25 that")
    say("  signature (+0.911 vs -0.036) meant one arm was computed on raw expression.")
}

say("")
say("=== GOLD (the within-run anchor) ===")
say("  ceiling %.2f (rank %g)   nolith %.2f (rank %g)   of %d",
    gC, m$medianRank_ceil[m$module=="gold"], gN, m$medianRank_nolith[m$module=="gold"], nrow(m))

say("")
say("=== THE QUESTION: do bp_nolith networks differ from held-out controls? ===")
say("  %d real modules, paired over identical gene sets", nrow(r))
say("")
wz <- suppressWarnings(wilcox.test(r$Zsummary_ceil, r$Zsummary_nolith, paired=TRUE))
wr <- suppressWarnings(wilcox.test(r$medianRank_ceil, r$medianRank_nolith, paired=TRUE))
say("  Zsummary   median %.2f -> %.2f   lower in nolith: %d/%d   paired Wilcoxon p = %.4g",
    median(r$Zsummary_ceil), median(r$Zsummary_nolith),
    sum(r$Zsummary_nolith < r$Zsummary_ceil), nrow(r), wz$p.value)
say("  medianRank median %.1f -> %.1f                       paired Wilcoxon p = %.4g",
    median(r$medianRank_ceil), median(r$medianRank_nolith), wr$p.value)
say("")
say("  medianRank is the statistic to believe: ordinal within a run, so immune to the")
say("  size inflation that makes gold score 21.74 while ranking last.")

r$d_rank <- r$medianRank_nolith - r$medianRank_ceil
r <- r[order(-r$d_rank), ]
say("")
say("  most degraded by rank (positive = worse in bp_nolith):")
print(head(r[,c("module","medianRank_ceil","medianRank_nolith","d_rank")], 8), row.names=FALSE)
write.csv(r, file.path(OUT,"ceiling_vs_nolith_matched.csv"), row.names=FALSE)
say("")
say("ANALYSIS_018_DONE")
