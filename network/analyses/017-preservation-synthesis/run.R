# What can and cannot be concluded from base-002-ceil / base-001 / inp-022.
#
# CANNOT: compare Zsummary across these arms. Each built its OWN reference
# network (36, 71 and 44 modules from four distinct ref_*.rds caches), and WGCNA
# names modules by size rank within a run, so matching on colour matches
# unrelated gene sets. An earlier version of this script did exactly that and
# produced p = 1e-4 and p = 3e-11 for comparisons that have no meaning.
#
# CAN: compare each arm's modules against ITS OWN gold, which is a random gene
# set drawn from that same network under that same instrument.
OUT <- "network/analyses/017-preservation-synthesis/results"
say <- function(...) cat(sprintf(...), "\n"); sd_ <- stats::sd

arms <- list(
  "ceiling (ctrl160 ref, test=74 held-out CONTROLS)" = "network/runs/base-002-ceil",
  "observed (ctrl234 ref, test=74 bp_lith)"          = "network/runs/base-001",
  "comp-adj (composition-residualised ref)"          = "network/runs/inp-022")

rows <- list()
for (nm in names(arms)) {
  f <- sort(Sys.glob(file.path(arms[[nm]], "lith_draw_*.rds")))
  if (!length(f)) { say("  %-48s MISSING", nm); next }
  d <- do.call(rbind, lapply(f, readRDS))
  g <- d[grepl("^gold$", d$module, ignore.case = TRUE), ]
  r <- d[!grepl("^gold$|^grey$", d$module, ignore.case = TRUE), ]
  nmod <- length(unique(r$module))
  # fraction of draws in which each module beats gold FROM THE SAME DRAW
  beat <- sapply(split(r, r$module), function(x) {
    gg <- g$Zsummary[match(x$draw, g$draw)]
    mean(x$Zsummary > gg, na.rm = TRUE) })
  rows[[nm]] <- data.frame(
    arm = nm, n_modules = nmod, draws = length(unique(d$draw)),
    gold_Z = mean(g$Zsummary), gold_Z_sd = sd_(g$Zsummary),
    gold_rank_pct = 100 * mean(g$medianRank) / (nmod + 1),
    beats_all = sum(beat == 1), beats_most = sum(beat >= 0.9), beats_few = sum(beat < 0.5))
}
T <- do.call(rbind, rows)

say("=== within-arm only: each arm's modules vs ITS OWN gold ===")
say("")
for (i in seq_len(nrow(T))) with(T[i,], {
  say("%s", arm)
  say("   modules %d over %d draws", n_modules, draws)
  say("   gold Zsummary %.2f (sd %.2f); gold sits at the %.0fth percentile of rank",
      gold_Z, gold_Z_sd, gold_rank_pct)
  say("   beat gold in ALL draws: %d/%d   in >=90%%: %d   in <50%%: %d",
      beats_all, n_modules, beats_most, beats_few)
  say("") })

say("=== what this does and does not settle ===")
say("  SETTLED: in every arm, real modules separate from random genes, and gold")
say("           sits near the bottom of the rank distribution. The decomposition")
say("           is real - which is what finding #5 denied on the raw matrix.")
say("")
say("  NOT SETTLED: whether control-vs-bp_lith preservation is lower than")
say("           control-vs-control. The ceiling arm builds its reference from 160")
say("           controls and the observed arm from 234, so their Zsummary values")
say("           come from different networks and are not on one scale.")
say("")
say("  NEEDED: ONE reference (the ctrl160 network) tested on BOTH the 74 held-out")
say("           controls AND 74 bp_lith subsamples, same perms, same n_test. Then")
say("           the two Zsummary columns ARE comparable and the question is a")
say("           paired test over identical gene sets.")
write.csv(T, file.path(OUT, "within_arm_gold_anchored.csv"), row.names = FALSE)
say("")
say("ANALYSIS_017_DONE")
