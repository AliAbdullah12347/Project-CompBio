#!/usr/bin/env Rscript
# ===========================================================================
# perm-059  The same permutation null as perm-054, on the LITHIUM contrast.
#
# analyses/019 scored BPD and lithium module enrichment against the same
# invalid null -- random size-matched gene sets, which are uncorrelated while
# module genes are co-expressed by construction. perm-054 redid the BPD side
# honestly and the null turned out to sit at ~16 modules, not ~0. Until the
# lithium side gets the same treatment the two are not comparable, and the
# headline contrast ("27 for BPD on 4 DEGs, 21 for lithium on 1,426") rests
# on a null that was wrong for both.
#
# This is not a copy. perm-054.R carries a CONTRASTS map keyed by row id and
# reads the id from its command line, so running THIS file with "perm-059"
# selects the lithium contrast and its own run directory. A duplicated
# implementation would drift from the one it is supposed to be identical to,
# which is the entire point of the comparison.
#
#   Rscript network/experiments/perm-059.R perm-059
# ===========================================================================
source("network/experiments/perm-054.R")
