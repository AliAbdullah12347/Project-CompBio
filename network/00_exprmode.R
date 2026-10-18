# Normalisation mode, in its own file so it can be sourced from anywhere without
# the cost of load_project(). Experiment scripts build CACHE_KEY near the top,
# often BEFORE they source 00_load.R, so the tag has to be available early.
#
# EXPR_MODE and EXPR_TAG must move together: the tag goes into every reference
# cache key, so a reference built one way can never be silently reused the other.
#
# Why the offset is removed: log2(TPM+1) keeps a per-sample global offset (TPM
# fixes each library's SUM, but the MEAN OF THE LOG still moves with how
# top-heavy the library is). On this data it put 62% of RANDOM gene pairs above
# r = 0.5 and differed by group (Welch p = 1.14e-09). Removing it took mean
# gene-gene r from 0.535 to 0.002, mean connectivity k from 1030.84 to 80.96,
# and the module count from 33 to 71.
EXPR_MODE <- "offset_removed"          # "offset_removed" (default) | "raw"
EXPR_TAG  <- if (identical(EXPR_MODE, "raw")) "" else "_offsetrm"
