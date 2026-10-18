# Raw-matrix artifacts, archived 20261018T211754Z

Everything here was computed from `expr_for()` BEFORE the per-sample library
offset was removed (see network/findings/2026-10-18-normalisation-rebuild-required.md).

The offset put 62% of RANDOM gene pairs above r = 0.5 and differed by group
(Welch p = 1.14e-09), so these module assignments, dissTOM and preservation
statistics describe a degenerate graph. Kept for provenance and for any
raw-vs-corrected comparison; **not valid as results**.

The `mod-000-resolution` grid here is the source of the "33 modules,
largest 11.68%" figure quoted throughout the pre-2026-10-18 findings.
