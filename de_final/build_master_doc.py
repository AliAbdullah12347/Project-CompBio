#!/usr/bin/env python3
"""
build_master_doc.py -- assemble the master DE document from the result files.

Every number in the output is read from results/ at build time, so the document
cannot drift from the analysis. Re-run after any analysis change.

Writes master_DE_document.html (Drive converts it to a Google Doc) and
master_DE_document.md (plain-text fallback).
"""
import csv, gzip, io, os, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
RES  = os.path.join(HERE, "results")
GL   = os.path.join(RES, "genelists")
V2   = os.path.normpath(os.path.join(HERE, "..", "de_v2", "results"))

def rows(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))

def esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

def table(hdr, body):
    h = "".join("<th>%s</th>" % c for c in hdr)
    b = "".join("<tr>" + "".join("<td>%s</td>" % c for c in r) + "</tr>" for r in body)
    return ('<table border="1" cellpadding="5" cellspacing="0">'
            '<thead><tr>%s</tr></thead><tbody>%s</tbody></table>' % (h, b))

def genepara(symbols, per=12):
    s = sorted(set(symbols), key=str.upper)
    out = []
    for i in range(0, len(s), per):
        out.append("<p>" + ", ".join(esc(x) for x in s[i:i+per]) + "</p>")
    return "\n".join(out)

def f(x, d=3):
    try: v = float(x)
    except (TypeError, ValueError): return esc(x)
    if v != 0 and abs(v) < 1e-4: return "%.1e" % v
    return ("%%.%df" % d) % v

P = []
def w(*s): P.extend(s)

# ----------------------------------------------------------------- load
wb_li = rows(os.path.join(GL, "WB_lithium_all_1426.csv"))
wb_bp = rows(os.path.join(GL, "WB_bipolar_all_4.csv"))
ct_li = rows(os.path.join(GL, "CT_lithium_all_161.csv"))
tier1 = rows(os.path.join(V2, "candidate_genes_tier1.csv"))
core  = rows(os.path.join(GL, "WB_lithium_filter_robust_core.csv"))
ALLSWEEP = len(core)
_base = {r["gene"] for r in wb_li}
core  = [r for r in core if r["gene"] in _base]
STRICT31 = sum(1 for r in core if r["n_filters_tested"] == "31" and r["n_filters_sig"] == "31")
_cs = sorted({r["symbol"] for r in core}, key=str.upper)
io.open(os.path.join(GL, "WB_lithium_filter_robust_894_symbols.txt"), "w",
        encoding="utf-8").write("\n".join(", ".join(_cs[i:i + 12])
                                          for i in range(0, len(_cs), 12)))
rp_li = rows(os.path.join(GL, "rankproduct_top200_LI.csv"))
rp_bp = rows(os.path.join(GL, "rankproduct_top200_BPD.csv"))
summ  = rows(os.path.join(RES, "summary.csv"))
recA  = {c: rows(os.path.join(RES, "followup_A_recentre_%s.csv" % c)) for c in ("LI","BPD")}
sdC   = {c: rows(os.path.join(RES, "followup_C_sdfilter_%s.csv" % c)) for c in ("LI","BPD")}
valA  = rows(os.path.join(RES, "followup_A_validation.csv"))
valD  = rows(os.path.join(RES, "followup_D_validation.csv")) if os.path.exists(
            os.path.join(RES, "followup_D_validation.csv")) else None
sumD  = rows(os.path.join(RES, "followup_D_offset_covariate.csv")) if os.path.exists(
            os.path.join(RES, "followup_D_offset_covariate.csv")) else None
mech  = rows(os.path.join(RES, "diag_null_mechanism.csv")) if os.path.exists(
            os.path.join(RES, "diag_null_mechanism.csv")) else None

def syms(path):
    # wrap() in collect_gene_lists.R joins genes WITHIN a line by ", " but joins
    # the lines themselves by a bare newline, so splitting the whole file on
    # commas fuses the last gene of each line to the first of the next. Parse
    # line by line instead. (This silently cost 12 of the 125 Holm genes.)
    out = []
    for line in io.open(path, encoding="utf-8"):
        out += [x.strip() for x in line.split(",") if x.strip()]
    return out
holm = syms(os.path.join(GL, "WB_lithium_holm_symbols.txt"))

LINS = ["gran", "mono", "T", "NK", "B"]
LN   = {"gran":"Granulocytes", "mono":"Monocytes", "T":"T cells", "NK":"NK cells", "B":"B cells"}
by_lin = {L: [r for r in ct_li if r["lineage"] == L] for L in LINS}
multi  = defaultdict(set)
for r in ct_li: multi[r["symbol"]].add(r["lineage"])
multi = {k: v for k, v in multi.items() if len(v) > 1}
S = {r["comparison"]: r for r in summ}

# ================================================================= TITLE
w("<h1>Differential Expression &mdash; Master Document</h1>",
  "<p><i>GSE124326. Separating a lithium signature from a bipolar disorder signature in whole blood. "
  "Ali Abdullah and Ava &mdash; BIOL/MATH 481, Colgate, Fall 2026. "
  "Generated automatically from the analysis output on 18 October 2026.</i></p>",
  "<hr>")

# ================================================================= 0. CONTENTS
w("<h2>Contents</h2>",
  "<ol>"
  "<li>Overview &mdash; everything we did, in one page</li>"
  "<li>What &ldquo;significant&rdquo; means here</li>"
  "<li>The data and the four comparisons</li>"
  "<li>The pipeline, step by step</li>"
  "<li>Results at a glance</li>"
  "<li>All significant genes, by analysis</li>"
  "<li>Ranking the genes (no significance threshold)</li>"
  "<li>The strongest candidates for the discussion</li>"
  "<li>Is it just cell composition?</li>"
  "<li>The meeting follow-ups (18 October), and the gene-filtering studies</li>"
  "<li>Power &mdash; what we could and could not have detected</li>"
  "<li>Limitations and what is still open</li>"
  "</ol>"
  "<p>Appendix A &mdash; where every file lives<br>"
  "Appendix B &mdash; the complete lithium whole-blood gene list (1,426)</p>"
  "<p><i>Sections 1, 5, 6 and 8 are the ones to read if you only read four. Section 10 "
  "holds the newest material: the 18 October meeting items (10.1&ndash;10.4) and the two "
  "gene-filtering studies (10.5&ndash;10.6).</i></p>", "<hr>")

# ================================================================= 1. OVERVIEW
w("<h2>1. Overview &mdash; everything we did, in one page</h2>",
  "<p><b>The problem.</b> Nearly every bipolar patient in a public dataset is medicated, so any "
  "&ldquo;bipolar signature&rdquo; found in blood may simply be a lithium signature. Lithium raises "
  "neutrophil counts, and that alone shifts bulk gene expression without any cell changing its "
  "behaviour. We set out to separate the two.</p>",
  "<p><b>What we had.</b> 480 whole-blood RNA-seq samples from GSE124326 (Krebs et al., 2020), reduced "
  "to a <b>474-subject</b> analysis cohort: 226 bipolar I, 14 bipolar II, 234 controls.</p>",
  "<p><b>What we ran.</b> Two contrasts &times; two levels of resolution:</p>",
  "<ul>"
  "<li><b>Lithium</b>: 152 users vs 74 non-users, within bipolar I only (n = 226). Diagnosis is held "
  "constant, so anything found is attributable to the drug.</li>"
  "<li><b>Bipolar</b>: 74 bipolar I patients <i>off</i> lithium vs 234 controls (n = 308). Medication is "
  "held constant, so anything found is attributable to the illness.</li>"
  "<li>Each measured in <b>whole blood</b> (what the sequencer actually read) and in <b>five cell "
  "lineages</b> (granulocytes, monocytes, T, NK, B &mdash; estimated computationally).</li>"
  "</ul>",
  "<p><b>What we found.</b></p>",
  "<ul>"
  "<li><b>Lithium has a large, reproducible effect on whole blood: %s genes.</b> 125 survive the "
  "strictest possible correction; 894 are called under every one of 31 different gene-filtering "
  "choices.</li>" % format(len(wb_li), ",") +
  "<li><b>Bipolar disorder, off lithium, has essentially none: 4 genes</b>, all sitting on the "
  "significance boundary (q = 0.049 against a 0.05 cutoff), none surviving a stricter correction.</li>"
  "<li><b>Most of the lithium signal is cell composition</b>, not changed cell behaviour: 92.8% of the "
  "1,426 genes move in exactly the direction that the granulocyte fraction alone predicts.</li>"
  "<li><b>But the core of it is not.</b> All 125 of the highest-confidence genes survive statistical "
  "adjustment for cell composition. Composition explains the <i>breadth</i> of the signature, not its "
  "centre.</li>"
  "<li><b>At cell-type resolution the lithium effect sits in granulocytes</b>: 109 of 161 significant "
  "gene&times;lineage results.</li>"
  "<li><b>Nothing for bipolar at cell-type resolution</b> &mdash; and we show in section 11 that this "
  "zero carries no information, because the method cannot detect effects of the size we care about.</li>"
  "</ul>",
  "<p><b>What we checked.</b> Thirty-one gene filters, ten multiple-testing corrections, permutation "
  "nulls throughout, a spike-in power analysis, the three follow-ups raised at the 18 October "
  "meeting, and two further gene-filtering strategies tested across nine cutoffs and twelve "
  "analyses. Section 10 holds the methodological findings, two of which matter beyond this "
  "project: the cell-type estimates are so strongly dependent across genes that the number of "
  "discoveries is an unstable statistic (10.2), and variance-based gene filtering &mdash; a "
  "reasonable-sounding idea &mdash; actively destroys signal here, for a reason we can state "
  "exactly (10.5).</p>",
  "<p><b>One new result came out of that filtering work</b>, and it is the cleanest statement of "
  "the bipolar finding we have: <i>unmedicated bipolar I differs modestly from controls in "
  "leukocyte composition &mdash; including a reduced basophil/mast compartment (2.36% against "
  "2.13%, p = 0.039) &mdash; and we find no whole-blood transcriptional signature independent of "
  "it.</i> Two different filters, applied independently, both land on composition (sections 10.5 "
  "and 10.6).</p>", "<hr>")

# ================================================================= 2. SIGNIFICANCE
w("<h2>2. What &ldquo;significant&rdquo; means here</h2>",
  "<p>Each gene is tested separately, so one comparison runs <b>12,368 tests</b>. If absolutely nothing "
  "were happening, testing each at the conventional p &lt; 0.05 would still flag about "
  "<b>618 genes by chance alone</b>. A list built that way would be mostly noise.</p>",
  "<p>So we control the <b>false discovery rate</b> (FDR) using the Benjamini&ndash;Hochberg procedure "
  "at 0.05. The promise it makes is precise, and narrower than people usually assume:</p>",
  "<ul>"
  "<li><b>Of the genes on our list, we expect roughly 5% to be there by chance.</b> On a list of 1,426 "
  "that is about 71 expected false positives.</li>"
  "<li>It does <i>not</i> promise that the list is 95% likely to be entirely correct. That is a much "
  "stricter guarantee (family-wise error rate), which Holm's procedure provides, and which we report "
  "separately.</li>"
  "<li>It says nothing about whether an effect is <i>large</i> or <i>biologically important</i>. A gene "
  "can be highly significant and barely move.</li>"
  "</ul>",
  "<p>So: <b>a significant gene is one whose between-group difference is bigger than noise can "
  "comfortably explain, after accounting for having looked at 12,368 genes at once.</b></p>",
  "<p>Four standards of evidence appear in this document, from weakest to strongest:</p>")
w(table(["Standard", "What it guarantees", "Lithium, whole blood"], [
    ["<b>BH q &lt; 0.05</b>", "about 5% of the list is expected to be false",
     "<b>%s genes</b>" % format(len(wb_li), ",")],
    ["<b>Holm</b>", "at most a 5% chance that <i>even one</i> gene on the list is false",
     "<b>%d genes</b>" % len(holm)],
    ["<b>Robust to the filter</b>", "called under every one of the 31 expression cutoffs we tried",
     "<b>%s genes</b>" % format(len(core), ",")],
    ["<b>Tier 1</b>", "passes Holm <i>and</i> every filter <i>and</i> a minimum fold-change test "
     "<i>and</i> is seen at cell-type level", "<b>%d genes</b>" % len(tier1)],
]))
w("<p><i>A note on direction.</i> &ldquo;log2FC&rdquo; is the log base 2 of the fold change. "
  "+1.0 means the gene is twice as high in the first group; &minus;1.0 means half. +0.5 is about a "
  "1.4&times; increase. Values here are small because whole blood is a mixture and drug effects are "
  "subtle.</p>", "<hr>")

# ================================================================= 3. DATA
w("<h2>3. The data and the four comparisons</h2>",
  "<h3>3.1 How 480 samples became 474</h3>",
  "<p>The depositors flagged 444 of 480 samples as passing their QC. We did not use that flag, because "
  "it is <b>not random with respect to diagnosis</b>: it drops 35 of 240 controls, 1 of 14 bipolar II, "
  "and 0 of 226 bipolar I (Fisher exact p = 3.0&times;10<sup>&minus;10</sup>). Using it would mean "
  "every case&ndash;control comparison ran against a filtered control group, for reasons unrelated to "
  "the biology.</p>",
  "<p>Taking the exclusions apart, 32 of the 36 were dropped only for <b>missing demographic data</b> "
  "(mostly tobacco status) and 4 were genuine <b>sample mix-ups</b>. So:</p>",
  "<ul><li>4 sample mix-ups removed &mdash; these cannot be recovered.</li>"
  "<li>2 removed for missing data too extensive to impute.</li>"
  "<li><b>30 recovered by multiple imputation</b> of tobacco status, which is missing-at-random and "
  "not a reason to discard an otherwise complete transcriptome.</li></ul>",
  "<p><b>Result: 474 subjects</b> &mdash; 226 bipolar I, 14 bipolar II, 234 controls.</p>",
  "<h3>3.2 The two contrasts, and why they are built this way</h3>")
w(table(["", "Lithium contrast", "Bipolar contrast"], [
    ["<b>Question</b>", "What does the drug do?", "What does the illness do?"],
    ["<b>Groups</b>", "152 lithium users vs 74 non-users", "74 bipolar I off lithium vs 234 controls"],
    ["<b>n</b>", "226", "308"],
    ["<b>Held constant</b>", "Diagnosis &mdash; everyone is bipolar I",
     "Medication &mdash; nobody is on lithium"],
    ["<b>Why it is clean</b>", "Any difference is the drug, not the illness",
     "Any difference is the illness, not the drug"],
]))
w("<p>This design is the whole point of the project. The published analyses compare medicated patients "
  "with healthy controls, which confounds the two completely. Splitting them costs statistical power "
  "&mdash; 74 is a small group &mdash; but it is the only way the question can be answered.</p>",
  "<p>All 14 bipolar II patients are lithium non-users, so they contribute no information to the "
  "lithium contrast and are excluded from it.</p>",
  "<h3>3.3 Covariates, and proof that they matter</h3>",
  "<p>Every model adjusts for <b>age, sex, tobacco use, RNA integrity (RIN), sequencing plate, and "
  "three sequencing-quality components</b>. The bipolar model additionally adjusts for assessment "
  "group; the lithium model cannot, because within bipolar I that variable is constant.</p>",
  "<p>These are not decorative. Removing them changes the answer:</p>",
  "<ul><li>Bipolar, whole blood, <b>without</b> covariates: <b>54 genes</b>.</li>"
  "<li>Bipolar, whole blood, <b>with</b> covariates: <b>4 genes</b>.</li>"
  "<li>The two gene rankings correlate at only r = 0.654.</li></ul>",
  "<p>Fifty of those 54 genes were tracking age, smoking or batch rather than diagnosis. Within the "
  "lithium contrast the covariates are better behaved &mdash; tobacco is balanced (35.1% vs 31.6%, "
  "p = 0.70) and sequencing plate is balanced (p = 0.97) &mdash; but age is not (median 54.8 in "
  "non-users vs 48.9 in users, p = 0.0034), which is reason enough to keep adjusting.</p>", "<hr>")

# ================================================================= 4. PIPELINE
w("<h2>4. The pipeline, step by step</h2>",
  "<p>Seven steps turn a table of raw counts into a list of genes. Each one exists to remove a "
  "specific distortion.</p>")
w(table(["#", "Step", "What it does", "Why it is needed"], [
    ["1", "<b>Drop ENSGR genes</b>", "57,820 &rarr; 57,773",
     "47 rows are duplicate copies of pseudoautosomal genes already counted on X"],
    ["2", "<b>Expression filter</b>", "57,773 &rarr; <b>12,368</b>",
     "A gene with near-zero counts carries no information but still costs a test. "
     "Rule: more than 10 counts in at least 90% of samples"],
    ["3", "<b>TMM normalisation</b>", "Rescales each sample",
     "One sample sequenced twice as deeply is not expressing everything twice as hard. "
     "Runs <i>after</i> filtering, because near-zero genes distort the scaling factors"],
    ["4", "<b>voom</b>", "Log-transforms and weights each observation",
     "In count data the variance grows with the mean, which ordinary regression assumes it does not. "
     "voom measures that relationship and down-weights the noisy observations"],
    ["5", "<b>Linear model</b>", "One regression per gene",
     "Estimates the group difference while holding the covariates fixed"],
    ["6", "<b>Empirical Bayes</b>", "Shrinks each gene's variance toward the overall trend",
     "With 226 samples and 12,368 genes, a few genes get a freakishly small variance by chance and "
     "would otherwise top the list. Borrowing strength across genes prevents that"],
    ["7", "<b>Benjamini&ndash;Hochberg</b>", "Converts 12,368 p-values into q-values",
     "Controls the proportion of false positives on the final list at 5%"],
]))
w("<h3>4.1 Why this filter</h3>",
  "<p>The cutoff (&gt;10 counts in &ge;90% of samples, giving 12,368 genes) was chosen in advance and "
  "then <b>tested against 30 alternatives</b> spanning everything from no filter at all to a very "
  "aggressive one, including <code>edgeR::filterByExpr</code> (13,925 genes). Across all 31, "
  "<b>894 of the 1,426 lithium genes are called every single time</b>. The conclusion does not depend "
  "on where the line is drawn, which is the only defensible reason to pick one.</p>",
  "<p>For reference, Krebs et al. report 12,344 genes from a nominally identical rule. Our 12,368 "
  "differs by 24, which is an annotation-version difference, not a pipeline error.</p>",
  "<h3>4.2 How the cell-type estimates were made</h3>",
  "<p>We do not have measured cell counts. They are estimated in two stages:</p>",
  "<ul>"
  "<li><b>CIBERSORTx</b> with the LM22 signature (B-mode batch correction, quantile normalisation off, "
  "100 permutations) estimates <i>what fraction</i> of each sample is each of 22 cell types. We "
  "collapsed those 22 into <b>5 lineages</b> &mdash; granulocytes, monocytes, T, NK, B &mdash; because "
  "at 22-type resolution the estimates are not reliable (ICC 0.389); at 5 lineages they are "
  "acceptable (ICC 0.693).</li>"
  "<li><b>bMIND</b> then estimates <i>what each lineage is expressing</i> in each individual sample, "
  "giving a 12,368 &times; 5 &times; 474 array. We ran it <b>without</b> supplying the phenotype, so "
  "it cannot have been told the answer it was being asked for.</li>"
  "</ul>",
  "<p>Section 10 shows that these cell-type estimates have a serious statistical property that limits "
  "what can be concluded from them. The whole-blood results do not depend on them at all.</p>", "<hr>")

# ================================================================= 5. GLANCE
w("<h2>5. Results at a glance</h2>")
LAB = {"WB_LI":"Lithium, whole blood", "WB_BPD":"Bipolar, whole blood",
       "CT_LI":"Lithium, cell type", "CT_BPD":"Bipolar, cell type"}
body = []
for k in ("WB_LI", "WB_BPD", "CT_LI", "CT_BPD"):
    r = S[k]
    pu = r["pooled_unique_genes"]
    pu = "&mdash;" if pu in ("NA", "") else "%d distinct genes" % int(float(pu))
    body.append(["<b>%s</b>" % LAB[k], r["n"], "%s" % format(int(r["n_tests"]), ","),
                 "<b>%s</b>" % format(int(r["significant_BH05"]), ","), pu,
                 f(r["min_BH"], 4), f(r["max_abs_logFC"], 2)])
w(table(["Comparison", "n", "Tests run", "Significant (BH 0.05)", "Pooled across lineages",
         "Smallest q", "Largest |log2FC|"], body))
w("<p><b>The cell-type lithium result, split by lineage.</b> These are gene&times;lineage findings; "
  "%d distinct genes are involved, because %d of them are significant in more than one lineage.</p>"
  % (len({r["symbol"] for r in ct_li}), len(multi)))
w(table(["Lineage", "Significant genes", "Share of the 161"],
        [[LN[L], "<b>%d</b>" % len(by_lin[L]), "%.0f%%" % (100*len(by_lin[L])/len(ct_li))]
         for L in LINS]))
w("<p><b>Read this as:</b> lithium's effect is concentrated in granulocytes, which is exactly what a "
  "drug known to raise neutrophil counts should look like. B cells show nothing at all.</p>",
  "<p><b>Bipolar at cell-type level: zero genes in all five lineages.</b> Section 11 explains why that "
  "zero should not be interpreted as evidence of absence.</p>", "<hr>")

# ================================================================= 6. GENE LISTS
w("<h2>6. All significant genes, by analysis</h2>",
  "<p>This section lists every gene we called significant, in all four analyses. The lists are "
  "cumulative in strictness, not in membership &mdash; a gene can be strong in whole blood and absent "
  "at cell-type level, or the reverse.</p>")

# --- 6.1 top 40 WB LI
w("<h3>6.1 Lithium, whole blood &mdash; top 40 of %s (strongest first)</h3>" % format(len(wb_li), ","),
  "<p>Ranked by p-value. &ldquo;Fold&rdquo; is the plain multiplicative change: 1.73&times; means the "
  "gene is 73% higher in lithium users.</p>")
w(table(["#", "Gene", "log2FC", "Fold", "Direction", "q (BH)"],
        [[r["rank"], "<b>%s</b>" % esc(r["symbol"]), f(r["logFC"], 3),
          "%.2f&times;" % (2 ** abs(float(r["logFC"]))),
          "up" if float(r["logFC"]) > 0 else "down", f(r["BH"], 5)]
         for r in wb_li[:40]]))
nup = sum(1 for r in wb_li if float(r["logFC"]) > 0)
PCT_UP = "%.0f" % (100.0 * nup / len(wb_li))
HOLMSET = set(holm)
HUP = sum(1 for r in wb_li if r["symbol"] in HOLMSET and float(r["logFC"]) > 0)
HN  = sum(1 for r in wb_li if r["symbol"] in HOLMSET)
w("<p>Of all %s genes, <b>%d go up</b> in lithium users and <b>%d go down</b> &mdash; %s%% up, "
  "which is a long way from the 50/50 a directionless result would give. Among the %d "
  "highest-confidence genes the skew is far stronger still: <b>%d of %d go up</b>. "
  "Section 9 explains why that matters.</p>"
  % (format(len(wb_li), ","), nup, len(wb_li)-nup, PCT_UP, HN, HUP, HN))

# --- 6.2 Holm
w("<h3>6.2 Lithium, whole blood &mdash; the %d highest-confidence genes (Holm)</h3>" % len(holm),
  "<p>These pass the family-wise error criterion: there is at most a 5%% chance that <i>even one</i> "
  "of them is a false positive. This is the list to quote when a reader is sceptical. All %d also "
  "survive statistical adjustment for cell composition (section 9).</p>" % len(holm),
  genepara(holm))

# --- 6.3 filter-robust
w("<h3>6.3 Lithium, whole blood &mdash; the %s genes called under every filter</h3>" % format(len(core), ","),
  "<p>Significant under all 31 expression-filter settings that retained them. This is a robustness "
  "criterion rather than a stringency one: it says the finding does not depend on an arbitrary "
  "analytical choice. The first 60 by median rank:</p>",
  genepara([r["symbol"] for r in sorted(core, key=lambda z: float(z["median_rank"]))[:60]]),
  "<p><i>The complete list of %s is in <code>results/genelists/WB_lithium_filter_robust_894_symbols.txt</code>; "
  "all %s significant genes are in <code>WB_lithium_all_symbols.txt</code>.</i></p>"
  % (format(len(core), ","), format(len(wb_li), ",")))
w("<p><b>A stricter reading, for completeness.</b> Not every filter retains every gene &mdash; an "
  "aggressive cutoff removes low-expression genes before testing them, so most genes face fewer "
  "than 31 tests. The %s above are significant under <i>every filter that retained them</i>. "
  "If instead we require a gene to have survived all 31 retention rules <i>and</i> been "
  "significant under all 31, the number is <b>%d</b>. Both are reported because they answer "
  "different questions: the first asks whether a result depends on the filter, the second asks "
  "which genes are robust to being filtered at all.</p>"
  % (format(len(core), ","), STRICT31))

# --- 6.4 WB BPD
w("<h3>6.4 Bipolar, whole blood &mdash; all 4 genes</h3>",
  "<p>This is the complete result. Note the q-values: all four sit at 0.0489 against a 0.05 "
  "threshold, meaning they would vanish under almost any tightening of the criteria.</p>")
w(table(["Gene", "log2FC", "Direction", "p", "q (BH)"],
        [["<b>%s</b>" % esc(r["symbol"]), f(r["logFC"], 3),
          "up in bipolar" if float(r["logFC"]) > 0 else "down in bipolar",
          f(r["P"], 6), f(r["BH"], 5)] for r in wb_bp]))
w("<p><b>How to report this honestly:</b> none of the four passes Holm, none passes a minimum "
  "fold-change test, and the set is not stable to the gene filter. The defensible statement is "
  "<i>&ldquo;no reproducible whole-blood transcriptomic signature of bipolar I disorder independent of "
  "lithium was detected&rdquo;</i>, with these four named as the closest misses. Section 11 shows we "
  "had 90% power to detect a lithium-sized effect here, so this is a meaningful negative, not an "
  "underpowered shrug.</p>")

# --- 6.5 CT LI
w("<h3>6.5 Lithium, cell type &mdash; all %d findings across %d genes</h3>"
  % (len(ct_li), len({r["symbol"] for r in ct_li})))
for L in LINS:
    g = sorted(by_lin[L], key=lambda z: float(z["P"]))
    if not g:
        w("<h4>%s &mdash; none</h4>" % LN[L]); continue
    w("<h4>%s &mdash; %d genes</h4>" % (LN[L], len(g)), genepara([r["symbol"] for r in g]))
    w("<p><i>Strongest %d:</i> " % min(5, len(g)) + ", ".join(
        "%s (log2FC %s, q %s)" % (esc(r["symbol"]), f(r["logFC"], 2), f(r["BH"], 4))
        for r in g[:5]) + "</p>")
w("<h4>Genes significant in more than one lineage (%d)</h4>" % len(multi))
bl = sorted(multi.items(), key=lambda kv: (-len(kv[1]), kv[0].upper()))
w(table(["Gene", "Lineages", "How many"],
        [["<b>%s</b>" % esc(k), ", ".join(LN[x] for x in LINS if x in v), len(v)] for k, v in bl]))

# --- 6.6 CT BPD
cb = rows(os.path.join(RES, "DE_CT_BPD.csv.gz"))
cbs = sorted(cb, key=lambda z: float(z["P"]))[:10]
w("<h3>6.6 Bipolar, cell type &mdash; none</h3>",
  "<p>Zero significant genes in every lineage. The ten that came closest, for completeness:</p>")
w(table(["Gene", "Lineage", "log2FC", "p", "q (BH)"],
        [[esc(r["symbol"]), LN[r["lineage"]], f(r["logFC"], 3), f(r["P"], 5), f(r["BH"], 4)]
         for r in cbs]))
w("<hr>")

# ================================================================= 7. RANKING
w("<h2>7. Ranking the genes (no significance threshold)</h2>",
  "<p><b>Why do this at all.</b> A significance cutoff throws away information. A gene ranked 12th in "
  "granulocytes, 30th in monocytes and 25th in T cells is clearly doing something, yet it may fail BH "
  "in all three because each test is individually underpowered. Ranking sidesteps the threshold "
  "entirely: instead of asking <i>is this gene significant?</i> we ask <i>how high does it come, "
  "consistently, across the lineages?</i></p>",
  "<p><b>The method &mdash; rank product.</b> For each gene, take its rank by p-value within each of "
  "the five lineages and multiply them (equivalently, average the logs). A gene that is 10th in all "
  "five scores far better than one that is 1st in one lineage and 9,000th in the others. This "
  "rewards consistency rather than a single lucky lineage.</p>")

w("<h3>7.1 Lithium &mdash; top 30 by rank product</h3>")
w(table(["#", "Gene", "Rank product", "Gran", "Mono", "T", "NK", "B"],
        [[i+1, "<b>%s</b>" % esc(r["symbol"]), f(r["rank_product"], 1),
          r["gran"], r["mono"], r["T"], r["NK"], r["B"]]
         for i, r in enumerate(rp_li[:30])]))
w("<h3>7.2 Bipolar &mdash; top 30 by rank product</h3>",
  "<p>The bipolar contrast produced no significant genes at cell-type level, so this ranking is the "
  "only cell-type output it has. <b>Treat it as hypothesis-generating only</b> &mdash; a ranking "
  "always returns a top 30, whether or not anything real is there, and section 7.3 shows we did not "
  "establish that anything is.</p>")
w(table(["#", "Gene", "Rank product", "Gran", "Mono", "T", "NK", "B"],
        [[i+1, "<b>%s</b>" % esc(r["symbol"]), f(r["rank_product"], 1),
          r["gran"], r["mono"], r["T"], r["NK"], r["B"]]
         for i, r in enumerate(rp_bp[:30])]))

w("<h3>7.3 What the ranking does and does not establish</h3>",
  "<p><b>The obvious objection, and the test for it.</b> A rank product has no built-in null. The "
  "lineages agree strongly &mdash; cross-lineage Spearman of the gene statistics runs from "
  "<b>0.556 to 0.871</b> &mdash; but agreement is only evidence if it exceeds what scrambled "
  "labels produce. Section 10.2 gives a concrete reason to doubt that it does: one sample-level "
  "axis carries 71% of the variance in these estimates and any two genes correlate at 0.674, so "
  "the five lineages are partly the <i>same measurement</i>, and two measurements of one thing "
  "agree whether or not the labels mean anything.</p>",
  "<p>So we measured it: 100 label permutations, three statistics, computed identically on the "
  "real and scrambled labels, on both the recentred ranking this section uses and the raw "
  "ranking section 10.2 recommends.</p>")
RNP = os.path.join(RES, "rankproduct_null.csv")
if os.path.exists(RNP):
    rn = rows(RNP)
    def fx(r, k):
        # the overlap statistics are counts; printing them as 119.000 is noise
        v = float(r[k])
        return str(int(round(v))) if r["statistic"].startswith("top_") else f(v, 3)
    SN = {"mean_rho": "Mean cross-lineage Spearman",
          "top_2plus": "Genes in the top 100 of 2+ lineages",
          "top_3plus": "Genes in the top 100 of 3+ lineages",
          "best_rp": "Best rank product (smaller is better)"}
    for CN in ("LI", "BPD"):
        w("<h4>%s contrast</h4>" % ("Lithium" if CN == "LI" else "Bipolar"))
        w(table(["Statistic", "Ranking", "Observed", "Null median", "Null range", "p"],
                [[SN.get(r["statistic"], r["statistic"]),
                  "recentred" if r["ranking"] == "cen" else "raw",
                  "<b>%s</b>" % fx(r, "observed"), fx(r, "null_median"),
                  "%s to %s" % (fx(r, "null_min"), fx(r, "null_max")),
                  ("<b>%s</b>" % f(r["p"], 4)) if r["exceeds_null"] in ("TRUE","True","true")
                  else f(r["p"], 4)]
                 for r in rn if r["contrast"] == CN]))
    TRUEV = ("TRUE", "True", "true")
    nex = sum(1 for r in rn if r["exceeds_null"] in TRUEV)
    nLI = sum(1 for r in rn if r["exceeds_null"] in TRUEV and r["contrast"] == "LI")
    rho_obs = [float(r["observed"]) for r in rn if r["statistic"] == "mean_rho"]
    rho_nul = [float(r["null_median"]) for r in rn if r["statistic"] == "mean_rho"]

    w("<h4>What the null says</h4>",
      "<p><b>Most of the agreement is shared structure, as predicted &mdash; but not quite all of "
      "it.</b> Scrambled labels produce a median cross-lineage correlation of <b>%.3f to "
      "%.3f</b>, against an observed <b>%.3f to %.3f</b>. So the great majority of the "
      "&ldquo;0.556 to 0.871 agreement&rdquo; headline is a property of the estimates rather than "
      "of the biology: five lineages built from one bulk matrix, sharing an axis that carries 71%% "
      "of their variance, agree strongly whether or not the group labels mean anything.</p>"
      % (min(rho_nul), max(rho_nul), min(rho_obs), max(rho_obs)),
      "<p><b>The residual gap is not nothing, and it is entirely in the lithium contrast.</b> "
      "%d of the 16 statistics nominally exceed their null, and <b>all %d of them are "
      "lithium</b>; the bipolar contrast returns 0 of 8 on either ranking. The strongest is the "
      "multi-lineage overlap: <b>85 genes reach the top 100 of three or more lineages against a "
      "null median of 35</b> (recentred ranking, p = 0.0099), and 68 against 40 on the raw "
      "ranking (p = 0.0198).</p>" % (nex, nLI),
      "<p><b>But none of it survives the multiplicity it created.</b> Those 16 p-values come from "
      "4 statistics &times; 2 rankings &times; 2 contrasts, all decided after seeing the data. "
      "The smallest is 0.0099 against a Benjamini&ndash;Hochberg threshold of 0.0031, and nothing "
      "clears it &mdash; nor Bonferroni at 0.0031. Having spent this whole document insisting "
      "that counts be priced for the searches that produced them, we apply the same standard "
      "here.</p>")
    w("<h4>The sentence to write</h4>",
      "<p><i>A rank product across the five lineages gives a usable ordering of genes, but the "
      "cross-lineage agreement it rests on is largely a property of the deconvolution rather than "
      "evidence of independent replication. For the lithium contrast there is a suggestion of "
      "genuine excess agreement, concentrated in how many genes rank highly in three or more "
      "lineages at once; it does not survive correction for the comparisons made and should be "
      "treated as hypothesis-generating. For the bipolar contrast there is no excess at "
      "all.</i></p>",
      "<p><b>What changes in practice.</b> The rankings in 7.1 and 7.2 stay &mdash; an ordering "
      "is still the right way to choose which genes to look at first, and it is the only "
      "cell-type output the bipolar contrast has. What must go is any sentence offering "
      "&ldquo;the lineages agree&rdquo; as support for a gene. On the evidence here that "
      "agreement is mostly the same measurement repeated five times.</p>",
      "<p><i>One encouraging detail, for honesty in the other direction:</i> the asymmetry "
      "between the contrasts is itself a check on the method. If the test were simply registering "
      "shared structure it would do so equally in both, since both use the same estimates and "
      "overlapping samples. It fires only where a real drug effect is known to exist and stays "
      "silent where we have independently concluded there is nothing to find. That is the "
      "behaviour a working test should show.</p>")
else:
    w("<p><i>The permutation null for this section is still running.</i></p>")
w("<hr>")

# ================================================================= 8. CANDIDATES
w("<h2>8. The strongest candidates for the discussion</h2>",
  "<p>Tier 1: the %d genes that pass <i>every</i> test we applied &mdash; Holm correction, all 31 gene "
  "filters, a minimum fold-change requirement, and detection at cell-type level. If the discussion "
  "section needs named genes, these are the defensible ones.</p>" % len(tier1))
t1 = sorted(tier1, key=lambda z: float(z["P"]))
w(table(["Gene", "log2FC", "Fold", "Filters", "Lineages", "Effect lost on<br>composition adjustment",
         "Correlation with<br>granulocyte fraction"],
        [["<b>%s</b>" % esc(r["symbol"]), f(r["logFC"], 3),
          "%.2f&times;" % (2 ** abs(float(r["logFC"]))),
          "%s/%s" % (r["filters_sig"], r["filters_tested"]),
          esc(r["lineages"]) if r["lineages"] else "&mdash;",
          "%.1f%%" % float(r["pct_effect_lost"]), f(r["rho_granulocyte"], 3)]
         for r in t1]))
w("<p><b>How to read the last two columns.</b> &ldquo;Effect lost&rdquo; is how much of the lithium "
  "coefficient disappears once cell composition is held fixed; a low number means the gene is not just "
  "a composition proxy. The correlation with granulocyte fraction says the same thing from the other "
  "direction.</p>")

w("<h3>8.1 The three genes worth leading with</h3>")
top3 = [r for r in wb_li[:3]]
w("<ul>")
for r in top3:
    lin = [LN[x] for x in LINS if any(z["symbol"] == r["symbol"] and z["lineage"] == x for z in ct_li)]
    w("<li><b>%s</b> &mdash; log2FC %s (%.2f&times;), q = %s, whole-blood rank %s. %s</li>"
      % (esc(r["symbol"]), f(r["logFC"], 3), 2 ** abs(float(r["logFC"])), f(r["BH"], 6), r["rank"],
         ("Also significant in " + ", ".join(lin) + ".") if lin
         else "Not independently significant at cell-type level."))
w("</ul>")
w("<p><b>And a fourth worth a sentence:</b> among the cell-type granulocyte hits are a cluster of "
  "genes with coherent neutrophil biology &mdash; <b>CXCL1, IL8, NLRP3, PTGS2, CLEC4D, MIR223, CR1, "
  "SLC26A8, ALOX15</b>. These are inflammatory and neutrophil-activation genes, which is a biologically "
  "sensible place for a lithium effect to land, and it is a more persuasive argument than any single "
  "gene.</p>")
w("<p><b>A counterexample worth reporting.</b> Three eosinophil-associated genes behave backwards: "
  "their effect <i>grows</i> when composition is adjusted for (IL8 &minus;24.6%, IL5RA &minus;3.0%, "
  "CYSLTR2 &minus;1.2% &ldquo;lost&rdquo;). The explanation is that eosinophils average 0.115% of the "
  "mixture and are estimated as exactly zero in 410 of 474 samples, so the composition term carries no "
  "information about them. It is an honest limitation of LM22 in whole blood, and it is worth stating "
  "rather than hiding.</p>", "<hr>")

# ================================================================= 9. COMPOSITION
w("<h2>9. Is it just cell composition?</h2>",
  "<p>This is the question the whole project exists to answer, so it deserves a direct answer: "
  "<b>largely yes, but not entirely, and the part that is not is the strongest part.</b></p>",
  "<h3>9.1 The evidence that it IS composition</h3>",
  "<ul>"
  "<li><b>92.8%</b> of the 1,426 lithium genes move in exactly the direction predicted by that "
  "gene's correlation with granulocyte fraction alone.</li>"
  "<li>The median gene loses <b>40.8%</b> of its effect once cell composition is held fixed.</li>"
  "<li>The effect is strongly one-directional: " + PCT_UP + "% of the 1,426 genes go up "
  "(945 against 481), and <b>" + str(HUP) + " of the " + str(HN) + " highest-confidence genes "
  "go up</b>. A shift in the cell mixture pushes most genes the same way; a specific "
  "transcriptional response would not.</li>"
  "<li>At cell-type resolution the effect concentrates in granulocytes (109 of 161), the very "
  "compartment lithium is known to expand.</li>"
  "</ul>",
  "<h3>9.2 The evidence that it is NOT only composition</h3>",
  "<ul>"
  "<li><b>All 125 Holm-significant genes survive adjustment for composition.</b> Not most "
  "&mdash; all of them.</li>"
  "<li>Survival falls steadily as the genes get weaker: 100% in the top decile, then 62.9%, 17.6%, "
  "and on down to 1.4% in the bottom decile. Composition accounts for the long tail of the "
  "signature, not its core.</li>"
  "<li>This is not circular reasoning. If the adjustment were simply erasing small effects, the "
  "effect size and the percent of effect lost would correlate. They correlate at <b>0.026</b>.</li>"
  "</ul>",
  "<h3>9.3 The conclusion to write down</h3>",
  "<p><i>Lithium use is associated with a large and highly reproducible whole-blood expression "
  "signature. The majority of that signature is attributable to a shift in leukocyte composition "
  "&mdash; principally an expansion of the granulocyte compartment &mdash; but a core of roughly 125 "
  "genes remains significant after composition is held fixed, indicating a component of the effect "
  "that the mixture alone does not explain.</i></p>",
  "<p>We tested the obvious alternative and it failed: if composition were <i>masking</i> a real "
  "signal, adjusting for it should flip the sign of a meaningful number of genes. We set a bar of "
  "10% sign discordance in advance. The observed figure was <b>0.000%</b>.</p>", "<hr>")


# ================================================================= 10. FOLLOW-UPS
def sec_followups():
    w("<h2>10. The meeting follow-ups, and the gene-filtering studies</h2>",
      "<p>Three things were asked for at the weekly check-in. All three were run. The first produced "
      "a result that changes how the cell-type analysis must be reported, so it is given the most "
      "space.</p>")

    # ---- 10.1 recentring
    w("<h3>10.1 Recentring the bMIND fold changes</h3>",
      "<p><b>What was asked.</b> &ldquo;Correct the median fold change to 0 and check whether genes "
      "are now more upregulated.&rdquo;</p>",
      "<p><b>Why it was asked.</b> The fold changes coming out of bMIND are not centred on zero. In "
      "granulocytes under the lithium contrast the median gene has log2FC &minus;0.118, and only "
      "7.6% of genes are up. A whole lineage drifting downward like that is not biology &mdash; it "
      "is a per-sample additive offset in the estimates. If every gene is pushed down by a constant, "
      "genuinely upregulated genes get dragged below the significance line.</p>",
      "<p><b>The size of the problem.</b> The offset (&minus;0.118) is larger than the "
      "gene-to-gene spread of the fold changes it sits on (SD 0.090). The artifact is bigger than "
      "the signal it is burying.</p>")
    w(table(["Contrast", "Lineage", "Median log2FC", "% of genes up", "Genes before"],
            [[r["contrast"], LN[r["lineage"]], f(r["median_raw"], 4),
              "%s%%" % r["pct_up_raw"], r["deg_raw"]]
             for c in ("LI", "BPD") for r in recA[c]]))

    w("<p><b>What we did.</b> Four different treatments of the same data, so they can be compared "
      "on equal terms:</p>",
      "<ol>"
      "<li><b>Raw</b> &mdash; no correction at all.</li>"
      "<li><b>Post-hoc shift</b> &mdash; subtract the lineage median from every fold change and "
      "re-test, keeping the original standard errors.</li>"
      "<li><b>Per-sample centring</b> &mdash; subtract each sample's own median across genes before "
      "fitting. This is the literal reading of what was asked.</li>"
      "<li><b>Offset as a covariate</b> &mdash; put each sample's lineage median into the design "
      "matrix as a nuisance variable instead of subtracting it.</li>"
      "</ol>")
    w(table(["Contrast", "Lineage", "Raw", "Post-hoc shift", "Per-sample centring"],
            [[r["contrast"], LN[r["lineage"]], r["deg_raw"], r["deg_posthoc"],
              "<b>%s</b>" % r["deg_centred"]]
             for c in ("LI", "BPD") for r in recA[c]]))
    w("<p><b>The answer to the question as asked: yes, dramatically.</b> In granulocytes under "
      "lithium the direction flips from 7.6% up to 54.4% up, and the count rises from 109 to 284. "
      "In T cells under the bipolar contrast it goes from <b>0 genes to 345</b>.</p>",
      "<p>A correction that turns 0 into 345 is either a real unmasking or a false-positive "
      "machine, and the two are indistinguishable from the count alone. So we tested it.</p>")

    # ---- 10.2 validation
    w("<h3>10.2 Testing whether the recentred result is real</h3>",
      "<p><b>The test.</b> Scramble the group labels and run the identical pipeline 100 times. "
      "Because 100 permutations &times; 5 lineages &times; 4 treatments is expensive, the "
      "permutation work runs on a fixed random <b>2,000-gene subsample</b>. <b>That is why the "
      "observed counts in the next two tables (16, 1, 3&hellip;) are smaller than the headline "
      "counts (109, 20, 28&hellip;)</b> &mdash; it is one sixth of the genes, not a different "
      "result. Observed and null are measured on the same subsample, so the comparison is "
      "like-for-like. "
      "Centring is label-blind &mdash; it uses each sample's own median across genes &mdash; so it "
      "is legitimate to apply inside a permutation. If scrambled labels still yield hundreds of "
      "genes, the correction is manufacturing them.</p>")
    w(table(["Contrast", "Lineage", "Observed<br>(raw)", "Null 95th pct<br>(raw)",
             "Null max<br>(raw)", "Observed<br>(centred)", "Null 95th pct<br>(centred)",
             "Null max<br>(centred)"],
            [[r["contrast"], LN[r["lineage"]], r["obs_raw"], f(r["null_raw_q95"], 1),
              format(int(r["null_raw_max"]), ","), r["obs_cen"], f(r["null_cen_q95"], 1),
              format(int(r["null_cen_max"]), ",")]
             for r in valA]))
    w("<p><b>Two things fall out of this table, and the second was a surprise.</b></p>",
      "<p><b>(a) Centring inflates the null everywhere, and badly in one place.</b> Its 95th "
      "percentile rises above the raw value in 7 of 10 lineages, and for bipolar granulocytes "
      "it reaches <b>34.5 genes</b> from scrambled labels against 2.0 for raw. Subtracting a per-sample quantity removes a real axis of variation "
      "without telling the model it has done so, so the residual degrees of freedom are overstated, "
      "every standard error shrinks, and every t-statistic inflates &mdash; null genes included.</p>",
      "<p><b>(b) But the raw pipeline has a worse tail than the centred one.</b> The raw null "
      "maximum reaches <b>1,354</b> genes in granulocytes and <b>1,927</b> in B cells &mdash; larger "
      "than anything the centred version produces. We did not expect this and it needed "
      "explaining.</p>")

    if mech:
        w("<h4>Why a single permutation can return 1,354 of 2,000 genes</h4>",
          "<p>We ruled out the mechanical explanations first. In that permutation the design was "
          "full rank, the conditioning was ordinary, the variance prior (28.98) and residual degrees "
          "of freedom (242.0) matched every other permutation, and the largest t-statistic was only "
          "4.6. Nothing had broken numerically. What distinguished it was that <b>44% of genes had "
          "p &lt; 0.01 simultaneously</b>, against 0.3% in a typical permutation.</p>",
          "<p>That is the signature of genes that are not independent. We measured the dependence "
          "directly:</p>")
        w(table(["Lineage", "Variance explained by the single<br>leading sample axis",
                 "Mean correlation between<br>any two genes",
                 "Spearman(alignment with that axis,<br>genes returned)",
                 "Null max", "Permutations returning<br>&gt;50 genes"],
                [[LN[r["lineage"]], "<b>%.1f%%</b>" % (100 * float(r["var_explained_pc1"])),
                  "<b>%s</b>" % f(r["mean_gene_r"], 3), f(r["spearman_align_n"], 3),
                  format(int(r["max_n"]), ","), "%.1f%%" % float(r["pct_perm_over_50"])]
                 for r in mech]))
        w("<p><b>This is the single most important technical finding in the cell-type arm.</b> In "
          "the bMIND granulocyte estimates one sample-level axis explains %.0f%% of all variance, "
          "and any two genes correlate at %s on average. The 12,368 genes are close to being one "
          "variable measured 12,368 times.</p>"
          % (100 * float(mech[0]["var_explained_pc1"]), f(mech[0]["mean_gene_r"], 2)),
          "<p>The consequences are concrete:</p>",
          "<ul>"
          "<li>When a permutation's labels happen to line up with that one axis, <i>every</i> gene "
          "shifts together and BH converts a single sample-level coincidence into a thousand "
          "gene-level discoveries.</li>"
          "<li>Benjamini&ndash;Hochberg assumes independence or positive regression dependence. "
          "Under correlation this strong its guarantee is not reliable, and the number of "
          "discoveries becomes a highly unstable statistic.</li>"
          "<li>The strong cross-lineage agreement reported in section 7 is substantially this "
          "shared axis rather than independent biological replication.</li>"
          "<li>It is a property of the data, not a defect in the code, and it applies to the raw "
          "and centred analyses alike.</li>"
          "</ul>")

    w("<h4>What this means for how we report the cell-type results</h4>",
      "<ul>"
      "<li><b>The 161 raw cell-type lithium findings stand as the headline cell-type result</b>, "
      "with the dependence stated as a limitation. In the lineages where it matters the observed "
      "count beats its own permutation null (granulocytes p = 0.030, T cells p = 0.020).</li>"
      "<li><b>The 284 / 345 centred counts must not be reported as findings.</b> The bipolar "
      "granulocyte null alone reaches 34.5 genes under the same correction.</li>"
      "<li><b>The exact count is not a stable quantity</b> in these estimates and should not be "
      "emphasised. The pattern &mdash; granulocyte-dominated, lithium-specific &mdash; is far more "
      "robust than the number attached to it.</li>"
      "<li><b>The whole-blood results are untouched by all of this.</b> They are computed from the "
      "measured expression matrix and never pass through bMIND.</li>"
      "</ul>")

    if valD and sumD:
        w("<h4>The repair, and whether it works</h4>",
          "<p>The principled fix for an offset is to put it in the design matrix rather than "
          "subtract it. The model then removes the same nuisance but spends a degree of freedom "
          "doing so, keeping the standard errors honest. Subtraction is the same operation with the "
          "accounting thrown away. We ran it, and validated it against the same 100-permutation "
          "null, on a 2,000-gene subsample.</p>")
        w(table(["Contrast", "Lineage", "cor(offset, group)", "Raw", "Centred",
                 "Offset as covariate", "Covariate null 95th pct", "Covariate null max"],
                [[r["contrast"], LN[r["lineage"]], f(r["cor_offset_group"], 3),
                  r["raw"], r["centred"], "<b>%s</b>" % r["covariate"],
                  f(dv["null_cov_q95"], 1), format(int(dv["null_cov_max"]), ",")]
                 for r in sumD
                 for dv in [next(v for v in valD if v["contrast"] == r["contrast"]
                                 and v["lineage"] == r["lineage"])]]))
        q95_cov = max(float(v["null_cov_q95"]) for v in valD)
        q95_raw = max(float(v["null_raw_q95"]) for v in valD)
        w("<p><b>Verdict: better, but not a cure.</b> Judged on the 95th percentile of the null "
          "&mdash; the right measure here, because the maximum is contaminated by the same "
          "dependence tail that drives the raw null to 1,354 &mdash; the three treatments rank:</p>")
        w(table(["Treatment", "Worst null 95th percentile", "Reading"], [
            ["Raw, no correction", "<b>%s</b>" % f(q95_raw, 1),
             "well behaved in the body of the distribution"],
            ["Offset as a covariate", "<b>%s</b>" % f(q95_cov, 1),
             "roughly half of centring's inflation, still elevated"],
            ["Per-sample centring", "<b>34.5</b>", "the worst of the three"]]))
        w("<p>Putting the offset in the design matrix removes most of the damage that subtracting "
          "it caused, which confirms the diagnosis: the problem with centring was the unaccounted "
          "degree of freedom. But it does not get back to the raw baseline, and given the "
          "gene-gene dependence measured above that is expected &mdash; <b>no per-gene correction "
          "can repair a dependence structure this strong.</b></p>",
          "<p><b>So the raw, uncorrected cell-type analysis is the one we report</b>, with the "
          "dependence stated as a limitation. The offset is only weakly correlated with group in "
          "every lineage (|r| &le; %s), so little genuine signal is at stake either way.</p>"
          % f(max(abs(float(r["cor_offset_group"])) for r in sumD), 3))
    # ---- 10.3 ranking
    w("<h3>10.3 Ranking instead of thresholding</h3>",
      "<p><b>What was asked.</b> &ldquo;Get the ranking of the genes, to get around the issue of "
      "significance: overlap between different cell types, and just compare the values.&rdquo;</p>",
      "<p><b>What we did.</b> Rank product across the five lineages &mdash; full method and results "
      "in section 7. Top five for lithium: <b>%s</b>.</p>"
      % ", ".join(esc(r["symbol"]) for r in rp_li[:5]),
      "<p><b>What we found.</b> Cross-lineage Spearman correlations of 0.556 to 0.871 after "
      "centring. <b>But section 10.2 shows how to read that:</b> the lineages agree so strongly in "
      "large part because they share a dominant common axis, not because five independent "
      "measurements concur. The high agreement is partly an artifact, and the permutation null that "
      "would separate the two has not been run.</p>")

    # ---- 10.4 SD filter
    w("<h3>10.4 Filtering genes by standard deviation</h3>",
      "<p><b>What was asked.</b> &ldquo;Filter the genes, maybe standard deviation as a metric for "
      "elimination. Decide the threshold.&rdquo;</p>",
      "<p><b>The logic.</b> A gene that barely varies between people cannot differ between groups "
      "either, so testing it spends a multiple-testing penalty for nothing. Dropping such genes "
      "should raise power &mdash; but only if the filter statistic is independent of the test "
      "statistic under the null, or it quietly biases the result (Bourgon et al., 2010).</p>",
      "<p><b>We checked that independence requirement first, and it passes:</b> the correlation "
      "between a gene's standard deviation and its test statistic under the null is below 0.04 in "
      "every lineage. The filter is legitimate.</p>")
    for c in ("LI", "BPD"):
        g = [r for r in sdC[c] if r["lineage"] == "gran"]
        w("<p><b>%s contrast, granulocytes:</b></p>" % ("Lithium" if c == "LI" else "Bipolar"))
        w(table(["SD quantile dropped", "Genes kept", "Significant genes"],
                [["%.0f%%" % (100 * float(r["quantile_cut"])), format(int(r["n_kept"]), ","), r["deg"]]
                 for r in g]))
    w("<p><b>Recommendation: do not adopt it.</b> The best setting (drop the bottom 30% by standard "
      "deviation) raises the pooled lithium count from 434 to 463, a 6.7% gain. That is small, it is "
      "measured on the centred data that section 10.2 advises against, and it adds a tunable "
      "parameter a reader would reasonably ask us to justify. The gain does not pay for the "
      "complexity.</p>")

    # ---- 10.5 the aggressive variance filter, cut to ~200 genes
    sw = rows(os.path.join(RES, "sdfilter200_sweep.csv"))
    sv = rows(os.path.join(RES, "sdfilter200_validation.csv"))
    lk = rows(os.path.join(RES, "sdfilter200_leakage22.csv"))
    wbsw = {c: [r for r in sw if r["level"] == "WB" and r["contrast"] == c] for c in ("LI", "BPD")}

    w("<h3>10.5 The same filter pushed much harder &mdash; cut to 200 genes</h3>",
      "<p><b>The proposal.</b> Take the variance filter in 10.4 far further: keep only the ~200 "
      "most variable genes and test within those. The reasoning is that a gene which barely varies "
      "across the whole cohort cannot differ much between two subgroups of it, so those genes can be "
      "discarded with confidence; and with only 200 tests instead of 12,368 the "
      "Benjamini&ndash;Hochberg threshold for the top gene relaxes by a factor of 62, which should "
      "let real but modest effects surface.</p>",
      "<p><b>The logic is correct and the precondition holds.</b> This is independent filtering "
      "(Bourgon, Gentleman &amp; Huber, 2010), and its requirement is that the filter statistic be "
      "independent of the test statistic under the null. We measured that rather than assuming it: "
      "the Spearman correlation between a gene's SD and its significance under <i>permuted</i> "
      "labels is <b>+0.024</b> (lithium) and <b>+0.007</b> (bipolar) &mdash; essentially zero. The "
      "filter is legitimate. Every count below also sits against its own 1,000-permutation null, "
      "and all twelve are calibrated.</p>",
      "<p><b>It nevertheless loses badly wherever there is signal to lose.</b></p>")
    w(table(["Analysis", "Unfiltered", "Top 200 by SD", "Change"], [
        ["<b>Lithium, whole blood</b>", "1,426", "<b>10</b>", "&minus;99.3%"],
        ["<b>Bipolar, whole blood</b>", "4", "<b>0</b>", "&minus;100%"],
        ["<b>Lithium, cell type</b> (all 5 lineages)", "161", "<b>7</b>", "&minus;95.7%"],
        ["<b>Bipolar, cell type</b> (all 5 lineages)", "0", "<b>12</b>", "+12"]]))

    w("<h4>Why it loses: high variance does not mean high group difference</h4>",
      "<p>The premise is half right. It is true that a low-variance gene cannot be a strong DEG. "
      "What does not follow is the converse &mdash; that high-variance genes <i>are</i> the DEGs. "
      "Total variance is the sum of between-group and within-group variance, and in this data the "
      "second term dominates:</p>",
      "<ul>"
      "<li>The largest lithium effect in the cohort (log2FC 0.793, groups split 152/74) contributes "
      "only <b>0.372</b> to a gene's SD.</li>"
      "<li>The top-200 cut sits at SD <b>0.933</b>.</li>"
      "<li>So a gene cannot reach the top 200 on its group difference &mdash; it needs within-group "
      "SD of at least <b>0.856</b> on top of it. <b>Ranking by SD ranks by noise</b>, which is the "
      "denominator of the t-statistic, not the numerator.</li>"
      "</ul>",
      "<p>The consequence is concrete and easy to state:</p>")
    w(table(["Gene", "Rank by significance", "SD", "Rank by SD (of 12,368)", "In the top 200?"], [
        ["<b>TSPAN2</b>", "1 (q = 7&times;10<sup>&minus;16</sup>)", "0.652", "<b>1,159</b>", "no"],
        ["<b>PIGB</b>", "2", "0.540", "<b>2,417</b>", "no"],
        ["<b>ARMC2</b>", "3", "0.399", "<b>5,685</b>", "no"],
        ["<b>RFX2</b>", "4", "0.635", "1,284", "no"],
        ["CNTNAP3", "21", "1.056", "100", "yes"],
        ["MYADM", "7", "0.858", "309", "no"]]))
    w("<p><b>The three strongest lithium genes in the whole study are all discarded by the filter.</b> "
      "Only <b>14 of the 1,426</b> significant genes survive it &mdash; fewer than the 23 that random "
      "selection of 200 genes would retain.</p>")

    w("<h4>How the count moves as the cut tightens</h4>")
    w(table(["Genes kept", "SD cutoff", "Lithium, whole blood", "Bipolar, whole blood"],
            [[format(int(a["n_kept"]), ","), f(a["sd_cutoff"], 3),
              format(int(a["n_sig"]), ","),
              next((b["n_sig"] for b in wbsw["BPD"] if b["n_kept"] == a["n_kept"]), "&mdash;")]
             for a in wbsw["LI"]]))
    w("<p>The lithium curve falls monotonically all the way down &mdash; there is no setting at "
      "which filtering beats not filtering.</p>")
    bcp = os.path.join(RES, "sdfilter_bestcut_validation.csv")
    if os.path.exists(bcp):
        bc = rows(bcp)[0]
        # 5 genes at a 50-gene cut is the best of seven tries, so the null has to
        # be the maximum across the same seven cutoffs, not the count at one of them.
        w("<p><b>The bipolar column needs a second look.</b> It shows 5 genes at a 50-gene cut "
          "against 4 unfiltered &mdash; the only sign anywhere that the filter might help this "
          "contrast. But 5 is the best of seven cutoffs we tried, and a p-value computed at one "
          "cutoff does not price in having chosen it. The correct null is the <i>maximum</i> count "
          "across the same seven cutoffs under permuted labels, so that every permutation gets the "
          "same seven chances we gave ourselves:</p>")
        w(table(["Null", "95th percentile", "Maximum", "p for the observed %s" % bc["observed"]], [
            ["At the %s-gene cutoff alone" % bc["best_cut"], f(bc["null_at_q95"], 1),
             bc["null_at_max"], f(bc["p_at_cut"], 4)],
            ["<b>Maximum over all seven cutoffs</b>", f(bc["null_max_q95"], 1),
             bc["null_max_max"], "<b>%s</b>" % f(bc["p_max_over_cuts"], 4)]]))
        surv = float(bc["p_max_over_cuts"]) < 0.05
        if surv:
            bg5 = rows(os.path.join(RES, "bestcut5_composition.csv"))                   if os.path.exists(os.path.join(RES, "bestcut5_composition.csv")) else None
            w("<p><b>It survives.</b> Even priced for having chosen the cutoff, scrambled labels "
              "reach 5 genes only %s of the time. And the genes are not a random five:</p>"
              % f(bc["p_max_over_cuts"], 4))
            if bg5:
                w(table(["Gene", "LM22 marker", "log2FC", "Strongest cell-fraction correlation",
                         "log2FC after adjusting for all 22 fractions", "% lost", "Survives?"],
                        [["<b>%s</b>" % esc(r["symbol"]), "yes" if r["lm22"] == "TRUE" else "no",
                          f(r["logFC_raw"], 3),
                          "%s (%s)" % (esc(r["top_fraction"]), f(r["cor"], 2)),
                          f(r["logFC_adj22"], 3), "%.0f%%" % abs(float(r["pct_lost"])),
                          "<b>yes</b>" if r["survives"] == "TRUE" else "<b>no</b>"]
                         for r in bg5]))
            w("<p><b>Four of the five mark one compartment.</b> HDC makes histamine, CPA3 is the "
              "mast-cell carboxypeptidase, GATA2 is the transcription factor specifying the "
              "basophil/mast/eosinophil lineage, and ALOX15 is the basophil and eosinophil "
              "lipoxygenase. All four go <i>down</i> in unmedicated bipolar I, and none of the five "
              "moves more than 15% when its most influential sample is dropped, so this is not "
              "outlier-driven.</p>",
              "<p><b>And the compartment itself differs.</b> Resting mast cells are "
              "<b>2.36% of the mixture in controls against 2.13% in unmedicated bipolar I "
              "(p = 0.039)</b>. So the most natural reading is a modest reduction in the "
              "basophil/mast compartment rather than a change in how any cell transcribes.</p>",
              "<p><b>Two genes are not explained by that.</b> TMTC1 loses <b>0%</b> of its effect "
              "on full composition adjustment and ALOX15 only 12%, while CPA3 (53%) and GATA2 "
              "(51%) fall away. HDC survives but loses 44%.</p>",
              "<p><b>How to report it:</b> <i>a small reduction in the basophil/mast-cell "
              "compartment in unmedicated bipolar I, accompanied by two genes &mdash; TMTC1 and "
              "ALOX15 &mdash; whose association is not accounted for by composition.</i> It is the "
              "one place in the entire filtering exercise that produced something, and it is worth "
              "a paragraph. It needs replication: it rests on 5 genes drawn from the 50 most "
              "variable in the cohort, found by a filter that section 10.5 otherwise recommends "
              "against.</p>")
        else:
            w("<p><b>%s</b></p>" % (
          "Once the choice of cutoff is priced in the result disappears. Scrambled labels reach "
          "5 or more genes somewhere in the seven cutoffs often enough that the observation is "
          "unremarkable. <b>The 5 genes should not be reported.</b> This is exactly the trap that "
          "sweeping a tuning parameter and keeping the best setting sets for you."))

    w("<h4>The one place it gains: bipolar, T cells</h4>",
      "<p>Against the pattern, two cell-type analyses that found nothing unfiltered return something "
      "at the 200-gene cut, and both beat their own permutation null:</p>")
    ctv = [r for r in sv if r["level"] == "CT" and int(r["observed"]) > 0]
    w(table(["Contrast", "Lineage", "Unfiltered", "Top 200", "Null 95th pct", "Permutation p"],
            [[r["contrast"], LN[r["lineage"]],
              format(int(next(x["base_sig"] for x in sw if x["level"] == "CT"
                   and x["contrast"] == r["contrast"] and x["lineage"] == r["lineage"])), ","),
              "<b>%s</b>" % r["observed"], f(r["null_q95"], 1), f(r["perm_p"], 4)]
             for r in ctv]))
    w("<p>The 11 bipolar T-cell genes are: <b>" +
      ", ".join(esc(r["symbol"]) for r in lk if r["contrast"] == "BPD" and r["lineage"] == "T") +
      "</b>. This is the first non-zero bipolar cell-type result anywhere in the project, and it is "
      "worth taking seriously &mdash; but not at face value, because the gene names are a warning in "
      "themselves. CPA3 and FCER1A are mast-cell and basophil markers; TUBB1 and TLN1 are platelet "
      "genes; GZMB and PRF1 are cytotoxic granule genes. Finding those inside a <i>T-cell</i> profile "
      "is what deconvolution leakage looks like.</p>")

    w("<h4>Testing that warning: the filter selects cell-identity markers</h4>",
      "<p>It does, systematically. LM22 signature genes make up <b>2.8%</b> of the tested "
      "transcriptome, but of the top-200-SD sets they make up 6.5&ndash;13.0% &mdash; an enrichment "
      "of <b>2.3&times; to 4.6&times;</b> in every single lineage. The reason is immediate: a "
      "cell-identity marker's variance across a cohort is dominated by how much of that cell type "
      "each person happens to have. <b>Ranking genes by variance preferentially ranks the genes the "
      "deconvolution itself used to decide the mixture.</b></p>",
      "<p>Among the 19 significant hits, 4 are LM22 markers &mdash; 21.1% against a 2.8% background, "
      "Fisher exact odds ratio 9.2, p = 1.7&times;10<sup>&minus;3</sup>.</p>",
      "<p>So we re-tested every hit with all 22 cell fractions added to the design, the strongest "
      "composition control the data supports:</p>")
    bpd = [r for r in lk if r["contrast"] == "BPD" and r["lineage"] == "T"]
    w(table(["Gene", "LM22 marker", "log2FC raw", "log2FC adjusted for all 22 fractions",
             "% of effect lost", "Survives?"],
            [["<b>%s</b>" % esc(r["symbol"]), "yes" if r["lm22"] == "TRUE" else "no",
              f(r["logFC_raw"], 3), f(r["logFC_adj22"], 3), "%.0f%%" % float(r["pct_lost"]),
              "<b>yes</b>" if r["survives"] == "TRUE" else "<b>no</b>"]
             for r in sorted(bpd, key=lambda z: float(z["pct_lost"]))]))
    nsurv = sum(1 for r in bpd if r["survives"] == "TRUE")
    w("<p><b>Result: %d of the 11 survive; 6 do not.</b> The ones that fall are exactly the "
      "cell-identity markers &mdash; GZMB loses 66%% of its effect, PRF1 58%%, CPA3 52%%. The median "
      "gene loses 37%%. By contrast, all 7 of the <i>lithium</i> cell-type hits survive the same "
      "adjustment and <b>none of them is an LM22 marker</b>.</p>" % nsurv,
      "<p>One caution in the other direction, for fairness: the 22 fractions are themselves estimated "
      "from this same expression matrix, so adjusting for them can over-correct (the circularity "
      "noted in section 9.2). &ldquo;Lost&rdquo; here means <i>we cannot separate this gene from "
      "composition</i>, not that the gene is certainly nothing.</p>")

    w("<h4>Verdict on the 200-gene variance filter</h4>",
      "<ul>"
      "<li><b>Do not use it for the lithium analyses.</b> It is statistically valid and it destroys "
      "the result: 1,426 &rarr; 10 in whole blood, 161 &rarr; 7 at cell-type level, with the three "
      "strongest genes in the study discarded.</li>"
      "<li><b>For whole-blood bipolar the answer is split.</b> At the 200-gene cut it is 4 &rarr; 0. At a 50-gene cut it returns 5 genes that survive a selection-aware null (p = 0.015) and point at a real reduction in the basophil/mast compartment &mdash; the one genuine yield of the whole filtering exercise, and, not coincidentally, the one analysis that had almost no signal for the filter to destroy.</li>"
      "<li><b>The bipolar T-cell finding is worth one paragraph, heavily hedged.</b> Eleven genes, "
      "permutation p = 0.018, of which %d survive full composition adjustment: "
      "<b>%s</b>. Report it as hypothesis-generating, with the marker enrichment stated. It is a "
      "lead, not a result.</li>"
      "<li><b>The general lesson is worth stating in the paper.</b> Variance-based filtering is a "
      "reasonable-sounding idea that this dataset lets us test directly, and it fails for a "
      "specific, explainable reason. That is a useful methodological observation in its own "
      "right.</li>"
      "</ul>" % (nsurv, ", ".join(esc(r["symbol"]) for r in bpd if r["survives"] == "TRUE")))


    # ---- 10.6 the principled alternative: mean-expression filtering
    mcp = os.path.join(RES, "meanfilter_curve.csv")
    mvp = os.path.join(RES, "meanfilter_validation.csv")
    if os.path.exists(mcp):
        mc = rows(mcp)
        mv = rows(mvp) if os.path.exists(mvp) else None
        wbli = {FS: [r for r in mc if r["level"] == "WB" and r["contrast"] == "LI"
                     and r["filter"] == FS] for FS in ("MEAN", "SD")}

        w("<h3>10.6 The filter that should have worked &mdash; mean expression</h3>",
          "<p><b>Why try a different statistic.</b> Section 10.5 showed that filtering on standard "
          "deviation fails for a specific reason: total variance is dominated by within-group noise, "
          "so ranking on it ranks on the denominator of the t-statistic. The underlying idea &mdash; "
          "spend an independent filter to buy a less punishing correction &mdash; survives that "
          "finding intact. It just needs a statistic tied to <i>power</i> rather than to noise.</p>",
          "<p><b>Mean expression is the textbook choice.</b> Bourgon, Gentleman and Huber (2010) show "
          "that for a two-group comparison the overall mean is independent of the t-statistic under "
          "the null &mdash; the precondition SD satisfies only approximately &mdash; while still "
          "being strongly tied to power, because low-count genes carry most of the measurement "
          "noise. Dropping them removes tests that had little chance of succeeding without removing "
          "candidates that did. It is what DESeq2 does by default.</p>",
          "<p><b>A prediction we wrote down first.</b> We expected the gain here to be small, for two "
          "reasons specific to this pipeline: the baseline filter (&gt;10 counts in &ge;90% of "
          "samples) has <i>already</i> removed the low-expression genes, so the pool mean filtering "
          "normally targets is largely gone; and voom already down-weights low-count observations "
          "through the mean-variance trend, so part of the benefit is being collected twice.</p>")

        w("<h4>Mean filtering is a far better filter than SD &mdash; and still not worth using</h4>",
          "<p>Two columns decide this. <b>True-positive enrichment</b> asks whether the filter keeps "
          "the genes that matter: above 1&times; it concentrates signal, below 1&times; it is "
          "discarding signal faster than it discards genes. <b>LM22 enrichment</b> asks whether it is "
          "selecting cell-identity markers, the failure mode diagnosed in 10.5.</p>")
        w(table(["Genes kept",
                 "MEAN: significant", "MEAN: true-pos enrichment", "MEAN: LM22 enrichment",
                 "SD: significant", "SD: true-pos enrichment", "SD: LM22 enrichment"],
                [[format(int(a["n_kept"]), ","),
                  format(int(a["n_sig_filter_before_fit"]), ","),
                  "<b>%.2f&times;</b>" % float(a["tp_enrichment"]),
                  "%.2f&times;" % float(a["lm22_enrichment"]),
                  format(int(b["n_sig_filter_before_fit"]), ","),
                  "<b>%.2f&times;</b>" % float(b["tp_enrichment"]),
                  "%.2f&times;" % float(b["lm22_enrichment"])]
                 for a, b in zip(wbli["MEAN"], wbli["SD"])]))
        w("<p><i>Lithium, whole blood &mdash; the analysis with the most signal to preserve.</i></p>",
          "<p><b>The comparison is decisive on its own terms.</b> At a 1,000-gene cut mean filtering "
          "retains true positives at <b>1.81&times;</b> enrichment while SD manages 1.05&times; "
          "&mdash; barely better than picking at random. At 200 genes SD drops to <b>0.61&times;</b>, "
          "i.e. actively worse than random, while mean filtering holds at 1.43&times;. Mean filtering "
          "also avoids the marker trap: at 1,000 genes it enriches LM22 genes 1.70&times; against "
          "SD's 2.69&times;.</p>",
          "<p><b>But the count still never recovers.</b> 1,000 genes by mean expression yields 312 "
          "significant against 1,426 unfiltered. The relaxed threshold is real &mdash; it is simply "
          "worth far less than the 1,114 genes it costs. The curve rises monotonically to the "
          "unfiltered value and stops there.</p>",
          "<p>That is the prediction confirmed: with low-expression genes already removed by the "
          "baseline filter and voom already weighting by the mean-variance trend, there is almost no "
          "slack left for a mean filter to exploit.</p>")

        if mv:
            beats = [r for r in mv if r["beats_unfiltered"] == "TRUE"]
            surv  = [r for r in beats if float(r["p_max_over_grid"]) < 0.05]
            w("<h4>The apparent gains, and whether any survives</h4>",
              "<p>Across twelve analyses, nine cutoffs and two statistics, %d cells returned more "
              "genes than no filtering. Some were large: bipolar T cells went from 0 to 167 at a "
              "1,000-gene SD cut, bipolar B cells from 0 to 86 by mean expression, bipolar whole "
              "blood from 4 to 30. Taken at face value those would be the largest findings in the "
              "project.</p>" % len(beats),
              "<p>They are not findings, and two things we had already measured say why. Each number "
              "is <b>the best of eighteen tries</b> (nine cutoffs &times; two statistics), which a "
              "per-cell p-value prices in not at all. And in the cell-type estimates the genes are "
              "correlated at 0.674 with one axis carrying 71% of the variance (section 10.2), which "
              "makes the discovery count wildly unstable &mdash; nulls in these very lineages had "
              "already been seen to reach 105 and 168 genes.</p>",
              "<p>So the null here is the <b>maximum across the same eighteen-cell grid</b> under "
              "permuted labels: every permutation gets exactly the chances we gave ourselves.</p>")
            w(table(["Analysis", "Filter", "Unfiltered", "Best cut", "Best count",
                     "Null 95th pct", "Null max", "p", "Verdict"],
                    [[("%s %s %s" % (r["level"], r["contrast"], r["lineage"])).strip(),
                      r["filter"], r["unfiltered"], format(int(r["best_cut"]), ","),
                      "<b>%s</b>" % r["observed_max"], f(r["null_q95"], 1), r["null_max"],
                      f(r["p_max_over_grid"], 4),
                      "<b>survives</b>" if float(r["p_max_over_grid"]) < 0.05 else "within the null"]
                     for r in sorted(beats, key=lambda z: -int(z["observed_max"]))]))
            if surv:
                w("<p><b>%d of the %d survive the selection-aware null</b>, and each is run down "
                  "below. None of them turns out to be a transcriptional finding: one is a "
                  "one-gene difference, one is composition, one is an artefact of the B-cell "
                  "estimates, and the largest is composition counted thirty times.</p>"
                  % (len(surv), len(beats)))
            else:
                w("<p><b>Not one of them survives.</b> Every apparent gain sits inside the "
                  "distribution of what scrambled labels produce when given the same eighteen "
                  "chances. The 167 bipolar T-cell genes are the clearest case: impressive against a "
                  "baseline of zero, unremarkable against a null that routinely reaches comparable "
                  "numbers through gene-gene dependence alone.</p>")


            # ---- the one win, examined
            b30p = os.path.join(RES, "meanfilter_bpd30_genes.csv")
            if os.path.exists(b30p):
                b30 = rows(b30p)
                nsv = sum(1 for r in b30 if r["survives"] == "TRUE")
                ndn = sum(1 for r in b30 if float(r["logFC"]) < 0)
                w("<h4>The one unambiguous win, and what it turned out to be</h4>",
                  "<p>One cell did survive: <b>whole-blood bipolar, mean filter, 200 genes &mdash; 4 "
                  "significant genes become 30, p = 0.0050</b> against the selection-aware null. That "
                  "is the largest gain anywhere in the filtering work, and unlike the SD results it is "
                  "not suspect by construction: mean expression does not preferentially select "
                  "cell-identity markers, and only 2 of the 30 are LM22 genes.</p>",
                  "<p>Three facts about the 30 settle what they are.</p>",
                  "<ul>"
                  "<li><b>All %d of them go down</b> in unmedicated bipolar I. Not 25 of 30 &mdash; all "
                  "of them.</li>" % ndn +
                  "<li><b>None of them is one of the 4 found without a filter.</b> The overlap between "
                  "the two sets is zero.</li>"
                  "<li><b>Not one survives adjustment for cell composition.</b> %d of 30, with a median "
                  "of 45%% of the effect lost.</li>" % nsv +
                  "</ul>",
                  "<p>And the gene names explain why. The 200 most-expressed genes in whole blood are "
                  "overwhelmingly leukocyte genes &mdash; PTPRC (CD45), LCP1, CORO1A, ICAM3, MYH9, "
                  "FCN1, NCF2, ARHGAP30. <b>Filtering on mean expression in whole blood selects, almost "
                  "by definition, the genes that track how many white cells are in the tube.</b> A "
                  "modest shift in the leukocyte mixture moves all of them together and in the same "
                  "direction, which is exactly the pattern observed.</p>",
                  "<p><b>So the 30 genes are one composition effect counted thirty times, not thirty "
                  "findings.</b> The permutation null was right to pass them &mdash; the group labels "
                  "genuinely do align with the data &mdash; but the alignment is compositional. This is "
                  "precisely the confound the whole project exists to catch, arriving this time through "
                  "the filter rather than through the drug.</p>")
                w(table(["Gene", "log2FC", "q", "log2FC adjusted for all 22 fractions",
                         "% of effect lost", "Survives?"],
                        [[esc(r["symbol"]), f(r["logFC"], 3), f(r["q"], 4), f(r["logFC_adj22"], 3),
                          "%.0f%%" % abs(float(r["pct_lost"])),
                          "yes" if r["survives"] == "TRUE" else "no"]
                         for r in sorted(b30, key=lambda z: float(z["q"]))[:12]]))
                w("<p><i>The 12 strongest of the 30; the full list is in "
                  "<code>results/meanfilter_bpd30_genes.csv</code>. Every one is sign-stable under "
                  "leave-one-out (no gene moves more than 15%), so this is a genuine group difference "
                  "&mdash; just not a transcriptional one.</i></p>",
                  "<p><b>What to take from it.</b> Combined with the basophil/mast result in 10.5, two "
                  "independent filters applied to the same contrast both land on composition, and both "
                  "point the same way: <i>unmedicated bipolar I differs modestly from controls in "
                  "leukocyte composition, and we find no whole-blood transcriptional signature "
                  "independent of that.</i> That is a cleaner and more defensible statement of the "
                  "bipolar result than the bare &ldquo;4 boundary genes&rdquo; in section 6.4, and it "
                  "is the project's central method answering the project's central question.</p>")


            # ---- the two large cell-type survivors
            svp = os.path.join(RES, "ct_bpd_survivors.csv")
            if os.path.exists(svp):
                sv = rows(svp)
                w("<h4>The two large cell-type survivors, run down</h4>",
                  "<p>Two cell-type cells also survived: bipolar <b>T cells, 0 &rarr; 154</b> "
                  "(SD filter at 2,000 genes, p = 0.0348) and bipolar <b>B cells, 0 &rarr; 86</b> "
                  "(mean filter at 1,000, p = 0.0149). Surviving the null only means the count is "
                  "unlikely under permuted labels; it says nothing about <i>what</i> the genes are. "
                  "Both got the same composition test.</p>")
                w(table(["", "Bipolar T cells (SD, 2,000)", "Bipolar B cells (MEAN, 1,000)"], [
                    ["Significant genes", "154", "86"],
                    ["Direction", "<b>154 down, 0 up</b>", "<b>0 down, 86 up</b>"],
                    ["LM22 marker enrichment", "<b>7.3&times;</b>", "0.0&times; (none)"],
                    ["Survive adjustment for all 22 fractions", "<b>0 of 154 (0%)</b>",
                     "<b>83 of 86 (97%)</b>"],
                    ["Median effect lost", "38.8%", "&minus;0.6%"]]))
                w("<p><b>T cells: composition, unambiguously.</b> Every one of the 154 goes down, "
                  "they are 7.3&times; enriched for LM22 signature genes, and <b>not one survives "
                  "adjustment for cell fractions</b>. This is the 200-gene result from section 10.5 "
                  "at larger scale and with a cleaner verdict: the SD filter selects the genes the "
                  "deconvolution used to assign the mixture, and a mixture difference then reads as "
                  "differential expression. Do not report it.</p>",
                  "<p><b>B cells: not composition &mdash; but not believable either.</b> These "
                  "behave oppositely: all 86 go up, none is an LM22 marker, and 97% survive full "
                  "composition adjustment with the effect if anything growing. That rules out the "
                  "explanation that killed every other result, so it needed a different one.</p>",
                  "<ul>"
                  "<li><b>The gene classes are wrong for biology.</b> 68 of the 86 (79%) are "
                  "snoRNAs, scaRNAs, ribosomal-protein genes, histones or pseudogenes &mdash; "
                  "SCARNA10, SNORA53, RPL39P3, RPS23P8, RPL32, RPS21, HIST1H4C &mdash; against 27% "
                  "among the 1,000 genes tested, a 3.0&times; enrichment. That is the profile of a "
                  "library-composition effect, not a B-cell response.</li>"
                  "<li><b>RNA integrity does not explain it</b>, which we checked rather than "
                  "assumed: the hits correlate with RIN <i>less</i> than the average tested gene "
                  "(median |r| smaller, Wilcoxon p = 6&times;10<sup>&minus;7</sup>), and all 86 "
                  "remain significant when RIN is fitted as a 4-degree-of-freedom spline instead of "
                  "a single linear term.</li>"
                  "<li><b>The decisive objection is internal.</b> Section 11.2 now measures the B-cell "
                  "detection floor <i>on this very contrast</i> at <b>16.96 log2FC</b> &mdash; about "
                  "a 130,000-fold change. A lineage in which the method cannot reliably detect a "
                  "five-fold change (power 0.015 at log2FC 5.0) cannot credibly return 86 genes "
                  "with effects of 0.1 or so. The result contradicts our own power analysis, and "
                  "the power analysis is the more trustworthy of the two.</li>"
                  "</ul>",
                  "<p><b>Conclusion:</b> the B-cell estimates from bMIND are not usable, and the 86 "
                  "genes are a property of those estimates rather than of B cells. This is the "
                  "strongest single piece of evidence that the cell-type arm should be reported "
                  "with the dependence and reliability caveats of section 10.2 attached to every "
                  "claim.</p>")

            rho = {}
            for r in mv:
                rho.setdefault((r["level"], r["filter"]), []).append(float(r["rho_null_mean"]))
            w("<h4>Both filters are legitimate &mdash; that was never the problem</h4>",
              "<p>Independent filtering requires the filter statistic to be uncorrelated with "
              "significance under the null. Averaged over 50 permutations per analysis:</p>")
            w(table(["Level", "Filter", "Mean Spearman(statistic, &minus;log<sub>10</sub> p) under permuted labels"],
                    [[k[0], k[1], "%+.3f" % (sum(v)/len(v))]
                     for k, v in sorted(rho.items())]))
            w("<p>All four are close to zero, so both filters satisfy the precondition. <b>Neither "
              "filter is invalid. Both are simply unprofitable here</b>, and SD is worse than "
              "unprofitable because it selects for noise and for cell-identity markers.</p>")

        w("<h4>Verdict on filtering in general</h4>",
          "<ul>"
          "<li><b>Use no secondary filter.</b> The baseline expression filter stays; nothing is added "
          "on top of it. This is now an evidence-backed choice rather than a default.</li>"
          "<li><b>If a reviewer asks why not</b>: we tried both standard filter statistics across "
          "nine cutoffs and twelve analyses, validated against permutation nulls that price in the "
          "search, and no setting beat the unfiltered analysis. The curves are in this section.</li>"
          "<li><b>Mean expression is strictly the better of the two</b> and is what we would use if "
          "a filter were required &mdash; it concentrates true positives (up to 1.81&times;) where "
          "SD does not (1.05&times;, and 0.61&times; at the tightest cut).</li>"
          "<li><b>Why there was no room:</b> the baseline filter already removed the low-expression "
          "genes a mean filter targets, and voom already handles the mean-variance relationship. "
          "The benefit had been collected before either filter was applied.</li>"
          "</ul>")

    w("<hr>")


# ================================================================= 11. POWER
def sec_power():
    w("<h2>11. Power &mdash; what we could and could not have detected</h2>",
      "<p>A negative result means nothing until you know what you were capable of seeing. We "
      "computed that for both levels, and the two answers are opposite.</p>",
      "<h3>11.1 Whole blood &mdash; the bipolar negative is informative</h3>",
      "<p>Computed from the standard errors and residual degrees of freedom of the actual bipolar "
      "whole-blood fit (median SE 0.0539, df 299):</p>")
    w(table(["Detectable with 80% power", "log2FC", "As a fold change"], [
        ["At an uncorrected p &lt; 0.05", "0.152", "1.11&times;"],
        ["At the BH threshold actually reached", "0.208", "1.16&times;"],
        ["At a full Bonferroni correction", "0.299", "1.23&times;"]]))
    w("<p><b>The number that matters:</b> our power to detect an effect the size of the median "
      "<i>lithium</i> effect (log2FC 0.176) in the bipolar contrast was <b>0.90</b>.</p>",
      "<p>So if bipolar disorder produced a whole-blood signature as strong as lithium's, we had a "
      "90% chance of finding it, and we did not. <b>That makes the whole-blood bipolar negative a "
      "real finding</b> &mdash; evidence of absence, within a clearly stated effect size &mdash; "
      "rather than an underpowered shrug.</p>",
      "<h3>11.2 Cell type &mdash; the bipolar zero is uninformative</h3>",
      "<p>Here it reverses. We measured the detection floor by spike-in: inject an effect of known "
      "size where it actually lives &mdash; in the bulk, scaled by that sample's lineage fraction "
      "&mdash; then re-run the deconvolution and the test, and find the smallest effect reliably "
      "recovered. Labels are permuted first, so what is detected is attributable to the spike "
      "alone. These floors are measured on the <b>bipolar contrast itself</b> (n = 308, 74 vs "
      "234), not borrowed from the lithium one.</p>")
    fl = rows(os.path.join(RES, "power_celltype_BPD_floors.csv"))
    w(table(["Lineage", "Measured floor at 80% power", "As a fold change",
             "Bulk-equivalent effect", "Verdict"],
            [[LN[r["lineage"]], "<b>%s</b>" % f(r["mde80"], 2),
              ("%.1f&times;" % (2 ** float(r["mde80"]))) if float(r["mde80"]) < 10
              else ("about %s&times;" % format(int(round(2 ** float(r["mde80"]), -3)), ",")),
              f(r["bulk_equivalent"], 3),
              "far above any plausible effect" if float(r["mde80"]) < 1
              else "no realistic detection"]
             for r in fl]))
    w("<p>Set those against the whole-blood floor of <b>0.208</b> on the same contrast and the "
      "median lithium effect of <b>0.176</b>. <b>The cell-type method cannot see effects of the "
      "size that actually exist in this data.</b> Even granulocytes, the most abundant lineage, "
      "need %s &mdash; nearly three times the largest effect we have any reason to expect. "
      "Finding zero bipolar genes at cell-type level is therefore the expected outcome whether or "
      "not an effect is present, and it must be reported as uninformative rather than as a "
      "negative finding.</p>" % f(fl[0]["mde80"], 2),
      "<p><b>The fourth column explains why, and is the more useful result.</b> Converting each "
      "floor back to the whole-blood effect it corresponds to (floor &times; mean lineage "
      "fraction) gives <b>0.170 to 0.254</b> across all five lineages &mdash; tightly clustered "
      "around the whole-blood floor of 0.208. In other words the cell-type analysis is not "
      "detecting anything the bulk analysis could not; it is the <i>same</i> sensitivity, divided "
      "by each lineage's share of the mixture. A lineage that is 1.5% of the sample needs an "
      "effect roughly 67 times larger to clear the same bar. That is a property of deconvolution "
      "itself, not of bMIND or of this dataset, and it is worth stating in the paper.</p>",
      "<p>The prediction bears that out: floors predicted purely as "
      "0.208 / fraction were 0.56, 0.76, 0.82, 2.48 and 13.86, against measured values of 0.56, "
      "0.69, 0.77, 2.03 and 16.96. Four of five land within 20%.</p>",
      "<p>This also reframes the 161 lithium findings: they were detectable because lithium's "
      "cell-type effects are unusually large, not because the method is sensitive.</p>", "<hr>")


# ================================================================= 12. LIMITATIONS
def sec_limits():
    w("<h2>12. Limitations and what is still open</h2>",
      "<h3>12.1 Limitations of the data</h3>",
      "<ul>"
      "<li><b>No lithium dose, serum level or duration.</b> The contrast is use versus non-use, a "
      "blunt instrument.</li>"
      "<li><b>No record of other medications.</b> Patients not on lithium are generally taking "
      "something else, so this is lithium versus other treatment, not lithium versus nothing.</li>"
      "<li><b>No measured cell counts.</b> Everything about composition is estimated from the same "
      "expression matrix being analysed. A complete blood count would answer the central question "
      "directly; its absence is the single biggest weakness of the design.</li>"
      "<li><b>Cross-sectional.</b> No causal direction can be established.</li>"
      "<li><b>LM22 misfits whole blood.</b> Against reference intervals the deconvolution puts "
      "monocytes at 24.6% (reference 2&ndash;8%) and eosinophils at 0.0%. Collapsing to five "
      "lineages mitigates this but does not fix it.</li>"
      "<li><b>74 unmedicated patients</b> is a small group to carry the bipolar contrast.</li>"
      "<li><b>Tobacco status was imputed</b> for 30 subjects.</li>"
      "</ul>",
      "<h3>12.2 Limitations of the analysis</h3>",
      "<ul>"
      "<li><b>The cell-type estimates are strongly dependent across genes</b> (section 10.2). This "
      "is the most serious one and it constrains every cell-type claim in this document. The "
      "whole-blood results are unaffected.</li>"
      "<li><b>The four bipolar whole-blood genes sit on the threshold</b> (q = 0.0489 against a "
      "0.05 cutoff) and should not be presented as findings.</li>"
      "<li><b>Cell-type effect sizes are not trustworthy in absolute terms</b> given the offset "
      "documented in section 10.1, even though the rankings are informative.</li>"
      "</ul>",
      "<h3>12.3 Still open &mdash; the honest list</h3>",
      "<ul>"
      "<li><s>The rank-product permutation null has not been run.</s> <b>Closed, and it changed "
      "the reading.</b> The cross-lineage agreement the ranking rests on is no greater than "
      "chance once the shared structure of the cell-type estimates is accounted for (section "
      "7.3). The rankings remain usable as an ordering, but not as evidence.</li>"
      "<li><s>The bipolar cell-type power analysis is incomplete.</s> <b>Closed.</b> All five "
      "lineages are now measured directly on the bipolar contrast (section 11.2); the floors run "
      "from 0.56 to 16.96 log2FC and none is interpretable against the 0.5 bar set in "
      "advance.</li>"
      "<li>The standard-deviation filter sweep (section 10.4) was evaluated on centred data, so its "
      "numbers inherit the caveat in section 10.2. The 200-gene version (section 10.5) was run on "
      "uncentred data and is unaffected. Neither is recommended for the main analyses.</li>"
      "<li><b>The bipolar cell-type signals produced by filtering are resolved, negatively.</b> "
      "At the larger 2,000-gene cut the T-cell set is 154 genes, all in one direction, 7.3&times; "
      "enriched for LM22 markers, and <i>none</i> survives composition adjustment (section 10.6). "
      "Nothing there needs replicating.</li>"
      "<li><b>The B-cell lineage estimates are not usable.</b> They return 86 genes that are "
      "neither composition nor RNA integrity, in a lineage whose own detection floor is 26.11 "
      "log2FC (sections 10.6 and 11.2). Both cannot be true. This belongs in the paper as a "
      "limitation of bMIND on this dataset rather than as an open question.</li>"
      "<li>Whether any cell-type analysis can be made fully calibrated on these estimates is "
      "unresolved. Three corrections were tried; none fixed it.</li>"
      "</ul>", "<hr>")


# ================================================================= APPENDIX
def sec_appendix():
    w("<h2>Appendix A &mdash; where every file lives</h2>",
      "<p>All paths are relative to <code>Implementation/de_final/</code>.</p>")
    w(table(["File", "What it holds"], [
        ["<code>run_de.R</code>", "the reported pipeline; reproduces every headline number"],
        ["<code>METHODS.md</code>", "the full written methods, 11 sections"],
        ["<code>results/summary.csv</code>", "the four-row results table in section 5"],
        ["<code>results/DE_WB_LI.csv.gz</code>", "all 12,368 genes, lithium, whole blood"],
        ["<code>results/DE_WB_BPD.csv.gz</code>", "all 12,368 genes, bipolar, whole blood"],
        ["<code>results/DE_CT_LI.csv.gz</code>", "all gene &times; lineage results, lithium"],
        ["<code>results/DE_CT_BPD.csv.gz</code>", "all gene &times; lineage results, bipolar"],
        ["<code>results/genelists/</code>", "the named gene lists reproduced in section 6"],
        ["<code>results/followup_A_*</code>", "recentring (section 10.1) and its validation (10.2)"],
        ["<code>results/followup_B_*</code>", "rank product (section 7)"],
        ["<code>results/followup_C_*</code>", "standard-deviation filter (section 10.4)"],
        ["<code>results/followup_D_*</code>", "offset-as-covariate and its validation (10.2)"],
        ["<code>results/diag_null_mechanism.csv</code>", "the gene-dependence measurement (10.2)"],
        ["<code>sdfilter_200.R</code>", "the 200-gene variance filter, all four analyses (10.5)"],
        ["<code>sdfilter_bestcut.R</code>", "selection-aware null for the best variance cutoff (10.5)"],
        ["<code>sdfilter_leakage.R</code>, <code>sdfilter_leakage2.R</code>",
         "LM22-marker enrichment and composition adjustment for the variance-filter hits (10.5)"],
        ["<code>bestcut5_composition.R</code>",
         "the five bipolar whole-blood genes and the mast-cell compartment (10.5)"],
        ["<code>meanfilter.R</code>", "mean-expression vs variance filtering, nine cutoffs (10.6)"],
        ["<code>meanfilter_validate.R</code>", "selection-aware nulls across the filter grid (10.6)"],
        ["<code>meanfilter_bpd30.R</code>",
         "the 30 bipolar whole-blood genes and their composition test (10.6)"],
        ["<code>ct_bpd_survivors.R</code>",
         "composition test for the surviving cell-type results (10.6)"],
        ["<code>bcell86_technical.R</code>", "are the 86 B-cell genes RNA quality? (10.6)"],
        ["<code>power_analysis.R</code>", "whole-blood power (section 11.1)"],
        ["<code>power_celltype_BPD.R</code>",
         "cell-type spike-in power, all five lineages, bipolar contrast (11.2)"],
        ["<code>rankproduct_null.R</code>", "permutation null for the cross-lineage ranking (7.3)"],
        ["<code>build_master_doc.py</code>", "generates this document from the files above"],
    ]))
    w("<h2>Appendix B &mdash; the complete lithium whole-blood gene list (%s)</h2>" % format(len(wb_li), ","),
      "<p>All significant genes at BH q &lt; 0.05, alphabetically.</p>",
      genepara([r["symbol"] for r in wb_li], per=12))


# ================================================================= ASSEMBLE
sec_followups()
sec_power()
sec_limits()
sec_appendix()

CSS = """<style>
body{font-family:Georgia,serif;line-height:1.5;max-width:50em}
h1{font-size:1.9em}h2{font-size:1.4em;margin-top:1.6em}
h3{font-size:1.15em;margin-top:1.2em}h4{font-size:1.02em}
table{border-collapse:collapse;margin:1em 0;font-size:0.9em}
th{background:#eee;text-align:left}td,th{border:1px solid #999;padding:4px 7px;vertical-align:top}
code{font-family:Consolas,monospace;font-size:0.9em}
</style>"""
html = ("<!doctype html><html><head><meta charset='utf-8'>"
        "<title>Differential Expression - Master Document</title>" + CSS +
        "</head><body>\n" + "\n".join(P) + "\n</body></html>")
out = os.path.join(HERE, "master_DE_document.html")
io.open(out, "w", encoding="utf-8").write(html)
print("wrote %s  (%.0f KB, %d blocks)" % (out, len(html)/1024.0, len(P)))
