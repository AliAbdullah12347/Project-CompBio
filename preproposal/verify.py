# -*- coding: utf-8 -*-
"""Check every comment from AliAva_PreProposal_Feedback.pdf against the revised text."""
import zipfile, xml.etree.ElementTree as ET

W = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
z = zipfile.ZipFile('Ali, Ava Project Preproposal (revised).docx')
root = ET.fromstring(z.read('word/document.xml'))
T = ''.join(t.text or '' for t in root.iter(W + 't'))

# (feedback comment, substring that must now be present)
CHECKS = [
 ("p1 no section headings",                 "1. Background and Motivation"),
 ("p1 heritability misdefined",             "variance in liability"),
 ("p1 cite Krebs neutrophil result",        "themselves reported elevated neutrophil"),
 ("p1 Boltz misstated",                     "continued to detect lithium-associated genes"),
 ("p1 what does your study add",            "What this project adds beyond Boltz"),
 ("p1 other psychotropics recorded?",       "no record of other medications"),
 ("p1 lithium dose known?",                 "no lithium dose or serum level"),
 ("p1 diagnosis signal change unclear",     "This is over-adjustment"),
 ("p2 hypothesis packs 3 claims",           "H1 (primary, prediction)"),
 ("p2 define significant proportion",       "exceeds 0.5"),
 ("p2 equal prediction != mediation",       "would not by itself demonstrate mediation"),
 ("p2 which result separates them",         "residualised on the ILR coordinates"),
 ("p2 moderated mediation power",           "H4 (exploratory, declared underpowered)"),
 ("p2 QC flag tracks diagnosis?",           "3.04"),
 ("p2 which fields missing / imputation",   "zero missingness post-QC"),
 ("p2 never impute exposure/outcome",       "never imputed"),
 ("p2 impute inside training folds",        "fitted inside training folds only"),
 ("p2 fractions from same matrix",          "Circularity: the mediator is estimated"),
 ("p2 rerun w/ different deconv methods",   "Cross-method sensitivity"),
 ("p2 report matrix dimensions",            "57,773 × 480"),
 ("p2 best-method choice in folds",         "made inside training folds rather than once"),
 ("p2 CIBERSORTx portal + account now",     "requires a registered, approved academic account"),
 ("p3 10% test set cannot show equivalence","repeated nested cross-validation"),
 ("p3 zero replacement rule",               "multiplicative simple replacement"),
 ("p3 does that choice change ratios",      "sensitivity analysis over the replacement rule"),
 ("p3 mediation unmeasured confounding",    "no unmeasured confounding"),
 ("p3 how choose E-value bound",            "declared in advance rather than tuned"),
 ("p3 is BMI recorded?",                    "no BMI or body-composition field"),
 ("p3 pre-specify primary mediation model", "pre-specify a single primary model"),
 ("p3 which pair isolates lithium",         "on lithium against bipolar I off lithium"),
 ("p3 PPIXpress is the wrong tool",         None),   # must be ABSENT
 ("p3 state the margin now",                "ΔAUC = 0.05"),
 ("p3 why 'only' myeloid/lymphoid",         "which was too strong"),
 ("p3 reconcile 12,353 / 12,343",           "12,353 is correct"),
 ("p4 preservation threshold",              "Zsummary"),
 ("p4 arms are not independent",            "The three arms are not independent"),
 ("p4 backup if deconvolution fails",       "collapse to coarse lineage resolution"),
 ("p4 module stability across draws",       "check module stability before comparing"),
 ("p4 power calculation for mediation",     "6.3 Statistical power"),
 ("p5 which reference for fractions",       "reliability and plausibility"),
 ("p5 how pick best without blood counts",  "no ground truth against which to rank"),
 ("p6 test exposure-by-mediator interaction","exposure-by-mediator interaction"),
 ("p6 define lineage balance trees",        "hematopoietic lineage tree"),
 ("p6 define product of coefficients",      "product of coefficients (Baron & Kenny"),
 ("p7 inner tuning loop",                   "inner 5-fold loop"),
 ("p7 figure numbers + captions",           "Figure 1. Preprocessing module."),
 ("p7 why Elastic Net and Random Forest",   "Elastic net"),
 ("p8 '500-100 draws' typo",                "the figure is 500"),
 ("p8 pool preservation before BH",         "pooled across draws before multiplicity"),
 ("p9 remove borderline-personality ref",   None),   # must be ABSENT
 ("p9 BayesPrism never cited",              "Chu, T., Wang, Z."),
]

ABSENT = {"p3 PPIXpress is the wrong tool": ["PPIXpress", "Thorsten", "Will & Helms"],
          "p9 remove borderline-personality ref": ["Arranz", "borderline personality"]}

ok = fail = 0
for label, needle in CHECKS:
    if needle is None:
        missing = [s for s in ABSENT[label] if s in T]
        good = not missing
        detail = '' if good else f'  <-- still present: {missing}'
    else:
        good = needle in T
        detail = '' if good else f'  <-- MISSING: {needle!r}'
    print(f"{'PASS' if good else 'FAIL'}  {label}{detail}")
    ok, fail = (ok + 1, fail) if good else (ok, fail + 1)

print(f"\n{ok} passed, {fail} failed, {len(CHECKS)} checks")
print(f"word count: ~{len(T.split()):,}")
