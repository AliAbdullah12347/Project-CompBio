# Scale-free topology fails at every power and both sample sizes

**Date:** 2026-10-09
**Analysis:** `network/runs/base-001` (sft curve) + Phase 0 probes
**Evidence:** 15/25 · **Interest:** 7/10 · **Δ:** +1.0

## Claim
Signed R² never reaches 0.80 at any power 1-20, in any of the three groups, at n=74 or
n=234. Max 0.547 at power 20 (n=234 control, all 12,368 genes). `powerEstimate` returns
1 - the function failing, not a threshold. At power 12: R² 0.162, mean connectivity 1,031
(8.3% of the genome).

Not a small-sample artifact: replicates at 3x the sample size.

E=4 R=2 B=3 D=4 C=2. Interest 7: WGCNA's central assumption fails and the function
reports a plausible-looking wrong number rather than failing loudly.
