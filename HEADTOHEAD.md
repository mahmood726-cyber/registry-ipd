# Head-to-head: registry-structured vs figure-digitization (real IPD)

*A self-cutting validation. The paper's value proposition includes "zero figure-digitization error."
This script tests that on real data and reports a result that partly tempers the claim — which is the
honest outcome and, as it turns out, the strongest argument for the policy brief. Reproduce:
`node validate/headtohead.js` → `realipd/headtohead_results.json` (25 open-IPD datasets).*

## Setup

For each open dataset (true patient-level IPD), the same true Kaplan–Meier curve is turned into two
posted-data scenarios and both are reconstructed through the **same** engine and scored against truth
(Cox HR, median, RMST):

- **REGISTRY** — exact KM survival at the posted anchor timepoints (what ClinicalTrials.gov/AACT
  exposes). No pixel error.
- **DIGITIZED** — the same curve read off a *plotted figure*: coordinates sampled along the curve with
  realistic pixel noise (survival σ = 1 percentage point, time σ = 0.5% of t_max), monotone-enforced.
  This is fair-to-generous to digitization: it traces the exact published curve.

Two comparisons: **(A) equal density** (K=8 points each — isolates pure digitization-error cost) and
**(B) realistic** (registry few-exact K=8 vs digitizer many-noisy K=25).

## Result

| Comparison | Registry (exact) median \|logHR err\| | Digitized (noisy) median \|logHR err\| | Registry ≤ digitized |
|---|---:|---:|---:|
| **(A) Equal density, K=8** | **0.148** | 0.139 | **17 / 25** |
| **(B) Realistic, 8 exact vs 25 noisy** | 0.148 | **0.041** | 4 / 25 |

**Two honest takeaways:**

1. **Digitization noise itself costs little.** At equal density (A), exact anchors and ~1-pp-noised
   anchors are a near-tie (registry wins the majority of head-to-heads, 17/25, but the medians are
   within 0.01 on the log-HR scale). The registry's "zero digitization error" advantage is **real but
   small**.

2. **Anchor density dominates.** When a digitizer reads many points off the figure (25) against a
   registry's few exact anchors (8), the **digitizer reconstructs the HR markedly better** (median
   |logHR| 0.041 vs 0.148). Curve-*shape* density beats per-point exactness for hazard-ratio recovery.

## What this means (and why it strengthens, not weakens, the case)

The registry path's binding weakness is **anchor sparsity, not pixel noise.** That is exactly the
finding the coverage census already quantified (median posted curve = 3–4 timepoints; only ~34% reach
≥6) and exactly what `POLICY.md` asks registries to fix. So the honest positioning is:

> Registry-native reconstruction eliminates digitization error and carries full provenance, but its
> accuracy is **anchor-density-limited**. It is competitive with — and cleaner than — figure
> digitization **only when enough timepoints are posted** (the ≥5–6 standard). Below that, a densely
> digitized figure recovers the HR better. This makes the case for the reporting recommendation
> directly: the value of the registry path is unlocked by posting more KM timepoints, not by the
> method alone.

It also means the two approaches are **complementary**: where a publication figure exists, digitize it
(dense shape); where only the registry record exists, reconstruct natively (exact, provenanced) — and
push registries toward the anchor density that makes the native path competitive.

## Robustness to the noise level (sensitivity sweep)

`validate/noise_sweep.js` varies the digitization noise σ from 0 to 5 percentage points (the dashboard
plots it). Registry-exact is the flat noise-free reference (median |logHR| 0.148). The two findings are
robust:

| σ (survival, pp) | registry exact | digitized K=8 (equal) | digitized K=25 |
|---:|---:|---:|---:|
| 0   | 0.148 | 0.148 | 0.056 |
| 1   | 0.148 | 0.139 | 0.041 |
| 2   | 0.148 | 0.168 | 0.062 |
| 3   | 0.148 | 0.194 | 0.038 |
| 5   | 0.148 | 0.147 | 0.099 |

Equal-density digitized *crosses above* registry once σ exceeds ~1 pp (digitization noise does cost at
equal density), but **K=25 digitized stays well below registry across the whole 0–5 pp range** — anchor
density dominates per-point noise for HR recovery. The head-to-head conclusion is not an artefact of
the chosen noise level.

*Caveat: the noise model is deterministic and conservative. RMST and median (curve-derived) are
recovered well by both paths; the HR — which depends on event/censoring timing — is where density
matters most. Small trials (e.g. `gehan`, N≈21/arm) reconstruct poorly under both paths, consistent
with the small-N boundary.*
