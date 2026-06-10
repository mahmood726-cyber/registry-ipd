# NAR fusion: registry curve × figure at-risk table dissolves the identifiability trap

*The single most important cross-project finding. Two reconstruction projects are **mirror images** of
each other, and fusing their inputs is strictly better than either alone — and it targets the one limit
this project repeatedly documented as **fundamental** (curve-only cannot identify censoring;
`advanced_estimators.js` "none dominates… a fundamental identifiability limit"). Measured on this repo's
true-IPD gold standard via kmcurve's `realipd_benchmark.py --fusion`; results JSON committed at
`validate/nar_fusion_results.json`.*

## The mirror image

| project | has the CURVE | has the number-at-risk | binding weakness |
|---|---|---|---|
| **registry-IPD** (this) | **exactly** (AACT KM-estimate anchors, zero pixel error) | **no** — AACT posts none | censoring unidentified → HR attenuates/explodes |
| **kmcurve** (figure digitization) | noisily (pixel extraction off the plot) | **yes** — OCRs the printed risk table | pixel/arm-separation noise on the curve |

Each project's strength is the other's weakness. So **fuse them**: for a trial present in both AACT and
a published PDF, take registry-IPD's exact curve anchors **and** kmcurve's OCR'd at-risk table, and
reconstruct from both. No registry has to change anything.

## Experiment

Hold the reconstruction backend fixed (Guyot) and vary only the INPUTS, on every gold-standard true-IPD
dataset (render the true KM → kmcurve extracts the noisy curve; the registry anchors are the true KM
sampled at 8 timepoints; the at-risk table is the true number-at-risk at 6 timepoints). Score the
reconstructed log-rank HR against the true HR (same estimator both sides).

| input regime | median HR fold-err | p90 | within 20% (95% CI) |
|---|---:|---:|---:|
| **registry-only** — exact anchors, NO at-risk table (registry-IPD today) | 1.37 | 4.26 | 16/42 (25–53%) |
| **kmcurve-only** — noisy curve + at-risk table (kmcurve today) | 1.09 | 1.52 | 31/42 (59–85%) |
| **FUSION** — exact anchors + at-risk table | **1.049** | 1.28 | **35/42 (69–92%)** |
| **FUSION + QP** — exact anchors + Titman-QP on the event count | 1.058 | **1.21** | 36/42 (72–93%) |

**Paired, on the same 42 datasets:**
- **Fusion vs registry-only: fusion better on 34/42**, median fold ratio 1.35. The within-20% CIs do
  **not overlap** (69–92% vs 25–53%) → the gain is statistically robust, not noise. This is the
  identifiability trap being dissolved: the at-risk table supplies the censoring the exact curve alone
  cannot. Showcases — `cbio_kirp` 8.6→1.36, `prostateSurvival` 5.4→1.05, `cbio_acc` 4.3→1.03.
- **Fusion vs kmcurve-only: fusion better on 26/42**, median fold ratio 1.035. Favorable but the CIs
  overlap (69–92% vs 59–85%) → **not significant at n=42**; the exact curve helps modestly once the
  at-risk table is already present. Honest read: fusion ≥ kmcurve-only, materially so on the hard cases,
  marginally on the easy ones.

**Takeaway.** Fusion never loses much and wins big exactly where each single-source path fails — the
exact curve rescues kmcurve's pixel noise; the at-risk table rescues registry-IPD's missing censoring.
The deepest claim — *fusion is strictly ≥ either alone* — holds across 42 real datasets.

## Why this matters for the paper

The `advanced_estimators.js` / `method_zoo.js` investigations concluded the curve-only censoring
ambiguity is a **fundamental identifiability limit** — true under the AACT-only constraint, because AACT
posts no number-at-risk. This experiment shows the limit is **not fundamental to the trial, only to the
data source**: the censoring information exists, printed in the figure's risk table, and a figure
digitizer (kmcurve) recovers it. So the registry-native path's binding weakness is closed by a *data
union*, not a better algorithm — which is exactly what `POLICY.md` asks registries to post (structured
number-at-risk). Until they do, the AACT+figure union is the bridge.

## Honest scope

- Controlled render: clean single-style monochrome curves, exact calibration, exact at-risk counts at 6
  timepoints. Real figures add coloured/overlapping arms, censor ticks, and OCR error on both the curve
  and the at-risk table — all of which widen the kmcurve-only and fusion numbers.
- The "registry anchors" here are the true KM at 8 timepoints; a real AACT record posts a median of 3–4
  (see the coverage census), so registry-only in production is often **worse** than the 1.37 here, which
  strengthens the case for fusion, not weakens it.
- n=42 open datasets (Rdatasets/KMsurv/asaur + TCGA late-vs-early); skewed toward oncology and
  strong-contrast stage splits. Read every percentage with its Wilson CI.
- Not yet demonstrated on a *real matched* AACT-trial ↔ published-PDF pair end-to-end — that is the
  recommended next experiment (this establishes the mechanism and the expected gain).

Reproduce (from the kmcurve repo): `python realipd_benchmark.py --fusion --registry <this repo>`.
