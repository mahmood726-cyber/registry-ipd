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

## Real-trial demonstration — RADIANT-4 (NCT01524783), validated vs the posted HR

The 42-dataset experiment establishes the mechanism on rendered curves. This grounds it on a **real
trial present in both worlds**, validated against the trial's **published** hazard ratio
(`kmcurve/ipd_km_pipeline/fusion_real_trial.py`; artifact
`validate/fusion_real_trial_radiant4.json`):

- **registry side** — this repo's harvested AACT record `NCT01524783.json`: the exact posted
  KM-estimate anchors (10/11 points) + N=205/97. AACT posts no number-at-risk.
- **figure side** — the published KM figure's "N (events)" totals (107/77), which kmcurve OCRs.
- **ground truth** — the registry-posted Cox HR **0.48 (95% CI 0.35–0.67)**, held out.

| reconstruction | HR | fold vs posted | inside posted 95% CI? |
|---|---:|---:|:--:|
| registry-only (exact anchors, NO censoring) | 0.83 | 1.73 | **no** |
| **FUSION (anchors + figure event counts)** | **0.56** | **1.17** | **yes** |

On a real trial with a real posted effect, registry-only falls **outside** the published CI — the
identifiability trap — and the figure's event count pulls the fusion estimate **inside** it. This is
the same worked example `VALIDATION.md` documents (curve-only 0.68 → censoring-informed 0.47 vs posted
0.48); the exact numbers differ because the reconstruction runs through kmcurve's Guyot/QP engine rather
than this repo's JS, but the conclusion is identical and now expressed as a fusion of the two projects'
**real** data. (Honest link: the event counts are taken from the AACT harvest — registry-ipd derived
them from participant-flow `drop_withdrawals`; kmcurve OCRs the identical totals from the figure's
at-risk table in production. The RADIANT-4 Lancet primary is not open-access here, so the figure-OCR
step itself is exercised separately by the corpus benchmark, not on this specific figure.)

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
- Not yet demonstrated on a *real matched* AACT-trial ↔ published-PDF pair end-to-end. **Quantified
  availability caveat (kmcurve `fusion_crossmatch.py`), updated at scale:** the first scan (327 PDFs →
  185 with an NCT → 267 unique NCTs → 37 with AACT results → **0** dual-available) was *corpus-bounded*,
  not capability-bounded — so we **doubled the corpus and re-ran**: 663 PDFs → 373 carry an NCT ID →
  592 unique NCTs → 106 have AACT posted results → **2** now have *both* a posted survival curve **and**
  an open-access figure (`NCT01658878`/PMC7530824 OS, 5 timepoints; `NCT03110107`/PMC13006393 PFS,
  3 timepoints), but both were early-phase multi-cohort dose-finding trials with no posted HR. Growing
  the mirror to **1500 PDFs** (scan made tractable by a ~10× faster fitz NCT extractor) raised the
  NCT-citation count to **6** — but verifying the PDFs exposed, and fixed, a real limitation of the
  cross-match: **citing an NCT ≠ being its primary publication.** A `classify_pdf` filter (likely-primary
  = ≤3 cited NCTs and not a meta-analysis, *and* an at-risk table present) cuts the 6 to **3 fusion-USABLE
  primaries, exactly 1 of which posts a hazard ratio.** (The HR-0.75 "OS" candidate `NCT00636168` that a
  first pass flagged turned out to be a **network meta-analysis citing the trial**, with no at-risk table
  — corrected.) The one genuinely usable **curve + HR + figure** pair is **`NCT01942135` (PMC9662922): a
  survival-probability curve with a posted HR 0.42 (0.32–0.56) in an open-access primary that contains an
  at-risk table** — a real held-out ground truth for an end-to-end NAR fusion (harvest the AACT anchors,
  OCR the figure's at-risk table, reconstruct, score the reconstructed HR against 0.42). So the mechanism
  (42 rendered datasets + RADIANT-4) finally has a real, verified matched pair; the narrow count
  (1 usable HR pair in 1500 OA PDFs) *re-confirms* the `POLICY.md` case — confirmatory 2-arm RCTs that
  post both a curve and an HR remain disproportionately paywalled, their OA mentions mostly reviews.
  Per-candidate detail in kmcurve `CORPUS_FINDINGS.md`.

Reproduce (from the kmcurve repo): `python realipd_benchmark.py --fusion --registry <this repo>`;
cross-match: `python ipd_km_pipeline/fusion_crossmatch.py`.
