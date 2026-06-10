# Registry-Native Pseudo-IPD Reconstructor

Reconstruct **pseudo** individual-patient time-to-event data from **ClinicalTrials.gov / AACT
summary data only**, fully offline in the browser. No R, no Shiny, no server, no CDN.

---

## Read this first — the honest scope

> **AACT contains no Kaplan–Meier curves and no true individual-patient data — only summary
> statistics.** Every existing reconstruction tool (IPDfromKM, KM-GPT, RESOLVE-IPD,
> Guyot-on-image) starts from a *published curve image* and digitizes it, injecting pixel error.
> Under a ct.gov-only constraint there is no image to digitize, so:
>
> - Output is always **pseudo-IPD**, never true IPD, and **never produced for Tier C**.
> - The edge over digitization tools is **exact registry anchors + zero digitization error +
>   full provenance** on the subset of trials with rich registry survival data — **not** universal
>   superiority. A dense digitized curve carries more between-anchor shape than sparse registry
>   anchors; where the registry is sparse, this tool degrades explicitly rather than fabricating.

The scoped claim we stand behind: *on trials where ct.gov posts a structured KM-estimate curve,
AACT-only reconstruction uses the exact registry anchors (zero digitization error on the anchors)
with full provenance; elsewhere it degrades explicitly via tiered verdicts.* See **`VALIDATION.md`**
for what is and isn't validated — in short: HR recovery ~12% median / ~83–94% within the registry CI
(n≈18, not significant between methods); RMST/median "fidelity" numbers are **round-trip
self-consistency, not external accuracy**; a real curve-vs-digitization head-to-head has not been run
(the only such figure is a self-graded synthetic benchmark).

## How it works (two pieces)

```
AACT snapshot ──(harvester, Python + aact-kit)──▶ trial.json ──(engine, offline JS)──▶ pseudo-IPD + badge
```

The split keeps **both** invariants true at once: the data is *only* from ct.gov/AACT (the
harvester reads nothing else), and the reconstructor is *fully offline* (the browser engine never
fetches — it reads the harvested JSON).

### Tiered reconstruction (richness-aware)

| Tier | Registry data available | Method | Output |
|------|-------------------------|--------|--------|
| **A** rich | KM-estimate ≥3 timepoints + number-at-risk ≥2 + total events | **best-of two methods** (see below), selected by min 1-Wasserstein to anchors | deterministic pseudo-IPD |
| **B** medium | median + hazard ratio + N + events | **parametric** (exponential), event-count conditioned, **bootstrap envelope** | pseudo-IPD + uncertainty |
| **C** sparse | hazard ratio only | **fail closed** — reconstruction refused | analytic summary only |

### Statistical methods (not just Guyot)

Tier A runs **two reconstructions and keeps the one that fits the registry curve best**, measured by
the **1-Wasserstein (L1-area) distance** between the reconstructed KM and the registry KM-estimate
step function:

- **Guyot (2012) inverse-KM** — classic iterative reconstruction assuming constant censoring within
  number-at-risk blocks. Strong when censoring is spread through follow-up.
- **Censoring-informed / anchor-exact** (RESOLVE-IPD CEN-KM style, 2025) — deaths taken directly
  from the curve with at-risk held constant within each interval and censoring placed only at
  NAR boundaries, so the reconstructed KM passes through every registry anchor **≈ exactly**.
  Correct when censoring is administrative (concentrated at the cutoff) — the case Guyot mishandles.

Wasserstein is used as *both* the model selector and the reported fidelity metric. Force a single
method with `reconstruct(t,{method:'guyot'|'anchor-exact'})`.

> **On the synthetic benchmark** (`validate/gen_benchmark.js`): the censoring-informed method cut
> mean anchor sup-error 0.127→0.011 and the "AACT-only beat simulated digitization 20/20" figure is
> a **self-graded methodology demonstration** — the "digitized" comparator is the same anchors plus
> injected noise, both reconstructed by the same engine and graded against the anchors, so AACT-only
> (given the exact anchors) wins *by construction*. It illustrates the mechanism; it is **not** a
> real-trial head-to-head (that comparison has not been run on actual digitized publications).

> Upgrade path (documented, not yet shipped): **Royston–Parmar flexible parametric (spline) models**
> for Tier B in place of the exponential, and a Wasserstein-barycenter ensemble across methods.

### Self-audit → Bronze / Silver / Gold

Reconstructed pseudo-IPD is checked **back against the registry anchors**: total-event match (C1),
anchor survival fidelity (C2), median within 5% (C3), reconstructed-HR vs registry HR via Cox (C4),
monotonicity (C5, hard), number-at-risk consistency (C6), population conservation (C7, hard),
follow-up sanity (C8), and **HR direction integrity (C9, hard — the HR is never inverted, only arm
labels are resolved)**. Any hard-check failure or Tier C ⇒ badge `none` and export is blocked.

## What AACT actually contains (measured, not assumed)

Coverage scan over the **2026-06-01 AACT snapshot** (76,067 trials with posted results):

| | count | meaning |
|---|---|---|
| **Tier A** — KM curve (survival measure at ≥3 parseable timepoints) | **288** | reconstructable from the curve + N |
| **Tier B** — median + hazard ratio | **3,283** | parametric reconstruction |
| **Tier C** — sparser | 72,496 | refused |
| trials with any HR/Cox | 5,567 | |
| trials with a median | 18,761 | |
| **structured number-at-risk rows in *all* of AACT** | **0** | — |

Two findings that shaped the design:
1. **AACT has zero structured number-at-risk.** So Guyot's classic NAR-driven reconstruction cannot
   run as-published; the engine uses N + total events with tail censoring (IPDfromKM's "no
   number-at-risk" mode), and **Tier A does not require NAR**. The "exact-anchor beats digitization"
   edge is real but applies to a **small subset** (288 trials strictly; ≈514–605 under a broader
   survival-curve net) — we say so plainly.
2. **KM data is often stored as cumulative *incidence* ("probability of progression/event"), not
   survival.** The harvester detects orientation **data-drivenly** (a curve starting near 0 and
   rising is incidence ⇒ `S = 1 − value`), the robust form of the sign-flip lesson. Verified on a
   live reconstruction of **NCT00725985** (Tier A, Silver) — bundled in the tool as the "REAL ct.gov"
   example.

**Anchor-density reality (the reuse bottleneck).** Of the trials that *do* post a curve, the median
posts only **3–4 KM timepoints**, and only **≈34%** post the **≥5–6** our validation shows are needed
for reliable reconstruction (HR fold-error 1.40 at K=3 → 1.15 at K=5 → 1.08 by K=12). Reproduce the
full-snapshot census with `harvest/census_full_aact.py`; the per-arm detail over 605 real harvested
trials with `validate/census_cohort.js`. This is the evidence base for the reporting recommendation in
**`POLICY.md`** (two structured fields + a timepoint threshold for CTTI/AACT + EU CTIS).

A real-data **head-to-head** (`HEADTOHEAD.md`, `validate/headtohead.js`) tests the "no digitisation
error" claim honestly: at equal anchor density, exact-vs-pixel-noised is a near-tie, and a *densely*
digitised figure actually recovers the HR better than few exact registry anchors. The registry path's
binding weakness is **anchor sparsity, not noise** — so the two methods are complementary, and the
native path's value is unlocked precisely by the anchor-density standard above.

**Where this method sits vs prior art.** Every established KM-IPD reconstruction tool (Guyot 2012,
IPDfromKM, ipdfc, RESOLVE-IPD, KM-GPT) digitises a *figure image*; the closest tabular-input neighbour,
Titman 2026, uses *journal-article* tables. To our knowledge this is the first reconstruction native to
ClinicalTrials.gov/AACT **structured registry data** — novelty on data *provenance*, not file format.
Full systematic prior-art search and hedge wording in **`NOVELTY.md`**; every citation is
PubMed-verified in **`CITATIONS.md`**.

## Usage

### Reconstruct (offline, no install)
Open **`index.html`** (or the single-file **`dist/registry-ipd.html`**) in any browser. Load a
bundled example or upload a harvested `trial.json`, click **Reconstruct**, export pseudo-IPD CSV.
For Tier-A two-arm trials a **Run advanced analysis** button surfaces the cutting-edge methods:
multiple-imputation **credible intervals** (HR/median), **time-varying HR** with a non-PH check,
and the **Aalen–Johansen CIF** when competing-event counts are present.

### Harvest a real trial (needs a local AACT snapshot)
```bash
# 1. get an AACT pipe-delimited snapshot: https://aact.ctti-clinicaltrials.org/snapshots
#    then point aact-kit at it, e.g.:
export AACT_ZIP=/path/to/YYYYMMDD_pipe-delimited-export.zip   # PowerShell: $env:AACT_ZIP=...

# 2. feasibility: how much of AACT is actually reconstructable?
python harvest/coverage_scan.py -o coverage_report.json

# 3. harvest one trial -> trial JSON for the browser tool
python harvest/harvest_trial.py NCT01234567 -o NCT01234567.json
```

## Develop / test

```bash
node test/gen_fixtures.js          # regenerate self-consistent fixtures
node --test test/engine.spec.js    # 10 engine tests (Tier A Guyot parity, Tier B envelope, fail-closed, edges)
python -m pytest harvest/ -q       # 16 harvester parsing tests (no snapshot needed)
python test/smoke_browser.py       # headless-Chrome smoke (offline, all tiers)
node examples/gen_examples.js      # regenerate bundled examples
node build.js                      # inline -> dist/registry-ipd.html (single file)
```

### Test / tolerance contract
- **Tier A** (deterministic): exact total-event recovery; integer death-schedule parity vs a
  hand-specified Guyot round-trip; anchor fidelity ≤ 1e-3 on dense anchors.
- **Tier B** (stochastic): registry event count matched exactly via order-statistic cutoff;
  Monte-Carlo (≈3σ) bound on the central-draw median; registry median/logHR inside the bootstrap
  envelope. Seeded PRNG (`seed = hash(nct_id)`) for reproducibility.

## Reuse / provenance
Ports the Guyot inversion idea from `wasserstein/improved_guyot_algorithm.py`; trial-JSON harvest
uses `aact-kit`; the offline shell follows the `allmeta/km-reconstructor` single-file pattern. The
AACT-only structured-anchor reconstruction, the harvester's `outcome_measurements` mining, the
tiered engine, and the self-audit checks are net-new.

## Status
- Engine, harvester, HTML shell, single-file build, validation harness: **done & green**
  (14 JS tests + 16 Python tests + headless-Chrome smoke).
- **Live AACT run: done.** Coverage measured (table above); a real trial (NCT00725985) harvested
  and reconstructed end-to-end (Tier A, Silver) and bundled as an in-tool example.
- Remaining (optional next steps): scan the 288 Tier-A trials for a higher-event oncology OS demo;
  run the full AACT-only-vs-digitization head-to-head on real trials that also have a published
  curve; ship Royston–Parmar for Tier B.
