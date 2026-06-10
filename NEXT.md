# NEXT — handoff for a fresh session

Repo is complete, validated, green, and live. This note lets a new session resume cleanly on the
three "elevate toward breakthrough" items without re-reading the whole history.

## Current state (as of commit 809e4d5)
- **Engine** (`src/engine.js`): curve-only · censoring-informed (N-matched) · HR-calibration ·
  max-entropy multiple-imputation ensemble (calibrated UQ) · competing-risks (Aalen–Johansen) ·
  Royston–Parmar · fractional-polynomial time-varying HR. 24 JS tests.
- **Harvester** (`harvest/`): AACT → trial JSON; orientation; N-matched flow mapping. 19 Python tests.
- **Validation** (`validate/`, full results in `VALIDATION.md`): registry-internal → published →
  **true IPD across 25 open datasets** (15 adequate; HR ~11% median, 12/15 within 20%; median ~3%;
  RMST ~2%) → uncertainty **CI covers true HR 14/14** → competing-risks gold standard (colon + aidssi,
  naive 1−KM bias up to 16pp corrected) → **anchor-density operating curve** (≥5–6 timepoints).
- **Artifacts**: offline tool + advanced panel (`index.html`/`dist/`), `validation-dashboard.html`,
  `METHODS.md`, `PAPER.md` (PubMed-verified core DOIs). Live: https://mahmood726-cyber.github.io/registry-ipd/
- **Data note**: `realipd/` (open datasets) and `cohort/` (AACT harvest) are gitignored; re-download
  command in `validate/goldstandard.js` header. AACT snapshot env: `AACT_ZIP`.

## The three next items (the user's "1, 2, 3")
1. ✅ **DONE (2026-06-10).** Systematic novelty search → `NOVELTY.md`. Verdict: registry-native
   AACT-sourced reconstruction is **defensible as first**, hedged against **Titman 2026** (Stat Med,
   PMID 41775249) — closest concurrent *tabular-input* method, but journal tables not registry. Every
   reconstruction tool in the literature is figure-image based. Related Work paragraph added to
   `PAPER.md` (claims novelty on data *provenance*, not file format).
2. ✅ **DONE (2026-06-10).** Registry-policy brief → `POLICY.md`. 1-pager for CTTI/AACT + EU CTIS:
   post structured number-at-risk + total events + ≥5–6 KM timepoints. Anchored on the measured
   0-NAR finding and the reliability threshold.
   - Plus: **coverage census** quantifying the binding limitation precisely (was "hundreds"):
     `harvest/census_full_aact.py` (full 76,067-trial snapshot) + `validate/census_cohort.js` (605
     real trials). 0 structured NAR rows; 288–~605 reconstructable curves (0.4–0.8%); only ~34%
     clear ≥6 timepoints. Wired into `PAPER.md` abstract + §2 and `README.md`.
   - Plus: **all citations PubMed-verified** → `CITATIONS.md` (8/8 indexable refs exact).
3. **Credentialed large-scale true-IPD validation** — STILL OPEN, *needs the user's credentialed
   access* (Vivli / Project Data Sphere / YODA). Blocked until a data export is provided; then plug
   CSVs into `validate/goldstandard.js` configs to extend the gold standard to dozens–hundreds of
   trials. This is the one remaining lever and the only non-autonomous item.

## Also done this session (beyond items 1+2)
- ✅ **Real head-to-head** `HEADTOHEAD.md` / `validate/headtohead.js`: registry-structured vs
  figure-digitization on 25 open-IPD datasets. Honest finding — at equal density digitization noise
  costs little; a dense digitized figure beats few exact registry anchors on HR; the registry path's
  binding weakness is **anchor sparsity, not noise**. Tempered into `PAPER.md` Discussion; reinforces
  the policy brief. (Refactored `goldstandard.js` to export `loadArms`/`coxHR`, output byte-identical.)
- ✅ **All citations PubMed-verified** (`CITATIONS.md`, 8/8 indexable exact).

- ✅ **Dashboard panels** (`validation-dashboard.html`): coverage census (anchor-density histogram +
  headline stats), head-to-head grouped bars, and a digitization-noise sensitivity line chart. All
  offline SVG; headless-smoke verified (3 panels render, 0 severe console errors).
- ✅ **Noise sweep** (`validate/noise_sweep.js` → `noise_sweep_results.json`): σ 0–5 pp; confirms the
  head-to-head ordering is robust (K=25 digitized beats registry across the whole range).

- ✅ **Gold standard expanded 25 → 31** (2026-06-10) via open packaged datasets (KMsurv/asaur):
  alloauto, larynx, burn, pneumon, bfeed, hepatoCellular. Median HR fold-error holds at 1.12;
  bfeed is the documented worst case (discrete-time). Configs in `goldstandard.js`; CSVs gitignored
  (download paths in the file header).

- ✅ **Uncertainty re-run on expanded set**: 14/14 → 16/17 (94%), bfeed the miss (2026-06-10).
- ✅ **cBioPortal/TCGA integration** (`harvest/fetch_cbioportal.js`): 8 real cancer OS cohorts pulled
  from the open API. Sex-split is near-null (true HR 0.7–1.2), so reported as a SEPARATE slice scored
  by absolute log-HR error (`goldstandard_cbio.js`, mean ~0.18), kept out of the headline aggregate +
  its own dashboard panel. Honest characterisation of the low-signal regime.

## What's autonomously left (optional, in priority order)
- **Strong-contrast cancer IPD**: TCGA pan-can stage/metastasis attributes are sparsely populated, so
  the cBioPortal cohorts ended up near-null. To get strong real oncology contrasts (large HR), use
  non-pan-can `*_tcga`/`*_tcga_pub` studies with clinical stage, or Dryad/Zenodo trial deposits
  (CC0/CC-BY). NOT credentialed repos (Vivli/YODA/PDS need a DUA — do not scrape).
- **Re-run competing-risks** on any new datasets that have competing-event structure (the current
  new additions don't).
- Surface the census/head-to-head panels inside the live tool (`index.html`), not just the dashboard.

## Honest framing (carry forward)
Novel registry-native data path + calibrated uncertainty + unusually rigorous validation; bounded by
**coverage** (a minority of trials post a structured curve) → it complements, not replaces,
figure-digitization. "Useful, validated, possibly-first method + policy-relevant findings", not a
universal IPD solution. Keep the no-overclaim discipline.

## Resume commands
```
node --test test/engine.spec.js validate/metrics.spec.js   # 24 JS tests
python -m pytest harvest/ -q                                # 19 Python tests
node validate/goldstandard.js                               # true-IPD gold standard (needs realipd/*.csv)
```
