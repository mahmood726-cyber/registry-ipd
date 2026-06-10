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

## What's autonomously left (smaller polish, optional)
- Surface the coverage census + head-to-head as panels in `validation-dashboard.html` / the live tool.
- Sensitivity of the head-to-head to the digitization noise level (currently fixed σ; the ordering is
  robust to doubling, but a sweep would document it).

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
