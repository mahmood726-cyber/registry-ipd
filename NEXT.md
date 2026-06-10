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
1. **Systematic novelty search** — confirm "first to reconstruct IPD from AACT *structured* data"
   (not figure-digitization). PubMed + arXiv + medRxiv + Google Scholar/grey lit. *Fully autonomous.*
   Mainstream check already done (IPDfromKM/KM-GPT/RESOLVE-IPD are all figure-based); needs a thorough
   systematic sweep before any "first" claim. Output: a short novelty memo + a related-work section.
2. **Registry-policy brief** — turn the "AACT has 0 structured number-at-risk" finding + the
   anchor-density ≥5–6-timepoint standard into a concrete reporting recommendation (1-pager aimed at
   CTTI/AACT + EU CTIS). *Autonomous to draft; adoption is external.*
3. **Credentialed large-scale true-IPD validation** — Vivli / Project Data Sphere / YODA exports to
   extend the gold standard to dozens–hundreds of trials. *Needs the user's credentialed access* —
   blocked until a data export is provided; then plug CSVs into `validate/goldstandard.js` configs.

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
