# Extending the gold standard with credentialed IPD (Vivli / YODA / Project Data Sphere)

**Roadmap item 3.** The open gold standard (45 datasets, `VALIDATION.md`) is bounded by what is openly
licensed. Credentialed repositories hold hundreds of real RCTs with patient-level time-to-event data.
This is the plug-in path: when you have a credentialed export, validating against it is **one command**
— no code edits. The adapter is `validate/ingest_ipd.js`; it runs the *same* pipeline as the open gold
standard (true HR/median/RMST → registry-style coarse summary → curve-only + Titman-QP censoring-informed
reconstruction → multiple-imputation uncertainty) and scores each trial against truth.

> **DUA discipline.** Credentialed IPD is data-use-agreement protected. It lives under
> `realipd/credentialed/`, which is **gitignored** — it must never be committed or pushed. Only the
> aggregate result numbers (median fold-error, coverage) are safe to quote; never commit per-patient
> rows or per-trial CSVs. Access the repositories through their own application process; do not scrape.

## 1. Folder layout

```
realipd/credentialed/
  <study-or-batch>/
    manifest.json        # describes each trial's CSV + column mapping (schema below)
    NCT01234567.csv      # one CSV per trial (the export you downloaded)
    NCT07654321.csv
    ...
```

## 2. Manifest schema (`manifest.json`)

```json
{
  "source": "vivli",                          // vivli | yoda | pds | mixed (label only)
  "note": "diabetes CV-outcome RCTs, 2026 export",
  "datasets": [
    {
      "id": "NCT01234567",
      "label": "DRUG-X CV outcomes",
      "source": "vivli",
      "format": "cdisc-adtte",                // CDISC ADaM time-to-event (most Vivli/YODA exports)
      "csv": "NCT01234567.csv",
      "arm": "TRTP", "exp": "Drug X", "ctl": "Placebo"
      // time defaults to AVAL, censor flag to CNSR — override with "time"/"censorFlag" if named differently
    },
    {
      "id": "PDS-9999",
      "label": "Trial from Project Data Sphere",
      "source": "pds",
      "format": "generic",                    // plain time/status columns
      "csv": "PDS-9999.csv",
      "time": "OS_DAYS", "status": "OS_EVENT", "eventVal": 1,
      "arm": "ARM", "exp": "1", "ctl": "0"
    }
  ]
}
```

**Per-dataset fields.** `format` is `cdisc-adtte` or `generic`.
- **`cdisc-adtte`** — `time` (default `AVAL`), `arm` (default `TRTP`/`TRT01P`), and the **censor flag**
  `censorFlag` (default `CNSR`). ⚠️ **CDISC `CNSR` is inverted: `1` = censored, `0` = event** — the
  adapter flips it for you. This is the single most common ingestion bug for Vivli/YODA exports; getting
  it wrong silently inverts every HR.
- **`generic`** — `time`, `status` (1 = event), optional `eventVal` (treat only this status value as the
  event), `arm`. Same shape as `validate/goldstandard.js` configs.
- Both: `exp` / `ctl` are the two arm-column values to compare (experimental vs comparator). Optional
  `K` (posted timepoints, default 8), `time_unit`.

The two arms must each have **≥20 rows** (≥100/arm to enter the headline aggregate and the uncertainty
coverage check).

## 3. Run

```
node validate/ingest_ipd.js realipd/credentialed/<batch>/manifest.json
```

It prints a summary and writes `ingest_results.json` next to the manifest:

```
curve_only_median_fold / within20      # the conservative (no event count) floor
censoring_informed_qp_median_fold      # the Titman-QP result when total_events is present
uncertainty_coverage                   # 95% credible interval covers the true HR, on >=100/arm trials
per_dataset[...]                        # true_HR, curve/QP fold-error, CI coverage per trial
```

These merge directly with `VALIDATION.md`: report the credentialed median fold-error and coverage
alongside the open 45-dataset figures (curve-only ~1.12, QP ~1.05, coverage 24/25).

## 4. Validated end-to-end (synthetic CDISC mock)

The pipeline is proven on a synthetic CDISC-ADTTE fixture (`realipd/credentialed/_mock/`, gitignored —
regenerate from the generator in the commit that added this file). Three mock trials with design HRs
2.0 / 3.5 / 1.6:

| trial | N exp/ctl | true HR | curve-only fold | QP fold | CI covers true |
|---|---|---|---|---|---|
| MOCK-001 | 256/260 | 1.78 | 1.11 | **1.03** | yes |
| MOCK-002 | 173/179 | 3.96 | 1.06 | **1.03** | yes |
| MOCK-003 | 399/399 | 1.67 | 1.07 | **1.00** | yes |

Aggregate: curve-only 1.067 (3/3 within 20%), **QP 1.025 (3/3)**, uncertainty **3/3**. The sensible true
HRs confirm the `CNSR` inversion is handled — a mishandled flag would invert or destroy every HR.

## 5. What you need to provide

Only the export + a manifest. Drop a credentialed batch under `realipd/credentialed/<batch>/`, write the
`manifest.json` mapping (copy the schema above), and run the one command. The results extend the gold
standard from 45 open datasets toward the dozens–hundreds the credentialed repositories hold — the only
non-autonomous lever in the roadmap.
