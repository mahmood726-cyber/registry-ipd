# Production gallery — the method on real ClinicalTrials.gov trials

*The reconstruction running on the **real AACT trials it is built for** (not the open-IPD
validation set). Each row: a real trial harvested from the 2026-06-01 AACT snapshot, reconstructed
with the shipped engine (Titman-QP default), its reconstructed Cox HR scored against the
**registry-reported HR** (a coarse held-out truth), with full NCT provenance. Reproduce:
`node validate/gallery.js`. (Numbers from `realipd/gallery_results.json`.)*

## Cohort-wide

- **250** Tier-A trials reconstructed end-to-end; **30** also report a Cox HR.
- Median reconstructed-vs-registry HR fold-error: **1.13**.
- Self-audit badges: 23 gold · 227 silver.

## Expanded held-out HR validation (sibling-outcome fix)

The census finds **112** validation-grade trials (curve **+** HR; see `POLICY.md`), but this gallery's
HR scoring originally covered only **30** — because the harvester scoped the HR lookup to the *curve's
own outcome* and silently dropped HRs posted in a **sibling survival outcome** (the common "OS rate over
time" curve + separate "Overall Survival" HR layout). Fixing that (`harvester.select_trial_hr`, with a
one-pass `harvest/backfill_validation_hr.py` over the validation-grade set) recovers those HRs and
roughly **doubles** the held-out validation, scored by HR source (`node validate/gallery_expanded.js`):

| HR source | trials scored | median fold | within 1.2× |
|---|---:|---:|---:|
| **curve outcome** (the original basis) | 25 | **1.119** | 17/25 |
| **same-endpoint sibling outcome** (newly recovered) | 32 | 1.190 | 18/32 |
| **all validation-grade scored** | **57** | **1.149** | 35/57 |

The sibling fallback is **endpoint-aware**: an OS curve is never scored against a PFS hazard ratio —
explicit endpoint mismatches (8 of them) are *dropped*, not validated, so the sibling rows are a genuine
like-for-like check (fold **1.19**, vs the contaminated 1.205 before endpoint-matching). The curve-sourced
subset reproduces the original **1.13** (consistency check). Either way the production held-out HR
evidence base goes from 30 to **57** real, endpoint-clean trials. Numbers from
`realipd/gallery_expanded.json`; the recovery's provenance counts (curve / same-endpoint sibling /
endpoint-mismatch dropped) are in `realipd/validation_hr_backfill.json`.

## Worked examples (diverse conditions, best-fit per condition)

| NCT | condition | N exp/ctl | events | anchors | badge | registry HR | reconstructed HR | fold |
|---|---|---|---|---|---|---|---|---|
| [NCT00878709](https://clinicaltrials.gov/study/NCT00878709) | Breast Cancer | 1420/1420 | 316 | 5 | gold | 0.952 | 0.959 | 1.01 |
| [NCT04305496](https://clinicaltrials.gov/study/NCT04305496) | Locally Advanced (Inoperable) or Metastatic Br | 355/353 | 411 | 3 | silver | 0.6 | 0.597 | 1.01 |
| [NCT00207142](https://clinicaltrials.gov/study/NCT00207142) | HIV Infections | 87/87 | 30 | 6 | silver | 0.97 | 1 | 1.03 |
| [NCT01151137](https://clinicaltrials.gov/study/NCT01151137) | Atrial Fibrillation | 1619/1617 | 41 | 6 | silver | 2.294 | 2.377 | 1.04 |
| [NCT01983683](https://clinicaltrials.gov/study/NCT01983683) | Clostridium Difficile Infection | 290/290 | 152 | 10 | silver | 1.04 | 1 | 1.04 |
| [NCT01286272](https://clinicaltrials.gov/study/NCT01286272) | Ann Arbor Stage III Grade 1 Follicular Lymphom | 66/66 | 56 | 4 | silver | 0.94 | 1 | 1.06 |
| [NCT00185211](https://clinicaltrials.gov/study/NCT00185211) | Multiple Sclerosis | 292/176 | 124 | 3 | silver | 0.764 | 0.827 | 1.08 |
| [NCT01335399](https://clinicaltrials.gov/study/NCT01335399) | Multiple Myeloma | 374/374 | 286 | 5 | silver | 0.93 | 1 | 1.08 |
| [NCT02485899](https://clinicaltrials.gov/study/NCT02485899) | Jansky-Bielschowsky Disease | 23/42 | 54 | 6 | silver | 0.14 | 0.151 | 1.08 |
| [NCT00261443](https://clinicaltrials.gov/study/NCT00261443) | Bipolar Disorder | 168/169 | 33 | 14 | silver | 0.348 | 0.31 | 1.12 |

*The HR is the hard estimand and the registry HR is itself coarse; read these as a triangulation
check, not bit-exact recovery. RMST/median reconstruct far more tightly (see `VALIDATION.md`).*
