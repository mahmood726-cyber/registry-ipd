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
