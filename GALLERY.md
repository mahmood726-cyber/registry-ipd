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

| HR source | trials scored | median fold | within 1.2× | **recon HR in registry 95% CI** |
|---|---:|---:|---:|---:|
| **curve outcome** (the original basis) | 25 | **1.119** | 17/25 | **24/25 (96%)** |
| **same-endpoint sibling outcome** (newly recovered) | 32 | 1.190 | 18/32 | 23/32 (72%) |
| **all validation-grade scored** | **57** | **1.149** | 35/57 | **47/57 (82%)** |

The sibling fallback is **endpoint-aware**: an OS curve is never scored against a PFS hazard ratio —
explicit endpoint mismatches (8 of them) are *dropped*, not validated, so the sibling rows are a genuine
like-for-like check (fold **1.19**, vs the contaminated 1.205 before endpoint-matching). The curve-sourced
subset reproduces the original **1.13** (consistency check).

The strongest column is the last one. Rather than treat the registry HR as exact, we ask whether the
reconstructed HR falls inside the registry's **own posted 95% CI** — i.e. whether reconstruction error is
smaller than the trial's sampling uncertainty. For curve-sourced trials it is **24/25 (96%)**; even
sibling-sourced (a different outcome record) is 72%, and overall **47/57 (82%)**. The production held-out
HR evidence base goes from 30 to **57** real, endpoint-clean trials. Numbers from
`realipd/gallery_expanded.json`; recovery provenance (curve / same-endpoint sibling / endpoint-mismatch
dropped) in `realipd/validation_hr_backfill.json`.

## Independent validation: reconstructed HR vs the *published* HR

The registry HR shares provenance with the posted curve (same sponsor submission). The trial's
**published HR**, parsed from its primary paper's PubMed abstract, is an **independent** source. For each
validation-grade trial we get its PMID from AACT `study_references`, efetch the abstract (NCBI
E-utilities), and extract the published HR with a deterministic, unit-tested parser
(`harvest/abstract_hr.py`: handles `[HR]`/`(HR)` labels, `95% CI`/`[95% CI, x to y]` forms, skips
prognostic/per-unit covariate HRs, and flags multi-HR abstracts). Restricting to **high-confidence**
extractions (a CI present, ≤2 HR candidates — multi-arm/multi-endpoint abstracts are too ambiguous to
auto-match the right pairwise HR) gives **16** independent comparisons:

| comparison | result |
|---|---|
| reconstructed vs **published** HR, median fold | **1.165** (≈ the 1.149 vs registry) |
| reconstructed HR inside the **published** 95% CI | **13/16** |
| **registry vs published** HR agree (both held-out sources) | **13/16** |

Two things matter here. First, the reconstruction agrees with the *independent* published HR about as
well as with the registry HR — it is not merely echoing the registry. Second, the two held-out "truths"
themselves agree only **13/16 (81%)**: registry and published HRs diverge on ~1 in 5 trials (different
analysis populations, ITT vs per-protocol, data cuts, or endpoint). So a residual ~1.15× fold is partly
irreducible *"which HR did you mean"* noise, not reconstruction error — a ceiling on what *any* method can
score against registry truth. Numbers from `realipd/pubmed_validation.json`; `python
harvest/pubmed_validation.py` to reproduce. Honest limit: the abstract HR is the first non-covariate
HR-with-CI; the high-confidence filter and per-row provenance keep this a triangulation, not a gold
standard (e.g. the one CheckMate-067-style 3-arm trial where the reconstruction picks the wrong arm pair
shows up correctly as an out-of-CI miss).

### Scaled to the whole reconstructed cohort (registry-independent)

The validation-grade set above needs both a curve *and* a registry HR. But the published HR is held-out
truth even when **AACT posts no HR at all** — so we extend the check to **every** reconstructed 2-arm
trial (`node validate/cohort_recon_export.js` → `python harvest/cohort_pubmed.py`): 250 reconstructed →
164 with a PMID → 155 abstracts → **20 high-confidence published HRs**, of which **12 have no registry
HR** (pure, registry-independent validations the gallery could not provide):

| full-cohort independent check | result |
|---|---|
| reconstructed vs published HR, median fold | **1.097** |
| reconstructed HR inside the published 95% CI | **17/20 (85%)** |
| …of these with **no** registry HR (registry-independent) | **12/20** |

More data tightened the central estimate (1.097 vs the validation-grade 1.165). The 3 out-of-CI cases are
all *genuine* (verified): two are hard multi-arm reconstructions where the engine picks the wrong arm pair
(`NCT01721772` = CheckMate-067), one is a real registry-vs-published source divergence — none are
extraction artifacts. Numbers from `realipd/cohort_pubmed_validation.json`.

### Does the reconstruction's *uncertainty* hold up on real curves?

The gold-standard check (`goldstandard_uncertainty.js`) showed the multiple-imputation 95% credible
interval covers the **true** HR **28/29 (97%)** — but on clean open IPD. The real-world question: on the
actual coarse, number-at-risk-less registry curves, does the reconstructed interval cover the
**independent published HR**? Running `reconstructEnsemble` (M=200) on each high-confidence trial
(`node validate/cohort_uncertainty_validation.js`):

| | result |
|---|---|
| published HR inside the reconstructed 95% credible interval | **17/20 (85%)** |
| median CI width ratio (hi/lo) | **2.56** (gold standard 2.46 — comparable, not inflated) |

So on real registry input the interval still contains the independent published effect **85%** of the
time, with a width essentially the same as on clean IPD. The 85% (vs the gold standard's 97%) is expected
and honest: it folds in the registry-vs-published effect divergence (~1 in 5 trials) on top of
reconstruction uncertainty — so it is a conservative real-world floor, not a clean coverage probability.
The 3 misses are the same verified-genuine cases (the CheckMate-067 multi-arm reconstruction failure, one
degenerate fit, and one registry-vs-published divergence where the interval correctly covers the registry
HR it was built from). Numbers from `realipd/cohort_uncertainty_validation.json`.

### …and vs the published *median* (the tightest estimand)

The HR is the reconstruction's hardest quantity; the **median** is its tightest (~3% on the open gold
standard). We check that against the same abstracts: extract the **published per-arm medians**
(`harvest/abstract_median.py` — handles two-arm "X vs Y months" forms, Lancet middle-dot decimals,
weeks/years, skips "improved by X months" differences) for the **endpoint matching the reconstructed
curve** (the curve's OS/PFS family from `validation_hr_backfill.json` — an OS abstract-median is never
scored against a PFS curve, exactly as for the HR), reconstruct the per-arm medians, and compare by
sorted magnitude. Across **5** endpoint-matched trials (**10** arm-medians): median fold **1.071**,
**7/10 within 20%**, with the OS and clean-PFS trials at ~3% (e.g. `NCT00861614` published 10.0/11.2 →
reconstructed 10.1/11.2). This is an estimand independent of *both* the registry and the HR, confirming
the median claim holds on real registry curves. Endpoint-matching was decisive: before it, an
OS-median-vs-PFS-curve mix inflated the fold to 1.21 (one trial 4× off); matching collapsed that to
1.071. Numbers from `realipd/pubmed_median_validation.json`. Honest limit: small n (abstracts that state
a clean two-arm same-endpoint median pair), and coarsely-sampled PFS curves run looser (~20–35% on two
trials) — shown per-row, not hidden.

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
