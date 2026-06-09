# Registry-Native Pseudo-IPD Reconstruction — Methods

*A consolidated, paper-style summary of the method set and its empirical validation. Companion to
`README.md` (usage) and `VALIDATION.md` (full result tables).*

## 1. Problem and scope

Reconstruct individual-patient time-to-event data (pseudo-IPD) for survival meta-analysis **using
only ClinicalTrials.gov / AACT data** — no published-figure digitization, fully offline. The central
constraint defines the contribution: every existing tool (IPDfromKM, KM-GPT, RESOLVE-IPD,
Guyot-on-image) digitizes a *published curve image*, incurring pixel error. AACT exposes *structured*
summary data instead — so where it posts Kaplan–Meier estimates, reconstruction uses **exact
registry anchors with zero digitization error and full provenance** to the NCT record. Output is
always **pseudo-IPD**, never true IPD, and is refused where the registry is too sparse.

## 2. Data source — what AACT actually contains (2026-06-01 snapshot)

Measured over 76,067 trials with posted results:

| | count |
|---|---|
| Tier A — survival curve at ≥3 parseable timepoints | 288 (broad regex: 595 incl. "probability of event"/"proportion") |
| Tier B — median + hazard ratio | 3,283 |
| Tier C — sparser | 72,496 |
| trials with any HR/Cox analysis | 5,567 |
| trials with a reported median | 18,761 |
| **structured number-at-risk rows (entire database)** | **0** |

Two facts shaped every design choice: (i) **AACT has no structured number-at-risk**, so classical
Guyot reconstruction cannot run as-published and Tier A must not require it; (ii) survival is often
posted as **cumulative incidence** ("probability of progression/event"), not survival.

## 3. Harvesting (AACT → trial JSON)

`outcome_measurements` (KM-estimate timepoints), `outcome_analyses` (HR/CI/method), `outcomes`
(endpoint/units), `result_groups`+`outcome_counts` (arms and per-arm N), `milestones`+
`drop_withdrawals` (participant flow). Key robustness measures, each forced by a real failure:

- **Group filtering by `outcome_id`** — `result_groups` repeats arms per outcome (147 rows ⇒ a
  spurious "147 arms" until filtered).
- **Data-driven survival/incidence orientation** — a curve starting near 0 and rising is cumulative
  incidence ⇒ `S = 1 − value` (robust form of the sign-flip lesson; verified on NCT00725985).
- **N per arm from `outcome_counts`**, not measurements.
- **N-matched flow↔outcome mapping** — map a participant-flow group to an outcome arm by its
  milestone `STARTED` count == the arm's N (robust to AACT's title/order inconsistencies, which made
  title/suffix mapping produce absurd HRs).

## 4. Tiered reconstruction

Classified per trial by richness (`tier = min over arms`):

- **Tier A** (KM curve ≥3 timepoints + N): reconstruct from the structured anchors. Two methods are
  run and the one with the lower **1-Wasserstein distance to the registry anchors** is selected:
  - **Guyot inverse-KM** (Guyot 2012): iterative event/censor recovery, constant censoring within
    number-at-risk blocks; here driven by N + total events when NAR is absent.
  - **Censoring-informed anchor-exact** (RESOLVE-IPD CEN-KM style, 2025): deaths taken directly from
    the curve with at-risk held constant within each interval; censoring only at boundaries ⇒ the
    reconstructed KM passes through the registry anchors ≈ exactly. Correct when censoring is
    administrative (which Guyot's constant-censoring assumption mishandles).
  Population is conserved by a forward capacity walk; reconciliation to `total_events` distributes
  the correction **proportional to the curve's death profile** (never piling events at t=0).
- **Tier B** (median + HR + N + events): parametric (exponential) per arm, event count matched
  exactly via an order-statistic cutoff, arms coupled by a proportional-hazards anchor; uncertainty
  by seeded bootstrap envelope (seed = hash(nct_id)).
- **Tier C**: fail closed — reconstruction refused (any IPD here would be fabrication).

## 5. Refinements (validated, opt-in)

- **Censoring-informed event counts** (`add_event_counts.py`): `drop_withdrawals` records *why*
  participants left; for a time-to-event endpoint the event-type reasons (progression/death/relapse
  for PFS; death for OS) give the per-arm event count, the rest is censoring. With N-matched group
  mapping this corrects the over-counting that otherwise attenuates the HR.
- **HR-calibration** (`{calibrateHR:true}`): 1-D bisection on the experimental arm's censoring level
  so the reconstructed Cox HR reproduces the registry HR, **preserving the anchors**. *Imposes* the
  published effect — the correct object for downstream IPD-MA, not a recovery claim.
- **Royston–Parmar flexible parametric** (`{smooth:'rp', extrapolateTo}`): OLS fit of `log(−log S)`
  on `log t` with a restricted cubic spline (no MLE — exact anchors available). Used for
  **extrapolation beyond observed follow-up** (HTA mean-survival), *not* within-data (the step curve
  is already near-exact there).

## 6. Self-audit (Bronze / Silver / Gold)

Reconstructed pseudo-IPD is checked back against the registry: total-event match (C1), anchor
survival fidelity (C2), median within 5% (C3), reconstructed-vs-registry HR via ridge-penalized Cox (C4),
monotonicity (C5, hard), number-at-risk consistency (C6, when present), population conservation
(C7, hard), follow-up sanity (C8), **HR direction integrity (C9, hard — the HR is never inverted,
only arm labels resolved)**. Hard-check failure or Tier C ⇒ badge `none` and export blocked.

## 7. Validation

**HR vs registry HR (held-out ground truth, paired on 30 two-arm trials; within-CI role-oriented,
no inverse clause; Wilson 95% CIs; ~18 direction-determinable):**

| | curve-only | N-matched censoring-informed |
|---|---|---|
| median HR fold-error | 1.12 (~12%) | 1.10 |
| within registry CI (oriented) | 83% (15/18, CI 61–94) | 94% (17/18, CI 74–99) |
| direction correct | 83% (15/18) | 89% (16/18) |

Paired McNemar: informed gains 2 within-CI, loses 0 — **directionally favorable, not significant at
n≈18.** Per-event-count strata omitted (informed shifts stratum membership → non-paired).

**Round-trip consistency (NOT external accuracy):** RMST recon-vs-anchor 0.19% median (n=1256);
median recon-vs-curve-0.5-crossing 0% (n=605). These are near-tautological — anchor-exact is built to
pass through the anchors — so they confirm the expand step preserves the curve, nothing more. The one
*external* check (recon median vs separately-reported registry median) is 48%, attributed (supported,
not proven) to endpoint mismatch. RADIANT-4 hand-verified illustration: recon median 12.0/4.0 mo vs
published PFS 11.0/3.9 mo; HR curve-only 0.68 → censoring-informed **0.47** (event counts set
manually from `drop_withdrawals`) vs registry **0.48**.

## 8. Honest conclusions

1. **RMST and median** are recovered with **near-zero round-trip error** — but that is internal
   self-consistency (the reconstruction is built from the curve), **not validated external accuracy**.
   A clean same-endpoint external RMST/median check is the key missing validation. They remain the
   theoretically better-behaved estimands (curve-derived; preferred under non-proportional hazards),
   but "primary source" status is **not yet earned by evidence**.
2. **HR** is the hard estimand (needs event timing/censoring AACT doesn't fully report). Curve-only
   recovers it to ~12% (median) and inside the registry CI ~83% of determinable trials; N-matched
   censoring-informed nudges that to ~94%, but at n≈18 the two are **statistically indistinguishable**.
   Useful as a triangulation input; HR-calibration imposes the HR exactly when consistency with the
   published effect is required (imposition, not recovery).
3. The recurring bottleneck was **AACT's group/endpoint mapping, not the statistics.** Every apparent
   method failure (147 arms, HR=20 explosions, 48% "median error") traced to a mapping/endpoint
   mismatch; fixing the mapping (N-matching) is what unlocked the gains.
4. **Coverage is the real limit**: registry-native reconstruction applies to the few hundred trials
   that post a structured curve, with the censoring-informed refinement gated on a clean flow↔outcome
   N-match. We report this explicitly rather than overclaiming scale.

## 9. Reproducibility

```
# coverage + cohort (needs a local AACT snapshot; set AACT_ZIP)
python harvest/coverage_scan.py            # Tier A/B/C universe
python harvest/harvest_cohort.py           # batch-harvest Tier-A trials -> cohort/
python harvest/add_event_counts.py cohort  # N-matched censoring-informed event counts
# validation
node   validate/validate_hr.js cohort      # HR vs registry HR
node   validate/validate_rmst.js cohort    # RMST + median fidelity
# tests / tool
node --test test/engine.spec.js            # 17 engine + method tests
python -m pytest harvest/ -q               # 16 harvester tests
node build.js                              # -> dist/registry-ipd.html (offline single file)
```

*References: Guyot P. et al. BMC Med Res Methodol 2012;12:9. Royston P, Parmar MKB. Stat Med
2002;21:2175. RESOLVE-IPD (arXiv:2511.01785, 2025). IPDfromKM (BMC MRM 2021).*
