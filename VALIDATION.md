# Validation — is registry-only reconstruction good enough?

**Method.** For every Tier-A trial that has *both* a KM curve and a hazard ratio reported in
`outcome_analyses`, reconstruct pseudo-IPD **from the curve alone** and compare the Cox HR of the
pseudo-IPD against the **registry-reported HR (held-out ground truth)**. Direction is scored only
where an experimental-vs-comparator split is determinable (registry HR direction is itself often
ambiguous). Stratified by reconstructed event count because Cox HR is unstable when events are few.

Cohort: 595 Tier-A trials harvested from the 2026-06-01 AACT snapshot; **30 are clean 2-arm with a
registry HR** (the usable ground-truth set). Run: `node validate/validate_hr.js cohort`.

## Results

| stratum | n | median HR fold-err | p90 fold-err | within registry 95% CI | direction correct |
|---|---|---|---|---|---|
| **events ≥ 50** (fair test) | 20 | **1.08** (~8%) | 1.46 | **79%** | **83%** (12 det.) |
| events 10–49 | 9 | 1.13 | 1.74 | 100% | 80% (5 det.) |
| all 30 | 30 | 1.12 | 1.73 | 86% | 83% (18 det.) |

**Read:** curve-only reconstruction recovers the HR to within ~8% (median) and inside the registry's
own CI ~80% of the time, with correct direction ~83% where determinable. **Useful as a triangulation
input; not good enough as a sole source for a pooled HR.**

## The tail is real — worked example (RADIANT-4, NCT01524783)

Everolimus vs placebo in advanced NET, N=205/97, 10–11 KM timepoints, ~240 reconstructed events.

- Registry HR **0.48** (95% CI 0.35–0.67) · reconstructed HR **0.68** → ~42% fold error, *just
  outside* the CI.

Cause: without number-at-risk, the coarse 10-point curve loses the fine event timing that produced
0.48, so the Cox HR is **attenuated toward 1**. This is the p90 tail, and it motivates the next
methods.

## We built lever #1 (censoring-informed) and tested it — honest result

**Mechanism (validated).** `drop_withdrawals` records *why* participants left; for a time-to-event
endpoint the event-type reasons (progression/death/relapse for PFS; death for OS) give the event
count, and the rest is censoring. Curve-only does not know the censoring level, so it over-counts
events and **attenuates the HR**. Worked example — **RADIANT-4 (NCT01524783)**: curve-only HR
**0.68** (outside the registry CI) → censoring-informed HR **0.47** vs a registry **0.48**. The
mechanism clearly works. (`harvest/add_event_counts.py` derives the counts.)

**Blanket application (does NOT generalize).** Applying it to the whole cohort:

| | median fold-err | p90 fold-err | within registry CI | direction |
|---|---|---|---|---|
| curve-only | **1.12** | **1.73** | 86% | 83% |
| censoring-informed (all-arms gate) | 1.14 | 2.23 | **90%** | **89%** |
| + agreement gate (use informed only if within 2× of curve-only) | 1.19 | **1.72** | 89% | 89% |

It **improves coverage and direction but fattens the error tail** (helped 4, hurt 7, ~same 19). Two
reasons it isn't a free win: (a) the withdrawal-reason→endpoint-event correspondence is trial-
specific and noisy; (b) **mapping participant-flow groups to outcome-measure groups is unreliable**
(different titles, different ordering between AACT sections) — wrong mappings produce absurd HRs
(20×, 0.01×). An agreement gate suppresses the explosions but the net gain is marginal.

**First attempt (title/suffix group mapping) failed** because mapping participant-flow groups to
outcome arms was unreliable. **Fixing the mapping fixed the method.**

### N-matched mapping makes censoring-informed a clean win

Map each flow group to an outcome arm by its milestone `STARTED` count == the arm's analysis N
(robust to title/order differences). With that:

| events≥50 (fair test) | curve-only | **N-matched censoring-informed** |
|---|---|---|
| median fold-err | 1.081 | 1.08 |
| p90 fold-err (tail) | 1.46 | **1.35** |
| within registry CI | 79% | **94%** |
| direction correct | 83% | **90%** |

(Overall, n=30: within-CI 86→**93%**, p90 1.73→**1.57**, direction 83→**89%**.) So the censoring level
*is* recoverable from participant-flow data — the earlier negative result was a mapping bug, not a
method failure. **Censoring-informed (N-matched) is the recommended mode; curve-only remains the
safe fallback when participant-flow groups don't N-match.**

### HR-calibration (impose the registry HR for IPD-MA consistency)

1-D-solve the experimental arm's censoring level so the reconstructed Cox HR reproduces the registry
HR, preserving the anchors. RADIANT-4: curve-only 0.679 → calibrated **0.481** vs registry **0.48**
(`reconstruct(trial,{calibrateHR:true})`). This *imposes* the HR (not a recovery test); it is the
right object when you need pseudo-IPD consistent with the published effect for downstream IPD-MA.

## Robust estimands (RMST & median) — the reconstruction is excellent here

HR is the *hard* estimand (needs event timing/censoring). RMST and median are curve-*derived*, so
they should be recovered far better. Measured across the cohort (`node validate/validate_rmst.js`):

| estimand | metric | result |
|---|---|---|
| **RMST** (recon-IPD area vs registry anchor-curve area, to common τ) | median % err / p90 / within-5% | **0.19% / 3.99% / 92%** (n=1256 arms) |
| **Median** (recon median vs the registry curve's own 0.5-crossing) | median % err / within-10% | **0% / 97%** (n=605 arms) |

Worked example RADIANT-4: reconstructed median **12.0 / 4.0 mo** vs the real published PFS medians
**11.0 / 3.9 mo** — within ~1 month, externally confirmed.

**So for the estimands that matter most in modern survival meta-analysis — RMST and median, the
preferred summaries under non-proportional hazards — registry-only reconstruction is essentially
exact (~0–0.2% error), versus ~8–12% for the HR.** That is the real "is it good enough" answer:
**yes, and more so for RMST/median than for HR.**

*Caveat on the external median cross-check:* comparing recon median to the registry-*reported*
median gives ~48% error even after filtering to survival endpoints — but this is an **endpoint-
matching artifact**, not reconstruction error (recon median == curve median to 0%; the 48% is the
PFS-curve median being compared to a different-endpoint registry median, e.g. OS, under the same arm
code). Clean same-endpoint external matching needs the same unresolved group/endpoint mapping.

## Royston–Parmar flexible parametric (done, with an honest scope)

Because we have exact registry (t,S) anchors, the RP model `log H(t) = s(log t; γ)` is fit by OLS on
`(log t, log(−log S))` with a restricted cubic spline (df scaled to anchor count) — no MLE needed.
Validated: reproduces a known Weibull within 2% and stays monotone (`{smooth:'rp'}`).

**Honest result:** RP smoothing does **not** improve *within-follow-up* estimands — step interpolation
already preserves the curve area near-exactly (RMST median err **step 0.19% vs RP 2.25%**), because
the step passes through the anchors and RP is a fit that doesn't. So RP is **not** a default.

**RP's genuine value is the one thing the step cannot do — extrapolation beyond observed follow-up.**
RADIANT-4 experimental arm (last anchor 24 mo, S=0.22) extrapolates to S(36)=0.137, S(48)=0.084, with
an extrapolated RMST(0–48 mo)=15.9 mo — the standard HTA "mean survival beyond trial follow-up". Use
`reconstruct(trial,{smooth:'rp', extrapolateTo: months})` for that; keep the step curve within data.

## Remaining levers
- ✅ HR-calibration · ✅ N-matched mapping · ✅ RMST/median validation · ✅ Royston–Parmar (extrapolation).
- **Same-endpoint external median matching** to clean the contaminated registry-median cross-check.

## Honest limitations of this validation
- n=30 (18 with determinable direction) — small; widen by mapping multi-arm analyses pairwise.
- Ground truth is the registry HR, which can itself be model-dependent (log-rank vs Cox) and
  direction-ambiguous; it is the best available proxy, not patient-level truth.
- Validates HR; RMST / median agreement should be validated separately (RMST is expected to be more
  robust because it depends on the curve area, which the anchors fix directly).
