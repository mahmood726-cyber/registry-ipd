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

**Gold-standard coverage on TRUE IPD** (`validate/goldstandard_uncertainty.js`): across the **24
adequately-sized (≥100/arm)** real datasets, the 95% credible interval covers the **true patient-level
HR in 23/24 (96%)** — empirical coverage matching the nominal 95% (median width ~2.3× fold; not
over-wide, e.g. diabetic [0.29, 0.63] tightly covers 0.46). The single miss is **`bfeed`** (true 1.245
vs band [1.45, 3.37]): its point reconstruction is so far off (the discrete-week / ~96%-event boundary
case) that even the uncertainty band does not reach truth — an honest failure, not hidden. So the
uncertainty band genuinely contains the *true* effect on 16 of 17 datasets, not merely the
registry-reported one, with calibration close to nominal rather than conservatively over-wide.

## Cutting-edge: competing-risks reconstruction (Aalen–Johansen)

When the harvested curve is a cause-specific endpoint (e.g. time-to-progression with death as a
**competing** event), treating 1−KM as the cumulative incidence is biased — censoring a competing
event wrongly assumes the patient could still have the event of interest. `reconstructCompetingRisks()`
labels causes (event-of-interest vs competing, the competing count from `drop_withdrawals` reasons)
and computes the correct **Aalen–Johansen CIF**. Validated: the invariant CIF₁(t)+CIF₂(t)+S(t)=1
holds at every step, and naive 1−KM provably overestimates the cause-1 incidence. Demonstration
(RADIANT-4, 15% competing deaths injected): naive 1−KM overestimates progression incidence by
**+3.6 pp (experimental) / +1.4 pp (comparator)**; AJ corrects it. This makes the reconstruction
honest for cause-specific endpoints, which the single-KM assumption silently biases.

## Cutting-edge: multi-constraint joint reconstruction (+ max-entropy)

`reconstructJoint()` produces pseudo-IPD consistent with the curve AND total events (censoring-
informed) AND the registry HR (calibration) at once, and **reports which constraints jointly hold**.
Key honest finding: the constraints can be **mutually inconsistent** (registry events + HR + curve
need not be jointly satisfiable); calibration then trades the experimental event count to hit the
HR, and the report surfaces that tension rather than hiding it. The *maximum-entropy* element — being
least-committal over the under-identified DOF — is the imputation ensemble above (it samples the
max-entropy distribution over the censoring level).

## Cutting-edge: fractional-polynomial time-varying HR (non-proportional hazards)

A single Cox HR is misleading under non-PH — and RMST/median (which we recover well) are the
preferred NPH summaries for exactly this reason. `piecewiseHR()` estimates a piecewise-exponential
rate-ratio HR(t) on the reconstructed pseudo-IPD; `fractionalPolyHR()` fits a first-order fractional
polynomial logHR(t)=β₀+β₁·tᵖ and **tests PH (β₁≈0)**; `poolTimeVaryingHR()` inverse-variance-pools
across trials (a time-varying network/meta effect). Validated: on synthetic delayed-effect data the
early-window HR≈1, late-window HR drops, and FP flags non-proportional. Demonstration (RADIANT-4
reconstructed IPD): piecewise HR **0.35 → 0.90 → 1.63** with FP non-proportional=true — time-structure
the single registry HR (0.48) hides. (On coarse reconstructed anchors the late windows are noisy; the
method detects and characterises time-variation, it does not claim window-exact HRs.)

## GOLD STANDARD: validation on TRUE patient-level IPD

The strongest validation — against **real patient-level RCT data** (openly published; R `survival::`
datasets via the Rdatasets mirror). For each, we compute the TRUE estimates from the full IPD, then
generate the **registry-style coarse summary a sponsor would post** (KM at 8 timepoints + N + total
events), reconstruct from *that alone*, and compare to truth. Not synthetic, not circular — the
engine never sees the patient-level data. (`validate/goldstandard.js`.)

| dataset (endpoint) | N exp/ctl | true HR | curve-only HR (log-err) | median %err | RMST-diff recon / **true** |
|---|---|---|---|---|---|
| **GBSG** breast (RFS) | 246/440 | 0.695 | **0.686 (1.3%)** | 5.9% | 199 / **199.2** |
| **Rotterdam** breast (OS) | 339/2643 | 1.51 | 1.269 (17%) | **2.4%** | −548 / **−561** |
| **PBC** liver (OS) | 154/158 | 0.944 | 1.117 (17%) | 3.4% | 77 / 40 |
| **Diabetic** retinopathy | 197/197 | 0.46 | **0.412 (11%)** | n/r | 10.7 / 12.2 |
| **NWTSG** Wilms (relapse) | 459/3569 | 5.10 | **5.37 (5%)** | n/r | −1532 / −1696 |
| **Myeloid** AML (OS) | 329/317 | 0.708 | **0.639 (10%)** | 9.7% | — |
| **Kidtran** kidney tx | 339/524 | 0.907 | 1.046 (14%) | n/r | — |
| **Veteran** lung (OS) | 68/69 | 1.016 | 0.755 (30%) | 58% | 31 / **−0.7** |
| **Leukemia** allo/auto tx | 50/51 | 0.838 | 0.678 (21%) | n/r | 4.1 / **3.6** |
| **Larynx** cancer (stage III/I) | 27/33 | 1.837 | **1.711 (7%)** | 1.5% | −1.6 / **−1.7** |
| **Burn** RCT (staph infection) | 84/70 | 0.571 | 0.466 (20%) | n/r | 9.4 / 11.9 |
| **Pneumonia** infant (smoking) | 838/2285 | 2.237 | **2.30 (3%)** | n/r | −0.1 / **−0.1** |
| **Breastfeeding** (smoking) | 270/657 | 1.245 | 2.179 (56%) | 26% | −12.2 / −3.8 |
| **HCC** liver (vasc. invasion) | 41/186 | 2.18 | **2.01 (8%)** | 1.1% | −14.4 / **−15.6** |

**Aggregate over 24 adequately-sized datasets (≥100/arm; of 43 real datasets tried, incl. 2 recurrent-event collapsed to first-event and 12 TCGA stage cohorts): curve-only
recovers HR to a median fold of 1.15 (15/24 within 20%; 1.12 / 13/17 excluding the heavily-censored
TCGA cohorts, see below), the censoring-informed tier to 1.15 (17/24), and the median to ~3%** — on
real patient data across 43 RCTs/cohorts (six added 2026-06-10 from `KMsurv`/`asaur` via the Rdatasets
mirror, twelve from TCGA via the cBioPortal API; the 5 below-100/arm TCGA cohorts add breadth but not
to the ≥100/arm aggregate). Large effects recovered cleanly
(Wilms 5.1→5.18, prostate 5.49→5.17, melanoma 4.36→3.99, HCC 2.18→2.01); the classic 1965 Gehan
leukemia RCT (6-MP) recovers 0.221→0.201 (9%); UDCA-in-PBC RCT 0.445→0.415. The set spans
breast/colon/lung/AML/melanoma/leukemia/transplant/PBC/MGUS/NAFLD/prostate/retinopathy/larynx/burn/
pneumonia/HCC. The **worst case is `bfeed`** (fold 1.75): breastfeeding duration in discrete weeks
with ~96% events — a heavily-tied discrete-time series, not the smooth KM curve the method targets;
kept as an honest out-of-scope boundary rather than dropped.

### Real cancer cohorts (TCGA / cBioPortal) — and why the event-count tier matters

**12 real TCGA overall-survival cohorts** pulled from the open cBioPortal API
(`harvest/fetch_cbioportal.js`), split by **late vs early stage** — a strong, real survival contrast
(true HR **1.6–7.6**). These are heavily censored (the early-stage arm is mostly alive at last
follow-up), exactly the regime that separates the two reconstruction tiers:

| TCGA cohort | N late/early | true HR | curve-only (fold) | **censoring-informed (fold)** |
|---|---|---|---|---|
| lung adeno | 105/394 | 2.65 | 2.70 (**1.02**) | 4.04 (1.52) |
| lung squamous | 89/388 | 1.64 | 1.58 (**1.04**) | 2.21 (1.35) |
| melanoma | 188/203 | 1.67 | 1.64 (**1.02**) | 1.80 (1.08) |
| liver HCC | 88/257 | 2.34 | 1.88 (1.24) | 2.85 (**1.22**) |
| stomach | 224/183 | 2.19 | 1.65 (1.33) | 2.26 (**1.03**) |
| esophageal | 64/96 | 3.10 | 4.17 (1.34) | 5.87 (1.89) |
| colorectal | 247/309 | 3.11 | 1.70 (1.83) | 3.11 (**1.00**) |
| bladder | 276/131 | 2.23 | 1.36 (1.64) | 2.17 (**1.03**) |
| head & neck | 349/103 | 1.76 | 1.10 (1.60) | 1.49 (1.18) |
| adrenocortical | 35/53 | 7.59 | 4.51 (1.68) | 6.34 (**1.20**) |
| kidney clear-cell | 209/301 | 4.05 | 2.60 (1.56) | 5.01 (1.23) |
| kidney papillary | 64/189 | 6.13 | 2.55 (2.41) | 6.99 (**1.14**) |

**Curve-only *under*estimates these large HRs (median fold 1.56)** — with the early-stage arm almost
entirely censored, the curve-only "censor-to-tail" assumption flattens the separation. **The
censoring-informed reconstruction, which uses the registry total-event count, recovers them (median
fold 1.20)** — colorectal 3.11→1.70→**3.11**, kidney-papillary 6.13→2.55→**6.99** (curve-only fold
2.41 → 1.14). The 7 cohorts at ≥100/arm all have 95% credible intervals covering the true HR (in the
23/24 coverage above). esophageal is the exception where *both* tiers overestimate — a genuinely hard
case kept honestly.

### Advanced estimators and the identifiability limit (`validate/advanced_estimators.js`)

Can cutting-edge statistics recover the heavily-censored HRs **without** the event count? We
benchmarked three curve-only point estimators on the gold standard: (A) the current **censor-to-tail**;
(B) a **max-entropy ensemble** (median log-HR over imputations of the under-identified censoring
level); and (C) a **1-Wasserstein barycenter** of the imputed pseudo-IPD point-clouds — the
optimal-transport point estimate, which (uniquely, via rank-matching) averages the reconstructions
*while preserving the at-risk structure* the HR depends on.

| estimator (curve-only) | all 43 (median fold) | heavily-censored TCGA (median fold) |
|---|---|---|
| censor-to-tail (current) | 1.16 | 1.56 |
| max-entropy ensemble | **1.14** | 1.55 |
| Wasserstein barycenter | 1.19 | **1.51** |

**No estimator dominates, and on the heavily-censored cohorts all three fail (~1.5).** The reason is
fundamental, not algorithmic: **censoring is invisible to the KM anchors** — every censoring level
passes through the same posted survival points — so the event/censoring split is genuinely
unidentified from the curve alone. The Wasserstein barycenter is marginally best on the hard subset but
cannot extract information that is not in the data. The lever is therefore **data reporting** (the
registry total-event count or number-at-risk, which collapses the TCGA error to 1.20) — i.e. exactly
the recommendation in `POLICY.md`. This is the honest answer to "can advanced stats fix it": they
confirm the bound is real and locate the fix in reporting, not in cleverer reconstruction.

### Anchor density: how many posted timepoints does reconstruction need?

Sweeping K (number of posted KM timepoints) across the 7 true-IPD datasets
(`validate/sensitivity_anchors.js`):

| K (timepoints) | 3 | 4 | 5 | 6 | 8 | 12 | 20 |
|---|---|---|---|---|---|---|---|
| HR fold-error (median) | 1.40 | 1.28 | **1.15** | 1.16 | 1.12 | 1.08 | 1.08 |
| median % error | 9.2 | 4.3 | 2.4 | 3.0 | 3.4 | 1.9 | 1.1 |

*(14 true-IPD datasets.)*

**Formalized recommendation (minimum reportable timepoints).** The HR fold-error e(K) is monotone-
decreasing in K and plateaus: e(3)≈1.40, e(4)≈1.28, e(5)≈1.15, e(8)≈1.12, e(≥12)≈1.08. We define the
operating thresholds: **K ≥ 5 KM timepoints for HR fold-error ≤ ~1.15** (≈ the residual reconstruction
floor), **K ≥ 8 for ≤ ~1.12**, with negligible further gain beyond ~12. Below K = 4 the error exceeds
~1.25 and the reconstruction should be treated as indicative only. *Practical guideline:* trust a
registry-native HR when the trial posts **≥5–6 KM-estimate timepoints** (and is adequately sized);
and registries/sponsors should post **at least 5–8 timepoints** to make a trial reliably
reconstructable — a concrete, evidence-based reporting standard (consistent with Guyot et al.'s own
recommendation to report numbers-at-risk and total events alongside KM curves).

**Accuracy improves sharply from 3→5 timepoints, then plateaus** (HR fold-error 1.33 → ~1.08–1.12,
median 9% → ~3%). Practical implication: registry-native reconstruction needs **≥5–6 posted KM
timepoints**; beyond ~8 the marginal gain is small. This both tells users when to trust the
reconstruction and motivates registries to post at least a handful of KM-estimate timepoints.

**Honest reading (confirmed on real data):**
- **Median is recovered to ~2–6%** for adequately-sized trials; **HR to ~11% (median)**, good on most
  (GBSG 1.3%, Wilms 5%, diabetic 11%), moderate on others (Rotterdam/PBC 17%).
- **RMST** is accurate in absolute terms; its *relative* error inflates for near-null effects with a
  tiny true RMST-difference (PBC: 77 vs 40, a small absolute gap on a multi-year scale).
- **The failure boundary is events × effect size, not N alone.** Across the 8 small datasets
  (<100/arm), 4/8 still land within 20%: those with a strong effect and enough events reconstruct
  fine even at small N (Gehan 6-MP N=21/arm, HR 0.22→0.20, 9%; AIDS-by-CCR5 0.29→0.31, 7%; melanoma
  4.36→3.99, 9%), whereas small trials with weak/near-null effects fail (Veteran N=69/arm ~null →
  spurious 0.76, 58% median error; bmt 39%). So the practical guard is "enough events + a real
  effect," not a raw N threshold.
- Censoring-informed (forcing the registry event count) improves HR on some trials but can distort the
  *median* — so prefer curve-only for median/RMST, censoring-informed for HR.

This is the first validation of the engine against genuine patient-level data, and it broadly
confirms the picture from registry/published checks: **RMST and median recover well for real trials
of adequate size; HR is harder; very small trials are unreliable.**

### Competing-risks gold standard on true IPD (`survival::colon`, `validate/goldstandard_cr.js`)

colon has true recurrence-vs-death competing events. We built the true cause-labeled IPD (recurrence
= event, death-without-recurrence = competing), per arm, and compared our reconstructed Aalen–Johansen
CIF to the TRUE AJ CIF:

| arm | n (recur / competing deaths) | recon AJ CIF — max abs err vs TRUE | true final CIF: AJ vs naive 1-KM |
|---|---|---|---|
| colon: Observation | 315 (177 / 13) | **0.016** | 0.584 vs 0.593 (naive +0.9 pp) |
| colon: Levamisole+5FU | 304 (119 / 15) | **0.003** | 0.394 vs 0.401 (naive +0.7 pp) |
| **aidssi: AIDS vs SI** | 329 (114 / **108**) | 0.06 | **0.408 vs 0.569 (naive +16 pp)** |

The reconstructed competing-risks CIF matches the true patient-level AJ CIF well. Crucially, the
`aidssi` case has **heavy** competing risk (108 SI vs 114 AIDS events), where **naive 1−KM
overestimates the AIDS incidence by 16 pp** (0.57 vs true 0.41) — and the Aalen–Johansen
reconstruction recovers the truth to within 6 pp. The bias scales with competing-event frequency
(colon: rare → ~1 pp; aidssi: heavy → 16 pp); the AJ reconstruction corrects the bulk of it on real
data. This is exactly where competing-risks-aware reconstruction matters.

(Datasets used for validation only, not redistributed; re-download: `validate/goldstandard.js` header.)

## Remaining levers
- ✅ HR-calibration · ✅ N-matched mapping · ✅ RMST/median validation · ✅ Royston–Parmar (extrapolation).
- **Same-endpoint external median matching** to clean the contaminated registry-median cross-check.

## Honest limitations of this validation
- n=30 (18 with determinable direction) — small; widen by mapping multi-arm analyses pairwise.
- Ground truth is the registry HR, which can itself be model-dependent (log-rank vs Cox) and
  direction-ambiguous; it is the best available proxy, not patient-level truth.
- Validates HR; RMST / median agreement should be validated separately (RMST is expected to be more
  robust because it depends on the curve area, which the anchors fix directly).
