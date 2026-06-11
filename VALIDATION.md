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

**On the engine default (re-run 2026-06-10).** The censoring-informed row is now the **Titman QP** (the
engine default when events are posted). On this registry-HR set it gives median fold **1.10 / within-CI
94% / 17 of 18** — essentially unchanged from before, and the McNemar is +3/−1. *This is expected, not
a let-down:* the QP's large gold-standard gains (1.15→1.04, §"GOLD STANDARD") are concentrated in
**heavily-censored, large-effect** trials, and the 30 registry-HR trials are mostly large,
lightly-censored oncology RCTs measured against a *coarse* registry-reported HR. A stratified analysis
(`validate/censoring_stratified.js`) confirms the nuance: the QP's advantage over curve-only is
**consistently positive (mean ~0.19 fold) but not predictable from summary features** — pooled censoring
barely correlates with curve-only error (Spearman −0.02); the weak drivers are censoring **asymmetry
between arms** (0.24) and **effect size** (0.21). So there is **no safe "curve-only is fine here" rule**;
the operational conclusion is to *always* prefer the event count (QP) where the registry posts it — and
a figure number-at-risk table where it does not (§"NAR fusion") — reserving bare curve-only for
triangulation with its credible interval.

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

**Gold-standard coverage on TRUE IPD** (`validate/goldstandard_uncertainty.js`): across the **29
adequately-sized (≥100/arm)** real datasets, the 95% credible interval covers the **true patient-level
HR in 28/29 (97%)** — empirical coverage matching the nominal 95% (median width ~2.3× fold; not
over-wide, e.g. diabetic [0.29, 0.63] tightly covers 0.46). The single miss is **`bfeed`** (true 1.245
vs band [1.45, 3.37]): its point reconstruction is so far off (the discrete-week / ~96%-event boundary
case) that even the uncertainty band does not reach truth — an honest failure, not hidden. So the
uncertainty band genuinely contains the *true* effect on 28 of 29 datasets, not merely the
registry-reported one.

**Multi-level calibration — the intervals are conservative, not over-confident**
(`validate/uncertainty_calibration.js`). Marginal 95% coverage is necessary but weak; an interval can
hit 95% while being miscalibrated elsewhere. Checking coverage at every nominal level (29 datasets,
M=400 imputations):

| nominal | 50% | 80% | 90% | 95% |
|---|---|---|---|---|
| empirical coverage | 83% | 90% | 97% | 97% |

Empirical coverage is **above** nominal at every level. The honest reading: the **95% interval is
well-calibrated** (97% ≈ 95%), but **narrower intervals over-cover** — the band is sized for the tails,
so it is *conservative* (wider than strictly needed) at the centre. This is the **safe** direction of
miscalibration: the intervals will not mislead by being too narrow, but a 50% interval should not be
read as a 50% bet (it covers ~83%). The likely cause is the deliberately wide censoring-level sampling
band ([0.55·E₀, E₀]); tightening it would sharpen the low-level intervals at some risk to tail coverage.
We keep the conservative band and report this honestly rather than tune to a calibration target.

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

**Aggregate over 29 adequately-sized datasets (≥100/arm; of 51 real datasets tried, incl. 2 recurrent-event collapsed to first-event and 12 TCGA stage cohorts): curve-only
recovers HR to a median fold of 1.15 (15/24 within 20%; 1.12 / 13/17 excluding the heavily-censored
TCGA cohorts, see below), and the censoring-informed Titman-QP tier (engine default when an event
count is posted) to 1.05 (23/24 within 20%), with the median to ~3%** — on
real patient data across 43 RCTs/cohorts (six added 2026-06-10 from `KMsurv`/`asaur` via the Rdatasets
mirror, fourteen from TCGA and six from non-TCGA cohorts (METABRIC/SU2C/MSK) via the cBioPortal API; the 5 below-100/arm TCGA cohorts add breadth but not
to the ≥100/arm aggregate). Large effects recovered cleanly
(Wilms 5.1→5.18, prostate 5.49→5.17, melanoma 4.36→3.99, HCC 2.18→2.01); the classic 1965 Gehan
leukemia RCT (6-MP) recovers 0.221→0.201 (9%); UDCA-in-PBC RCT 0.445→0.415. The set spans
breast/colon/lung/AML/melanoma/leukemia/transplant/PBC/MGUS/NAFLD/prostate/retinopathy/larynx/burn/
pneumonia/HCC. The **worst case is `bfeed`** (fold 1.75): breastfeeding duration in discrete weeks
with ~96% events — a heavily-tied discrete-time series, not the smooth KM curve the method targets;
kept as an honest out-of-scope boundary rather than dropped.

### Real cancer cohorts (TCGA / cBioPortal) — and why the event-count tier matters

**14 real TCGA overall-survival cohorts** pulled from the open cBioPortal API
(`harvest/fetch_cbioportal.js`), split by **late vs early stage** — a strong, real survival contrast
(true HR **1.6–7.6**). These are heavily censored (the early-stage arm is mostly alive at last
follow-up), exactly the regime that separates the two reconstruction tiers:

The **censoring-informed** column is the **Titman-2026 quadratic program** (the engine default when a
total-event count is posted; see the next subsection):

| TCGA cohort | N late/early | true HR | curve-only (fold) | **Titman-QP censoring-informed (fold)** |
|---|---|---|---|---|
| lung adeno | 105/394 | 2.65 | 2.70 (**1.02**) | 2.67 (**1.01**) |
| lung squamous | 89/388 | 1.64 | 1.58 (**1.04**) | 1.66 (**1.01**) |
| melanoma | 188/203 | 1.67 | 1.64 (**1.02**) | 1.82 (**1.09**) |
| liver HCC | 88/257 | 2.34 | 1.88 (1.24) | 2.45 (**1.05**) |
| stomach | 224/183 | 2.19 | 1.65 (1.33) | 2.31 (**1.05**) |
| esophageal | 64/96 | 3.10 | 4.17 (1.34) | 3.09 (**1.00**) |
| colorectal | 247/309 | 3.11 | 1.70 (1.83) | 3.28 (**1.05**) |
| bladder | 276/131 | 2.23 | 1.36 (1.64) | 2.27 (**1.02**) |
| head & neck | 349/103 | 1.76 | 1.10 (1.60) | 1.67 (**1.06**) |
| breast | 264/788 | 2.70 | 1.66 (1.63) | 2.74 (**1.02**) |
| kidney clear-cell | 209/301 | 4.05 | 2.60 (1.56) | 4.28 (**1.06**) |
| mesothelioma* | 60/25 | 0.99 | 1.20 (1.21) | 1.12 (**1.13**) |
| adrenocortical | 35/53 | 7.59 | 4.51 (1.68) | 5.66 (1.34) |
| kidney papillary | 64/189 | 6.13 | 2.55 (2.41) | 15.60 (2.54) |

**Curve-only *under*estimates these large HRs (median fold 1.56)** — with the early-stage arm almost
entirely censored, the curve-only "censor-to-tail" assumption flattens the separation. **The Titman-QP
censoring-informed reconstruction, using the registry total-event count, recovers them (median fold
1.05)** — colorectal 3.11→1.70→**3.28**, kidney-clear-cell 4.05→2.60→**4.28**, esophageal 3.10→4.17→
**3.09**. The **7 cohorts at ≥100/arm all land within 20%** (and their 95% credible intervals cover the
true HR, in the 23/24 above). The two honest outliers are both <100/arm and at extreme HR: adrenocortical
(QP 5.66 vs true 7.59) undershoots, and **kidney-papillary overshoots (15.6 vs 6.13)** — the QP can
overcorrect when the true effect is very large and the late-stage arm is tiny and heavily censored.

### Beyond TCGA — non-TCGA published cohorts and credentialed repositories

To diversify beyond TCGA, we added three contrasts from **METABRIC** (the 1,980-patient published breast
cohort, via the same open cBioPortal API; `harvest/fetch_cbio_cohort.js`): grade 3 vs 1 (true HR 1.63,
curve-only 1.51 → **QP 1.03**), tumour-stage 2 vs 1 (1.81, 1.38 → **1.02**), and ER-negative vs positive
(1.16, recovered to **1.00**). These are real patient-level survival from a *different source and cancer
series* than TCGA, and the QP recovers all three to ≤3%. (They are three clinical axes of one cohort, so
not three independent trials — noted for honesty.) cBioPortal exposes **412 non-TCGA studies**; this is a
sampler, and the fetcher generalises to any of them.

The genuinely *credentialed* repositories — **Vivli, YODA, ClinicalStudyDataRequest, Project Data
Sphere** — are gated behind committee-reviewed data-use agreements (and PDS has no open API), so they
cannot be pulled autonomously; that is by design for patient privacy. The ingestion path
(`validate/ingest_ipd.js`, `CREDENTIALED.md`) makes any such export a one-command drop-in (CDISC ADTTE
`CNSR` handled), so when a credentialed export is obtained it extends this gold standard immediately.

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
registry total-event count or number-at-risk, which collapses the TCGA error to **1.05** (24/25 uncertainty coverage) via the QP
below) — i.e. exactly the recommendation in `POLICY.md`. This is the honest answer to "can advanced
stats fix it": they confirm the bound is real and locate the fix in reporting, not in cleverer
reconstruction.

### Titman-2026 quadratic program — the censoring-informed default (`validate/titman_qp.js`)

When the registry *does* post a total-event count, advanced statistics pays off directly. We implement
a **Titman-2026-style quadratic program**: on the cumulative-hazard scale the posted curve fixes the
per-interval hazards `h_k`, so events are `d_k = h_k n_k` and the at-risk recursion
`n_{k+1} = n_k(1−h_k) − c_k` is **linear** in the unknown censoring counts; the total-event count is a
**linear** constraint `E = Σ h_k n_k`; and the leftover censoring degree-of-freedom is resolved by the
convex QP `min ½‖c‖²` s.t. `E(c)=E, c≥0`, with closed-form minimum-norm solution `c_k = max(0, λ A_k)`.
Events are then *spread* within each interval (not piled at the anchor), which is what makes the
at-risk sets — and hence the Cox HR — correct.

| reconstruction (with posted event count) | all 43 | ≥100/arm (24) | TCGA (12) |
|---|---|---|---|
| anchor-exact (previous default) | 1.15 (32/51) | 1.15 (19/29) | 1.20 (8/14) |
| **Titman QP (new default)** | **1.06 (46/51)** | **1.04 (27/29)** | **1.11 (12/14)** |

The QP is now the engine default whenever `total_events` is posted (it cannot be chosen by the
anchor-Wasserstein best-of, because censoring is invisible to the anchors, so it is selected by
data-availability instead). It lifts within-20% from 32→46 of 51 datasets and is unit-tested
(`test/engine.spec.js`). This is the constructive half of the identifiability story: cleverer
reconstruction cannot manufacture the missing event count, but *given* it, the QP extracts the HR
near-exactly.

### Is the QP near-optimal? A 12-method benchmark (`validate/method_zoo.js`)

To check whether the QP is leaving accuracy on the table, we benchmarked **12 advanced
reconstruction/HR methods** on the full gold standard (HR fold-error vs true; arms subsampled to ≤400
for tractability, so true HRs shift slightly from the headline):

| # | method | all-median fold | within-20% | worst |
|---|---|---|---|---|
| 1 | QP, roughness-penalised censoring | **1.034** | 42/45 | 2.38 |
| 2 | QP, max-entropy censoring | 1.037 | **43/45** | 2.69 |
| 3 | QP + ridge-Cox | 1.037 | 42/45 | **2.03** |
| 4 | **QP, L2 min-norm (current default)** | 1.043 | 41/45 | 2.54 |
| 6 | QP + Firth-penalised Cox | 1.045 | 41/45 | 2.42 |
| 7 | max-entropy imputation ensemble | 1.064 | 38/45 | 1.76 |
| 8 | Rubin-pooled log-HR | 1.069 | 38/45 | 1.78 |
| 9 | QP + cumulative-hazard-ratio HR | 1.083 | 31/45 | 2.10 |
| 10 | anchor-exact (pre-QP default) | 1.137 | 31/45 | 1.90 |
| 11 | Guyot | 1.141 | 30/45 | 2.05 |
| 12 | 1-Wasserstein barycenter | 1.194 | 24/45 | 17.97 |

**Findings.** (i) The **QP family dominates** — every QP variant (~1.03–1.04) is far ahead of the
classical Guyot / anchor-exact (~1.14). (ii) The censoring **regulariser is second-order**: roughness
(1.034) and max-entropy (1.037) edge the L2 min-norm (1.043) on the subsampled set, but on the *full*
gold standard the gap closes to within Monte-Carlo noise (L2 1.050 vs max-entropy 1.045) and L2 has the
better worst-case — so we keep the closed-form L2 as the default. (iii) **Ridge-Cox gives the best
worst-case** (2.03) by shrinking extreme overshoots, at the cost of a small global bias. (iv) **Nothing
fixes the two extreme small-arm outliers** (kidney-papillary true 6.13, adrenocortical 7.59): ridge
only pulls kirp 15.6→12.5, and the Wasserstein barycenter *blows up* on it (17.97). Those are a genuine
small-sample limit, not a method deficiency. The honest conclusion: **the current QP is near-optimal —
the remaining error is irreducible from the posted summary, not a better-algorithm problem.**

### NAR fusion — the figure's risk table substitutes for the missing event count (`validate/nar_fusion.js`)

The identifiability limit above is "irreducible *from the posted summary*". But a *figure* of the same
trial prints a **numbers-at-risk table** the sibling `kmcurve` project OCRs — exactly the field AACT
lacks. Fusing the registry-exact curve with that NAR (sparse ~4-column table, ±3% OCR noise, **and no
registry event count**) via a NAR-aware QP reconstruction recovers the HR to QP level:

| set | curve-only (no NAR/events) | NAR via anchor-exact | **fusion: curve + figure NAR (NAR-aware QP)** | QP: curve + registry events |
|---|---|---|---|---|
| all 51 | 1.18 (26/51) | 1.14 (35/51) | **1.05 (46/51)** | 1.05 (47/51) |
| ≥100/arm | 1.15 | 1.13 | **1.03 (28/29)** | 1.04 (28/29) |
| heavily-censored TCGA | 1.38 | 1.19 | **1.06 (20/20)** | 1.05 (18/20) |

**The figure's at-risk table fully substitutes for the missing registry event count** — the fusion
matches the event-count QP and is *20/20 within 20%* on the heavily-censored cohorts. Two notes: (a) the
gain needs the **QP's event-spreading together with the NAR at-risk path** — feeding NAR to the
anchor-exact method alone only reaches 1.14 (it piles events at the anchors); (b) denser/cleaner NAR did
not beat the sparse OCR'd table, so a real 4-column risk table suffices. This is the concrete
**best-of-both registry+figure reconstruction**: where a trial is in both AACT and a publication, the
exact registry curve + the figure's NAR dissolves the identifiability limit *without* the registry
changing its reporting. (`KMCURVE-SYNERGY.md` idea 2, now validated.)

### Abstract event-count lever — the in-scope version of the NAR fusion (`harvest/abstract_events.py`)

The NAR fusion above imports a *figure*, which is **out of this project's production data scope** (AACT +
PubMed abstracts only — a figure is validation-only). The censoring lever the QP actually needs is just a
per-arm total-event count, and the **PubMed abstract** — which IS in scope — routinely prints it. So the
production-legal analogue of the figure's at-risk table is a deterministic abstract event-count
extractor. The full chain is grounded on **real** data, in two halves (the two halves fall on different
trials because the triple intersection — real curve + posted HR + abstract event counts for one trial —
is empty in local data; see the honest caveat below).

**Half 1 — abstract → per-arm event count (real, fresh, non-cache).** On the DAPA-HF primary publication
(McMurray et al., *N Engl J Med* 2019;381:1995-2008, [DOI](https://doi.org/10.1056/NEJMoa1911303),
PubMed PMID 31535829), `abstract_events` extracts the per-arm composite-outcome counts directly from the
sentence *"the primary outcome occurred in 386 of 2373 patients … and in 502 of 2371 patients …"* →
`events=[386,502]`, `ns=[2373,2371]`, matched to arms by N; `abstract_hr` → 0.74 (directionally
consistent, 16.3% vs 21.2%). *Per PubMed attribution requirements, this article is cited with its DOI.*
On the 161-abstract local cache the guards give **100% precision** (1 true positive — a mortality count —
0 false positives); recall is honestly low because abstracts post per-arm "X of N" counts less often than
a median or HR. A flagship counter-example: RADIANT-4's Lancet abstract (Yao et al., *Lancet*
2016;387:968-977, [DOI](https://doi.org/10.1016/S0140-6736(15)00817-X), PMID 26703889) reports only
*adverse-event* "X of N" counts, which `abstract_events` correctly returns **None** for (no fabricated
efficacy count), while `abstract_hr` recovers the posted 0.48.

**Half 2 — event count → QP → recovered HR vs truth (real curve + truth; `validate/abstract_lever_realtrial.js`).**
On RADIANT-4 (NCT01524783), holding the exact AACT curve and the posted PFS HR 0.48 (95% CI 0.35-0.67):

| reconstruction | HR | fold vs posted | inside posted 95% CI? |
|---|---:|---:|:--:|
| curve-only (no lever) | 0.679 | 1.42 | **no** (the identifiability trap) |
| censoring-informed (107/77 → QP) | 0.577 | 1.20 | **yes** |
| **abstract-HR calibrated (0.48 → calibrateHR)** | **0.481** | **1.00** | **yes** (solves 112 exp events) |

Two in-scope levers, both pulling the estimate from **outside** the posted CI to **inside** it: (a) the
per-arm **event count** fed to the QP (here from AACT participant-flow, since RADIANT-4's abstract posts no
efficacy count) → 0.577; and (b) the **abstract HR** (0.48, extracted by `abstract_hr` from the real
Lancet abstract) fed to `calibrateHR`, which solves the experimental arm's event count (→112, near the
true 107) so the Cox HR reproduces the target almost exactly → 0.481. The second is the **high-coverage**
lever (the abstract supplies an HR in 15% of trials where AACT posts none; see the coverage table below),
the first the rare-but-exact one. Locked by `test/abstract_lever.spec.js`.

**Honest caveat.** No single local trial has *all three* of {exact curve, posted HR, abstract event
counts}: RADIANT-4 has the curve + HR but its abstract gives only AE counts; DAPA-HF's abstract gives the
counts but its curve is not in local AACT data. This restates, for the abstract path, the same scarcity
the figure path found ("1 usable pair in 1500 OA PDFs", `FUSION.md`): the data sources exist and each
half is validated on real data, but the trials that publish *every* piece openly are rare. (`KMCURVE-SYNERGY.md` idea 5.)

**Real cohort coverage of the abstract levers (`harvest/abstract_events_coverage.py`).** Measured on the
155 harvested trials that have a cached PubMed abstract (no network, no AACT re-harvest), the abstract's
three in-scope levers cover very different shares:

| lever | available | share |
|---|---:|---:|
| HR → calibration / cross-check | 34 | **22%** |
| published median → cross-check | 17 | 11% |
| per-arm event count → QP lever | 1 | 1% |
| **any lever** | **42** | **27%** |

So **the PubMed abstract enriches ~27% of reconstructable trials**, but the value is overwhelmingly the
**HR** (22%) and median (11%); the per-arm "X of N" *event count* is genuinely rare (1%, and the lone hit
did not N-match an arm → 0 marginal gain here). This is the honest yield: the event-count lever is precise
and validated but seldom available, which is exactly why the unified enrichment (`abstract_enrich.py`)
captures all three — the HR carries the coverage, the event count is the occasional high-value bonus when
a trial happens to print it.

**The HR lever's *marginal* value is the real headline.** AACT posts a usable HR for only **25/155** of
these trials, while the abstract yields one for 34 — and in **23/155 (15%)** the registry posts **no HR at
all but the abstract supplies one**. That is the production payoff: for ~1 in 7 reconstructable trials the
abstract hands the reconstruction an HR (for calibration / independent cross-check) it otherwise could not
have — versus **0** trials where the abstract event count filled a gap here. The lesson the coverage run
makes quantitative: within the AACT + PubMed-abstract contract, the **abstract HR is the high-coverage
lever and the event count is the rare-but-exact one**; `abstract_enrich.py` is right to promote a confident
HR to `trial.hr` when AACT has none.

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

## Tier B (median + HR, no KM curve) — validated, with an honest shape limit (`validate/tierb_validation.js`)

Tier B fires when a trial posts a **median + HR + N + events but no Kaplan–Meier curve** (≈3,263 AACT
trials); the engine reconstructs each arm as an **exponential** parametric model. This tier had never
been tested against truth. On the 30 gold-standard datasets where the median is defined (HR and median
are *inputs*, so recovering them is circular — the honest test is **RMST**, which depends on the
survival *shape*):

- **On average the exponential is adequate**: per-arm RMST median error **~7%** (23/30 within 20%);
  RMST-difference median absolute error ~8 time-units.
- **But it has a real failure mode on strongly non-exponential survival** (early-heavy hazards): RMST
  error reaches **40–58%** on `bmt`, `veteran`, `pharmacoSmoking`, and several TCGA cancers — the
  constant-hazard assumption cannot bend.

A closed-form **2-parameter Weibull** fit from the *same* inputs (median + events/N pin the shape:
`k = ln(R/L)/ln(t_max/median)`) fixes exactly those cases — `bmt` 58%→**14%**, `pharmacoSmoking`
38%→**3%**, `veteran` 39%→**16%**, `tongue` 25%→**12%** — but is a **wash on average** (per-arm RMST
7.6% vs 6.6%; it adds variance on the easy near-exponential arms) and complicates the imposed-HR
coupling (Weibull PH needs a shared shape, `eMed = cMed/HR^{1/k}`). So **we keep the exponential as the
stable default and document the failure mode** rather than ship a change that regresses the average and
the bootstrap envelope; the Weibull is a *targeted* improvement worth adopting only with a
shape-confidence gate. Net: Tier B is a usable last resort for RMST/median when no curve is posted, with
a flagged caveat for non-exponential survival.

## The use case: IPD meta-analysis fidelity (`validate/ipd_meta_fidelity.js`)

Per-trial HR recovery is the means; the **end** is IPD meta-analysis. Does reconstruction error survive
pooling, or wash out? We pooled the **14 TCGA stage cohorts** (a real meta-analytic question — the
pooled prognostic effect of advanced vs early stage across cancers, with genuine between-cancer
heterogeneity) with a proper random-effects model — **REML τ², HKSJ confidence interval (floored),
prediction interval on t_{k−1}, log-HR pooling** — on the **true IPD** and the **QP-reconstructed**
pseudo-IPD:

| | pooled HR | 95% CI (HKSJ) | 95% prediction interval | τ² | I² |
|---|---|---|---|---|---|
| true IPD | **2.50** | [1.89, 3.30] | [1.03, 6.04] | 0.150 | 76% |
| QP-reconstructed | **2.64** | [1.91, 3.65] | [0.95, 7.40] | 0.204 | 79% |

**The pooled estimate is recovered within 6%** (2.50 vs 2.64), the confidence and prediction intervals
overlap, and reconstruction **slightly inflates the heterogeneity** (τ² 0.15→0.20, I² 76→79%). That
inflation is expected and in the *safe* direction: per-trial reconstruction error adds a little apparent
between-study variance, widening the prediction interval rather than producing a falsely precise pooled
effect. So a meta-analysis built on reconstructed pseudo-IPD reaches essentially the same conclusion as
one built on true IPD — the method is **fit for its actual purpose**, with a conservative heterogeneity
bias to disclose. (Reconstruction noise does not bias the *central* pooled effect — consistent with the
per-trial errors being roughly symmetric on the log-HR scale.)

## Remaining levers
- ✅ HR-calibration · ✅ N-matched mapping · ✅ RMST/median validation · ✅ Royston–Parmar (extrapolation).
- **Same-endpoint external median matching** to clean the contaminated registry-median cross-check.

## Systematic published-literature validation (registry-independent) — `GALLERY.md`, `harvest/pubmed_validation.py`

The registry HR shares provenance with the posted curve. To validate against a genuinely independent
source we parse the trial's **published** HR and median from its PubMed abstract (PMIDs via AACT
`study_references`; abstracts via NCBI E-utilities; deterministic endpoint-matched, covariate-guarded
extractors `harvest/abstract_hr.py` + `abstract_median.py`, with unit tests):

- **Published HR (full reconstructed cohort, endpoint-matched, high-confidence):** 20 trials, **12 with
  no registry HR at all** (pure registry-independent). Reconstructed-vs-published median fold **1.10**;
  reconstructed HR inside the published 95% CI **17/20 (85%)**. The 3 misses are verified genuine (a
  3-arm CheckMate-067 reconstruction picking the wrong arm pair; one degenerate fit; one registry-vs-
  published divergence where the interval covers the registry HR it was built from).
- **Published median (endpoint-matched):** 5 trials / 10 arm-medians, median fold **1.071**, OS and
  clean-PFS to ~3% (e.g. `NCT00861614` published 10.0/11.2 → reconstructed 10.1/11.2).
- **Uncertainty on real curves:** the reconstruction's own 95% credible interval covers the independent
  published HR **17/20 (85%)** (gold standard 28/29 on clean IPD; comparable interval width 2.56 vs 2.46).
- **Held-out-truth divergence:** the registry HR and the published HR *themselves* agree only **~80%**
  (13/16), so a residual ~1.1–1.15 fold is partly irreducible analysis-population/data-cut disagreement —
  a ceiling on what any method can score against registry truth.

Reproduce: `node validate/cohort_recon_export.js && python harvest/cohort_pubmed.py`,
`python harvest/pubmed_medians.py && node validate/pubmed_median_validation.js`,
`node validate/cohort_uncertainty_validation.js`.

## Honest limitations of this validation
- n=30 (18 with determinable direction) — small; widen by mapping multi-arm analyses pairwise.
- The AACT-internal check uses the registry HR, which can itself be model-dependent (log-rank vs Cox) and
  direction-ambiguous; it is the best *internal* proxy, not patient-level truth. The systematic
  published-literature section above is the registry-independent answer (and shows the registry HR itself
  diverges from the published HR ~20% of the time).
- Validates HR; RMST / median agreement should be validated separately (RMST is expected to be more
  robust because it depends on the curve area, which the anchors fix directly).
