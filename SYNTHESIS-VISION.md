# Where this fits in evidence synthesis — and a new way to see the space

*A positioning + development note. The claim: registry-native pseudo-IPD reconstruction with **calibrated
uncertainty** is the missing substrate that lets a single synthesis pool trials of **mixed data
granularity** without manufacturing false precision. This reframes evidence synthesis from a binary
(aggregate-data vs IPD) into a continuous, partially-identified missing-data problem — and gives this tool
a load-bearing role in it.*

## 1. The space today, and why it's stuck

Evidence synthesis is organised by the granularity of the data each trial shares:

| input a trial shares | what you can identify | method |
|---|---|---|
| 2×2 table / a single effect (logHR, SE) | pooled effect, τ², PI | aggregate-data (AD) meta-analysis |
| a Kaplan–Meier **curve** | + median, RMST, S(t); HR only *partially* (censoring unobserved) | pseudo-IPD reconstruction (Guyot 2012) |
| curve + **numbers-at-risk / events** | + HR identified, non-PH, time-varying | reconstruction → IPD-like |
| full **individual patient data** | everything + within-trial interactions, deconfounded subgroups | IPD meta-analysis (the gold standard) |

The field treats this as a **binary**: a review is "AD" or "IPD". IPD is unobtainable for most trials
(governance, DUAs, months–years, cost), so almost every review collapses to the **lowest common
denominator** — an AD meta-analysis that throws away everything the richer trials posted (curves, medians,
time-resolved structure). The trials that posted more get pooled as if they posted less.

## 2. The reframing: one missing-data problem on a granularity manifold

Drop the binary. Every trial sits on a **data-completeness manifold**: it publishes some statistics and
leaves the patient-level data **latent**. What's missing (the event/censoring split, the event times, the
covariates) is **partially identified** by what's published — and *which estimand is identified depends on
which statistics were posted*:

- curve alone → median, RMST, S(t) **point-identified**; **HR partially identified** (an identification
  region, because every censoring level passes through the same survival points — proven structural, see
  `validate/advanced_estimators.js`);
- curve + total events → **HR point-identified** (the Titman-QP result, fold ≈1.05);
- curve + number-at-risk → fully identified (Guyot);
- abstract HR alone → logHR point-identified, nothing time-resolved.

So **reconstruction is not "faking IPD". It is principled multiple imputation of the unobserved
patient-level data conditional on the posted statistics**, carrying an *honest identification region* for
whatever the posted statistics leave open. This tool already does the hard part: it produces the
reconstruction **and a calibrated credible interval for the under-identified censoring level** (coverage
of the true HR 28/29 on real IPD; `uncertainty_calibration.js`).

That single capability — *reconstruction with calibrated uncertainty* — is what lets a curve-only or
HR-only trial join an **IPD-level model** as a **partially-identified IPD trial**, contributing exactly as
much information as its posted statistics support, no more.

## 3. The powerful role: "one model, many granularities"

The target is a single synthesis in which **each trial contributes a likelihood at its own granularity**:

- IPD trials → the full patient-level likelihood (gold standard, where available);
- curve trials → the reconstructed pseudo-IPD likelihood, with the reconstruction (censoring-level)
  **uncertainty propagated**, not discarded;
- HR-/abstract-only trials → a logHR likelihood (the in-scope abstract HR lever already supplies this for
  ~15% of trials the registry omits — `harvest/abstract_enrich.py`).

The estimand (a flexible-parametric or fractional-polynomial **survival NMA**, RMST contrasts, a
time-varying HR) is then fit across **all** trials at once, each weighted by what it actually identifies.
This dissolves the AD/IPD binary and **expands the evidence base without inflating false precision** —
which is the whole point, and the part that must be earned statistically (§4).

In the methods landscape this is a concrete extension of three established lines:
- **Phillippo's multilevel network meta-regression (ML-NMR)** combines IPD + AD in one model by
  integrating over covariate distributions — but for *binary/continuous* outcomes. Reconstruction-with-UQ
  supplies the missing piece to extend ML-NMR to **time-to-event**: it turns an AD survival trial into a
  partially-identified IPD survival trial.
- **Jansen's fractional-polynomial survival NMA** needs patient-level-like data and today uses
  digitised pseudo-IPD **as if it were exact** — ignoring reconstruction uncertainty and so
  over-confident. Propagating the calibrated interval is the fix.
- **Manski-style partial identification**: each trial contributes an *identified set*, not a point; the
  synthesis is over sets-with-priors. This is the rigorous name for what the calibrated interval is.

## 4. The linchpin, proven: propagate the reconstruction uncertainty or distort the synthesis

Adding reconstructed trials to a pool is only honest if the reconstruction variance is carried through.
A reconstructed logHR has **two** variance components: the trial's sampling variance `s²` **and** the
reconstruction (censoring-level) variance `r²` that the credible interval already quantifies. Pooling on
`s²` alone — treating pseudo-IPD as exact IPD — is the same error the lab's own advanced-stats rules flag
for multiverse pooling (*"CIs collapse below truth"*).

`validate/honest_pooling_sim.js` (seeded Monte-Carlo, 4000 reps, k=12, true τ²=0.05) proves the effect and
the fix — and the result is more subtle than "the CI gets too narrow":

| pooling | mean τ² (truth 0.05) | pooled-mean coverage | PI width |
|---|---:|---:|---:|
| true IPD (gold) | 0.045 | 0.95 | 0.94 |
| **NAIVE** (ignore `r²`) | **0.089 (≈2×)** | 0.94 | **1.34** |
| **HONEST** (Rubin `s²+r²`) | **0.044** | 0.95 | 0.92 |

**The pooled mean survives either way** — REML τ² silently *absorbs* the ignored reconstruction noise — but
that is exactly the danger: NAIVE pooling **mis-reads reconstruction noise as between-trial heterogeneity**,
roughly **doubling τ²** and **inflating the prediction interval ~40%**. You would conclude the trials
disagree far more than they do — corrupting the PI, I², subgroup and meta-regression decisions (the lab's
PI-Atlas territory). Propagating `r²` via Rubin's rules recovers the true τ² and a calibrated PI. The
residual honest cost — reconstructed point estimates are simply noisier, slightly degrading PI coverage
(0.83 vs 0.89) — is now *visible* instead of hidden inside a falsely-large τ². Locked by
`test/honest_pooling.spec.js`.

**This is the constructive counterpart to the false-robustness work** (`spec-collapse-atlas`,
`evidence-observatory`): those projects *detect* synthesis that is more robust than the data warrant; this
shows how to *add* evidence (curve- and HR-only trials) **without** committing that sin — additive
evidence with honest uncertainty.

## 4b. Phase 2 — the linchpin on REAL reconstructions, and a sharper finding

`validate/phase2_real_pooling.js` repeats the three-way pool on the **14 TCGA stage cohorts where we hold
the true IPD**, replacing the simulation's `r²` with each trial's *real* reconstruction variance, measured
from the engine's multiple-imputation ensemble (`reconstructEnsemble`, M=200): `r_i² = ((ln hi − ln lo)/
(2·1.96))²` from the HR credible interval.

| pool | pooled HR (true 2.49) | τ² (true 0.131) | PI width (log) |
|---|---:|---:|---:|
| true IPD | 2.49 | 0.131 | 1.66 |
| **NAIVE** (ignore `r²`) | 2.63 | **0.176 (overstates)** | 1.92 |
| **HONEST** (Rubin `s²+r²`) | **2.51** | **0.077 (understates)** | 1.35 |

Real data sharpens the lesson. The **honest pooled HR is recovered** (2.51 vs 2.49; naive drifts to 2.63),
and **naive overstates heterogeneity** exactly as the simulation predicted. But here reconstruction
uncertainty is **huge** — mean per-trial variance inflation **5.4×**, because these are heavily-censored
*curve-only* cohorts — and at that magnitude **full Rubin propagation over-corrects: it absorbs the genuine
between-cancer heterogeneity into the large within-trial variance, *under*-stating τ²**. So both extremes
are biased when `r²` dominates `s²`: naive reads reconstruction noise as heterogeneity; honest reads
heterogeneity as reconstruction noise. The truth needs `r²` to be *moderate*.

**This unifies the whole project.** The value of the censoring lever — the Titman-QP event count, the
abstract event count, the NAR fusion — is not only per-trial HR accuracy; it is that **each lever shrinks
`r²`**, moving a curve-only trial out of the "uninformative about heterogeneity" regime into the range
where it can honestly contribute to a synthesis's τ²/PI. A curve with no censoring information is nearly
weightless for heterogeneity (`r²` swamps the signal); the same curve plus an event count or at-risk table
becomes a real synthesis participant. The reconstruction-uncertainty magnitude is the currency, and the
lever is how you earn it. Locked by `test/phase2_pooling.spec.js`.

## 4c. Phase 2b — the lever shrinks r², and the honest corrective

`validate/phase2b_lever_shrinks_r2.js` measures `r²` on the same 14 cohorts in two regimes from the
ensemble: **curve-only** (no event count → full censoring band) vs **event-pinned** (the engine's new
`reconstructEnsemble({pinEvents:true})`, which samples a tight ±5% band when a total-event count is posted).

**The mechanism is confirmed: the event count shrinks the reconstruction SD 5.2× (0.219 → 0.042)**, for
≥80% of cohorts individually. That is the lever turning a near-weightless curve-only trial into an
informative one.

But pooling honest under each `r²` surfaced a corrective that *deepens* the framing rather than confirming
the tidy story:

| honest pool | τ² (true 0.131) | why |
|---|---:|---|
| curve-only `r²` (large) | 0.122 | large within-trial variance **absorbs** the spread |
| event-pinned `r²` (tiny) | 0.173 | tiny within-trial variance → residual spread → **τ² overstated** |

**Shrinking the variance is not enough.** Even with the event count, the reconstructed point estimates
retain a residual ~5–6% per-trial error (pooled HR 2.63 vs true 2.49 — the documented QP fold ≈1.05), and
when `r²` is driven small that residual **bias** is mis-read as heterogeneity, *over*-stating τ². The
ensemble's variance captures the censoring *uncertainty* but not the censoring *bias*. So the honest
synthesis needs both: the lever to shrink `r²`, **and** the reconstruction treated as a partially-identified
**set with a possibly-off centre** (the Manski point in §3), not a point-with-variance. This is the real,
non-obvious lesson — and it sets the agenda: model the reconstruction bias as an identified-set offset, not
just inflate the variance. Locked by `test/phase2b_lever.spec.js`.

## 4d. Phase 2c — reconstruction bias is partially-identified heterogeneity

Phase 2b's corrective has a precise statistical name. A **deterministic per-trial reconstruction bias `b_i`
is observationally identical to between-trial heterogeneity** — "this trial reconstructed 8% high" and
"this trial truly differs by 8%" leave the same footprint in the data. So **τ² is not point-identified from
reconstructed trials; it is partially identified** (Manski). The honest object is an *identified set*,
calibrated against a true-IPD gold standard (which the project has).

`validate/phase2c_bias_offset.js` does exactly this, all **leave-one-out** (each held-out cohort calibrated
from the other 13, so it is out-of-sample): the systematic offset `β = mean(e)` is identifiable, so we
**de-bias** it; the residual per-trial bias is bounded `|b_i| ≤ 1.64·SD(e)` and the pooled effect and τ²
range over all configurations inside those boxes.

| quantity | true IPD (held out) | naive event-pinned POINT | de-biased POINT | **identified SET** | set ∋ truth? |
|---|---:|---:|---:|---:|:--:|
| pooled HR | 2.49 | 2.63 | **2.44** | **[1.59, 3.75]** | **yes** |
| τ² | 0.131 | 0.173 | 0.182 | **[0, 0.765]** | **yes** |

Two clean results. **De-biasing the LOO systematic offset (0.075 log-HR ≈ 7.8%) recovers the pooled HR**
(2.63 → 2.44 vs true 2.49). And **the calibrated identified set brackets the truth on both axes**, where
the naive point is a precise-but-wrong single number. The τ² set is wide ([0, 0.765]) — that width is the
*honest price* of residual reconstruction bias: it is exactly how much the synthesis cannot resolve, and it
is what the lever (Phase 2b, shrinking the per-trial error) and a better reconstruction narrow. Reporting
the set, not a false point, is the correct treatment. Locked by `test/phase2c_set.spec.js`.

**The arc closes.** Phase 1: propagate the reconstruction variance (Rubin) — necessary. Phase 2: on real
data, when `r²` is large it absorbs τ²; the lever shrinks `r²`. Phase 2b: the lever shrinks `r²` 5.2×, but
variance ≠ bias. Phase 2c: the residual bias is partially-identified heterogeneity, handled by a
gold-standard-calibrated identified set that brackets the truth. A reconstructed trial joins a synthesis as
a **de-biased point with a calibrated identification interval** — not a fake-exact IPD row, and not a
discarded curve.

## 4e. Phase 3 — granularity-mixed synthesis, and the evidence-completeness curve

The capstone (`validate/phase3_granularity_mixed.js`): one pooled survival contrast over a corpus whose
trials sit at **different granularities**, each contributing what it identifies — IPD and HR-only trials as
**identified points** `(logHR, s²)`, reconstructed-curve trials as the **Phase-2c object** (de-biased
point + inflated variance + residual-bias half-width `δ`). The pool is REML/HKSJ; because curve trials
carry identification intervals, the output is an **identified set** (the pooled HR and τ² range over the
curve trials' bias boxes; IPD/HR-only are fixed). All curve calibration is leave-one-out.

**A realistic mixed corpus (5 IPD · 5 curve · 4 HR-only)** pools to HR **2.57, identified set [2.28, 2.92]**
— which **brackets the all-IPD truth (2.49)** — with a τ² set [0.07, 0.51]. A synthesis that *mixes*
granularities gives an honest answer with an interval that reflects exactly how much the low-granularity
trials cost.

And the new artefact — the **evidence-completeness curve**: sweep the fraction of trials that are curve-only
(rest IPD) and the identified-set width grows **smoothly from a point to its widest**, while the central HR
stays put:

| curve-only fraction | pooled HR | HR identified set | set width |
|---:|---:|---:|---:|
| 0.00 (all IPD) | 2.49 | [2.49, 2.49] | **0** |
| 0.25 | 2.46 | [2.13, 2.86] | 0.74 |
| 0.50 | 2.42 | [1.93, 3.12] | 1.19 |
| 0.75 | 2.51 | [1.80, 3.57] | 1.77 |
| 1.00 (all curve) | 2.44 | [1.59, 3.75] | **2.15** |

This is the reframing made quantitative: you **can** add trials that posted only a curve (or only an HR) to
a synthesis, and the model **honestly widens the interval** to price their reduced information instead of
faking precision — the central estimate is stable, the *uncertainty* tracks the evidence granularity. The
curve also reads as guidance: it says how much precision a review buys by obtaining IPD / a posted HR vs
working from curves, and (with Phase 2b) how much each censoring lever recovers. Locked by
`test/phase3_mixed.spec.js`.

This realises the §3 vision — "one model, many granularities" — end to end on real reconstructions: a
reconstructed trial is neither a discarded curve nor a fake-exact IPD row, but a **de-biased point with a
calibrated identification interval**, pooled alongside IPD and HR-only trials in a single honest synthesis.

## 5. Development roadmap

- **Phase 1 — honest pooling (DONE).** Rubin's-rules propagation; Monte-Carlo proof it recovers τ²/PI where
  naive pooling distorts them (`honest_pooling_sim.js`).
- **Phase 2 — real reconstructions (DONE, §4b).** Same three-way pool on the 14 true-IPD TCGA cohorts with
  real ensemble `r²`. Finding: honest recovers the pooled HR and beats naive, but when `r²` dominates `s²`
  (heavily-censored curve-only trials) τ² is under-identified from either side — the censoring lever is
  what shrinks `r²` into the usable range.
- **Phase 2b — lever-shrinks-r² (DONE, §4c).** The event count pins the censoring, shrinking reconstruction
  SD **5.2×** — but it surfaced an honest corrective: variance propagation is necessary, not sufficient.
- **Phase 2c — reconstruction bias as a partially-identified set (DONE, §4d).** De-bias the LOO systematic
  offset (recovers the pooled HR) and report the τ² / pooled-effect identified set over residual-bias boxes;
  validated to bracket the held-out truth where a naive point does not.
- **Phase 3 — granularity-mixed synthesis (DONE, §4e).** One pooled contrast over {IPD, reconstructed-curve
  as a de-biased point + identification interval, HR-only}; the mixed corpus brackets the all-IPD truth, and
  the evidence-completeness curve quantifies the precision↔granularity trade-off.
- **Phase 3b — extend to a survival NMA** (multiple treatments, indirect comparisons): carry the §4e
  per-trial identified-interval contribution through a network model with consistency checks (the Jansen
  fractional-polynomial NMA extended with propagated UQ + partial identification). Reuse the lab's existing
  NMA machinery rather than re-implement (see `allmeta`/capsule node-split + design-by-treatment).
- **Phase 4 — an evidence-completeness atlas.** For a real review question, harvest every trial, classify
  each by posted-statistics granularity, compute a per-trial *identification/information score* for the
  target estimand, and map what fraction of the evidence is point- vs partially-identified. A new artefact:
  the synthesis "information map" of a question, before any pooling.
- **Phase 5 — position formally** vs ML-NMR (time-to-event extension), Jansen survival NMA (uncertainty
  propagation), and Manski partial identification (the identified-set formalism).

## 5b. Reuse map — the portfolio already has the synthesis engines

Phase 3b should **wire into existing, verified machinery in the portfolio**, not re-implement it. Surveyed
and path-verified:

| need | existing engine | path / entry point | how registry-IPD plugs in |
|---|---|---|---|
| honest pooling of under-identified specs | **spec-collapse-atlas** (Py) | `spec_collapse/aggregators.py::weighted_likelihood(specs, cl, weights)` | the canonical tool my Phase 1 (Rubin) / Phase 2c (set) approximate. Map each trial to a spec `{theta: logHR, var: s²+r²+(δ/z)², k}`; its mixture-CDF interval is never narrower than a single spec (no false robustness) |
| random-effects NMA + consistency | **advanced-nma-pooling** (Py) | `src/nma_pool/models/core_ad.py::ADNMAPooler.fit`; `validation/inconsistency.py::design_by_treatment_test`, `node_splitting_diagnostics` | feed reconstructed-curve trials as contrast-level inputs with the Phase-2c inflated variance |
| **ML-NMR** (mixed IPD + AD) | **advanced-nma-pooling** (Py) | `src/nma_pool/models/ml_nmr.py::MLNMRPooler.fit` | the §3 "extend ML-NMR to time-to-event" is largely *already built*; the new input is registry-IPD's reconstructed pseudo-IPD-with-UQ |
| **non-PH survival NMA** | **advanced-nma-pooling** (Py) | `src/nma_pool/models/survival_nph.py::SurvivalNPHPooler.fit` (piecewise-exponential, interval-specific effects) | reconstructed pseudo-IPD → per-interval events/at-risk as a survival-AD trial |
| **fractional-polynomial survival NMA** (Jansen) | **allmeta** (JS) | `HTA/src/engine/fpNMA.js::fitFPCoefficients(times, logHRs, weights, powers)` | encode the partial-ID variance in the weight `w = 1/(var+δ²)` |
| browser IPD+AgD NMA, forest/league/network UI | **IPDNMA** (JS) | `IPDNMA/index.html` | validation / visualization reference |

The key realisation: the vision's "extends ML-NMR (Phillippo) and Jansen survival NMA" is not aspirational —
**those engines exist in the portfolio**, and registry-IPD supplies the missing input they were never given:
a curve- or HR-only trial as a *de-biased point with a calibrated identification interval*. Phase 3b is an
integration, and `weighted_likelihood` is the honest aggregator to route the mixed-granularity contributions
through (a cross-check against this repo's `metaRE` is the natural first step, despite the Py↔JS boundary).

## 6. Honest risks (carry forward)

- Reconstruction **bias** (not just variance) is not always imputed away; some estimands stay
  *partially identified* (regions, not points). Report sets, never false points.
- Coverage / false-robustness of the mixed model must be **validated, not assumed** — the same discipline
  the linchpin sim applies, on every new estimand.
- Coverage is still the binding limit for full-IPD estimands (§ `PAPER.md` "usable for meta-analysis?"):
  this widens the *contributable* set and makes the pooling honest, but it does not conjure trials that
  posted nothing.
