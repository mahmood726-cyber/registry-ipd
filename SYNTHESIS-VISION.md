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

## 5. Development roadmap

- **Phase 1 — honest pooling (DONE, this note).** Rubin's-rules propagation of the reconstruction credible
  interval into the meta-analytic variance; Monte-Carlo proof that it recovers τ²/PI where naive pooling
  distorts them. Next: run it on the 14-cohort true-IPD set (`ipd_meta_fidelity.js`) using the real
  per-trial ensemble intervals, not a simulation.
- **Phase 2 — granularity-mixed synthesis.** One random-effects model ingesting {IPD, reconstructed-pseudo-
  IPD-with-UQ, HR-only} trials, each contributing its identified information; start pairwise, then a
  survival NMA with propagated reconstruction uncertainty (the Jansen extension).
- **Phase 3 — an evidence-completeness atlas.** For a real review question, harvest every trial, classify
  each by posted-statistics granularity, compute a per-trial *identification/information score* for the
  target estimand, and map what fraction of the evidence is point- vs partially-identified. A new artefact:
  the synthesis "information map" of a question, before any pooling.
- **Phase 4 — position formally** vs ML-NMR (time-to-event extension), Jansen survival NMA (uncertainty
  propagation), and Manski partial identification (the identified-set formalism).

## 6. Honest risks (carry forward)

- Reconstruction **bias** (not just variance) is not always imputed away; some estimands stay
  *partially identified* (regions, not points). Report sets, never false points.
- Coverage / false-robustness of the mixed model must be **validated, not assumed** — the same discipline
  the linchpin sim applies, on every new estimand.
- Coverage is still the binding limit for full-IPD estimands (§ `PAPER.md` "usable for meta-analysis?"):
  this widens the *contributable* set and makes the pooling honest, but it does not conjure trials that
  posted nothing.
