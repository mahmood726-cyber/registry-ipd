# Registry-native reconstruction of survival pseudo-IPD from ClinicalTrials.gov, with calibrated uncertainty

*Working manuscript draft. Companion artifacts: `METHODS.md` (method detail), `VALIDATION.md` (full
results), live tool + dashboard at https://mahmood726-cyber.github.io/registry-ipd/. Citations below
should be DOI-resolved programmatically before submission.*

## Abstract

**Background.** Reconstructing individual-patient time-to-event data (IPD) from published Kaplan–Meier
(KM) curves underpins much of survival evidence synthesis, but every established tool (Guyot's
iterative method, IPDfromKM, and recent automated pipelines KM-GPT and RESOLVE-IPD) begins by
*digitising a figure image*, incurring pixel error and requiring access to the publication. We ask
whether usable pseudo-IPD can instead be reconstructed from **ClinicalTrials.gov / AACT structured
summary data alone**, fully offline, and we quantify — against true patient-level data — exactly which
estimands are trustworthy.

**Methods.** From the 2026-06 AACT snapshot we harvest, per arm, the structured KM-estimate timepoints,
number-at-risk, total events (and, via participant-flow `drop_withdrawals`, event/censoring/competing
counts), and the reported hazard ratio (HR). Reconstruction is tiered by data richness and, for the
rich tier, selects per trial between a faithful Guyot inverse-KM and a censoring-informed anchor-exact
method by minimum 1-Wasserstein distance to the registry anchors. We add: (i) **multiple-imputation
uncertainty** that samples the under-identified censoring level to produce calibrated credible
intervals; (ii) **competing-risks** reconstruction with the Aalen–Johansen cumulative-incidence
estimator; (iii) **HR-calibration** that imposes the reported HR for downstream IPD meta-analysis;
(iv) **fractional-polynomial time-varying-HR** analysis for non-proportional hazards. Validation
proceeds up a ladder of increasing independence: AACT-internal HR, primary publications, and finally
**true patient-level IPD** from eight open RCT datasets (R `survival::` and related).

**Results.** AACT contains **zero structured number-at-risk** rows; ~hundreds of trials post a KM-estimate
curve. Against true patient-level data across 7 adequately-sized RCTs (≥100/arm), curve-only
reconstruction recovers the HR to a median fold-error of **1.12 (~11%)** and the median to **~3–6%**;
RMST is recovered to **~2%**. The multiple-imputation 95% credible interval covers the **true HR in
7/7** datasets (median width ~2.2×). Reconstructed Aalen–Johansen CIFs match the true CIFs to
**0.3–1.6 percentage points** (`survival::colon`). Accuracy rises sharply with the number of posted KM
timepoints and **plateaus at ≥5–6** (HR fold-error 1.33 at K=3 → ~1.09 at K=5). Very small trials
(N≈137) do not reconstruct reliably.

**Conclusions.** Registry-native, image-free reconstruction is **good enough for RMST/median-based
survival synthesis** of adequately-sized trials and provides a useful, uncertainty-quantified HR
triangulation input — without digitisation error and with full provenance to the registry record. It
is bounded honestly by registry coverage and by the under-identified censoring level, which our
calibrated intervals make explicit.

## 1. Introduction

Pseudo-IPD reconstruction from KM curves (Guyot et al. 2012) is standard in HTA and meta-analysis;
IPDfromKM (Liu et al. 2021) and automated pipelines (KM-GPT 2025; RESOLVE-IPD 2025) all operate on the
*figure image*. This injects digitisation error and presumes access to the publication's plot. We
invert the premise: ClinicalTrials.gov, via the AACT relational mirror, exposes structured
KM-estimate measurements, participant flow, and reported effect sizes for tens of thousands of trials.
Reconstructing from these **exact registry values** removes digitisation error and yields provenance to
the NCT record — but raises a different question we answer empirically: *which reconstructed estimands
are trustworthy, and how trustworthy?*

## 2. Data

The AACT pipe-delimited snapshot (2026-06-01) holds 76,067 trials with posted results. Survival-relevant
tables: `outcome_measurements` (KM-estimate timepoints, median), `outcome_analyses` (HR/CI/method),
`outcome_counts` (per-arm N), `result_groups`/`milestones`/`drop_withdrawals` (arms, randomised counts,
withdrawal reasons). Two empirical facts drove the design: **(i)** AACT contains *no* structured
number-at-risk; **(ii)** survival is frequently posted as cumulative *incidence* ("probability of
event"), requiring data-driven orientation to a survival scale.

## 3. Methods

**Harvesting.** A Python harvester maps AACT → a per-trial JSON. Robustness measures forced by real
failures: filtering result-groups by `outcome_id`; data-driven survival/incidence orientation;
N-matched mapping of participant-flow groups to outcome arms via milestone `STARTED` counts.

**Tiered reconstruction.** Tier A (KM curve ≥3 timepoints + N) runs two estimators and keeps the lower
1-Wasserstein fit: a faithful port of Guyot's iterative inverse-KM, and a censoring-informed
anchor-exact estimator (RESOLVE-IPD CEN-KM style) that holds at-risk constant within intervals so the
reconstructed curve passes through the registry anchors. Tier B (median + HR + N) uses an exponential
parametric model with a seeded bootstrap envelope; Tier C fails closed.

**Self-audit.** Nine checks (event count, anchor fidelity, median, Cox-derived HR with ridge penalty,
monotonicity, number-at-risk consistency, population conservation, follow-up sanity, and a hard HR
**direction-integrity** check) yield a Bronze/Silver/Gold badge; the HR is never inverted.

**Uncertainty (novel).** The registry curve does not uniquely determine the IPD — the censoring level
is under-identified. We sample that degree of freedom (a maximum-entropy stance) plus method and
anchor-rounding to produce an imputation ensemble and percentile **credible intervals** for HR,
median and RMST.

**Competing risks.** With participant-flow reason counts, cause-labeled pseudo-IPD is reconstructed and
the **Aalen–Johansen** cause-specific cumulative incidence computed, in place of the biased 1−KM.

**HR-calibration and time-varying HR.** A 1-D solve on the experimental censoring level imposes the
reported HR (for IPD-MA consistency); a piecewise-exponential rate-ratio with a first-order
fractional-polynomial fit (Royston/Jansen style) estimates HR(t) and tests proportional hazards.

## 4. Validation

A ladder of increasing independence (full numbers in `VALIDATION.md`):

1. **AACT-internal** (registry HR as held-out truth, 30 two-arm trials): HR within the registry CI
   83→94% (curve-only→censoring-informed), median fold-error ~1.1.
2. **Primary publication** (RADIANT-4, Yao et al. *Lancet* 2016): reconstructed HR 0.47–0.48 vs
   published 0.48; median 11/4 vs 11.0/3.9 months.
3. **True patient-level IPD**, 8 open RCTs (GBSG, Rotterdam, PBC, diabetic retinopathy, NWTSG Wilms,
   myeloid AML, kidney-transplant, Veterans lung). For the 7 adequately-sized trials: HR median
   fold-error **1.12**, median **~3–6%**, RMST **~2%**; large effects recovered cleanly (Wilms 5.1→5.4).
4. **Uncertainty coverage**: the 95% credible interval covers the **true HR 7/7**.
5. **Competing-risks gold standard** (`survival::colon`): reconstructed AJ CIF within **0.3–1.6 pp** of
   the true CIF.
6. **Anchor-density operating curve**: accuracy plateaus at **≥5–6** posted timepoints.

All reconstruction code is unit-tested (deterministic estimators to 1e-6; stochastic to seeded
Monte-Carlo tolerances) with a headless-browser smoke test of the offline tool.

## 5. Discussion

The reconstruction is **accurate where theory predicts**: RMST and median are curve-derived and recover
to single-digit percent; the HR depends on event/censoring timing the registry under-reports, and is a
~11% triangulation input best used with the calibrated interval. The recurring practical bottleneck was
not the statistics but **registry group/endpoint mapping**; N-matching resolved the main case. Two
implications follow: for **users**, trust the reconstruction when ≥5–6 KM timepoints are posted and the
trial is adequately sized, and read the credible interval; for **registries**, posting a handful of
KM-estimate timepoints makes a trial reconstructable. Unlike digitisation tools, the method has zero
anchor digitisation error and full registry provenance, at the cost of applying only to the subset of
trials that post structured survival data.

## 6. Limitations

Registry coverage is the binding limit (hundreds of trials, not all). The censoring level is
under-identified; we surface this as honest interval width rather than a false point. Very small trials
(N≈137) do not reconstruct. The external-median check is sensitive to endpoint matching. True-IPD
validation used eight open datasets; credentialed repositories (Vivli, Project Data Sphere) would
extend it. Tier B is exponential-only.

## 7. Availability

Code, tests, harvester, validation scripts, offline tool and dashboard: MIT-licensed at
https://github.com/mahmood726-cyber/registry-ipd (live: https://mahmood726-cyber.github.io/registry-ipd/).

## References (DOI-resolve before submission)

1. Guyot P, Ades AE, Ouwens MJNM, Welton NJ. *BMC Med Res Methodol.* 2012;12:9. doi:10.1186/1471-2288-12-9 *(PubMed-verified; its own recommendation that RCTs report numbers-at-risk + total events alongside KM curves directly supports our censoring-informed method and anchor-density finding.)*
2. Liu N, Zhou Y, Lee JJ. IPDfromKM. *BMC Med Res Methodol.* 2021;21(1):111. doi:10.1186/s12874-021-01308-8 *(PubMed-verified.)*
3. Royston P, Parmar MKB. Flexible parametric proportional-hazards and proportional-odds models for censored survival data. *Stat Med.* 2002;21(15):2175–2197. doi:10.1002/sim.1203 *(PubMed-verified.)*
4. Aalen OO, Johansen S. An empirical transition matrix for non-homogeneous Markov chains based on censored observations. *Scand J Stat.* 1978;5(3):141–150. *(Foundational competing-risks paper; not PubMed-indexed.)*
5. Jansen JP. Network meta-analysis of survival data with fractional polynomials. *BMC Med Res Methodol.* 2011;11:61. doi:10.1186/1471-2288-11-61 *(PubMed-verified.)*
6. Rubin DB. *Multiple Imputation for Nonresponse in Surveys.* Wiley, 1987. *(Book; no DOI.)*
7. Yao JC, et al. RADIANT-4. *Lancet.* 2016;387:968–977. doi:10.1016/S0140-6736(15)00817-X
8. RESOLVE-IPD. arXiv:2511.01785 (2025). · KM-GPT. arXiv:2509.18141 (2025).
9. AACT / Clinical Trials Transformation Initiative, ClinicalTrials.gov. Validation datasets: R
   `survival` package (Therneau) and Rdatasets mirror.
