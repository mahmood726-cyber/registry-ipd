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
**true patient-level IPD** from 31 open RCT/cohort datasets (R `survival::`, `KMsurv`, `asaur` and related).

**Results.** Of the 76,067 AACT trials with posted results, **zero** contain a structured
number-at-risk row, and only **288–~600** (0.4–0.8%, depending on detection strictness) post a
reconstructable structured KM curve — the binding coverage limit, quantified by census
(`census_full_aact.py`). Against true patient-level data across **17 adequately-sized RCTs/cohorts** (≥100/arm; of 31
real datasets), curve-only reconstruction recovers the HR to a median fold-error of **1.12 (~12%;
13/17 within 20%)** and the median to **~3%**; RMST to **~2%**. The multiple-imputation 95% credible
interval covers the **true HR in 16/17 (94%)** of the uncertainty-validation datasets (median width
~2.3×; the single miss is `bfeed`, the discrete-time outlier) — empirical coverage matching the
nominal 95%. Reconstructed Aalen–Johansen
CIFs match the true CIFs even under heavy competing risk (`aidssi`: naive 1−KM overstates the AIDS
incidence by 16 pp, AJ recovers truth within 6 pp). Accuracy rises sharply with posted KM timepoints
and **plateaus at ≥5–6** (HR fold-error 1.40 at K=3 → 1.15 at K=5 → 1.08 by K=12). Very small trials
(N≈137) do not reconstruct reliably.

**Conclusions.** Registry-native, image-free reconstruction is **good enough for RMST/median-based
survival synthesis** of adequately-sized trials and provides a useful, uncertainty-quantified HR
triangulation input — without digitisation error and with full provenance to the registry record. It
is bounded honestly by registry coverage and by the under-identified censoring level, which our
calibrated intervals make explicit.

## 1. Introduction

Pseudo-IPD reconstruction from KM curves (Guyot et al. 2012) is standard in HTA and meta-analysis;
IPDfromKM (Liu et al. 2021), ipdfc (Wei & Royston 2017) and the recent automated pipelines KM-GPT
(2025) and RESOLVE-IPD (2025) all operate on the *figure image*. This injects digitisation error and
presumes access to the publication's plot. We invert the premise: ClinicalTrials.gov, via the AACT
relational mirror, exposes structured KM-estimate measurements, participant flow, and reported effect
sizes for tens of thousands of trials. Reconstructing from these **exact registry values** removes
digitisation error and yields provenance to the NCT record — but raises a different question we answer
empirically: *which reconstructed estimands are trustworthy, and how trustworthy?*

**Related work and positioning.** Since Guyot et al. (2012), survival-IPD reconstruction has relied on
digitising the plotted KM curve — reading (time, survival) coordinates off the figure with software
such as DigitizeIt and inverting the KM equations using numbers-at-risk transcribed from the figure's
risk table. This figure-digitisation paradigm underlies IPDfromKM and ipdfc, the automated pipelines
RESOLVE-IPD and KM-GPT, and the only reconstruction route in HTA methods guidance (NICE DSU TSD 19).
Independently, the ClinicalTrials.gov ecosystem has been harvested at scale for *aggregate* evidence
synthesis (the AACT database, Tasneem et al. 2012; CT.gov results knowledge graphs) — but these efforts
extract efficacy/safety fields and never reconstruct time-to-event IPD. The closest tabular-input
neighbour is Titman (2026), who reconstructs pseudo-IPD — including competing-risks data from
cumulative-incidence functions — from *published* numbers-at-risk and marked censoring times via
quadratic programming, but draws those tables from journal articles rather than registry-posted
structured results. **To our knowledge this is the first method to reconstruct survival IPD natively
from the structured survival tables ClinicalTrials.gov/AACT exposes, bypassing figure digitisation**;
we cite Titman (2026) as concurrent tabular-input work and claim novelty on data *provenance*
(registry-posted tables), not merely on accepting tabular input. A systematic prior-art search
supporting this positioning is recorded in `NOVELTY.md`.

## 2. Data

The AACT pipe-delimited snapshot (2026-06-01) holds 76,067 trials with posted results. Survival-relevant
tables: `outcome_measurements` (KM-estimate timepoints, median), `outcome_analyses` (HR/CI/method),
`outcome_counts` (per-arm N), `result_groups`/`milestones`/`drop_withdrawals` (arms, randomised counts,
withdrawal reasons). Two empirical facts drove the design: **(i)** AACT contains *no* structured
number-at-risk (**0 rows** across all 76,067 results-trials); **(ii)** survival is frequently posted
as cumulative *incidence* ("probability of event"), requiring data-driven orientation to a survival
scale. A full-snapshot census (`census_full_aact.py`) sizes the reconstructable population: **288**
trials post a survival curve under strict detection (Kaplan-Meier/survival/PFS/EFS at ≥3 timepoints)
and **≈514–605** under a broader net (adding disease-free survival, cumulative incidence); a further
**3,263** are parametrically reconstructable from a posted median + HR (Tier B). The median
curve-posting trial reports only **3–4** timepoints, and only **≈34%** post the **≥5–6** needed for
reliable reconstruction (see §4.6) — motivating the registry-reporting recommendation in `POLICY.md`.

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
3. **True patient-level IPD**, 31 open datasets (breast/colon/lung/AML/melanoma/leukemia/transplant/
   PBC/MGUS/NAFLD/prostate/retinopathy/AIDS/larynx/burn/pneumonia/HCC; R `survival::`, `KMsurv`,
   `asaur`). For the 17 adequately-sized (≥100/arm): HR median fold-error **1.12 (13/17 within 20%)**,
   median **~3%**, RMST **~2%**; large effects clean (Wilms 5.1→5.2, melanoma 4.4→4.0), classic Gehan
   6-MP RCT 0.22→0.20. The worst case is `bfeed` (fold 1.75) — breastfeeding duration in discrete
   weeks with ~96% events, a heavily-tied discrete-time series rather than the smooth KM curve the
   method targets; retained as an honest out-of-favour boundary.
4. **Uncertainty coverage**: the 95% credible interval covers the **true HR 16/17 (94%)** (median width
   2.3×); empirical coverage matches the nominal 95%, with `bfeed` the sole miss.
5. **Competing-risks gold standard**: reconstructed AJ CIF within ~1 pp of truth where competing risk
   is rare (`survival::colon`) and recovers truth within 6 pp where it is heavy (`aidssi`: naive 1−KM
   overstates AIDS incidence by 16 pp).
6. **Anchor-density operating curve** (formalized): HR fold-error e(K) is monotone-decreasing and
   plateaus — e(3)≈1.40, e(5)≈1.15, e(≥12)≈1.08 — giving a concrete reporting standard: **post ≥5–6
   KM timepoints** (ideally ≥8) for a reliably reconstructable trial.

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

A head-to-head on the same real datasets (`HEADTOHEAD.md`) sharpens — and honestly tempers — the
"no-digitisation-error" advantage. Reconstructing the same true curve from exact registry anchors
versus a pixel-noised figure shows that at *equal* anchor density the digitisation noise itself costs
little (median |log HR| 0.148 exact vs 0.139 noised; exact wins 17/25), while a densely digitised
figure (25 points) recovers the HR *better* than few exact registry anchors (8 points) — median |log
HR| 0.041 vs 0.148. The registry path's binding weakness is therefore **anchor sparsity, not pixel
error**: it is competitive with and cleaner than digitisation *only when enough timepoints are posted*
(the ≥5–6 standard), and the two methods are complementary — digitise a published figure for dense
shape; reconstruct natively from the registry for exactness and provenance. This directly motivates the
reporting recommendation in `POLICY.md`: the native path's value is unlocked by anchor density.

## 6. Limitations

Registry coverage is the binding limit (hundreds of trials, not all). The censoring level is
under-identified; we surface this as honest interval width rather than a false point. Very small trials
(N≈137) do not reconstruct. The external-median check is sensitive to endpoint matching. True-IPD
validation used 31 open datasets (R `survival`/`KMsurv`/`asaur`); credentialed repositories (Vivli,
Project Data Sphere, YODA) would extend it further to dozens–hundreds of trials. Tier B is exponential-only.

## 7. Availability

Code, tests, harvester, validation scripts, offline tool and dashboard: MIT-licensed at
https://github.com/mahmood726-cyber/registry-ipd (live: https://mahmood726-cyber.github.io/registry-ipd/).

## References

*All PubMed-indexed references below were programmatically DOI→PMID resolved and field-matched
(author/year/journal/volume/pages) against PubMed on 2026-06-10 (verification log:
`CITATIONS.md`). arXiv preprints are flagged as such; the two non-indexed classics (Aalen–Johansen,
Rubin) are a journal article and a book respectively.*

1. Guyot P, Ades AE, Ouwens MJNM, Welton NJ. *BMC Med Res Methodol.* 2012;12:9. PMID 22297116.
   doi:10.1186/1471-2288-12-9 *(Verified. The abstract's own recommendation — "all RCTs should
   report information on numbers at risk and total number of events alongside KM curves" — is direct
   prior-authority support for our censoring-informed method and the ≥5–6-timepoint anchor-density
   finding.)*
2. Liu N, Zhou Y, Lee JJ. IPDfromKM. *BMC Med Res Methodol.* 2021;21(1):111. PMID 34074267.
   doi:10.1186/s12874-021-01308-8 *(Verified; abstract confirms figure-image input — "extract raw
   data coordinates from published K-M curves".)*
3. Royston P, Parmar MKB. Flexible parametric proportional-hazards and proportional-odds models for
   censored survival data. *Stat Med.* 2002;21(15):2175–2197. PMID 12210632. doi:10.1002/sim.1203
   *(Verified.)*
4. Aalen OO, Johansen S. An empirical transition matrix for non-homogeneous Markov chains based on
   censored observations. *Scand J Stat.* 1978;5(3):141–150. *(Foundational competing-risks paper;
   not PubMed-indexed.)*
5. Jansen JP. Network meta-analysis of survival data with fractional polynomials. *BMC Med Res
   Methodol.* 2011;11:61. PMID 21548941. doi:10.1186/1471-2288-11-61 *(Verified.)*
6. Rubin DB. *Multiple Imputation for Nonresponse in Surveys.* Wiley, 1987. *(Book; no DOI.)*
7. Yao JC, Fazio N, Singh S, et al. Everolimus for advanced, non-functional neuroendocrine tumours
   of the lung or gastrointestinal tract (RADIANT-4). *Lancet.* 2016;387(10022):968–977. PMID
   26703889. doi:10.1016/S0140-6736(15)00817-X *(Verified; print 2016, online 2015-12-17. Published
   PFS HR 0.48 [95% CI 0.35–0.67] matches our reconstruction.)*
8. RESOLVE-IPD. arXiv:2511.01785 (2025). · KM-GPT. arXiv:2509.18141 (2025). *(arXiv preprints —
   figure-image reconstruction pipelines; identifiers carried verbatim, not PubMed-indexed.)*
9. Titman AC. Using quadratic programming to reconstruct data from published survival and competing
   risks analyses. *Stat Med.* 2026;45(6-7):e70474. PMID 41775249. doi:10.1002/sim.70474 *(Verified;
   closest concurrent prior art — tabular numbers-at-risk input via QP, but from journal articles, not
   ClinicalTrials.gov/AACT. See `NOVELTY.md`.)*
10. Wei Y, Royston P. Reconstructing time-to-event data from published Kaplan-Meier curves (ipdfc).
    *Stata J.* 2017;17(4):786–802. PMID 29398980; PMCID PMC5796634. *(Figure-digitisation
    implementation.)*
11. Tasneem A, et al. The database for Aggregate Analysis of ClinicalTrials.gov (AACT). *PLoS One.*
    2012;7(3):e33677. PMCID PMC3306288. doi:10.1371/journal.pone.0033677 *(The AACT substrate; aggregate
    harvesting, no IPD reconstruction.)*
12. AACT / Clinical Trials Transformation Initiative, ClinicalTrials.gov. Validation datasets: R
    `survival` package (Therneau) and Rdatasets mirror.
