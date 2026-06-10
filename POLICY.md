# Two small reporting changes would make trial survival results reusable

**A reporting recommendation for ClinicalTrials.gov / AACT (CTTI) and EU CTIS**

*From the registry-IPD reconstruction study (`PAPER.md`, `VALIDATION.md`). Evidence base: the
2026-06-01 AACT snapshot (76,067 trials with posted results) plus validation against true
patient-level data from 25 open RCT/cohort datasets. Reproduce with `harvest/census_full_aact.py`
and `validate/census_cohort.js`.*

---

## The problem, in two numbers

1. **Zero.** Across all **76,067** ClinicalTrials.gov trials with posted results, the AACT relational
   mirror contains **0 structured "number at risk" rows.** Number-at-risk — the single field that
   most improves the reconstruction of a survival (Kaplan–Meier) curve into reusable patient-level
   data — is *never* posted as structured data, even when it appears in the figure of the matching
   publication.

2. **Too few timepoints.** Only **288–~600 trials** (0.4%–0.8% of the 76,067) post a structured
   survival curve at all — and of those, the median trial posts just **3–4 KM-estimate timepoints.**
   Our validation against true patient-level data shows reconstruction is weak at 3–4 timepoints and
   only becomes reliable at **≥5–6** (hazard-ratio fold-error falls from ≈1.40 at 3 points to ≈1.15
   at 5 and ≈1.08 by 12). Only about **one third (≈34%)** of curve-posting trials currently clear
   that bar.

### Evidence (2026-06-01 AACT snapshot)

| Quantity | Count | % of 76,067 |
|---|---:|---:|
| Trials with posted results | 76,067 | 100% |
| Structured *number-at-risk* rows | **0** | 0% |
| Posts a reconstructable KM curve (≥3 timepoints) | 288–~600¹ | 0.4–0.8% |
| — of those, ≥6 timepoints (reliably reconstructable) | ≈174 | ≈34% of curve-posters |
| — of those, only 3–4 timepoints (weak) | ≈281 | ≈55% of curve-posters |
| Reconstructable at *some* tier (curve, or median+HR) | ≈3,800 | ≈5.0% |

¹ 288 when a "survival curve" is detected strictly (title says Kaplan-Meier / survival / progression-
or event-free); ≈514–605 under a broader net (adds disease-free survival, "probability of event",
cumulative incidence). Either way it is well under 1% of posted-results trials. Reproduce:
`harvest/census_full_aact.py` (full snapshot) and `validate/census_cohort.js` (per-arm detail).

The result: a survival result that an agency *required* the sponsor to post is, in most cases, not
accurately reusable for the meta-analysis and cost-effectiveness work it was meant to inform — not
because the data are secret, but because a few structured fields are missing. This is a free,
already-collected evidence asset being lost at the last step.

This is not a new ask. Guyot et al. (2012), whose algorithm underpins essentially all KM
reconstruction, concluded their validation paper with exactly this recommendation: *"all RCTs should
report information on numbers at risk and total number of events alongside KM curves"*
([DOI](https://doi.org/10.1186/1471-2288-12-9), PMID 22297116). Fourteen years later, the structured
registries still do not capture it.

## The recommendation (two fields, one threshold)

For any trial posting a time-to-event outcome, the results schema should capture, **as structured
data** (not only inside an uploaded figure):

1. **Number at risk per arm at each posted timepoint.** This is the highest-value single change. It
   removes the principal source of under-identification in reconstruction (the censoring pattern) and
   is information the sponsor already has and already plots.

2. **Total events per arm** (and, where applicable, competing-event counts via the
   existing participant-flow / withdrawal-reason structure).

3. **A minimum anchor density of ≥5–6 KM-estimate timepoints** per curve (ideally ≥8), reasonably
   spaced across follow-up. This is the empirically-derived threshold above which reconstruction
   becomes reliable; below it, a posted curve is decorative rather than reusable.

None of this asks for new data collection — only that fields already in the figure be entered in the
structured results form. For **AACT/CTTI**, items 1–2 map onto the existing `outcome_measurements`
and participant-flow tables (a new `number_at_risk` classification, and an events count). For **EU
CTIS**, the same three items belong in the structured summary-results template.

## Why it is worth doing

- **Reuse without re-collection.** Hundreds of already-posted survival results become accurately
  reconstructable for evidence synthesis and HTA, with full provenance to the NCT/EU record and
  *zero* figure-digitisation error.
- **Cheap and bounded.** Two fields and a timepoint-count guideline. The data already exist in the
  sponsor's own figure.
- **Auditable.** A registry can measure its own compliance directly: the same census scripts that
  produced the numbers above (anchor-density distribution, NAR-coverage rate) can run as a quarterly
  data-quality metric.

## What we are *not* claiming

We are not claiming reconstruction replaces individual-patient-data sharing (Vivli, YODA, Project Data
Sphere) — true IPD remains the gold standard. We are claiming that a small, free improvement to
*structured results reporting* would make a large, currently-wasted share of registry survival data
reusable, and that the ≥5–6-timepoint threshold gives registries a concrete, measurable target.

---

*Contact / reproducibility: all figures regenerate from the cited scripts against a public AACT
snapshot; method and validation in `PAPER.md`, `METHODS.md`, `VALIDATION.md`.*
