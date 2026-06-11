# Synergy with `kmcurve` — the figure-extraction complement

Ideas from reviewing the sibling project **`kmcurve`** (`C:/Projects/kmcurve`), a PDF→figure→IPD
pipeline (PyMuPDF raster, k-medoids colour clustering, OCR axis calibration, **OCR numbers-at-risk**,
Guyot reconstruction). `kmcurve` is the *figure-digitisation* half; registry-IPD is the
*structured-registry* half. They are natural complements, and reviewing `kmcurve` surfaces four
concrete cross-pollination ideas — one of which directly attacks registry-IPD's binding limitation.

## The two projects, side by side

| | registry-IPD (this repo) | kmcurve |
|---|---|---|
| Input | AACT structured KM-estimate anchors | published PDF figure (rasterised) |
| Curve fidelity | **exact** (registry values) | pixel/OCR error |
| Numbers-at-risk | **absent** (0 rows in all of AACT) | **OCR'd from the figure's risk table** |
| Provenance | NCT record | publication |
| Coverage | ~600 trials posting a structured curve | any trial with a published KM figure |
| Reconstruction | Tiered + **Titman-2026 QP** (validated 1.05) | **Guyot** (`guyot.py`, needs NAR) |

The complementarity is exact: registry-IPD has the *curve* precisely but lacks *numbers-at-risk*;
kmcurve recovers *numbers-at-risk* (by OCR) but the *curve* is noisy.

## Idea 1 — the Titman QP is a drop-in upgrade for kmcurve's Guyot backend (proven here)

kmcurve reconstructs with Guyot (`ipd_km_pipeline/guyot.py`). Our 12-method benchmark
(`validate/method_zoo.js`, `VALIDATION.md`) shows the **Titman-2026 QP beats Guyot decisively on the
true-IPD gold standard (median HR fold-error ~1.05 vs ~1.14)**. Because kmcurve's OCR gives it a
numbers-at-risk table (hence the per-interval event count), the QP applies directly — and would extract
the HR more accurately than Guyot from the *same* digitised curve. Actionable: vendor the QP
(`reconstructArmQP`) as an alternative backend in kmcurve, or expose it as a shared library.

## Idea 2 — NAR fusion resolves registry-IPD's #1 limitation ★ — **VALIDATED (2026-06-10)**

**Result (`validate/nar_fusion.js`, 51-dataset gold standard):** fusing the registry-exact curve with a
figure's numbers-at-risk table (sparse, OCR-noise-modelled), via a NAR-aware QP reconstruction, recovers
the HR to **median fold 1.05 *without* a registry event count** — matching the Titman QP that *uses* the
registry event count, and lifting the heavily-censored TCGA cohorts to **20/20 within 20%**:

| | curve-only (no NAR/events) | **fusion: curve + figure NAR** | QP: curve + registry events |
|---|---|---|---|
| all 51 | 1.18 | **1.05** | 1.05 |
| ≥100/arm | 1.15 | **1.03** | 1.04 |
| heavily-censored TCGA | 1.38 | **1.06 (20/20)** | 1.05 (18/20) |

So **the figure's at-risk table substitutes for the missing registry event count** — dissolving the
binding identifiability limit this project documents. (Feeding NAR to the *anchor-exact* method only
gets to 1.14; the gain needs the QP's event-spreading + the NAR at-risk path together — see
`nar_fusion.js`.) This is the concrete best-of-both registry+figure reconstruction, validated on real
true-IPD. Below is the original motivation.



registry-IPD's binding accuracy limit is **the absent number-at-risk** (Section "identifiability limit"
in `VALIDATION.md`: censoring is unidentified from the curve alone, so curve-only underestimates large
HRs ~1.5-fold). kmcurve **OCRs the numbers-at-risk table** off the published figure. For a trial that is
*both* in AACT (exact curve) *and* published (figure with a risk table), fuse them:

> **exact registry curve anchors (registry-IPD) + OCR'd numbers-at-risk (kmcurve) → reconstruction with
> both the precise curve AND the censoring identified.**

This is strictly better than either alone: registry-IPD avoids the pixel error kmcurve has on the
curve, and kmcurve supplies the NAR registry-IPD lacks. It would lift the curve-only regime out of the
~1.5 identifiability trap *without* needing the registry to change its reporting — a concrete experiment
worth running on a handful of dual-available trials.

## Idea 3 — a *real* head-to-head replaces our simulated digitisation

`HEADTOHEAD.md` compares registry-structured vs **simulated** figure-digitisation (Gaussian pixel
noise). kmcurve is a **real** digitisation engine with a corpus of real trial PDFs
(`ipd_km_pipeline/corpus/`, `acquire_corpus.py`). For trials in both AACT and the corpus, run the
genuine comparison — registry-structured-QP vs kmcurve-figure-extraction, both scored against the
registry-reported HR (or true IPD where available). That upgrades the head-to-head from a noise model
to an empirical result.

## Idea 4 — cross-validate both engines on this gold standard

registry-IPD has a 51-dataset **true-IPD gold standard** (R packages + 14 TCGA + 6 non-TCGA cBioPortal).
kmcurve has competitor benchmarks (`benchmark_competitors.py`, `EXCEEDING_SURVDIGITIZER.md`) but
validates on digitised curves. Running kmcurve's Guyot backend and registry-IPD's QP on the *same*
true-IPD coarse summaries (already built here) gives a clean, shared, apples-to-apples engine
comparison — and the gold-standard harness (`validate/goldstandard.js`) is the natural host.

## Idea 5 — the abstract event-count lever: NAR fusion's gain, INSIDE the data-scope contract ★★ — **BUILT (2026-06-11)**

Idea 2 (NAR fusion) is scientifically the strongest, but it imports a *figure's* at-risk table — and a
published figure is **out of this project's production data scope** (AACT + PubMed abstracts only; figures
are validation-only). The key realisation: **the censoring lever the QP needs is a per-arm total-event
count, and the PubMed abstract — which IS in scope — routinely prints exactly that**:

> "death occurred in **107 of 205** patients in the everolimus group versus **77 of 97** ..."
> "the mortality rate was 20.2% (**19 of 94**) for sabizabulin versus 45.1% (**23 of 51**) ..."

So the production-legal analogue of kmcurve's number-at-risk OCR is a deterministic **abstract
event-count extractor** — same lever, in-scope source, no figure. Built and validated:

- **`harvest/abstract_events.py`** — bounded-regex per-arm "X of N" event-fraction extractor +
  `match_to_arms` (aligns a count to the right arm by registry N, label-independent) +
  `enrich_trial_events` (fills a trial's missing `total_events`, stamps `events_source='pubmed_abstract'`,
  **never overwrites** an AACT participant-flow count). Mirrors the existing `abstract_hr.py` /
  `abstract_median.py` design.
- **High-precision by construction** (a wrong count corrupts the QP, so over-inclusion is the enemy):
  clause-scoped guards reject the dominant false-positive classes — **adverse-event / safety counts**
  ("serious adverse events ... 22 of 218"), **enrolment / response** fractions, **negation** ("no
  deaths"), **count > N**, and **drug-name slashes** ("CDK4/6 inhibitor" → not 4/6). 19 unit tests.
- **Validated on the real abstract cache** (`.pubmed_cache.json`, 161 abstracts): after the guards,
  **1 hit, a true positive** (mortality 19/94 vs 23/51) and **0 false positives** — 100% precision.
  Recall is honestly low: abstracts report per-arm "X of N" event counts far less often than they report
  a median or an HR, so this is **one lever among several** (HR→`calibrateHR`, median→cross-check), not a
  universal fix. When an abstract does post the fractions, the QP gets its censoring lever with no figure.

This closes the constraint-respecting version of the identifiability fix: registry-IPD's binding limit is
the missing event count; the **abstract supplies it within the AACT+PubMed-abstract contract**, and the
figure-NAR fusion (Idea 2) remains the *validation-only* upper bound on what an out-of-scope figure could
add.

**Grounded end-to-end on real PubMed data (2026-06-11).** Both halves of the chain are now validated on
real trials (different trials, because no local trial publishes all three pieces openly — the same
scarcity the figure path found): (1) **abstract → count** on DAPA-HF (PMID 31535829, NEJM 2019) —
`abstract_events` pulls `386 of 2373` vs `502 of 2371` straight from "the primary outcome occurred in …",
matched to arms by N; (2) **count → QP → truth** on RADIANT-4 (NCT01524783) — the 107/77 event count fed
to the QP moves the HR from **0.679 (outside the posted 95% CI) to 0.577 (inside it)** vs truth 0.48.
RADIANT-4's own abstract posts only adverse-event counts (correctly rejected → None) + the HR 0.48, so
its in-scope event source is AACT participant-flow. See `VALIDATION.md` "Abstract event-count lever";
locked by `test/abstract_lever.spec.js`.

## Recommended next step

**Idea 2 (NAR fusion)** is the scientifically strongest: it directly dissolves the identifiability
limit this project documents, using a data source (the figure's risk table) that kmcurve already
extracts. Prototype: take 5–10 trials present in both AACT and a published PDF, OCR their NAR via
kmcurve, feed `{registry curve anchors, kmcurve NAR}` into the QP, and measure the HR fold-error against
truth versus registry-curve-only. If it lands near the censoring-informed 1.05 *without* a registry
event count, that is a publishable "best-of-both" registry+figure reconstruction.
