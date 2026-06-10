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

## Idea 2 — NAR fusion resolves registry-IPD's #1 limitation ★ (highest value)

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

## Recommended next step

**Idea 2 (NAR fusion)** is the scientifically strongest: it directly dissolves the identifiability
limit this project documents, using a data source (the figure's risk table) that kmcurve already
extracts. Prototype: take 5–10 trials present in both AACT and a published PDF, OCR their NAR via
kmcurve, feed `{registry curve anchors, kmcurve NAR}` into the QP, and measure the HR fold-error against
truth versus registry-curve-only. If it lands near the censoring-informed 1.05 *without* a registry
event count, that is a publishable "best-of-both" registry+figure reconstruction.
