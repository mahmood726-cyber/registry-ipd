# Novelty assessment — registry-native survival IPD reconstruction

*Systematic prior-art search to test the paper's central novelty claim. Method: a fan-out multi-source
search (PubMed, arXiv, CRAN/PyPI, HTA grey literature, ClinicalTrials.gov-mining literature) with
3-vote adversarial verification of each claim; 103 search/verify agents, 23 confirmed claims. The
single pivotal new reference (Titman 2026) was additionally DOI→PMID resolved and abstract-matched
against PubMed by hand. Date: 2026-06-10.*

## The claim under test

> First method to reconstruct individual-patient-level survival data (pseudo-IPD) **natively from
> ClinicalTrials.gov / AACT structured summary data** — posted KM-estimate timepoints, number-at-risk
> tables, participant-flow event/censoring counts, and reported hazard ratios — **without digitizing a
> Kaplan–Meier figure image**.

## Verdict — DEFENSIBLE, with one required hedge (confidence: high)

No verified prior work reconstructs survival IPD from **registry-posted structured tables**. Every
established reconstruction method takes a **digitized figure image** as input. The distinguishing axis
is therefore data **provenance** (registry-posted survival tables vs figure-digitized coordinates),
**not** merely input file format — and the claim must be phrased on provenance to be airtight.

One concurrent 2026 method, **Titman (Stat Med)**, reconstructs pseudo-IPD from *published tabular*
numbers-at-risk via quadratic programming — the nearest neighbour on the "tabular input" axis — but
draws its tables from **journal articles, not ClinicalTrials.gov/AACT**. It does not falsify the
registry-native claim, but it must be cited explicitly as concurrent prior art to pre-empt a reviewer
challenge, and it also independently covers the competing-risks-from-CIF angle.

**Safe wording:** "first to reconstruct survival IPD natively from ClinicalTrials.gov/AACT *structured
registry data*, bypassing figure digitization" / "first registry-native (AACT-sourced) KM-IPD
reconstruction."
**Unsafe wording (do not use):** bare "first to reconstruct pseudo-IPD from tabular data" — Titman 2026
has tabular input.

## Closest prior art

| Work | Input modality | Reconstructs survival IPD? | Uses CT.gov/AACT structured tables? | Relation to our claim |
|---|---|---|---|---|
| **Guyot et al. 2012** (BMC MRM; PMID 22297116) | Digitized KM **figure** (DigitizeIt) + at-risk transcribed from figure | Yes (the canonical algorithm) | No | Foundational figure-digitization method we port and invert the *input* of. |
| **IPDfromKM** — Liu et al. 2021 (BMC MRM; PMID 34074267) | Mouse-click coordinates on a KM **bitmap**; accepts a tabular *coordinate* file, but coordinates are figure-read | Yes | No | Most-used implementation; "tabular" option is still figure-provenance. |
| **ipdfc** — Wei & Royston 2017 (Stata J; PMC5796634) | DigitizeIt curve coordinates off the **figure** | Yes | No | Stata sibling of Guyot; figure input. |
| **RESOLVE-IPD 2025** (arXiv:2511.01785) | Vector-graphics/raster KM **figures**; MAPLE uses published HR/median from papers | Yes | **No (states it does not use CT.gov/AACT)** | Recent automated pipeline; figure provenance. |
| **KM-GPT 2025** (arXiv:2509.18141) | KM **figure** images (LLM-assisted digitization) | Yes | No | Recent automated pipeline; figure provenance. |
| **Titman 2026** (Stat Med 45(6-7):e70474; PMID 41775249) | **Published tabular** numbers-at-risk + marked censoring times (and CIFs) via quadratic programming | Yes | **No (journal-article tables, not registry)** | **Closest neighbour.** Same "tabular input" idea, concurrent; different provenance. Cite explicitly. |
| **AACT / Tasneem 2012** (PLoS ONE; PMC3306288) | CT.gov structured results (aggregate) | **No** | n/a (is the registry mirror) | Provides the data substrate we mine; never attempts reconstruction. |
| **CT.gov results knowledge-graph mining** (PMC10771511) | CT.gov structured results | **No** (extracts p-values/CIs/AE counts) | Reads them, no survival reconstruction | Shows registry-results mining exists but is disjoint from KM/IPD reconstruction. |
| **NICE DSU TSD 19 (2017)** | Published KM **curves** (Guyot + Hoyle & Henley) | Yes (guidance) | No | Authoritative HTA guidance; only describes figure-based reconstruction. |
| **Coordinate-extraction tool reviews 2025** (CurveSnap, ScanIt, SurvdigitizeR) | KM **figure** images | Yes (digitizers) | No | Entire tool family is image-input. |

## Related Work paragraph (paste-ready for PAPER.md)

> Reconstruction of individual-patient survival data from published trials has, since Guyot et al.
> (2012), relied on digitizing the plotted Kaplan–Meier curve: an analyst reads (time, survival)
> coordinates off the figure image with software such as DigitizeIt, then inverts the KM equations
> using auxiliary numbers-at-risk transcribed from the risk table. This figure-digitization paradigm
> underlies the widely used IPDfromKM (Liu et al., 2021) and ipdfc (Wei & Royston, 2017)
> implementations, the recent automated pipelines RESOLVE-IPD (2025) and KM-GPT (2025), and is the
> only reconstruction route described in HTA methods guidance (NICE DSU TSD 19). Independently, the
> ClinicalTrials.gov ecosystem has been harvested at scale for aggregate evidence synthesis (the AACT
> database; CT.gov results knowledge graphs), but these efforts extract aggregate efficacy/safety
> fields and never reconstruct time-to-event IPD. The closest tabular-input neighbour is Titman
> (2026), who reconstructs pseudo-IPD — including competing-risks data from cumulative-incidence
> functions — from published numbers-at-risk and marked censoring times via quadratic programming, but
> draws those tables from journal articles rather than registry-posted structured results. To our
> knowledge, no prior method reconstructs survival IPD natively from the structured survival tables
> that ClinicalTrials.gov/AACT exposes (posted KM-estimate timepoints, number-at-risk tables,
> participant-flow event/censoring counts, and reported hazard ratios) without digitizing a figure.

## Sources

Verified against, among others: Guyot 2012 (PMC3313891), IPDfromKM (PMC8168323 + CRAN), Wei & Royston
(PMC5796634), RESOLVE-IPD (arXiv:2511.01785), Titman 2026 (doi:10.1002/sim.70474, PMID 41775249), AACT
(PMC3306288), CT.gov KG mining (PMC10771511), NICE DSU TSD 19 (Sheffield), parametric curve-fitting
(PMC11719069). Full verification log in the workflow transcript.
