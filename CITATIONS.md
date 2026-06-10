# Citation verification log

Programmatic DOI→PMID resolution and field-match of every PubMed-indexable reference in `PAPER.md`,
per the project's citation-integrity discipline (LLM-drafted reference lists carry a ~4% baseline
misattribution rate, so eyeball inspection is not sufficient — every citation is machine-resolved).

**Method.** Each DOI was resolved to a PMID via NCBI ID Converter (PubMed MCP `convert_article_ids`),
then PubMed `get_article_metadata` / `lookup_article_by_citation` was used to confirm the returned
title, first author, year, journal, volume and pages match the citation text.

**Date.** 2026-06-10. **Source.** PubMed (NCBI E-utilities).

| # | Citation (PAPER.md) | DOI | PMID | Match | Notes |
|---|---------------------|-----|------|-------|-------|
| 1 | Guyot et al. 2012, *BMC Med Res Methodol* 12:9 | 10.1186/1471-2288-12-9 | 22297116 | ✅ exact | Title "Enhanced secondary analysis of survival data: reconstructing the data from published Kaplan-Meier survival curves." Abstract explicitly recommends RCTs report numbers-at-risk + total events alongside KM curves — supports our anchor-density finding. |
| 2 | Liu, Zhou, Lee 2021 (IPDfromKM), *BMC Med Res Methodol* 21(1):111 | 10.1186/s12874-021-01308-8 | 34074267 | ✅ exact | Title confirms figure-image input ("reconstruct individual patient data from published Kaplan-Meier survival curves"; "extract raw data coordinates"). |
| 3 | Royston & Parmar 2002, *Stat Med* 21(15):2175–2197 | 10.1002/sim.1203 | 12210632 | ✅ exact | Pages 2175-97. Flexible parametric PH/PO spline models. Did not resolve through idconv DOI lookup (Wiley) but matched cleanly via citation lookup. |
| 4 | Aalen & Johansen 1978, *Scand J Stat* 5(3):141–150 | — | — | n/a | Not PubMed-indexed (1978 statistics journal). Foundational competing-risks transition-matrix paper; cited from the primary source. |
| 5 | Jansen 2011, *BMC Med Res Methodol* 11:61 | 10.1186/1471-2288-11-61 | 21548941 | ✅ exact | "Network meta-analysis of survival data with fractional polynomials." Supports the time-varying-HR / fractional-polynomial method. |
| 6 | Rubin 1987, *Multiple Imputation for Nonresponse in Surveys*, Wiley | — | — | n/a | Book; no DOI. Underpins the multiple-imputation uncertainty method. |
| 7 | Yao et al. 2016 (RADIANT-4), *Lancet* 387(10022):968–977 | 10.1016/S0140-6736(15)00817-X | 26703889 | ✅ exact | Print 2016, online 2015-12-17. Published PFS HR 0.48 (95% CI 0.35–0.67); matches our external-publication reconstruction check (0.47–0.48). |
| 8 | RESOLVE-IPD arXiv:2511.01785; KM-GPT arXiv:2509.18141 | — | — | arXiv | 2025 preprints, figure-image reconstruction pipelines. Not PubMed-indexed; arXiv identifiers carried verbatim. (Cross-checked by the novelty search; see `NOVELTY.md`.) |
| 9 | AACT / CTTI, ClinicalTrials.gov; R `survival` (Therneau) + Rdatasets | — | — | data | Data sources, not bibliographic claims. |

**Result.** 5/5 PubMed-indexable references verified exact, 0 misattributions, 0 fabrications. The 4
remaining entries are a non-indexed 1978 journal article, a book, two arXiv preprints, and the data
sources — all correctly categorised. The two arXiv preprints are additionally scrutinised by the
prior-art / novelty search (`NOVELTY.md`).
