#!/usr/bin/env python
"""Make the PubMed abstract a first-class ENRICHMENT source for a harvested trial (in scope: AACT +
PubMed abstracts only).

Until now the abstract extractors had two separate roles:
  - `abstract_events.py`  ENRICHES the reconstruction input (fills per-arm `total_events` -> QP lever);
  - `abstract_hr.py` / `abstract_median.py`  only VALIDATE the reconstruction after the fact.

This unifies them into one entry point that applies every in-scope abstract lever to a trial dict, each
with explicit provenance and conservative precision gating, and NEVER overwrites an AACT-posted value:

  1. per-arm event counts  -> `arm.total_events`           (the QP censoring lever)
  2. a confident primary HR -> `trial.hr` if AACT posted none, else `trial.hr_abstract` (cross-check)
  3. published median(s)    -> `trial.median_abstract`     (cross-check only; never drives the curve)

The HR is gated hard: only a CI-backed HR that is either the abstract's *only* HR or endpoint-matched is
trusted as the trial HR — a multi-HR abstract (primary + secondary + subgroup) is attached as a flagged
cross-check, never promoted, because a wrong HR would mis-calibrate the reconstruction.

Pure (no network): the caller fetches the abstract. Returns a provenance summary.
"""
from __future__ import annotations
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import abstract_events as _ev      # noqa: E402
import abstract_hr as _hr          # noqa: E402
import abstract_median as _med     # noqa: E402


def _has_usable_hr(trial: dict) -> bool:
    hr = trial.get("hr")
    return bool(hr and hr.get("value") is not None)


def _hr_is_confident(h: dict) -> bool:
    """A CI-backed HR that is unambiguous (the only HR) or matched to the requested endpoint."""
    return bool(h and h.get("has_ci") and
                (h.get("n_hr_candidates", 99) <= 1 or h.get("endpoint_matched") is True))


def enrich_from_abstract(trial: dict, abstract: str, endpoint: str = None) -> dict:
    """Apply every in-scope abstract lever to `trial` in place. Returns a provenance summary."""
    summary = {"events": {"patched": 0}, "hr": {"set_as": None, "confident": False},
               "median": {"attached": False}}
    if not abstract:
        return summary

    # 1. event counts -> total_events (delegates to the precision-guarded extractor)
    ev = _ev.enrich_trial_events(trial, abstract, endpoint=endpoint)
    summary["events"] = {"patched": ev["patched"],
                         "source": "pubmed_abstract" if ev["patched"] else None}

    # 2. HR -> trial.hr (only if AACT posted none AND the abstract HR is confident), else cross-check
    h = _hr.extract_hr(abstract, endpoint=endpoint)
    if h:
        confident = _hr_is_confident(h)
        rec = {"value": h["value"], "ci_low": h["ci_low"], "ci_high": h["ci_high"],
               "method": "PubMed abstract", "source": "pubmed_abstract",
               "endpoint": endpoint, "n_hr_candidates": h["n_hr_candidates"]}
        if confident and not _has_usable_hr(trial):
            trial["hr"] = rec
            summary["hr"] = {"set_as": "trial.hr", "confident": True, "value": h["value"]}
        else:
            trial["hr_abstract"] = rec                 # independent cross-check, never overrides AACT
            summary["hr"] = {"set_as": "trial.hr_abstract", "confident": confident, "value": h["value"]}

    # 3. median(s) -> cross-check only (never feeds the reconstruction)
    m = _med.extract_medians(abstract, endpoint=endpoint)
    if m and (m.get("medians") or m.get("not_reached")):
        trial["median_abstract"] = {"medians": m["medians"], "not_reached": m["not_reached"],
                                    "source": "pubmed_abstract", "endpoint": endpoint}
        summary["median"] = {"attached": True, "medians": m["medians"]}

    return summary


def enrich_trial_with_fetcher(trial: dict, pmid, fetch_abstract, endpoint: str = None) -> dict:
    """Resolve the abstract via the injected `fetch_abstract(pmid) -> str|None` and enrich the trial.

    The PMID->abstract source is injected so the production path (cached PubMed efetch) and tests (a stub)
    share one code path. Fail-soft: a missing PMID or a fetch that returns nothing is a no-op (enrichment
    is additive — it must never block or corrupt a harvest). Returns the enrichment summary + the pmid.
    """
    base = {"pmid": pmid, "events": {"patched": 0}, "hr": {"set_as": None, "confident": False},
            "median": {"attached": False}}
    if not pmid:
        return base
    try:
        abstract = fetch_abstract(pmid) or ""
    except Exception as e:  # noqa: BLE001  -- a fetch failure must not break the harvest
        base["error"] = f"{type(e).__name__}: {e}"
        return base
    summary = enrich_from_abstract(trial, abstract, endpoint=endpoint)
    summary["pmid"] = pmid
    return summary


if __name__ == "__main__":
    import json
    trial = json.load(open(sys.argv[1], encoding="utf-8"))
    abstract = open(sys.argv[2], encoding="utf-8").read() if len(sys.argv) > 2 else ""
    print(json.dumps(enrich_from_abstract(trial, abstract), indent=2))
