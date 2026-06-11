#!/usr/bin/env python
"""Apply the unified PubMed-abstract enrichment across the whole harvested cohort (in scope: AACT +
PubMed abstracts).

This is the "use it for real" step after `abstract_events_coverage.py` measured the yield: for every
cohort trial that has an AACT-resolved results-publication PMID (from `cohort_pubmed_validation.json`,
which read AACT `study_references`) and a cached abstract, run `abstract_enrich.enrich_from_abstract` to
fill in `total_events` / a confident HR / a median — each with provenance, never overwriting an AACT value.

Writes the enriched trials to `cohort_enriched/<nct>.json` (gitignored data) and a committable manifest
`realipd/cohort_enrichment_manifest.json` summarising exactly what each trial gained. Pure-local: uses the
already-cached abstracts (no network) and the existing harvested cohort JSONs (no AACT re-harvest).

Run from repo root:  python harvest/enrich_cohort.py
"""
from __future__ import annotations
import io
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import abstract_enrich as AE  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(__file__))
RIPD = os.path.join(ROOT, "realipd")
COHORT = os.path.join(ROOT, "cohort")
OUT = os.path.join(ROOT, "cohort_enriched")


def _trial_has_hr(t):
    hr = t.get("hr")
    return bool(hr and hr.get("value") is not None)


def enrich_rows(rows, abstracts, trials_by_nct):
    """Pure: enrich each trial in place from its abstract; return (manifest, {nct: enriched_trial}).
    `rows`: [{nct, pmid, curve_endpoint}]; `abstracts`: {pmid: text}; `trials_by_nct`: {nct: trial}."""
    manifest = {"considered": 0, "enriched": 0, "events_filled": 0, "hr_promoted": 0,
                "hr_crosscheck": 0, "median_attached": 0, "trials": []}
    enriched = {}
    for r in rows:
        nct, pmid = r.get("nct"), r.get("pmid")
        ab = abstracts.get(pmid) if pmid else None
        trial = trials_by_nct.get(nct)
        if not ab or trial is None:
            continue
        manifest["considered"] += 1
        had_hr = _trial_has_hr(trial)
        summary = AE.enrich_from_abstract(trial, ab, endpoint=r.get("curve_endpoint") or None)
        gained = []
        if summary["events"]["patched"]:
            gained.append(f"events(+{summary['events']['patched']})")
            manifest["events_filled"] += 1
        if summary["hr"]["set_as"] == "trial.hr":
            gained.append(f"hr={summary['hr']['value']}")
            manifest["hr_promoted"] += 1
        elif summary["hr"]["set_as"] == "trial.hr_abstract":
            gained.append(f"hr_xcheck={summary['hr']['value']}")   # independent HR for triangulation
            manifest["hr_crosscheck"] += 1
        if summary["median"]["attached"]:
            gained.append("median")
            manifest["median_attached"] += 1
        if gained:
            manifest["enriched"] += 1
            enriched[nct] = trial
            manifest["trials"].append({"nct": nct, "pmid": pmid, "had_aact_hr": had_hr, "gained": gained})
    return manifest, enriched


def run():
    rows = json.load(open(os.path.join(RIPD, "cohort_pubmed_validation.json"), encoding="utf-8")).get("rows", [])
    abstracts = json.load(open(os.path.join(ROOT, ".pubmed_cache.json"), encoding="utf-8"))
    trials = {}
    for r in rows:
        p = os.path.join(COHORT, f"{r.get('nct')}.json")
        if os.path.exists(p):
            try:
                trials[r["nct"]] = json.load(open(p, encoding="utf-8"))
            except Exception:
                pass
    manifest, enriched = enrich_rows(rows, abstracts, trials)
    os.makedirs(OUT, exist_ok=True)
    for nct, trial in enriched.items():
        json.dump(trial, open(os.path.join(OUT, f"{nct}.json"), "w"), indent=2)
    json.dump(manifest, open(os.path.join(RIPD, "cohort_enrichment_manifest.json"), "w"), indent=2)
    return manifest


def main():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    m = run()
    print("=== cohort enrichment from PubMed abstracts (in scope) ===")
    print(f"  considered (PMID + cached abstract + cohort JSON): {m['considered']}")
    print(f"  trials enriched (gained >=1 in-scope field)      : {m['enriched']}")
    print(f"    HR promoted to trial.hr (AACT had none)        : {m['hr_promoted']}")
    print(f"    HR attached as cross-check (AACT/ambiguous)    : {m['hr_crosscheck']}")
    print(f"    per-arm event counts filled                    : {m['events_filled']}")
    print(f"    medians attached                               : {m['median_attached']}")
    promoted = [t for t in m["trials"] if any(g.startswith("hr=") for g in t["gained"])]
    if promoted:
        print(f"\n  trials that gained a usable HR they otherwise lacked ({len(promoted)}):")
        for t in promoted:
            hrv = next(g for g in t["gained"] if g.startswith("hr="))
            print(f"    {t['nct']} (PMID {t['pmid']}): {hrv}")
    print(f"\n  wrote cohort_enriched/*.json + realipd/cohort_enrichment_manifest.json")
    return m


if __name__ == "__main__":
    main()
