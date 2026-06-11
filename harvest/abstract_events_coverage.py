#!/usr/bin/env python
"""Measure the REAL marginal coverage of the abstract event-count lever on the harvested cohort.

Question the lever has to answer honestly: across real reconstructable trials, how often does the PubMed
abstract actually supply a per-arm total-event count that AACT did NOT post? That marginal gain is the
lever's whole value (the QP needs the count; AACT participant-flow frequently lacks it).

Pure analysis, fully local (no network, no AACT mount): pairs each trial↔PMID row with its cached abstract
and the harvested cohort JSON, runs `abstract_events`, and tallies where the abstract fills a gap. Run
from repo root:  python harvest/abstract_events_coverage.py  ->  realipd/abstract_events_coverage.json
"""
from __future__ import annotations
import io
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import abstract_events as E   # noqa: E402
import abstract_hr as H       # noqa: E402
import abstract_median as M   # noqa: E402

ROOT = os.path.dirname(os.path.dirname(__file__))
RIPD = os.path.join(ROOT, "realipd")
COHORT = os.path.join(ROOT, "cohort")


def _arm_ns(trial):
    seen, ns = set(), []
    for a in trial.get("arms", []):
        if a.get("label") in seen or not a.get("N"):
            continue
        seen.add(a.get("label"))
        ns.append(int(a["N"]))
    return ns


def _aact_has_events(trial):
    return any(a.get("total_events") is not None for a in trial.get("arms", []))


def analyze(rows, abstracts, trials_by_nct):
    """rows: [{nct, pmid, curve_endpoint}]; abstracts: {pmid: text}; trials_by_nct: {nct: trial dict}.
    Returns a coverage summary + per-trial gain records (where the abstract fills a missing event count)."""
    stat = {"rows": len(rows), "with_abstract": 0,
            # the three in-scope abstract levers, measured side by side
            "events_available": 0, "hr_available": 0, "median_available": 0, "any_lever": 0,
            # HR marginal value: the registry posts no HR but the abstract supplies one
            "registry_hr_present": 0, "hr_marginal_gain": 0,
            # event-count specifics (the QP's censoring lever)
            "two_arm_fractions": 0, "aact_already_has_events": 0,
            "marginal_gain": 0, "n_matched_both_arms": 0}
    gains = []
    for r in rows:
        nct, pmid = r.get("nct"), r.get("pmid")
        ab = abstracts.get(pmid) if pmid else None
        trial = trials_by_nct.get(nct)
        if not ab or trial is None:
            continue
        stat["with_abstract"] += 1
        ep = r.get("curve_endpoint") or None
        ev = E.extract_events(ab, endpoint=ep)
        hr = H.extract_hr(ab, endpoint=ep)
        md = M.extract_medians(ab, endpoint=ep)
        reg_hr = r.get("registry_HR")
        if reg_hr is not None:
            stat["registry_hr_present"] += 1
        if hr:
            stat["hr_available"] += 1
            if reg_hr is None:                         # abstract supplies an HR the registry lacks
                stat["hr_marginal_gain"] += 1
        if md and (md.get("medians") or md.get("not_reached")):
            stat["median_available"] += 1
        if ev or hr or md:
            stat["any_lever"] += 1
        if not ev:
            continue
        stat["events_available"] += 1
        if ev["n_fractions"] >= 2:
            stat["two_arm_fractions"] += 1
        aact = _aact_has_events(trial)
        if aact:
            stat["aact_already_has_events"] += 1
        ns = _arm_ns(trial)
        mapping = E.match_to_arms(ev, ns)
        if len(mapping) >= 2:
            stat["n_matched_both_arms"] += 1
        if not aact and mapping:                       # abstract fills a gap AACT left
            stat["marginal_gain"] += 1
            gains.append({"nct": nct, "pmid": pmid, "events": ev["events"], "ns": ev["ns"],
                          "arm_ns": ns, "mapped_arms": len(mapping), "context": ev["context"][:120]})
    return {"summary": stat, "gains": gains}


def main():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    val = json.load(open(os.path.join(RIPD, "cohort_pubmed_validation.json"), encoding="utf-8"))
    rows = val.get("rows", [])
    abstracts = json.load(open(os.path.join(ROOT, ".pubmed_cache.json"), encoding="utf-8"))
    trials = {}
    for r in rows:
        nct = r.get("nct")
        p = os.path.join(COHORT, f"{nct}.json")
        if os.path.exists(p):
            try:
                trials[nct] = json.load(open(p, encoding="utf-8"))
            except Exception:
                pass
    out = analyze(rows, abstracts, trials)
    os.makedirs(RIPD, exist_ok=True)
    json.dump(out, open(os.path.join(RIPD, "abstract_events_coverage.json"), "w"), indent=2)
    s = out["summary"]
    n = max(s["with_abstract"], 1)
    pct = lambda x: f"{100 * x / n:.0f}%"
    print("=== abstract enrichment levers - real cohort coverage ===")
    print(f"  trial<->PMID rows                : {s['rows']}")
    print(f"  with a cached abstract + cohort  : {s['with_abstract']}")
    print(f"  -- in-scope levers (share of {s['with_abstract']} abstracts) --")
    print(f"  HR available                     : {s['hr_available']}  ({pct(s['hr_available'])})")
    print(f"    registry posts an HR            : {s['registry_hr_present']}")
    print(f"    HR MARGINAL GAIN (registry lacks one, abstract supplies): {s['hr_marginal_gain']}  ({pct(s['hr_marginal_gain'])})")
    print(f"  median available                 : {s['median_available']}  ({pct(s['median_available'])})")
    print(f"  per-arm event count available    : {s['events_available']}  ({pct(s['events_available'])})")
    print(f"  any lever                        : {s['any_lever']}  ({pct(s['any_lever'])})")
    print(f"  -- event-count specifics (the QP censoring lever) --")
    print(f"    two-arm fractions              : {s['two_arm_fractions']}")
    print(f"    N-matched to both arms         : {s['n_matched_both_arms']}")
    print(f"    AACT already posted events     : {s['aact_already_has_events']}")
    print(f"    MARGINAL GAIN (fills AACT gap) : {s['marginal_gain']}")
    for g in out["gains"]:
        print(f"      + {g['nct']} (PMID {g['pmid']}): events {g['events']} of {g['ns']} -> arms {g['arm_ns']}")
    print(f"\n  wrote realipd/abstract_events_coverage.json")
    return out


if __name__ == "__main__":
    main()
