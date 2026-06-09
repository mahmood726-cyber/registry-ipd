#!/usr/bin/env python
"""Phase 0 feasibility/coverage scan: how much of AACT is actually reconstructable?

Classifies the universe of trials with posted survival results into:
  Tier A  rich     : KM-estimate points (>=3) + number-at-risk (>=2) + total events
  Tier B  medium   : median + hazard ratio + N + events (parametric)
  Tier C  sparse   : hazard ratio only (or less) -> not reconstructable

This SIZES THE WHOLE PREMISE. If Tier A is vanishingly small, the headline shifts to
Tier B and the scoped claim narrows accordingly (documented, not hidden).

Requires a resolvable AACT snapshot. Fails closed otherwise.
Usage: python harvest/coverage_scan.py [--limit N] [-o coverage_report.json]
"""
from __future__ import annotations
import argparse
import json
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(__file__))
import harvester as H  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None, help="cap trials scanned (sampling)")
    ap.add_argument("-o", "--out", default="coverage_report.json")
    args = ap.parse_args(argv)

    try:
        from aact_kit import load_table, resolve_aact_location
        loc = resolve_aact_location()
    except Exception as e:  # noqa: BLE001
        print("ERROR: no AACT snapshot resolvable. Set AACT_ZIP/AACT_SQLITE/AACT_TSV_DIR.\n",
              file=sys.stderr)
        print(str(e), file=sys.stderr)
        return 2

    print(f"AACT location: {loc}")
    # 1) trials with a hazard-ratio analysis
    print("loading outcome_analyses ...")
    analyses = load_table("outcome_analyses", location=loc,
                          columns=["nct_id", "outcome_id", "param_type", "param_value", "method"])
    hr_mask = analyses["param_type"].fillna("").str.contains(H._HR_PARAM) | \
        analyses.get("method", "").fillna("").str.contains("cox", case=False)
    hr_trials = set(analyses[hr_mask]["nct_id"])
    print(f"  trials with HR/Cox analysis: {len(hr_trials)}")

    # 2) measurement-level signals. NOTE (verified on the 2026-06 snapshot): AACT has ZERO
    #    structured "number at risk" rows, so Tier A must NOT require NAR. A KM curve is a
    #    survival-ish measure reported at >=3 classifications that PARSE AS TIMEPOINTS
    #    (excludes category classifications like "mortality"/"progression"). Reconstruction then
    #    uses N + total events with censoring at the tail (IPDfromKM "no number-at-risk" mode).
    print("loading outcome_measurements (large) ...")
    meas = load_table("outcome_measurements", location=loc,
                      columns=["nct_id", "outcome_id", "title", "param_type", "classification"])
    blob = (meas["param_type"].fillna("") + " " + meas["title"].fillna("")).str.lower()
    surv = meas[blob.str.contains("kaplan|survival|progression-free|event-free", regex=True)].copy()
    med_trials = set(meas[blob.str.contains("median")]["nct_id"])
    nar_available = int(blob.str.contains("number at risk", regex=False).sum())

    # KM curve = survival-ish measure with >=3 distinct PARSEABLE timepoints
    km_pts = defaultdict(set)
    for nct, oid, cls, title in zip(surv["nct_id"], surv["outcome_id"],
                                    surv["classification"].fillna(""), surv["title"].fillna("")):
        t = H.parse_timepoint_to_months(cls) or H.parse_timepoint_to_months(title)
        if t is not None:
            km_pts[(nct, oid)].add(round(t, 3))

    tierA = set(k[0] for k, pts in km_pts.items() if len(pts) >= 3)
    tierB = (med_trials & hr_trials) - tierA
    universe = set(meas["nct_id"]) | hr_trials
    if args.limit:
        universe = set(list(universe)[:args.limit])
    tierC = universe - tierA - tierB

    report = {
        "aact_location": str(loc),
        "universe_trials_with_results": len(universe),
        "tier_A_km_curve": len(tierA),
        "tier_B_medium": len(tierB),
        "tier_C_sparse": len(tierC),
        "trials_with_hr": len(hr_trials),
        "trials_with_median": len(med_trials),
        "structured_number_at_risk_rows": nar_available,
        "note": "Tier A = survival measure at >=3 PARSEABLE timepoints (NAR not required — AACT has "
                "~0 structured NAR; reconstruction uses N+events, censoring to tail). B = median+HR. "
                "C = rest. This is the corrected, reality-matched definition.",
        "sample_tierA_nct": list(tierA)[:30],
    }
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    print(json.dumps({k: v for k, v in report.items() if k != "sample_tierA_nct"}, indent=2))
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
