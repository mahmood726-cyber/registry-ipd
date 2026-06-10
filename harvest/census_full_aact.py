#!/usr/bin/env python
"""Definitive full-AACT reconstructability census + anchor-density distribution.

Reconciles the strict coverage_scan count (288) with the harvested-cohort count (605) by reporting
BOTH detection criteria over the ENTIRE 76k-trial results universe, and adds the anchor-density
(KM-timepoint) histogram across the whole Tier-A population — the evidence base for the policy brief
and the paper's binding-coverage limitation.

  STRICT survival curve : title/param matches kaplan|survival|progression-free|event-free
  BROAD  survival curve : harvester._SURV_RE (adds disease-free, probability-of, cumulative-incidence,
                          proportion) — the same net the cohort harvest used.
Tier A (either) = a survival-ish (nct, outcome) posted at >=3 distinct PARSEABLE timepoints.

Requires a resolvable AACT snapshot (AACT_ZIP). Fails closed otherwise.
Usage: python harvest/census_full_aact.py -o realipd/census_full_aact.json
"""
from __future__ import annotations
import argparse
import json
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(__file__))
import harvester as H  # noqa: E402

STRICT_RE = re.compile(r"kaplan|survival|progression-free|event-free", re.I)


def km_timepoints(meas_subset):
    """(nct,outcome) -> set of distinct parseable timepoints."""
    pts = defaultdict(set)
    for nct, oid, cls, title in zip(meas_subset["nct_id"], meas_subset["outcome_id"],
                                    meas_subset["classification"].fillna(""),
                                    meas_subset["title"].fillna("")):
        t = H.parse_timepoint_to_months(cls) or H.parse_timepoint_to_months(title)
        if t is not None:
            pts[(nct, oid)].add(round(t, 3))
    return pts


def tierA_and_hist(pts):
    """trials with >=3 timepoints on some outcome; binding K = max timepoints over that trial's
    qualifying outcomes (best reconstructable endpoint), histogrammed."""
    best = {}
    for (nct, oid), ts in pts.items():
        k = len(ts)
        if k >= 3:
            best[nct] = max(best.get(nct, 0), k)
    hist = defaultdict(int)
    for nct, k in best.items():
        hist['12+' if k >= 12 else str(k)] += 1
    return set(best), dict(sorted(hist.items(), key=lambda kv: (kv[0] == '12+', int(kv[0].rstrip('+'))))), best


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", default="realipd/census_full_aact.json")
    args = ap.parse_args(argv)

    try:
        from aact_kit import load_table, resolve_aact_location
        loc = resolve_aact_location()
    except Exception as e:  # noqa: BLE001
        print("ERROR: no AACT snapshot resolvable. Set AACT_ZIP/AACT_SQLITE/AACT_TSV_DIR.",
              file=sys.stderr)
        print(str(e), file=sys.stderr)
        return 2

    print(f"AACT location: {loc}")
    print("loading outcome_analyses ...")
    analyses = load_table("outcome_analyses", location=loc,
                          columns=["nct_id", "param_type", "method"])
    hr_mask = analyses["param_type"].fillna("").str.contains(H._HR_PARAM) | \
        analyses["method"].fillna("").str.contains("cox", case=False)
    hr_trials = set(analyses[hr_mask]["nct_id"])

    print("loading outcome_measurements (large, one pass) ...")
    meas = load_table("outcome_measurements", location=loc,
                      columns=["nct_id", "outcome_id", "title", "param_type", "classification"])
    blob = (meas["param_type"].fillna("") + " " + meas["title"].fillna("")).str.lower()

    universe = set(meas["nct_id"]) | hr_trials
    med_trials = set(meas[blob.str.contains("median")]["nct_id"])
    nar_rows = int(blob.str.contains("number at risk", regex=False).sum())

    strict = meas[blob.str.contains(STRICT_RE)]
    broad = meas[blob.str.contains(H._SURV_RE)]

    A_strict, hist_strict, _ = tierA_and_hist(km_timepoints(strict))
    A_broad, hist_broad, best_broad = tierA_and_hist(km_timepoints(broad))

    def band(best, lo, hi=None):
        return sum(1 for k in best.values() if (k >= lo and (hi is None or k <= hi)))

    tierB = (med_trials & hr_trials) - A_broad
    report = {
        "aact_location": str(loc),
        "universe_trials_with_results": len(universe),
        "structured_number_at_risk_rows": nar_rows,
        "trials_with_hr": len(hr_trials),
        "trials_with_median": len(med_trials),
        "tierA_strict_kaplan_survival_pfs_efs": len(A_strict),
        "tierA_broad_harvester_surv_re": len(A_broad),
        "tierB_median_plus_hr": len(tierB),
        "pct_of_universe": {
            "tierA_strict": round(100 * len(A_strict) / len(universe), 3),
            "tierA_broad": round(100 * len(A_broad) / len(universe), 3),
            "tierA_broad_plus_tierB": round(100 * (len(A_broad) + len(tierB)) / len(universe), 3),
        },
        "anchor_density_hist_broad": hist_broad,
        "anchor_density_hist_strict": hist_strict,
        "reliability_bands_broad": {
            "reliable_K_ge_6": band(best_broad, 6),
            "borderline_K_eq_5": band(best_broad, 5, 5),
            "weak_K_3_to_4": band(best_broad, 3, 4),
        },
        "note": "Binding K here = MAX timepoints over a trial's qualifying survival outcomes (best "
                "reconstructable endpoint), so this histogram is an UPPER bound on per-trial anchor "
                "density; the cohort census (validate/census_cohort.js) uses the stricter binding "
                "MIN-across-arms. Strict vs broad brackets how liberally a 'survival curve' is "
                "detected. 0 structured NAR rows across all of AACT is the exact zero-NAR finding.",
    }
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    print(json.dumps(report, indent=2))
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
