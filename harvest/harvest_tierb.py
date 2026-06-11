#!/usr/bin/env python
"""Harvest the TIER-B population from AACT: trials that post a median + hazard ratio + N but NO
Kaplan-Meier curve (the ~3,263 trials the coverage census counts). This is the scale-up beyond the
Tier-A curve cohort (which is already saturated: 595 harvested >= 514 broad-census Tier-A).

For each trial with a survival-type MEDIAN reported in >=2 arms AND a hazard-ratio analysis, we build a
Tier-B trial JSON: two arms (median + N + total_events + follow_up) and the trial HR. These reconstruct
via the engine's parametric (exponential) Tier-B path -- lower fidelity than curve-based Tier A
(RMST ~7%, HR imposed), but it extends the reconstructable population several-fold.

Writes tierb_cohort/<nct>.json + tierb_cohort/manifest.json. Requires the AACT snapshot.
Run: python harvest/harvest_tierb.py [--limit N]
"""
from __future__ import annotations
import argparse, json, os, sys
from collections import defaultdict

os.environ.setdefault("AACT_ZIP", r"C:\Users\mahmo\AACT\20260601_pipe-delimited-export.zip")
sys.path.insert(0, r"C:\Projects\aact-kit\src")
sys.path.insert(0, os.path.dirname(__file__))
import harvester as H  # noqa: E402
from aact_kit import load_table  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(__file__))
OUT = os.path.join(ROOT, "tierb_cohort")
MEDIAN_RE = __import__("re").compile(r"median", __import__("re").I)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--limit", type=int, default=None); args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    print("loading outcome_analyses (HR) ...")
    analyses = load_table("outcome_analyses", columns=["nct_id", "outcome_id", "param_type", "param_value",
                                                       "ci_lower_limit", "ci_upper_limit", "method"])
    hr_by_nct = {}
    for nct, sub in analyses[analyses["param_type"].fillna("").str.contains(H._HR_PARAM)].groupby("nct_id"):
        hr = H.parse_hazard_ratio(sub)
        if hr and hr.get("value") is not None:
            hr_by_nct[nct] = hr
    print(f"  trials with a parseable HR: {len(hr_by_nct)}")

    print("loading outcome_measurements (medians, one pass) ...")
    meas = load_table("outcome_measurements", columns=["nct_id", "outcome_id", "ctgov_group_code",
                                                       "title", "param_type", "param_value", "units"])
    blob = (meas["param_type"].fillna("") + " " + meas["title"].fillna("")).str.lower()
    med = meas[blob.str.contains("median") & blob.str.contains(H._SURV_RE)]
    # median per (nct, outcome, group): value in months
    medians = defaultdict(dict)            # nct -> {group_code: median_months}
    for nct, oid, grp, pv, units, title in zip(med["nct_id"], med["outcome_id"], med["ctgov_group_code"].fillna(""),
                                               med["param_value"], med["units"].fillna(""), med["title"].fillna("")):
        if nct not in hr_by_nct:
            continue
        m = None
        try:
            v = float(str(pv).replace(",", ""))
            m = v * H._TIME_UNIT_TO_MONTHS.get(units.lower(), 1.0) if units.lower() in H._TIME_UNIT_TO_MONTHS else v
        except Exception:
            m = H.parse_timepoint_to_months(str(pv), None)
        if m is not None and m > 0 and grp:
            medians[nct].setdefault(grp, m)

    print("loading outcome_counts (N per arm) ...")
    counts = load_table("outcome_counts", columns=["nct_id", "ctgov_group_code", "count"])
    N_by = defaultdict(dict)
    for nct, grp, cnt in zip(counts["nct_id"], counts["ctgov_group_code"].fillna(""), counts["count"]):
        n = H.parse_count(cnt)
        if n and grp:
            N_by[nct][grp] = max(N_by[nct].get(grp, 0), n)

    written = 0; manifest = []
    for nct, gm in medians.items():
        if len(gm) < 2:
            continue
        groups = sorted(gm.items(), key=lambda kv: kv[1])     # by median asc
        (gL, mL), (gH, mH) = groups[0], groups[-1]            # most-separated pair
        NL, NH = N_by.get(nct, {}).get(gL), N_by.get(nct, {}).get(gH)
        if not NL or not NH or NL < 20 or NH < 20:
            continue
        hr = hr_by_nct[nct]
        trial = {
            "nct_id": nct, "source_url": f"https://clinicaltrials.gov/study/{nct}", "tier_hint": "B",
            "time_unit": "months",
            "arms": [
                {"arm_id": gH, "role": "experimental", "N": int(NH), "total_events": None,
                 "follow_up_max": round(mH * 3, 1), "median": {"value": round(mH, 2)}, "km_points": [], "nar_points": []},
                {"arm_id": gL, "role": "comparator", "N": int(NL), "total_events": None,
                 "follow_up_max": round(mL * 3, 1), "median": {"value": round(mL, 2)}, "km_points": [], "nar_points": []},
            ],
            "hr": {"value": hr["value"], "ci_low": hr.get("ci_low"), "ci_high": hr.get("ci_high"),
                   "favors_arm_id": gH if hr["value"] < 1 else gL},
        }
        with open(os.path.join(OUT, f"{nct}.json"), "w", encoding="utf-8") as f:
            json.dump(trial, f)
        manifest.append({"nct": nct, "median_exp": round(mH, 1), "median_ctl": round(mL, 1),
                         "N_exp": int(NH), "N_ctl": int(NL), "hr": hr["value"]})
        written += 1
        if args.limit and written >= args.limit:
            break

    with open(os.path.join(OUT, "manifest.json"), "w", encoding="utf-8") as f:
        json.dump({"n_trials": written, "trials": manifest}, f, indent=1)
    print(f"\nwrote {written} Tier-B trials to tierb_cohort/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
