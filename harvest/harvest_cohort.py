#!/usr/bin/env python
"""Batch-harvest ALL Tier-A trials in ONE pass (ZIP backend re-reads the whole table per query,
so we read each table once and filter in memory). Writes cohort/<nct>.json + cohort/manifest.json.

Each trial JSON also carries condition text (for the oncology-OS demo pick) and the registry HR
(held-out ground truth for validation). Usage: python harvest/harvest_cohort.py [--limit N]
"""
from __future__ import annotations
import argparse
import json
import os
import sys

os.environ.setdefault("AACT_ZIP", r"C:\Users\mahmo\AACT\20260601_pipe-delimited-export.zip")
sys.path.insert(0, r"C:\Projects\aact-kit\src")
sys.path.insert(0, os.path.dirname(__file__))
import harvester as H  # noqa: E402
from aact_kit import load_table  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(__file__))
OUT = os.path.join(ROOT, "cohort")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)

    print("reading outcome_measurements (one full pass) ...")
    meas = load_table("outcome_measurements",
                      columns=["nct_id", "outcome_id", "result_group_id", "ctgov_group_code",
                               "classification", "title", "param_type", "param_value", "units"])
    blob = (meas["param_type"].fillna("") + " " + meas["title"].fillna("") + " "
            + meas["units"].fillna("")).str.lower()
    surv = meas[blob.str.contains(H._SURV_RE)]

    # Tier-A (nct, outcome) = >=3 parseable timepoints; choose best outcome per trial
    from collections import defaultdict
    tp = defaultdict(set)
    for nct, oid, cls, title in zip(surv["nct_id"], surv["outcome_id"],
                                    surv["classification"].fillna(""), surv["title"].fillna("")):
        t = H.parse_timepoint_to_months(cls) or H.parse_timepoint_to_months(title)
        if t is not None:
            tp[(nct, oid)].add(round(t, 3))
    best_outcome = {}
    for (nct, oid), pts in tp.items():
        if len(pts) >= 3 and len(pts) > best_outcome.get(nct, (0, None))[0]:
            best_outcome[nct] = (len(pts), oid)
    ncts = list(best_outcome)
    if args.limit:
        ncts = ncts[:args.limit]
    print(f"Tier-A trials: {len(best_outcome)} (processing {len(ncts)})")

    print("reading supporting tables (one pass each) ...")
    nctset = set(ncts)
    def slice_to(tbl, cols=None):
        df = load_table(tbl, columns=cols)
        return df[df["nct_id"].isin(nctset)]
    outcomes = slice_to("outcomes")
    groups = load_table("result_groups")
    groups = groups[groups["nct_id"].isin(nctset)]
    if "result_type" in groups:
        groups = groups[groups["result_type"].astype(str).str.lower().str.contains("outcome")]
    counts = slice_to("outcome_counts")
    analyses = slice_to("outcome_analyses")
    conds = slice_to("conditions", ["nct_id", "name"]) if True else None

    meas_t = meas[meas["nct_id"].isin(nctset)]
    cond_by_nct = conds.groupby("nct_id")["name"].apply(lambda s: "; ".join(map(str, s))).to_dict()

    manifest = []
    for i, nct in enumerate(ncts):
        n_tp, oid = best_outcome[nct]
        o = outcomes[outcomes["nct_id"] == nct]
        m = meas_t[meas_t["nct_id"] == nct]
        g = groups[groups["nct_id"] == nct]
        c = counts[counts["nct_id"] == nct]
        a = analyses[analyses["nct_id"] == nct]
        try:
            arms = H.assemble_arms(o, m, g, c, oid)
            hr = H.parse_hazard_ratio(a[a["outcome_id"] == oid]) if "outcome_id" in a else H.parse_hazard_ratio(a)
            hr = H.resolve_hr_direction(hr, arms)
        except Exception as e:  # noqa: BLE001
            print(f"  skip {nct}: {type(e).__name__}: {e}")
            continue
        trial = {
            "nct_id": nct, "source_url": f"https://clinicaltrials.gov/study/{nct}",
            "outcome_id": int(oid), "time_unit": "months", "condition": cond_by_nct.get(nct, ""),
            "arms": arms, "hr": hr,
        }
        with open(os.path.join(OUT, f"{nct}.json"), "w", encoding="utf-8") as f:
            json.dump(trial, f)
        manifest.append({
            "nct": nct, "outcome_id": int(oid), "n_arms": len(arms),
            "n_timepoints": n_tp, "has_hr": hr is not None,
            "N_total": sum(a.get("N") or 0 for a in arms),
            "min_final_S": min([(arm["km_points"][-1]["S"] if arm["km_points"] else 1.0) for arm in arms] or [1.0]),
            "condition": cond_by_nct.get(nct, "")[:120],
        })
        if (i + 1) % 25 == 0:
            print(f"  {i+1}/{len(ncts)}")
    with open(os.path.join(OUT, "manifest.json"), "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    n_hr = sum(1 for r in manifest if r["has_hr"])
    n_2arm = sum(1 for r in manifest if r["n_arms"] == 2)
    print(f"\nwrote {len(manifest)} trial JSONs to {OUT}")
    print(f"  with registry HR: {n_hr}   2-arm: {n_2arm}   2-arm AND HR: "
          f"{sum(1 for r in manifest if r['has_hr'] and r['n_arms']==2)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
