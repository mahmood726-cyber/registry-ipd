#!/usr/bin/env python
"""Harvest the AACT-posted per-arm MEDIAN survival for the reconstructed cohort (registry-side cross-check).

The independent PubMed median validation is strong but small (~5 trials with a clean two-arm abstract
pair). A larger-n cross-check: the registry ITSELF posts median survival for many trials -- but, like the
hazard ratio, usually in a SIBLING survival outcome, not the curve outcome (so the harvester's
curve-scoped extraction left arms[].median empty for every Tier-A cohort trial). This recovers it in one
snapshot pass: for each reconstructed trial, find the median rows in a survival outcome whose endpoint
MATCHES the reconstructed curve (OS curve -> OS median, never a PFS median), take the two per-arm
medians, convert to months. Emits realipd/registry_medians.json.

This is registry-provenance (not independent like PubMed), but the reconstruction does NOT use the posted
median -- Tier A reconstructs from the KM timepoints -- so "reconstructed median vs separately-posted
median" is a genuine recovery check at cohort scale. validate/registry_median_validation.js compares.

Run from repo root: python harvest/registry_medians.py
"""
from __future__ import annotations
import json
import os
import re
import sys

os.environ.setdefault("AACT_ZIP", r"C:\Users\mahmo\AACT\20260601_pipe-delimited-export.zip")
sys.path.insert(0, r"C:\Projects\aact-kit\src")
sys.path.insert(0, os.path.dirname(__file__))
import harvester as H            # noqa: E402
from aact_kit import load_table  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(__file__))
RIPD = os.path.join(ROOT, "realipd")
COHORT = os.path.join(ROOT, "cohort")
_MEDIAN = re.compile(r"median", re.IGNORECASE)


def main():
    recon = json.load(open(os.path.join(RIPD, "cohort_recon.json"), encoding="utf-8"))
    ncts = {r["nct"] for r in recon}
    tte = {}
    for nct in ncts:
        try:
            tte[nct] = json.load(open(os.path.join(COHORT, f"{nct}.json"), encoding="utf-8")).get("outcome_id")
        except Exception:
            tte[nct] = None
    print(f"reconstructed trials: {len(ncts)}")

    print("loading outcome_measurements (medians) ...")
    me = load_table("outcome_measurements", columns=["nct_id", "outcome_id", "ctgov_group_code",
                                                     "title", "param_type", "param_value", "units"])
    me = me[me["nct_id"].isin(ncts)]
    blob = (me["param_type"].fillna("") + " " + me["title"].fillna("")).str.lower()
    surv = me[blob.str.contains(H._SURV_RE)]
    med = surv[(surv["param_type"].fillna("") + " " + surv["title"].fillna("")).str.contains(_MEDIAN)]

    def to_months(v, units):
        try:
            x = float(str(v).replace(",", ""))
        except (TypeError, ValueError):
            return None
        u = (units or "").lower()
        return x * H._TIME_UNIT_TO_MONTHS.get(u, 1.0) if u in H._TIME_UNIT_TO_MONTHS else x

    out = {}
    n_curve_ep = 0
    for nct in ncts:
        sub = med[med["nct_id"] == nct]
        if sub.empty:
            continue
        # endpoint of the reconstructed curve, from its outcome title
        oid = tte.get(nct)
        ctitles = me[(me["nct_id"] == nct) & (me["outcome_id"] == oid)]["title"].dropna().astype(str)
        curve_fam = H.endpoint_family(" ".join(ctitles.tolist()))
        if not curve_fam:
            continue
        # endpoint family of each median-bearing outcome; keep outcomes matching the curve endpoint
        cand = {}
        for oid2, g in sub.groupby("outcome_id"):
            otitles = me[(me["nct_id"] == nct) & (me["outcome_id"] == oid2)]["title"].dropna().astype(str)
            if H.endpoint_family(" ".join(otitles.tolist())) != curve_fam:
                continue
            vals = {}
            for grp, pv, un in zip(g["ctgov_group_code"].fillna(""), g["param_value"], g["units"].fillna("")):
                m = to_months(pv, un)
                if m is not None and m > 0 and grp:
                    vals.setdefault(grp, m)
            if len(vals) >= 2:
                cand[oid2] = sorted(vals.values())
        if not cand:
            continue
        # prefer the curve outcome's own medians, else any same-endpoint outcome with a 2-arm pair
        meds = cand.get(oid) or next(iter(cand.values()))
        out[nct] = {"endpoint": curve_fam, "medians": [round(meds[0], 2), round(meds[-1], 2)]}
        n_curve_ep += 1

    json.dump(out, open(os.path.join(RIPD, "registry_medians.json"), "w"), indent=1)
    print(f"recovered AACT-posted same-endpoint two-arm medians for {len(out)} trials")
    print("wrote realipd/registry_medians.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
