#!/usr/bin/env python
"""Harvest registry-reported MEDIAN survival per arm for the cohort, as external ground truth for
the median-accuracy validation (the cohort KM-curve outcome rarely co-reports a median, so we scan
ALL outcomes for a median-time measurement and key it by ctgov_group_code).

Writes cohort/registry_medians.json: { nct: { "OG000": median_months, ... } }
Run: python harvest/harvest_medians.py
"""
from __future__ import annotations
import json
import os
import sys
from collections import defaultdict

os.environ.setdefault("AACT_ZIP", r"C:\Users\mahmo\AACT\20260601_pipe-delimited-export.zip")
sys.path.insert(0, r"C:\Projects\aact-kit\src")
sys.path.insert(0, os.path.dirname(__file__))
import harvester as H  # noqa: E402
from aact_kit import load_table  # noqa: E402

COHORT = os.path.join(os.path.dirname(__file__), "..", "cohort")
TIME_UNITS = set(H._TIME_UNIT_TO_MONTHS)


def main():
    ncts = set(f[:-5] for f in os.listdir(COHORT) if f.endswith(".json")
               and not f.startswith(("manifest", "validation", "registry")))
    print(f"cohort: {len(ncts)} trials; reading outcome_measurements for median-time rows ...")
    m = load_table("outcome_measurements",
                   columns=["nct_id", "ctgov_group_code", "title", "param_type", "param_value", "units"])
    m = m[m["nct_id"].isin(ncts)]
    blob = (m["param_type"].fillna("") + " " + m["title"].fillna("")).str.lower()
    units = m["units"].fillna("").str.lower()
    is_median = blob.str.contains("median")
    is_time = units.apply(lambda u: any(tu in u for tu in TIME_UNITS))
    # MUST be a survival/time-to-EVENT median, NOT median follow-up/duration/age/exposure (those
    # also match "median"+time units and contaminated the first pass).
    SURV = ("survival", "progression-free", "progression free", "pfs", "overall survival",
            "disease-free", "event-free", "relapse-free", "recurrence-free", "time to progression",
            "time to event", "time to death", "time to relapse", "time to recurrence")
    EXCLUDE = ("follow-up", "follow up", "duration", "exposure", "age", "treatment", "on study",
               "on-study", "dosing", "time on", "response", "hospital", "number of", "adherence")
    is_surv = blob.apply(lambda b: any(k in b for k in SURV) and not any(k in b for k in EXCLUDE))
    med = m[is_median & is_time & is_surv]
    print(f"survival-median-time rows (filtered): {len(med):,}")

    out = defaultdict(dict)
    for _, r in med.iterrows():
        try:
            v = float(r["param_value"])
        except (TypeError, ValueError):
            continue
        if v <= 0 or v != v:  # nonpositive or NaN
            continue
        u = str(r["units"]).lower()
        factor = next((H._TIME_UNIT_TO_MONTHS[k] for k in H._TIME_UNIT_TO_MONTHS if k in u), 1.0)
        code = str(r.get("ctgov_group_code", "") or "")
        if code and code not in out[r["nct_id"]]:   # first reported median per arm
            out[r["nct_id"]][code] = round(v * factor, 4)

    out = {k: v for k, v in out.items() if v}
    with open(os.path.join(COHORT, "registry_medians.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    print(f"wrote registry medians for {len(out)} trials -> cohort/registry_medians.json")


if __name__ == "__main__":
    raise SystemExit(main())
