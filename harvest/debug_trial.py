#!/usr/bin/env python
"""Debug the real AACT structure for one trial so the harvester can extract its KM curve.
Usage: python harvest/debug_trial.py NCT00725985
"""
import os
import sys

os.environ.setdefault("AACT_ZIP", r"C:\Users\mahmo\AACT\20260601_pipe-delimited-export.zip")
sys.path.insert(0, r"C:\Projects\aact-kit\src")
from aact_kit import load_table  # noqa: E402

nct = sys.argv[1] if len(sys.argv) > 1 else "NCT00725985"
w = {"nct_id": nct}

outcomes = load_table("outcomes", where=w)
print(f"=== outcomes ({len(outcomes)}) ===")
print(outcomes[[c for c in ["id", "outcome_type", "title", "param_type", "units"] if c in outcomes]].to_string())

groups = load_table("result_groups", where=w)
print(f"\n=== result_groups ({len(groups)}) — columns: {list(groups.columns)} ===")
if "result_type" in groups:
    print("result_type counts:", groups["result_type"].value_counts().to_dict())
print(groups[[c for c in ["id", "ctgov_group_code", "result_type", "title", "outcome_id"] if c in groups]].head(12).to_string())

meas = load_table("outcome_measurements", where=w)
print(f"\n=== outcome_measurements ({len(meas)}) — columns: {list(meas.columns)} ===")
# find the survival-ish outcome with most parseable timepoints
sys.path.insert(0, os.path.dirname(__file__))
import harvester as H  # noqa: E402
blob = (meas["param_type"].fillna("") + " " + meas["title"].fillna("")).str.lower()
surv = meas[blob.str.contains("kaplan|survival|progression-free|event-free", regex=True)]
print(f"survival-ish measurement rows: {len(surv)}")
if len(surv):
    oid = surv["outcome_id"].value_counts().index[0]
    sub = meas[meas["outcome_id"] == oid]
    print(f"\n--- measurements for survival outcome_id={oid} ({len(sub)} rows) ---")
    cols = [c for c in ["result_group_id", "classification", "category", "title", "param_type", "param_value", "units"] if c in sub]
    print(sub[cols].head(30).to_string())
    print("\ndistinct result_group_id on this outcome:", sorted(sub["result_group_id"].unique())[:10])
    print("distinct classification:", sub["classification"].fillna("(none)").unique()[:15])
