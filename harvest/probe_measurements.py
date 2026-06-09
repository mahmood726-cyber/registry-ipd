#!/usr/bin/env python
"""Diligence probe: how does AACT actually store survival / KM / number-at-risk data?
Determines whether Tier A (~0) is the true registry or a detector miss.
Run: python harvest/probe_measurements.py
"""
import os
import sys

os.environ.setdefault("AACT_ZIP", r"C:\Users\mahmo\AACT\20260601_pipe-delimited-export.zip")
sys.path.insert(0, r"C:\Projects\aact-kit\src")
from aact_kit import load_table  # noqa: E402

print("loading outcome_measurements (param_type/title/classification/units) ...")
m = load_table("outcome_measurements",
               columns=["nct_id", "outcome_id", "title", "param_type", "classification", "units"])
print(f"rows: {len(m):,}")

print("\n=== top 25 param_type values ===")
print(m["param_type"].fillna("(none)").value_counts().head(25).to_string())

print("\n=== top 20 units values ===")
print(m["units"].fillna("(none)").value_counts().head(20).to_string())

blob = (m["param_type"].fillna("") + " | " + m["title"].fillna("") + " | " + m["units"].fillna("")).str.lower()
for kw in ["at risk", "number at risk", "kaplan", "survival", "median", "hazard", "event-free", "progression-free", "percent"]:
    print(f"  measurements mentioning '{kw}': {blob.str.contains(kw, regex=False).sum():,}")

# how many (nct,outcome) have >=3 distinct classification values among survival-ish measures?
surv = m[blob.str.contains("kaplan|survival|progression-free|event-free", regex=True)]
print(f"\nsurvival-ish measurements: {len(surv):,}")
g = surv.groupby(["nct_id", "outcome_id"])["classification"].nunique()
print(f"  (nct,outcome) with >=3 distinct classification timepoints: {(g >= 3).sum():,}")
print(f"  distinct trials among those: {surv[surv.set_index(['nct_id','outcome_id']).index.isin(g[g>=3].index)]['nct_id'].nunique() if (g>=3).any() else 0:,}")

print("\n=== sample classification values on survival-ish rows (do they encode timepoints?) ===")
print(surv["classification"].fillna("(none)").value_counts().head(20).to_string())

print("\n=== example: one survival outcome's measurement rows ===")
if len(surv):
    key = g[g >= 3].index[0] if (g >= 3).any() else surv.set_index(["nct_id", "outcome_id"]).index[0]
    ex = m[(m["nct_id"] == key[0]) & (m["outcome_id"] == key[1])]
    print(f"trial {key[0]} outcome {key[1]}:")
    print(ex[["title", "param_type", "classification", "units"]].head(20).to_string())
