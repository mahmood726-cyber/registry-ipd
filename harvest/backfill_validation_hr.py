#!/usr/bin/env python
"""Backfill the held-out registry HR for the VALIDATION-GRADE population using the sibling-outcome fix.

The full-snapshot census finds 112 trials posting both a reconstructable curve AND a hazard ratio, but
the production gallery only validated 30 -- because the harvester scoped the HR lookup to the curve's own
outcome and dropped HRs posted in a sibling survival outcome (now fixed: harvester.select_trial_hr).

This applies that fix in ONE snapshot pass over just the validation-grade NCTs (cheap vs re-harvesting
the whole cohort): for each, load its outcome_analyses + outcome_measurements, pick the curve-outcome HR
or a survival-sibling HR, and resolve direction from the EXISTING cohort JSON's arm roles. Emits
realipd/validation_hr_backfill.json -- consumed by validate/gallery_expanded.js to re-score the
reconstructed HR against this enlarged held-out set.

Honest: sibling-sourced HRs are flagged (`from_sibling_outcome`) because within survival outcomes the
fallback cannot distinguish OS from PFS; gallery_expanded reports them as a separate, caveated tier.

Run: python harvest/backfill_validation_hr.py   (requires the AACT snapshot)
"""
from __future__ import annotations
import json, os, sys

os.environ.setdefault("AACT_ZIP", r"C:\Users\mahmo\AACT\20260601_pipe-delimited-export.zip")
sys.path.insert(0, r"C:\Projects\aact-kit\src")
sys.path.insert(0, os.path.dirname(__file__))
import harvester as H  # noqa: E402
from aact_kit import load_table  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(__file__))
COHORT = os.path.join(ROOT, "cohort")
CENSUS = os.path.join(ROOT, "realipd", "census_full_aact.json")
OUT = os.path.join(ROOT, "realipd", "validation_hr_backfill.json")


def main():
    cen = json.load(open(CENSUS, encoding="utf-8"))
    vg = cen.get("validation_grade_ncts_broad") or []
    if not vg:
        print("no validation_grade_ncts_broad in census; re-run harvest/census_full_aact.py", file=sys.stderr)
        return 2
    vgset = set(vg)
    print(f"validation-grade NCTs: {len(vgset)}")

    print("loading outcome_analyses (filtered) ...")
    an = load_table("outcome_analyses", columns=["nct_id", "outcome_id", "param_type", "param_value",
                                                 "ci_lower_limit", "ci_upper_limit", "ci_percent",
                                                 "p_value", "method", "ci_n_sides"])
    an = an[an["nct_id"].isin(vgset)]
    print("loading outcome_measurements (filtered) ...")
    me = load_table("outcome_measurements", columns=["nct_id", "outcome_id", "title", "param_type"])
    me = me[me["nct_id"].isin(vgset)]

    out = {}
    n_curve = n_sibling = n_nodir = n_nohr = n_nocohort = 0
    for nct in sorted(vgset):
        fp = os.path.join(COHORT, f"{nct}.json")
        if not os.path.isfile(fp):
            n_nocohort += 1
            continue
        trial = json.load(open(fp, encoding="utf-8"))
        tte = trial.get("outcome_id")
        arms = trial.get("arms") or []
        sub_an = an[an["nct_id"] == nct]
        sub_me = me[me["nct_id"] == nct]
        surv_ids = H._survival_outcome_ids(sub_me)
        hr, from_sib = H.select_trial_hr(sub_an, tte, surv_ids)
        if hr is None:
            n_nohr += 1
            continue
        hr = H.resolve_hr_direction(hr, arms)
        if from_sib:
            hr["from_sibling_outcome"] = True
            n_sibling += 1
        else:
            n_curve += 1
        if hr.get("favors_arm_id") is None:
            n_nodir += 1
        out[nct] = {"value": hr["value"], "ci_low": hr.get("ci_low"), "ci_high": hr.get("ci_high"),
                    "favors_arm_id": hr.get("favors_arm_id"),
                    "from_sibling_outcome": bool(hr.get("from_sibling_outcome"))}

    report = {
        "n_validation_grade": len(vgset),
        "n_hr_recovered": len(out),
        "from_curve_outcome": n_curve,
        "from_sibling_outcome": n_sibling,
        "direction_unresolved": n_nodir,
        "no_hr_found": n_nohr,
        "not_in_cohort": n_nocohort,
        "hr": out,
    }
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=1)
    print(f"\nrecovered HR for {len(out)}/{len(vgset)} validation-grade trials")
    print(f"  from curve outcome: {n_curve} | from survival sibling: {n_sibling} | "
          f"direction unresolved: {n_nodir} | no HR: {n_nohr} | not in cohort: {n_nocohort}")
    print(f"  wrote {os.path.relpath(OUT, ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
