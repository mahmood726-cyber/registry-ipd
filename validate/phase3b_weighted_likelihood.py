#!/usr/bin/env python
"""PHASE 3b (step 1): pool each trial's reconstruction imputations with spec-collapse's weighted_likelihood.

The M imputations of one curve are correlated specs of ONE dataset. Inverse-variance pooling them collapses
the variance by ~M (false precision); spec-collapse-atlas's weighted_likelihood (mixture, never narrower
than one draw) is the honest within-trial aggregator. This loads the JS-exported imputations
(phase3b_imputations.json), pools each cohort BOTH ways, and checks whether the pooled interval covers the
held-out true log-HR. The honest reading: the within-trial reconstruction uncertainty must be aggregated by
weighted_likelihood, then carried (as one (theta, var) per trial) into the across-trial metaRE pool.

Cross-project reuse: imports C:\\Projects\\spec-collapse-atlas. Run from repo root:
  node validate/phase3b_export_imputations.js && python validate/phase3b_weighted_likelihood.py
"""
import io
import json
import math
import os
import sys

SPEC_COLLAPSE = r"C:\Projects\spec-collapse-atlas"
HERE = os.path.dirname(__file__)


def _load_aggregators():
    """Import spec-collapse-atlas lazily (kept out of module scope so importing this file is side-effect
    free — no sys.stdout reassignment, no hard import that breaks pytest capture)."""
    sys.path.insert(0, SPEC_COLLAPSE)
    from spec_collapse.aggregators import weighted_likelihood, naive_ivre_pool  # noqa: E402
    return weighted_likelihood, naive_ivre_pool


def covers(res, x):
    return bool(float(res["ci_low"]) <= x <= float(res["ci_high"]))


def main():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    weighted_likelihood, naive_ivre_pool = _load_aggregators()
    data = json.load(open(os.path.join(HERE, "phase3b_imputations.json"), encoding="utf-8"))
    rows = []
    for ds, d in data.items():
        if "specs" not in d:
            continue
        specs, true_log = d["specs"], d["true_logHR"]
        wl = weighted_likelihood(specs)            # honest mixture
        iv = naive_ivre_pool(specs)                # the cardinal sin (variance collapses by ~M)
        rows.append({
            "ds": d["ds"], "true_HR": d["true_HR"], "M": len(specs),
            "wl": {"HR": round(math.exp(wl["theta"]), 2),
                   "CI": [round(math.exp(wl["ci_low"]), 2), round(math.exp(wl["ci_high"]), 2)],
                   "width_log": round(wl["ci_high"] - wl["ci_low"], 3),
                   "var": round(wl["var"], 4), "covers_true": covers(wl, true_log)},
            "iv": {"HR": round(math.exp(iv["theta"]), 2),
                   "CI": [round(math.exp(iv["ci_low"]), 2), round(math.exp(iv["ci_high"]), 2)],
                   "width_log": round(iv["ci_high"] - iv["ci_low"], 3),
                   "var": round(iv["var"], 4), "covers_true": covers(iv, true_log)},
            "collapse_ratio": round((wl["ci_high"] - wl["ci_low"]) / (iv["ci_high"] - iv["ci_low"]), 1),
        })
    summary = {
        "n_cohorts": len(rows),
        "weighted_likelihood_covers_true": sum(1 for r in rows if r["wl"]["covers_true"]),
        "naive_ivre_covers_true": sum(1 for r in rows if r["iv"]["covers_true"]),
        "mean_collapse_ratio": round(sum(r["collapse_ratio"] for r in rows) / max(len(rows), 1), 1),
        "cohorts": rows,
    }
    json.dump(summary, open(os.path.join(HERE, "phase3b_weighted_likelihood_results.json"), "w"), indent=2)

    print("=== Phase 3b: within-trial imputation pooling — weighted_likelihood vs naive IV ===\n")
    print(f"  reused spec-collapse-atlas weighted_likelihood (cross-project)\n")
    print("  cohort     true HR   naive-IV CI (covers?)        weighted-likelihood CI (covers?)   widen")
    for r in rows:
        iv, wl = r["iv"], r["wl"]
        print(f"  {r['ds']:<9}  {r['true_HR']:<7}  [{iv['CI'][0]}, {iv['CI'][1]}] "
              f"{'Y' if iv['covers_true'] else 'N':<2}                  "
              f"[{wl['CI'][0]}, {wl['CI'][1]}] {'Y' if wl['covers_true'] else 'N':<2}              {r['collapse_ratio']}x")
    print(f"\n  covers true effect:  naive-IV {summary['naive_ivre_covers_true']}/{len(rows)}   "
          f"weighted-likelihood {summary['weighted_likelihood_covers_true']}/{len(rows)}")
    print(f"  mean interval widening (honest / naive): {summary['mean_collapse_ratio']}x")
    print("\n  Reading: IV-pooling the M correlated imputations collapses the variance (false precision);")
    print("  weighted_likelihood keeps the honest reconstruction interval. This is the within-trial")
    print("  aggregator; its (theta, var) per trial is what the across-trial metaRE pool then consumes.")
    print("\n  wrote validate/phase3b_weighted_likelihood_results.json")
    return summary


if __name__ == "__main__":
    main()
