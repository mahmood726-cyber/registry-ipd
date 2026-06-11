#!/usr/bin/env python
"""Full-cohort independent validation: reconstructed HR + median vs the PUBLISHED HR + median (PubMed).

Extends the validation-grade PubMed check (112 trials) to EVERY reconstructed 2-arm Tier-A cohort trial
(~250). Crucially this includes the ~220 trials where AACT posts NO hazard ratio -- there the published
HR (from the trial's abstract) is the ONLY held-out truth, so this is a pure, registry-independent
validation at much larger n.

Pipeline (one snapshot pass): realipd/cohort_recon.json (node export of reconstructed HR + per-arm
medians) -> AACT study_references for PMIDs + outcome_measurements for the curve endpoint -> efetch
abstracts (cached, shared with pubmed_validation) -> endpoint-aware abstract_hr / abstract_median ->
compare. Emits realipd/cohort_pubmed_validation.json.

Honest limits carried over: published HR = primary non-covariate HR-with-CI (high-confidence = CI + <=2
candidates); median matched to the curve endpoint; both are triangulation, not a gold standard.

Run from repo root: python harvest/cohort_pubmed.py
"""
from __future__ import annotations
import json
import math
import os
import sys

os.environ.setdefault("AACT_ZIP", r"C:\Users\mahmo\AACT\20260601_pipe-delimited-export.zip")
sys.path.insert(0, r"C:\Projects\aact-kit\src")
sys.path.insert(0, os.path.dirname(__file__))
import harvester as H            # noqa: E402
import abstract_hr as AH         # noqa: E402
import abstract_median as AM     # noqa: E402
from pubmed_validation import fetch_abstracts, CACHE  # reuse efetch + shared cache  # noqa: E402
from aact_kit import load_table  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(__file__))
RIPD = os.path.join(ROOT, "realipd")
COHORT = os.path.join(ROOT, "cohort")


def _clean(p):
    p = str(p).strip()
    return p[:-2] if p.endswith(".0") else p


def _fold(a, b):
    return round(math.exp(abs(math.log(a) - math.log(b))), 3)


def main():
    recon = json.load(open(os.path.join(RIPD, "cohort_recon.json"), encoding="utf-8"))
    by_nct = {r["nct"]: r for r in recon}
    ncts = set(by_nct)
    print(f"reconstructed 2-arm trials: {len(ncts)}")

    # curve outcome_id per trial (from the cohort JSON) -> classify the curve's endpoint
    tte = {}
    for nct in ncts:
        try:
            tte[nct] = json.load(open(os.path.join(COHORT, f"{nct}.json"), encoding="utf-8")).get("outcome_id")
        except Exception:
            tte[nct] = None

    print("loading study_references (PMIDs) ...")
    sr = load_table("study_references", columns=["nct_id", "pmid", "reference_type"])
    sr = sr[sr["nct_id"].isin(ncts)]
    pmid_by = {}
    for nct, sub in sr.groupby("nct_id"):
        res = sub[sub["reference_type"].fillna("").str.lower().str.contains("result")]
        pick = res if len(res) else sub
        pm = _clean(pick.iloc[0]["pmid"]) if len(pick) else None
        if pm and pm != "nan":
            pmid_by[nct] = pm

    print("loading outcome_measurements (curve endpoint) ...")
    me = load_table("outcome_measurements", columns=["nct_id", "outcome_id", "title"])
    me = me[me["nct_id"].isin(ncts)]
    endpoint_by = {}
    for nct in ncts:
        oid = tte.get(nct)
        if oid is None:
            continue
        titles = me[(me["nct_id"] == nct) & (me["outcome_id"] == oid)]["title"].dropna().astype(str)
        fam = H.endpoint_family(" ".join(titles.tolist()))
        if fam:
            endpoint_by[nct] = fam

    cache = json.load(open(CACHE, encoding="utf-8")) if os.path.exists(CACHE) else {}
    cache = fetch_abstracts(sorted(set(pmid_by.values())), cache)

    rows = []
    for nct in ncts:
        pmid = pmid_by.get(nct)
        ab = cache.get(pmid, "") if pmid else ""
        if not ab:
            continue
        ep = endpoint_by.get(nct)
        r = by_nct[nct]
        pub_hr = AH.extract_hr(ab, endpoint=ep)
        pub_med = AM.extract_medians(ab, endpoint=ep) if ep else None
        row = {"nct": nct, "pmid": pmid, "condition": r.get("condition"), "curve_endpoint": ep,
               "recon_HR": r["recon_HR"], "registry_HR": r.get("registry_HR")}
        if pub_hr:
            conf = pub_hr["has_ci"] and pub_hr["n_hr_candidates"] <= 2
            row.update({"published_HR": pub_hr["value"], "hr_high_confidence": conf,
                        "n_hr_candidates": pub_hr["n_hr_candidates"],
                        "recon_vs_published_HR_fold": _fold(r["recon_HR"], pub_hr["value"])})
            if pub_hr["has_ci"]:
                row["recon_in_published_HR_CI"] = bool(pub_hr["ci_low"] <= r["recon_HR"] <= pub_hr["ci_high"])
        if pub_med and len(pub_med["medians"]) == 2 and pub_med["n_numbers"] == 2 and r.get("recon_medians"):
            ps = sorted(pub_med["medians"]); rs = sorted(r["recon_medians"])
            row.update({"published_medians": ps, "recon_medians": rs,
                        "median_arm_folds": [_fold(rs[0], ps[0]), _fold(rs[1], ps[1])]})
        rows.append(row)

    def med(xs):
        s = sorted(xs)
        return round(s[len(s) // 2], 3) if s else None

    hc = [r for r in rows if r.get("hr_high_confidence")]
    hc_ci = [r for r in hc if "recon_in_published_HR_CI" in r]
    reg_hc = [r for r in hc if r.get("registry_HR") is not None]
    med_rows = [r for r in rows if "median_arm_folds" in r]
    arm_folds = [f for r in med_rows for f in r["median_arm_folds"]]
    no_registry_hr = [r for r in hc if r.get("registry_HR") is None]

    summary = {
        "reconstructed_trials": len(ncts),
        "with_pmid": len(pmid_by),
        "with_abstract": sum(1 for r in rows),
        "published_HR_high_confidence": len(hc),
        "of_which_no_registry_HR": len(no_registry_hr),     # pure registry-independent validations
        "HR": {
            "recon_vs_published_median_fold": med([r["recon_vs_published_HR_fold"] for r in hc]),
            "recon_within_published_95CI": f"{sum(1 for r in hc_ci if r['recon_in_published_HR_CI'])}/{len(hc_ci)}",
            "registry_vs_published_agree": f"{sum(1 for r in reg_hc if (r['registry_HR'] - 1) * (r['published_HR'] - 1) >= 0 and _fold(r['registry_HR'], r['published_HR']) < 1.25)}/{len(reg_hc)}",
        },
        "median": {
            "trials": len(med_rows), "arm_medians": len(arm_folds),
            "median_arm_fold": med(arm_folds),
            "within_20pct": f"{sum(1 for f in arm_folds if f < 1.2)}/{len(arm_folds)}",
        },
        "note": "Full reconstructed cohort vs the PUBLISHED HR + median (PubMed abstracts), endpoint-matched. "
                "Most trials have NO registry HR, so the published HR is the only held-out truth -- a pure "
                "registry-independent validation. High-confidence published HR = CI present + <=2 candidates; "
                "median matched to the curve endpoint. Triangulation, not a gold standard.",
    }
    json.dump({"summary": summary, "rows": rows},
              open(os.path.join(RIPD, "cohort_pubmed_validation.json"), "w"), indent=1)
    print("\n" + json.dumps(summary, indent=1))
    print("\nwrote realipd/cohort_pubmed_validation.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
