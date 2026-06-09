#!/usr/bin/env python
"""Harvest registry-reported MEDIAN survival per arm, matched to the SAME ENDPOINT as the
reconstructed KM curve — the genuinely external (non-circular) median ground truth.

The earlier version keyed medians by arm code across ANY outcome, so a PFS curve got compared to an
OS / follow-up median (=> spurious ~48% error). Here we match the median to the trial's *curve
outcome*:
  Tier 1 — a median reported in the SAME outcome_id as the curve (definitionally same endpoint);
  Tier 2 — a median in another outcome with the SAME endpoint class (PFS/OS) AND a similar title
           (Jaccard token overlap), keyed to the same arm (ctgov_group_code).

Writes cohort/registry_medians.json: { nct: { "OG000": median_months, ... } } (same-endpoint only),
plus cohort/registry_medians_meta.json documenting each match (tier, titles) for auditability.
Run: python harvest/harvest_medians.py
"""
from __future__ import annotations
import json
import os
import re
import sys
from collections import defaultdict

os.environ.setdefault("AACT_ZIP", r"C:\Users\mahmo\AACT\20260601_pipe-delimited-export.zip")
sys.path.insert(0, r"C:\Projects\aact-kit\src")
sys.path.insert(0, os.path.dirname(__file__))
import harvester as H  # noqa: E402
from aact_kit import load_table  # noqa: E402

COHORT = os.path.join(os.path.dirname(__file__), "..", "cohort")
TIME_UNITS = set(H._TIME_UNIT_TO_MONTHS)
_STOP = {"median", "the", "of", "for", "with", "in", "to", "by", "at", "per", "and", "time",
         "rate", "estimate", "estimates", "kaplan", "meier", "probability", "percentage",
         "participants", "number", "from", "months", "days", "years", "weeks", "a", "an"}


def endpoint_kind(title):
    t = (title or "").lower()
    if any(k in t for k in ("progression-free", "progression free", "pfs", "time to progression",
                            "disease-free", "event-free", "relapse-free", "recurrence-free",
                            "probability of disability", "probability of progression", "time to")):
        return "PFS"
    if "overall survival" in t or re.search(r"\bos\b", t) or "survival" in t:
        return "OS"
    return "PFS"


def toks(s):
    return {w for w in re.findall(r"[a-z0-9-]+", str(s or "").lower()) if len(w) > 2 and w not in _STOP}


def jaccard(a, b):
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def to_months(v, units):
    u = str(units).lower()
    return v * next((H._TIME_UNIT_TO_MONTHS[k] for k in H._TIME_UNIT_TO_MONTHS if k in u), 1.0)


def main():
    ncts = {}
    for f in os.listdir(COHORT):
        if f.endswith(".json") and not f.startswith(("manifest", "validation", "registry")):
            t = json.load(open(os.path.join(COHORT, f), encoding="utf-8"))
            ncts[t["nct_id"]] = t.get("outcome_id")
    print(f"cohort: {len(ncts)} trials; reading outcomes + outcome_measurements ...")
    outcomes = load_table("outcomes", columns=["id", "nct_id", "title"])
    outcomes = outcomes[outcomes["nct_id"].isin(ncts)]
    title_by_oid = dict(zip(outcomes["id"], outcomes["title"]))

    m = load_table("outcome_measurements",
                   columns=["nct_id", "outcome_id", "ctgov_group_code", "title", "param_type", "param_value", "units"])
    m = m[m["nct_id"].isin(ncts)]
    blob = (m["param_type"].fillna("") + " " + m["title"].fillna("")).str.lower()
    units = m["units"].fillna("").str.lower()
    is_median = blob.str.contains("median")
    is_time = units.apply(lambda u: any(tu in u for tu in TIME_UNITS))
    EXCLUDE = ("follow-up", "follow up", "duration", "exposure", "age", "treatment", "on study",
               "on-study", "dosing", "time on", "response", "hospital", "number of", "adherence")
    SURV = ("survival", "progression-free", "progression free", "pfs", "overall survival",
            "disease-free", "event-free", "relapse-free", "recurrence-free", "time to")
    is_surv = blob.apply(lambda b: any(k in b for k in SURV) and not any(k in b for k in EXCLUDE))
    med = m[is_median & is_time & is_surv]
    print(f"survival-median-time rows: {len(med):,}")

    out, meta = defaultdict(dict), defaultdict(dict)
    for nct, curve_oid in ncts.items():
        curve_title = title_by_oid.get(curve_oid, "")
        curve_kind = endpoint_kind(curve_title)
        ctoks = toks(curve_title)
        rows = med[med["nct_id"] == nct]
        if rows.empty:
            continue
        # candidates per arm with a match score (tier 1 = same outcome; tier 2 = same kind + title sim)
        for code, grp in rows.groupby("ctgov_group_code"):
            code = str(code or "")
            if not code:
                continue
            best = None  # (tier, score, value_months, med_title)
            for _, r in grp.iterrows():
                try:
                    v = float(r["param_value"])
                except (TypeError, ValueError):
                    continue
                if v <= 0 or v != v:
                    continue
                vm = round(to_months(v, r["units"]), 4)
                mtitle = r["title"]
                if r["outcome_id"] == curve_oid:
                    cand = (1, 1.0, vm, mtitle)
                else:
                    if endpoint_kind(mtitle) != curve_kind:
                        continue
                    sim = jaccard(ctoks, toks(mtitle))
                    if sim < 0.34:
                        continue
                    cand = (2, sim, vm, mtitle)
                if best is None or (cand[0], cand[1]) > (best[0], best[1]):
                    best = cand
            if best:
                out[nct][code] = best[2]
                meta[nct][code] = {"tier": best[0], "score": round(best[1], 3),
                                   "median_months": best[2], "curve_title": curve_title[:90],
                                   "median_title": str(best[3])[:90]}
    out = {k: v for k, v in out.items() if v}
    json.dump(out, open(os.path.join(COHORT, "registry_medians.json"), "w", encoding="utf-8"), indent=1)
    json.dump(meta, open(os.path.join(COHORT, "registry_medians_meta.json"), "w", encoding="utf-8"), indent=1)
    t1 = sum(1 for n in meta for c in meta[n] if meta[n][c]["tier"] == 1)
    t2 = sum(1 for n in meta for c in meta[n] if meta[n][c]["tier"] == 2)
    print(f"same-endpoint medians: {len(out)} trials, {t1} arms tier-1 (same outcome), {t2} arms tier-2 (matched endpoint)")


if __name__ == "__main__":
    raise SystemExit(main())
