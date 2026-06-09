#!/usr/bin/env python
"""Censoring-informed enhancement: derive per-arm total events from `drop_withdrawals`
(participant flow) and patch them into the harvested cohort JSONs.

Curve-only reconstruction does not know the censoring level, so it over-counts events and
ATTENUATES the hazard ratio. The registry reports *why* participants left the study; for a
time-to-event endpoint the event-type reasons (progression/death/relapse for PFS; death for OS)
give the event count, and the rest is censoring. Proof of concept on RADIANT-4 moved the
reconstructed HR 0.68 -> 0.47 vs a registry 0.48.

Maps flow-group (FGxxx) -> outcome-group (OGxxx) by numeric suffix, validated against arm N.
Writes total_events back into cohort/<nct>.json (only when derivable). Run:
    python harvest/add_event_counts.py [cohort_dir]
"""
from __future__ import annotations
import json
import os
import re
import sys
from collections import defaultdict

os.environ.setdefault("AACT_ZIP", r"C:\Users\mahmo\AACT\20260601_pipe-delimited-export.zip")
sys.path.insert(0, r"C:\Projects\aact-kit\src")
from aact_kit import load_table  # noqa: E402

COHORT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "..", "cohort")

_EVENT_PFS = ("progression", "relapse", "recurrence", "death", "died", "disease")
_EVENT_OS = ("death", "died")
_SUFFIX = re.compile(r"(\d+)\s*$")


def endpoint_kind(title: str) -> str:
    t = (title or "").lower()
    if any(k in t for k in ("progression-free", "progression free", "pfs", "time to progression",
                            "disease-free", "event-free", "relapse-free", "probability of disability",
                            "probability of progression", "time to")):
        return "PFS"
    if "overall survival" in t or re.search(r"\bos\b", t) or "survival" in t:
        return "OS"
    return "PFS"  # default: count progression+death (most ct.gov KM curves are PFS-like)


def main():
    files = [f for f in os.listdir(COHORT) if f.endswith(".json")
             and not f.startswith(("manifest", "validation", "registry"))]
    ncts = set(f[:-5] for f in files)
    print(f"cohort trials: {len(ncts)}; reading drop_withdrawals, outcomes, result_groups ...")
    dw = load_table("drop_withdrawals")
    dw = dw[dw["nct_id"].isin(ncts)]
    outcomes = load_table("outcomes", columns=["id", "nct_id", "title"])
    outcomes = outcomes[outcomes["nct_id"].isin(ncts)]
    title_by_oid = dict(zip(outcomes["id"], outcomes["title"]))
    # N-MATCHED MAPPING: map a participant-flow group to an outcome arm by its initial randomized
    # count (milestone STARTED) == the arm's analysis N. Robust to the title/ordering mismatches
    # that made title/suffix mapping unreliable.
    mil = load_table("milestones")
    mil = mil[mil["nct_id"].isin(ncts)]
    started = mil[mil["title"].astype(str).str.upper().str.strip() == "STARTED"]
    started_n = defaultdict(dict)   # nct -> { FGcode: max STARTED count }
    for _, r in started.iterrows():
        code = str(r.get("ctgov_group_code", "") or "")
        try:
            c = int(float(r.get("count")))
        except (TypeError, ValueError):
            continue
        if code:
            started_n[r["nct_id"]][code] = max(started_n[r["nct_id"]].get(code, 0), c)

    patched = skipped = 0
    for f in files:
        nct = f[:-5]
        p = os.path.join(COHORT, f)
        trial = json.load(open(p, encoding="utf-8"))
        kind = endpoint_kind(title_by_oid.get(trial.get("outcome_id"), ""))
        ev_kw = _EVENT_OS if kind == "OS" else _EVENT_PFS
        sub = dw[dw["nct_id"] == nct]
        # idempotent: clear any prior derived value first
        for arm in trial["arms"]:
            arm["total_events"] = None
        if sub.empty:
            json.dump(trial, open(p, "w", encoding="utf-8")); skipped += 1; continue
        # events per flow-group code (endpoint-aware)
        ev_by_fg = defaultdict(int)
        for _, r in sub.iterrows():
            reason = str(r.get("reason", "")).lower()
            try:
                c = int(float(r.get("count")))
            except (TypeError, ValueError):
                continue
            if any(k in reason for k in ev_kw):
                ev_by_fg[str(r.get("ctgov_group_code", "") or "")] += c
        # N-match each outcome arm to a flow group whose STARTED count == arm N
        starts = started_n.get(nct, {})
        derived = {}
        for arm in trial["arms"]:
            N = arm.get("N")
            if not N:
                continue
            fg = next((code for code, n in starts.items() if n == N), None)  # exact N match
            if fg is None:                                                    # nearest within 1%
                cand = [(abs(n - N), code) for code, n in starts.items()]
                if cand:
                    d, code = min(cand)
                    if d <= max(1, 0.01 * N):
                        fg = code
            if fg is not None and fg in ev_by_fg and 0 < ev_by_fg[fg] <= N:
                derived[arm["arm_id"]] = int(ev_by_fg[fg])
        # require ALL arms mapped (mixing informed + curve-only arms produces spurious HRs)
        if len(derived) == len(trial["arms"]) and len(set(derived.values())) >= 1:
            for arm in trial["arms"]:
                arm["total_events"] = derived[arm["arm_id"]]
            patched += 1
        else:
            skipped += 1
        json.dump(trial, open(p, "w", encoding="utf-8"))
    print(f"patched total_events into {patched} trials; skipped {skipped} (no/unmappable withdrawals)")


if __name__ == "__main__":
    raise SystemExit(main())
