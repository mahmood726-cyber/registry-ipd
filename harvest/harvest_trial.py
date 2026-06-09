#!/usr/bin/env python
"""CLI: harvest one trial's survival anchors from AACT into a trial JSON.

Usage:
    python harvest/harvest_trial.py NCT01234567 [-o out.json] [--outcome-id N]

Requires a resolvable AACT snapshot (set AACT_ZIP / AACT_SQLITE / AACT_TSV_DIR, or
place one under ~/AACT/YYYY-MM-DD/). Fails closed with an actionable message otherwise.
"""
from __future__ import annotations
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import harvester as H  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description="Harvest AACT survival anchors -> trial JSON")
    ap.add_argument("nct_id")
    ap.add_argument("-o", "--out", default=None, help="output path (default: <nct>.json)")
    ap.add_argument("--outcome-id", type=int, default=None, help="force a specific AACT outcome id")
    args = ap.parse_args(argv)

    try:
        from aact_kit import resolve_aact_location
        resolve_aact_location()
    except Exception as e:  # noqa: BLE001
        print("ERROR: no AACT snapshot resolvable.\n", file=sys.stderr)
        print(str(e), file=sys.stderr)
        return 2

    try:
        trial = H.harvest_trial(args.nct_id, outcome_id=args.outcome_id)
    except Exception as e:  # noqa: BLE001
        print(f"ERROR harvesting {args.nct_id}: {type(e).__name__}: {e}", file=sys.stderr)
        return 1

    out = args.out or f"{args.nct_id}.json"
    with open(out, "w", encoding="utf-8") as f:
        json.dump(trial, f, indent=2)
    # quick tier hint
    arms = trial["arms"]
    tier = "C"
    if arms and all(len(a["km_points"]) >= 3 and a["N"] is not None for a in arms):
        tier = "A"
    elif arms and all(a.get("median") and a["median"].get("value") is not None and a["N"] is not None for a in arms) and trial.get("hr"):
        tier = "B"
    print(f"wrote {out}  (arms={len(arms)}, likely tier {tier})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
