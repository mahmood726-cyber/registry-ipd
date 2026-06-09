#!/usr/bin/env python
"""Live orchestration: once the AACT snapshot is downloaded, run the gated phases end-to-end.

Steps (fails closed at each gate):
  1. verify the snapshot zip is complete + valid
  2. point aact-kit at it (AACT_ZIP) and confirm resolution
  3. Phase 0 coverage scan -> Tier A/B/C counts
  4. harvest the first Tier-A trial it finds -> trial JSON
  5. reconstruct it via the JS engine and print tier/badge/method/Wasserstein

Usage: python harvest/run_live.py [--zip PATH]
"""
from __future__ import annotations
import argparse
import json
import os
import subprocess
import sys
import zipfile

HERE = os.path.dirname(__file__)
ROOT = os.path.dirname(HERE)
DEFAULT_ZIP = r"C:\Users\mahmo\AACT\20260601_pipe-delimited-export.zip"
EXPECTED_BYTES = 2452984720


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--zip", default=DEFAULT_ZIP)
    args = ap.parse_args()

    # 1. verify zip
    if not os.path.exists(args.zip):
        print(f"GATE 1 FAIL: snapshot not found at {args.zip}", file=sys.stderr); return 2
    size = os.path.getsize(args.zip)
    print(f"[1] zip size {size} bytes ({size/EXPECTED_BYTES*100:.1f}% of expected)")
    if size < EXPECTED_BYTES:
        print("GATE 1 FAIL: download incomplete — wait for it to finish.", file=sys.stderr); return 2
    try:
        with zipfile.ZipFile(args.zip) as z:
            bad = z.testzip()
            names = z.namelist()
        if bad:
            print(f"GATE 1 FAIL: corrupt entry {bad}", file=sys.stderr); return 2
        print(f"[1] zip OK — {len(names)} entries (e.g. {', '.join(n for n in names[:4])})")
    except zipfile.BadZipFile:
        print("GATE 1 FAIL: not a valid zip yet (still downloading?)", file=sys.stderr); return 2

    # 2. point aact-kit at it
    os.environ["AACT_ZIP"] = args.zip
    sys.path.insert(0, r"C:\Projects\aact-kit\src")
    try:
        from aact_kit import resolve_aact_location
        loc = resolve_aact_location()
        print(f"[2] aact-kit resolved: {loc}")
    except Exception as e:  # noqa: BLE001
        print(f"GATE 2 FAIL: {e}", file=sys.stderr); return 2

    env = dict(os.environ)

    # 3. coverage scan
    print("[3] running coverage scan (this reads large tables; may take a few minutes) ...")
    cov = os.path.join(ROOT, "coverage_report.json")
    r = subprocess.run([sys.executable, os.path.join(HERE, "coverage_scan.py"), "-o", cov],
                       env=env, cwd=ROOT)
    if r.returncode != 0 or not os.path.exists(cov):
        print("GATE 3 FAIL: coverage scan did not complete", file=sys.stderr); return 1
    report = json.load(open(cov, encoding="utf-8"))
    print(f"[3] coverage: A(km-curve)={report['tier_A_km_curve']} B={report['tier_B_medium']} "
          f"C={report['tier_C_sparse']} (universe {report['universe_trials_with_results']}); "
          f"structured NAR rows in all of AACT = {report['structured_number_at_risk_rows']}")
    tierA = report.get("sample_tierA_nct", [])

    # 4. harvest a real Tier-A trial
    if not tierA:
        print("[4] NOTE: zero Tier-A trials found — premise narrows to Tier B (parametric). "
              "Harvesting a Tier-B/HR trial instead is the next step.")
        return 0
    nct = tierA[0]
    print(f"[4] harvesting first Tier-A trial: {nct}")
    out = os.path.join(ROOT, f"{nct}.json")
    r = subprocess.run([sys.executable, os.path.join(HERE, "harvest_trial.py"), nct, "-o", out],
                       env=env, cwd=ROOT)
    if r.returncode != 0 or not os.path.exists(out):
        print("GATE 4 FAIL: harvest did not produce JSON", file=sys.stderr); return 1

    # 5. reconstruct via node
    print(f"[5] reconstructing {nct} via the engine ...")
    node_snippet = (
        "const RIPD=require('./src/engine.js');const fs=require('fs');"
        f"const t=JSON.parse(fs.readFileSync({json.dumps(out)},'utf8'));"
        "const r=RIPD.reconstruct(t);"
        "console.log(JSON.stringify({nct:r.nct_id,tier:r.tier,badge:r.audit.badge,method:r.method,"
        "wasserstein:r.wasserstein_to_anchors,exportable:r.exportable,flags:r.flags},null,2));"
    )
    subprocess.run(["node", "-e", node_snippet], cwd=ROOT)
    print("\nLIVE RUN COMPLETE.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
