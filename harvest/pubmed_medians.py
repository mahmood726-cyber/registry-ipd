#!/usr/bin/env python
"""Extract published per-arm median survival for the validation-grade trials (independent of the registry).

Reuses the abstracts already cached by pubmed_validation.py. For each validation-grade trial with an
abstract, pulls the (up to two) point medians of the first median-survival statement and writes
realipd/pubmed_medians.json. validate/pubmed_median_validation.js then reconstructs the per-arm medians
and compares them by sorted magnitude (arm labelling never matters) -- an independent check on the
reconstruction's TIGHTEST quantity (median ~3% on the open gold standard).

Run from repo root: python harvest/pubmed_medians.py  ->  realipd/pubmed_medians.json
"""
from __future__ import annotations
import json
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import abstract_median as M  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(__file__))
RIPD = os.path.join(ROOT, "realipd")
CACHE = os.path.join(ROOT, ".pubmed_cache.json")


def _clean(p):
    p = str(p).strip()
    return p[:-2] if p.endswith(".0") else p


def main():
    if not os.path.exists(CACHE):
        print("no .pubmed_cache.json -- run harvest/pubmed_validation.py first", file=sys.stderr)
        return 2
    cache = json.load(open(CACHE, encoding="utf-8"))
    pmids = json.load(open(os.path.join(RIPD, "vg_pmids.json"), encoding="utf-8"))
    # curve endpoint per trial (OS/PFS/...) from the HR backfill -> match the abstract median to the
    # reconstructed curve's endpoint (never score an OS abstract-median against a PFS curve).
    try:
        curve_ep = json.load(open(os.path.join(RIPD, "validation_hr_backfill.json"),
                                  encoding="utf-8")).get("curve_endpoint", {})
    except FileNotFoundError:
        curve_ep = {}

    out = {}
    n_two = 0
    n_skip_unknown_ep = 0
    for nct, pmid in pmids.items():
        ab = cache.get(_clean(pmid), "")
        ep = curve_ep.get(nct)
        if not ep:
            # curve endpoint unclassifiable -> can't guarantee the abstract median is the same endpoint;
            # drop rather than risk an OS-vs-PFS comparison (the median analogue of the HR endpoint guard).
            n_skip_unknown_ep += 1
            continue
        r = M.extract_medians(ab, endpoint=ep) if ab else None
        if not r or not r["medians"]:
            continue
        out[nct] = {"pmid": _clean(pmid), "endpoint": ep, "medians": r["medians"],
                    "n_numbers": r["n_numbers"], "not_reached": r["not_reached"]}
        if len(r["medians"]) == 2 and r["n_numbers"] == 2:
            n_two += 1
    json.dump(out, open(os.path.join(RIPD, "pubmed_medians.json"), "w"), indent=1)
    print(f"extracted published medians for {len(out)} trials "
          f"({n_two} with a clean two-arm pair; {n_skip_unknown_ep} skipped: curve endpoint unknown)")
    print("wrote realipd/pubmed_medians.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
