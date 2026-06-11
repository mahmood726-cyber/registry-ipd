#!/usr/bin/env python
"""Independent held-out validation: reconstructed HR vs the PUBLISHED HR (from the PubMed abstract).

The production gallery scores the reconstructed HR against the REGISTRY HR -- but the registry HR shares
provenance with the posted curve (same sponsor submission). The trial's PUBLISHED HR, parsed from its
primary paper's abstract, is an INDEPENDENT source. This triangulates three numbers per validation-grade
trial:
    reconstructed HR  (from our pseudo-IPD)
    registry HR       (AACT, the gallery's held-out truth)
    published HR      (PubMed abstract, independent)
and reports (a) reconstructed-vs-published agreement and (b) whether the two held-out sources
(registry vs published) even agree -- a check on the held-out truth itself.

Pipeline: AACT study_references gave nct -> PMID (realipd/vg_pmids.json); here we efetch each abstract
(NCBI E-utilities, cached + rate-limited) and run the deterministic abstract_hr extractor. Honest limits:
the abstract HR is the FIRST-reported HR-with-CI (usually but not always the primary-endpoint, two-arm
effect); multi-HR abstracts are flagged (n_hr_candidates) and reported as a noisy triangulation, never a
gold standard. No patient data; abstracts + HRs are public.

Run from repo root: python harvest/pubmed_validation.py  ->  realipd/pubmed_validation.json
"""
from __future__ import annotations
import json
import os
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.dirname(__file__))
import abstract_hr as A  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(__file__))
RIPD = os.path.join(ROOT, "realipd")
CACHE = os.path.join(ROOT, ".pubmed_cache.json")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
_UA = "registry-ipd-pubmed-validation (research; mailto:research@example.org)"


def _clean_pmid(p):
    p = str(p).strip()
    return p[:-2] if p.endswith(".0") else p


def _get(url, timeout=30, retries=4):
    req = urllib.request.Request(url, headers={"User-Agent": _UA})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception:
            if attempt == retries - 1:
                return None
            time.sleep(2 ** attempt)
    return None


def fetch_abstracts(pmids, cache):
    """efetch abstract XML for any uncached PMIDs (batched 40/req, 3 req/s). Caches pmid -> abstract."""
    todo = [p for p in pmids if p not in cache]
    for i in range(0, len(todo), 40):
        chunk = todo[i:i + 40]
        url = EUTILS + "?" + urllib.parse.urlencode(
            {"db": "pubmed", "id": ",".join(chunk), "retmode": "xml", "rettype": "abstract"})
        data = _get(url)
        time.sleep(0.34)
        got = set()
        if data:
            try:
                root = ET.fromstring(data)
            except ET.ParseError:
                root = None
            if root is not None:
                for art in root.iter("PubmedArticle"):
                    pid = art.findtext(".//MedlineCitation/PMID")
                    texts = [(e.text or "") for e in art.iter("AbstractText")]
                    if pid:
                        cache[pid] = " ".join(t for t in texts if t).strip()
                        got.add(pid)
        for p in chunk:                       # cache misses as empty so we don't refetch forever
            if p not in got:
                cache.setdefault(p, "")
        json.dump(cache, open(CACHE, "w", encoding="utf-8"))
        print(f"  fetched {min(i + 40, len(todo))}/{len(todo)} new abstracts...", flush=True)
    return cache


def _fold(a, b):
    import math
    return round(math.exp(abs(math.log(a) - math.log(b))), 3)


def main():
    pmids_map = json.load(open(os.path.join(RIPD, "vg_pmids.json"), encoding="utf-8"))
    gx = json.load(open(os.path.join(RIPD, "gallery_expanded.json"), encoding="utf-8"))
    by_nct = {r["nct"]: r for r in gx["rows"]}

    cache = json.load(open(CACHE, encoding="utf-8")) if os.path.exists(CACHE) else {}
    # curve endpoint per trial (from the HR backfill) -> match the published HR to the curve's endpoint.
    try:
        curve_ep = json.load(open(os.path.join(RIPD, "validation_hr_backfill.json"),
                                  encoding="utf-8")).get("curve_endpoint", {})
    except FileNotFoundError:
        curve_ep = {}
    want = {nct: _clean_pmid(p) for nct, p in pmids_map.items() if nct in by_nct}
    print(f"scored validation-grade trials with a PMID: {len(want)}")
    cache = fetch_abstracts(sorted(set(want.values())), cache)

    rows, n_pub = [], 0
    for nct, pmid in want.items():
        ab = cache.get(pmid, "")
        ep = curve_ep.get(nct)
        pub = A.extract_hr(ab, endpoint=ep) if ab else None
        r = by_nct[nct]
        recon, reg = r["recon_HR"], r["registry_HR"]
        # high-confidence extraction: a CI present AND few HR candidates (low ambiguity). Endpoint
        # matching (when the curve endpoint is known) is used to PICK the right HR among candidates, not
        # to gate confidence -- a single primary HR is usually the curve's endpoint even if the abstract
        # doesn't repeat the word next to it.
        conf = bool(pub) and pub["has_ci"] and pub["n_hr_candidates"] <= 2
        row = {"nct": nct, "pmid": pmid, "condition": r.get("condition"), "curve_endpoint": ep,
               "recon_HR": recon, "registry_HR": reg,
               "published_HR": pub["value"] if pub else None,
               "published_CI": ([pub["ci_low"], pub["ci_high"]] if pub and pub["has_ci"] else None),
               "n_hr_candidates": pub["n_hr_candidates"] if pub else None,
               "endpoint_matched": pub.get("endpoint_matched") if pub else None,
               "high_confidence": conf}
        if pub:
            n_pub += 1
            row["recon_vs_published_fold"] = _fold(recon, pub["value"])
            if pub["has_ci"]:
                row["recon_in_published_CI"] = bool(pub["ci_low"] <= recon <= pub["ci_high"])
            if reg is not None:
                row["registry_vs_published_fold"] = _fold(reg, pub["value"])
                row["registry_published_agree"] = bool(
                    (reg - 1) * (pub["value"] - 1) >= 0 and _fold(reg, pub["value"]) < 1.25)
        rows.append(row)

    def med(xs):
        if not xs:
            return None
        s = sorted(xs)
        return round(s[len(s) // 2], 3)

    def frac(pred, base):
        b = [r for r in rows if base(r)]
        return f"{sum(1 for r in b if pred(r))}/{len(b)}"

    hc = [r for r in rows if r["high_confidence"]]
    summary = {
        "scored_trials_with_pmid": len(want),
        "published_HR_extracted": n_pub,
        "high_confidence_extractions": len(hc),
        "high_confidence": {
            "recon_vs_published_median_fold": med([r["recon_vs_published_fold"] for r in hc]),
            "recon_within_published_95CI": frac(lambda r: r.get("recon_in_published_CI"),
                                                lambda r: r["high_confidence"]),
            "registry_vs_published_agree": frac(lambda r: r.get("registry_published_agree"),
                                                lambda r: r["high_confidence"]),
        },
        "all_extractions": {
            "recon_vs_published_median_fold": med([r["recon_vs_published_fold"] for r in rows if r.get("recon_vs_published_fold")]),
            "recon_within_published_95CI": frac(lambda r: r.get("recon_in_published_CI"),
                                                lambda r: r.get("recon_in_published_CI") is not None),
        },
        "note": "published HR = primary non-covariate HR-with-CI in the abstract (deterministic abstract_hr "
                "extractor; prognostic/per-unit HRs skipped, parenthetical (HR) + bracket [95% CI] forms "
                "handled; when the curve endpoint is known the HR is matched to that endpoint to PICK among "
                "candidates -- though for single-HR abstracts the endpoint is rarely adjacent, so matching "
                "mainly disambiguates multi-HR abstracts). HIGH-CONFIDENCE = has CI and <=2 HR candidates "
                "(low ambiguity); multi-arm / multi-HR abstracts are excluded from the headline because the "
                "right pairwise HR can't be auto-matched. registry-vs-published agreement is a check on the "
                "held-out truth itself (two INDEPENDENT sources). Abstracts via NCBI E-utilities; cached.",
    }
    json.dump({"summary": summary, "rows": rows}, open(os.path.join(RIPD, "pubmed_validation.json"), "w"),
              indent=1)
    print("\n" + json.dumps(summary, indent=1))
    print(f"\nwrote realipd/pubmed_validation.json ({n_pub} published HRs extracted)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
