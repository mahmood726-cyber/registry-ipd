#!/usr/bin/env python
"""Extract a published hazard ratio + 95% CI from a PubMed abstract (deterministic, bounded regex).

Independent held-out truth: the trial's PUBLISHED HR (from its primary paper's abstract) is a source
independent of the registry HR (which shares provenance with the posted curve). This pulls the FIRST
HR-with-CI reported in an abstract -- usually the primary-endpoint effect -- and flags when several HRs
are present (ambiguity), so downstream code can treat multi-HR abstracts cautiously.

Honest limits: abstracts report multiple HRs (primary/secondary/subgroup/adjusted); the first-with-CI
heuristic is right for most two-arm primary-endpoint papers but is NOT a gold standard. Always pair the
extracted value with `n_hr_candidates` and treat agreement as triangulation, not proof. No network here
(pure text in -> dict out); the caller fetches abstracts.
"""
from __future__ import annotations
import html
import re

_DASH = r"[-‐‑‒–—]"          # hyphen + unicode dashes
# trigger ("hazard ratio" / "HR"), optional "(HR)"/"[HR]"/"of"/sep, then the value
_HR = re.compile(
    r"(?:hazard\s+ratio|\bHR\b)\s*(?:[\[(]\s*HR\s*[\])])?\s*(?:of\s+)?[=:,]?\s*"
    r"(\d\.\d{1,3}|\d{1,2}\.\d{1,3})"                  # group 1: HR value (bounded)
    r"(?:"                                             # optional CI block (bounded lookahead window)
    r"[ ,;(\[]{1,4}(?:adjusted\s+)?(?:9[05]\s*%\s*)?(?:CI|confidence\s+interval)?\s*"
    r"(?:[\[(]\s*CI\s*[\])])?\s*[:,]?\s*"
    r"(\d{1,2}\.\d{1,3})\s*(?:" + _DASH + r"|to)\s*(\d{1,2}\.\d{1,3})"   # groups 2,3: CI low/high
    r")?",
    re.IGNORECASE,
)

# looser pattern used ONLY to COUNT distinct HR mentions (for the ambiguity flag) -- it also matches the
# "hazard ratio for <endpoint>, 0.61" labeled form the strict extractor intentionally skips, so an
# abstract that reports several endpoint-labeled HRs is correctly flagged multi-HR (and excluded from the
# high-confidence tier) even though only one of them parses cleanly into a value+CI.
_HR_MENTION = re.compile(r"(?:hazard\s+ratio|\bHR\b)[^.;]{0,45}?(\d{1,2}\.\d{1,3})", re.IGNORECASE)

# context markers that mean an HR is a COVARIATE / prognostic effect (per-unit, per-year, biomarker
# predictor, multivariable covariate) rather than the treatment-arm comparison we want to validate.
_COVARIATE = re.compile(
    r"predict|prognost|independent\s+(?:adverse|risk|predictor|prognostic)|multivariab|multivariate|"
    r"per\s+(?:year|month|unit|mg|kg|point|increase|decrease|\d)|each\s+\d|"
    r"per[-\s]?(?:one|1|10|100)\b|increment",
    re.IGNORECASE,
)


def extract_hr(abstract: str):
    """Return {value, ci_low, ci_high, n_hr_candidates, context} for the first plausible
    HR-with-CI in the abstract (falls back to the first bare HR if none carry a CI), or None."""
    if not abstract:
        return None
    text = re.sub(r"\s+", " ", html.unescape(abstract))
    with_ci, bare = [], []
    for m in _HR.finditer(text):
        try:
            val = float(m.group(1))
        except (TypeError, ValueError):
            continue
        if not (0.05 <= val <= 20):                    # HRs outside this range are almost never real
            continue
        # skip covariate / prognostic HRs (per-unit, biomarker predictor, multivariable) -- the ~45
        # chars before the trigger reveal "an independent adverse predictor (HR 3.1)" vs a trial-arm HR.
        if _COVARIATE.search(text[max(0, m.start() - 45):m.start() + 3]):
            continue
        lo = float(m.group(2)) if m.group(2) else None
        hi = float(m.group(3)) if m.group(3) else None
        rec = {"value": val, "ci_low": lo, "ci_high": hi, "pos": m.start(),
               "context": text[max(0, m.start() - 10):m.start() + 60]}
        if lo is not None and hi is not None and lo < hi and lo <= val <= hi:
            with_ci.append(rec)
        else:
            bare.append(rec)
    pool = with_ci or bare
    if not pool:
        return None
    # honest ambiguity count: distinct non-covariate HR mentions (incl. endpoint-labeled forms the strict
    # pattern skips), so multi-endpoint abstracts are flagged even when only one value parses cleanly.
    n_mentions = 0
    for m in _HR_MENTION.finditer(text):
        try:
            v = float(m.group(1))
        except (TypeError, ValueError):
            continue
        if 0.05 <= v <= 20 and not _COVARIATE.search(text[max(0, m.start() - 45):m.start() + 3]):
            n_mentions += 1
    best = pool[0]                                      # first reported = usually primary endpoint
    return {"value": best["value"], "ci_low": best["ci_low"], "ci_high": best["ci_high"],
            "n_hr_candidates": max(len(with_ci) + len(bare), n_mentions), "has_ci": bool(with_ci),
            "context": best["context"].strip()}


if __name__ == "__main__":                             # smoke: run against a saved metadata JSON
    import json, sys
    path = sys.argv[1]
    d = json.load(open(path, encoding="utf-8"))
    for a in d.get("articles", []):
        ids = a.get("identifiers") or {}
        pmid = ids.get("pubmed") or ids.get("pmid") or "?"
        hr = extract_hr(a.get("abstract") or "")
        print(pmid, "->", hr["value"] if hr else None,
              (f"[{hr['ci_low']},{hr['ci_high']}] n={hr['n_hr_candidates']}" if hr else ""),
              ("| " + hr["context"]) if hr else "")
