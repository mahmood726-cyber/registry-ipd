#!/usr/bin/env python
"""Extract published median survival(s) from a PubMed abstract (deterministic, bounded).

Independent validation of the reconstruction's TIGHTEST quantity: abstracts routinely report
"median overall survival was 8.4 months ... versus 9.1 months", giving the per-arm medians the curve
reconstruction reproduces to ~3% (see VALIDATION.md). This pulls the (up to two) point medians for the
FIRST median-survival statement, in months, skipping the CI bounds inside parentheses.

Robustness: normalises mojibaked Lancet middle-dot decimals ("8�4" / "8·4" -> 8.4), converts
weeks/years/days to months, and detects "not reached" (common for the better oncology arm). Returns the
point medians only; the caller matches them to reconstructed arm medians by sorted magnitude (so arm
labelling never matters). Honest limit: an abstract may state medians for several endpoints/subgroups --
we take the first median-survival window and flag when more numbers than two arms appear.
"""
from __future__ import annotations
import html
import re

_TO_MONTHS = {"month": 1.0, "months": 1.0, "mo": 1.0, "week": 1 / 4.345, "weeks": 1 / 4.345,
              "wk": 1 / 4.345, "year": 12.0, "years": 12.0, "yr": 12.0, "day": 1 / 30.437,
              "days": 1 / 30.437}
_SURV = r"(?:overall\s+survival|progression[-\s]*free\s+survival|event[-\s]*free\s+survival|" \
        r"disease[-\s]*free\s+survival|recurrence[-\s]*free\s+survival|\bOS\b|\bPFS\b|\bEFS\b|" \
        r"\bDFS\b|\bRFS\b|survival)"
_MEDIAN_WIN = re.compile(r"median\s+" + _SURV + r".{0,220}", re.IGNORECASE)
# per-endpoint median windows -- used to match the abstract median to the reconstructed curve's endpoint
# (avoid scoring an OS-median against a PFS curve, the median analogue of the HR sibling-endpoint guard).
_ENDPOINT_SURV = {
    "OS": r"overall\s+survival|\bOS\b",
    "PFS": r"progression[-\s]*free\s+survival|\bPFS\b",
    "EFS": r"event[-\s]*free\s+survival|\bEFS\b",
    "DFS": r"disease[-\s]*free\s+survival|\bDFS\b",
    "RFS": r"recurrence[-\s]*free\s+survival|relapse[-\s]*free\s+survival|\bRFS\b",
}
_MEDIAN_WIN_BY_EP = {ep: re.compile(r"median\s+(?:" + rx + r").{0,220}", re.IGNORECASE)
                     for ep, rx in _ENDPOINT_SURV.items()}
_NUM_UNIT = re.compile(r"(\d{1,3}(?:\.\d{1,2})?)\s*(months?|mo|weeks?|wk|years?|yr|days?)\b", re.IGNORECASE)
_NOT_REACHED = re.compile(r"not\s+reached|\bNR\b|not\s+estimable|\bNE\b", re.IGNORECASE)


def _normalise(text: str) -> str:
    text = html.unescape(text)
    text = re.sub(r"(\d)[·�•](\d)", r"\1.\2", text)   # middle-dot / mojibake decimal
    return re.sub(r"\s+", " ", text)


def extract_medians(abstract: str, endpoint: str = None):
    """Return {medians:[months,...], not_reached:bool, n_numbers:int, context} for the first
    median-survival statement, or None. Point medians only (parenthetical CI bounds removed).

    If `endpoint` (OS/PFS/EFS/DFS/RFS) is given, only a median statement for THAT endpoint is read, so the
    published median is matched to the reconstructed curve's endpoint (an OS abstract-median is never
    scored against a PFS curve). With endpoint=None, falls back to the first median-survival statement."""
    if not abstract:
        return None
    text = _normalise(abstract)
    rx = _MEDIAN_WIN_BY_EP.get(endpoint) if endpoint else None
    win = (rx.search(text) if rx else None) or (None if endpoint else _MEDIAN_WIN.search(text))
    if not win:
        return None
    window = win.group(0)
    not_reached = bool(_NOT_REACHED.search(window))
    # drop parenthetical/bracket groups (they hold the CIs) before reading point medians
    bare = re.sub(r"[\(\[][^\)\]]*[\)\]]", " ", window)
    meds = []
    for m in _NUM_UNIT.finditer(bare):
        try:
            val = float(m.group(1)) * _TO_MONTHS.get(m.group(2).lower(), 1.0)
        except (TypeError, ValueError):
            continue
        # skip a DIFFERENCE/improvement amount ("improved ... by 10.7 months"), not an arm median
        if re.search(r"\bby\s*$", bare[max(0, m.start() - 5):m.start()], re.IGNORECASE):
            continue
        if 0.1 <= val <= 200:                       # plausible median survival in months
            meds.append(round(val, 2))
    if not meds and not not_reached:
        return None
    return {"medians": meds[:2], "not_reached": not_reached, "n_numbers": len(meds),
            "context": window[:120].strip()}


if __name__ == "__main__":
    import json, sys
    d = json.load(open(sys.argv[1], encoding="utf-8"))
    for a in d.get("articles", []):
        pmid = (a.get("identifiers") or {}).get("pubmed", "?")
        print(pmid, "->", extract_medians(a.get("abstract") or ""))
