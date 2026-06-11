#!/usr/bin/env python
"""Extract per-arm total EVENT COUNTS from a PubMed abstract (deterministic, bounded regex).

WHY THIS EXISTS — the censoring lever, sourced legally.
  registry-IPD's binding accuracy limit is the *unidentified censoring* in a curve-only reconstruction
  (VALIDATION.md "identifiability limit"): without a per-arm total-event count the QP engine
  (`reconstructArmQP`) cannot split deaths from censoring, so large HRs attenuate ~1.5-fold. The sibling
  `kmcurve` project dissolves this with the figure's at-risk / "N (events)" table — but a published
  *figure* is OUT of this project's production data scope (AACT + PubMed abstracts only). The PubMed
  abstract, however, is IN scope and routinely prints the same number:

      "death occurred in 107 of 205 patients in the everolimus group versus 77 of 97 ..."
      "disease progression or death occurred in 51 of 234 ... and 73 of 232 ..."

  This pulls those per-arm "X of N" event fractions deterministically, giving the QP its censoring lever
  WITHOUT any figure — the production-legal analogue of kmcurve's number-at-risk OCR.

Returns the point counts (and the N each is out of, for an N-match sanity check against the AACT arm N).
The caller matches an extracted (events, n) pair to a reconstructed arm by N — so arm labelling never
matters, exactly as `abstract_median.py` matches medians by sorted magnitude.

Honest limits (read every hit with its context):
  - Abstracts report fractions for many things (enrolment, response, adverse events, subgroups). We
    require an EVENT keyword (death/progression/recurrence/relapse/event) in the local window and skip
    NEGATED ("no deaths in ...") and ENROLMENT ("randomly assigned 205 of ...") contexts. This is a
    high-precision / lower-recall heuristic by design — a wrong event count corrupts the QP, so an
    over-inclusion is worse than a miss (cf. the dual-LLM screening lesson). Pair every count with
    `n_fractions` and treat it as triangulation, not gold truth.
  - No network here (pure text in -> dict out); the caller fetches the abstract.
"""
from __future__ import annotations
import html
import re

# "X of N" / "X/N" / "X out of N" — the count of an event out of an arm total.
# leading guard rejects a number glued to a token ("CDK4/6", "PI3K/AKT") or a decimal tail; trailing
# guard rejects a decimal continuation (dot+digit) or a unit letter but allows a sentence period ("118.").
_FRACTION = re.compile(
    r"(?<![\d.A-Za-z])(\d{1,4})\s*(?:of|/|out\s+of)\s*(\d{1,5})(?!\d)(?!\.\d)(?![A-Za-z])",
    re.IGNORECASE,
)

# an EVENT (the thing being counted) must be nearby, or it is not an event fraction.
_EVENT_KW = re.compile(
    r"death|died|dying|fatal|mortalit|\bevents?\b|progress(?:ion|ed|ive|ing)|"
    r"recurren(?:ce|t)|relaps(?:e|ed|ing)|metastas|\bPD\b|disease\s+progression",
    re.IGNORECASE,
)
# negation anywhere in the left window -> "no deaths occurred in 0 of 50" — skip the fraction.
_NEG = re.compile(r"\b(?:no|not|never|without|nor|neither)\b", re.IGNORECASE)
# a comparator joining two arms' fractions ("107 of 205 ... versus 77 of 97 ...") lets the second
# fraction inherit the event context the abstract states only once.
_COMPARATOR = re.compile(r"\b(?:versus|vs\.?|compared\s+(?:with|to)|and)\b|,", re.IGNORECASE)
# enrolment / baseline / response context -> the fraction is patients assigned/female/responding, not events.
_ENROLL = re.compile(
    r"random|enroll|recruit|assign|allocat|includ|eligib|female|male|\baged\b|baseline|"
    r"screen|response|respond|remission|\bORR\b|complete\s+response|partial\s+response",
    re.IGNORECASE,
)
# SAFETY / adverse-event context -> these "X of N had an event" counts are toxicity, NOT the efficacy
# time-to-event we reconstruct. This is the dominant false-positive class (bare "events" matches "adverse
# events") — skip whenever a safety marker is in the window. (lessons.md: generic "events" keyword trap.)
_ADVERSE = re.compile(
    r"adverse|toxicit|side[-\s]effect|adverse\s+reaction|treatment[-\s]emergent|serious\s+adverse|"
    r"\bSAE\b|\bAEs?\b|grade\s*[3-5]|\bsafety\b|hospitali[sz]|discontinu|infection|bleed",
    re.IGNORECASE,
)
# per-endpoint tokens (optional endpoint-scoped extraction, mirrors abstract_median's endpoint guard).
_ENDPOINT = {
    "OS": r"overall\s+survival|\bOS\b|death|died|mortalit",
    "PFS": r"progression[-\s]*free\s+survival|\bPFS\b|progress|disease\s+progression",
    "EFS": r"event[-\s]*free\s+survival|\bEFS\b",
    "DFS": r"disease[-\s]*free\s+survival|\bDFS\b|recurren|relaps",
    "RFS": r"recurrence[-\s]*free\s+survival|relapse[-\s]*free\s+survival|\bRFS\b|recurren|relaps",
}
_ENDPOINT_RX = {ep: re.compile(rx, re.IGNORECASE) for ep, rx in _ENDPOINT.items()}

_WIN = 45  # bounded context window each side of the fraction (no unbounded scans -> no ReDoS)


_BOUND = re.compile(r"[.;]\s")


def _normalise(text: str) -> str:
    text = html.unescape(text)
    text = re.sub(r"(\d)[·•](\d)", r"\1.\2", text)   # middle-dot decimal (rare in counts, harmless)
    return re.sub(r"\s+", " ", text)


def _clause(text: str, s: int, e: int):
    """The bounded sentence/clause (split on '. ' or '; ') containing the fraction at [s, e).
    Abstracts state the event TYPE once per clause ("Serious adverse events ... in 22 of 218 ..."),
    so guards (adverse/enrolment/negation/event-keyword) are evaluated against the whole clause, not a
    fixed character window — which is what lets the safety guard catch an 'adverse' qualifier 58 chars
    upstream. Returns (clause_text, clause_start)."""
    lo = max(0, s - 200)
    bounds = list(_BOUND.finditer(text[lo:s]))
    start = lo + bounds[-1].end() if bounds else lo
    hi = min(len(text), e + 200)
    after = _BOUND.search(text[e:hi])
    end = e + after.start() if after else hi
    return text[start:end], start


def extract_events(abstract: str, endpoint: str = None):
    """Return {events:[e1,e2], ns:[n1,n2], n_fractions:int, context} for the per-arm event fractions
    in the abstract, or None.

    Each (events[i], ns[i]) is an "X of N" event count for one arm. The caller matches the pair to a
    reconstructed arm by N (label-independent). If `endpoint` (OS/PFS/EFS/DFS/RFS) is given, only
    fractions whose local window mentions that endpoint's event are kept, so an OS death count is never
    fed to a PFS curve. With endpoint=None, any in-scope event fraction qualifies.
    """
    if not abstract:
        return None
    text = _normalise(abstract)
    ep_rx = _ENDPOINT_RX.get(endpoint) if endpoint else None
    # collect every plausible "X of N" fraction (count <= N), in order
    fracs = []
    for m in _FRACTION.finditer(text):
        count, n = int(m.group(1)), int(m.group(2))
        if count > n or n == 0:                       # an event count cannot exceed its arm total
            continue
        fracs.append((m.start(), m.end(), count, n, m.group(0)))

    events, ns, contexts = [], [], []
    last_q_end = None                                 # end pos of the last QUALIFIED fraction (anchor)
    for s, e, count, n, g in fracs:
        clause, cstart = _clause(text, s, e)
        pre = text[cstart:s]                          # clause text before the count (for negation)
        if _NEG.search(pre) or _ENROLL.search(clause) or _ADVERSE.search(clause):
            continue                                   # "no deaths" / enrolment-or-response / safety count
        direct = bool(_EVENT_KW.search(clause))
        # companion: a fraction with no event word of its own inherits from the previous qualified
        # fraction when the connecting text is short and carries a comparator (versus / and / ,).
        companion = False
        if not direct and last_q_end is not None and endpoint is None:
            gap = text[last_q_end:s]
            if 0 < len(gap) <= 70 and _COMPARATOR.search(gap):
                companion = True
        if not (direct or companion):
            continue
        if ep_rx is not None:                         # endpoint scope on a TIGHT window (avoid overlap)
            tight = text[max(0, s - 28):s] + " " + g + " " + text[e:e + 28]
            if not ep_rx.search(tight):
                continue
        events.append(count)
        ns.append(n)
        contexts.append(clause.strip()[:140])
        last_q_end = e
    if not events:
        return None
    return {"events": events[:2], "ns": ns[:2], "n_fractions": len(events),
            "context": contexts[0]}


def match_to_arms(extracted: dict, arm_ns):
    """Best-effort: align extracted (events, ns) fractions to arms given the arms' registry N list.
    Returns {arm_index: total_events} for arms whose N matches an extracted fraction's N within 2% (or
    exactly), so the abstract event count is attached to the RIGHT arm regardless of order/labelling.
    Returns {} when the N's do not line up (then the caller should not trust the pairing).
    """
    if not extracted or not extracted.get("events"):
        return {}
    out = {}
    used = set()
    for ai, an in enumerate(arm_ns):
        if an is None:
            continue
        for fi, (e, fn) in enumerate(zip(extracted["events"], extracted["ns"])):
            if fi in used:
                continue
            if fn == an or (an > 0 and abs(fn - an) / an <= 0.02):
                out[ai] = e
                used.add(fi)
                break
    return out


def enrich_trial_events(trial: dict, abstract: str, endpoint: str = None) -> dict:
    """Fill per-arm `total_events` from the PubMed abstract when AACT did not post it — the
    production-legal censoring lever for the QP engine (`reconstructArmQP`).

    Pure (no network): caller supplies the already-fetched abstract. Only patches an arm whose
    `total_events` is missing AND whose registry N matches an extracted "X of N" fraction's N (via
    `match_to_arms`), so the count is attached to the correct arm regardless of label order. Stamps
    `events_source='pubmed_abstract'` on every arm it fills, so provenance is explicit and the AACT
    participant-flow path (`events_source` absent / 'aact_flow') is never silently overwritten.
    Returns {'patched': n_arms_filled, 'extracted': <extract_events result or None>}.
    """
    arms = [a for a in trial.get("arms", []) if a.get("km_points") and a.get("N")]
    # de-dup repeated arms by label (harvested JSON repeats them), preserving order
    uniq, seen = [], set()
    for a in arms:
        if a.get("label") in seen:
            continue
        seen.add(a.get("label"))
        uniq.append(a)
    extracted = extract_events(abstract, endpoint=endpoint)
    if not extracted:
        return {"patched": 0, "extracted": None}
    mapping = match_to_arms(extracted, [int(a["N"]) for a in uniq])
    patched = 0
    for ai, ev in mapping.items():
        arm = uniq[ai]
        if arm.get("total_events") is None:           # never overwrite an AACT-sourced count
            arm["total_events"] = ev
            arm["events_source"] = "pubmed_abstract"
            patched += 1
    return {"patched": patched, "extracted": extracted}


if __name__ == "__main__":
    import json
    import sys
    d = json.load(open(sys.argv[1], encoding="utf-8"))
    for a in d.get("articles", []):
        pmid = (a.get("identifiers") or {}).get("pubmed", "?")
        print(pmid, "->", extract_events(a.get("abstract") or ""))
