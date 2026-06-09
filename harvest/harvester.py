"""
Registry-native survival harvester: ClinicalTrials.gov / AACT -> trial JSON.

The trial JSON is the ONLY thing the offline browser engine ever sees, so provenance
stays 100% registry. There is NO Kaplan-Meier curve and NO patient-level data in AACT;
we extract exact *summary* anchors (KM-estimate points, number-at-risk, median, HR).

This module is split into PURE parsing functions (DataFrame -> dict, fully unit-testable
without a live snapshot) and a thin I/O shell (`harvest_trial`) that calls aact-kit.

Defensive parsing follows the portfolio lessons:
  - negation guard on count extraction ("Not Randomized 1,807" must not become N)
  - percent-vs-proportion guard on KM estimates
  - never invert HR; only resolve/swap arm labels
"""
from __future__ import annotations
import re
import math
from typing import Optional

import pandas as pd

# ----------------------------------------------------------------- time parsing

_TIME_UNIT_TO_MONTHS = {
    "day": 1.0 / 30.4375, "days": 1.0 / 30.4375,
    "week": 7.0 / 30.4375, "weeks": 7.0 / 30.4375,
    "month": 1.0, "months": 1.0,
    "year": 12.0, "years": 12.0,
}
_NEGATION = re.compile(r"\b(not|non|never|excluded|ineligible)\b", re.I)


def parse_timepoint_to_months(text: str, default_unit: Optional[str] = None) -> Optional[float]:
    """Parse 'Month 12', '6 Months', 'Year 1', 'Week 24', 'Day 28', or a bare number
    (interpreted with default_unit) into months. Returns None if unparseable."""
    if text is None:
        return None
    s = str(text).strip().lower()
    if not s:
        return None
    # "<unit> <n>"  e.g. "month 12"
    m = re.search(r"\b(day|week|month|year)s?\b[^\d]*(\d+(?:\.\d+)?)", s)
    if m:
        return float(m.group(2)) * _TIME_UNIT_TO_MONTHS[m.group(1)]
    # "<n> <unit>" e.g. "12 months"
    m = re.search(r"(\d+(?:\.\d+)?)\s*(day|week|month|year)s?\b", s)
    if m:
        return float(m.group(1)) * _TIME_UNIT_TO_MONTHS[m.group(2)]
    # bare number with a default unit
    m = re.fullmatch(r"\d+(?:\.\d+)?", s)
    if m and default_unit:
        u = default_unit.strip().lower()
        if u in _TIME_UNIT_TO_MONTHS:
            return float(s) * _TIME_UNIT_TO_MONTHS[u]
    return None


def km_value_to_survival(param_value, param_type: str = "", units: str = "") -> Optional[float]:
    """Normalise a Kaplan-Meier estimate to a survival probability in [0,1].
    Handles percent (0-100) vs proportion (0-1). Returns None if not interpretable."""
    if param_value is None:
        return None
    try:
        v = float(param_value)
    except (TypeError, ValueError):
        return None
    if math.isnan(v):
        return None
    blob = f"{param_type} {units}".lower()
    is_percent = "%" in blob or "percent" in blob or v > 1.5
    if is_percent:
        v = v / 100.0
    if v < 0 or v > 1.0001:
        return None
    return min(1.0, max(0.0, v))


def parse_count(param_value, context: str = "") -> Optional[int]:
    """Parse an integer count, guarding against negated descriptors in `context`
    (e.g. 'Not Randomized 1807' must NOT be read as the randomized N)."""
    if param_value is None:
        return None
    if context and _NEGATION.search(str(context)):
        return None
    try:
        v = float(str(param_value).replace(",", ""))
    except (TypeError, ValueError):
        return None
    if math.isnan(v) or v < 0:
        return None
    return int(round(v))


# ----------------------------------------------------------------- HR parsing

_HR_PARAM = re.compile(r"hazard\s*ratio|^hr$|\bhr\b", re.I)


def parse_hazard_ratio(analyses: pd.DataFrame) -> Optional[dict]:
    """Extract the hazard ratio + CI from an outcome_analyses subframe.
    NEVER inverts the HR. Returns dict or None."""
    if analyses is None or analyses.empty:
        return None
    rows = analyses[analyses["param_type"].fillna("").str.contains(_HR_PARAM)]
    if rows.empty:
        # some sponsors record method='Cox...' without naming HR in param_type
        rows = analyses[analyses.get("method", pd.Series([""] * len(analyses))).fillna("").str.contains("cox", case=False)]
    if rows.empty:
        return None
    r = rows.iloc[0]
    try:
        value = float(r["param_value"])
    except (TypeError, ValueError, KeyError):
        return None
    if not math.isfinite(value) or value <= 0:
        return None

    def _f(col):
        try:
            x = float(r[col]); return x if math.isfinite(x) else None
        except (TypeError, ValueError, KeyError):
            return None
    ci_sides = r.get("ci_n_sides", None)
    return {
        "value": value,
        "ci_low": _f("ci_lower_limit"),
        "ci_high": _f("ci_upper_limit"),
        "ci_percent": _f("ci_percent"),
        "p_value": _f("p_value"),
        "method": (str(r.get("method", "")) or None),
        "one_sided": (str(ci_sides).strip() in ("1", "1-Sided", "1-sided")),
        "favors_arm_id": None,  # resolved downstream from group roles; never inferred by inverting
    }


# ----------------------------------------------------------------- arm assembly

def _group_role(title: str) -> Optional[str]:
    """Heuristic mapping of a result-group title to experimental/comparator."""
    if not title:
        return None
    t = title.lower()
    if any(k in t for k in ("placebo", "control", "comparator", "standard of care", "soc", "usual care")):
        return "comparator"
    return "experimental"


_SURV_RE = re.compile(r"kaplan|survival|progression-free|event-free|disease-free|"
                      r"probability of|cumulative incidence|proportion", re.I)


def orient_to_survival(raw: list[tuple]) -> list[dict]:
    """Convert a per-arm series of (time_months, raw_value) into survival probabilities S(t).

    DATA-DRIVEN orientation (robust version of the SROC/incidence sign-flip lesson): a true KM
    survival curve starts near 1 and is non-increasing; a cumulative-incidence / "probability of
    [event]" curve starts near 0 and rises. So if the earliest value is low and the series trends
    up, it is an event/incidence curve and S = 1 - value. Percent (>1.5) is divided by 100.
    """
    if not raw:
        return []
    norm = sorted(((t, (v / 100.0 if v > 1.5 else v)) for t, v in raw), key=lambda p: p[0])
    first, last = norm[0][1], norm[-1][1]
    incidence = (first <= 0.5) and (last >= first - 1e-9)
    merged = {}
    for t, vv in norm:
        s = (1.0 - vv) if incidence else vv
        merged[round(t, 4)] = min(1.0, max(0.0, s))
    return [{"t": t, "S": round(merged[t], 6)} for t in sorted(merged)]


def assemble_arms(outcomes: pd.DataFrame, measurements: pd.DataFrame,
                  groups: pd.DataFrame, counts: pd.DataFrame,
                  outcome_id) -> list[dict]:
    """Build per-arm anchors for ONE time-to-event outcome (real-AACT shape).

    groups should already be the result_groups for THIS outcome; counts is outcome_counts
    (result_group_id -> N). Falls back gracefully when those columns are absent (unit tests).
    """
    out_row = outcomes[outcomes["id"] == outcome_id]
    default_unit = (out_row["units"].iloc[0] if not out_row.empty and "units" in out_row else None)
    meas = measurements[measurements["outcome_id"] == outcome_id].copy()

    # N per arm from outcome_counts (scope 'Measure' preferred), keyed by result_group_id
    cmap = {}
    if counts is not None and not counts.empty:
        cc = counts[counts["outcome_id"] == outcome_id] if "outcome_id" in counts else counts
        if "scope" in cc and (cc["scope"].astype(str).str.lower() == "measure").any():
            cc = cc[cc["scope"].astype(str).str.lower() == "measure"]
        for _, r in cc.iterrows():
            cmap[r["result_group_id"]] = parse_count(r.get("count"))

    grp = groups[groups["outcome_id"] == outcome_id] if "outcome_id" in groups else groups
    if grp.empty:
        grp = groups[groups["id"].isin(meas["result_group_id"].unique())]

    arms = []
    for _, g in grp.iterrows():
        gid = g["id"]
        gmeas = meas[meas["result_group_id"] == gid]
        title = str(g.get("title", "") or "")
        arm = {
            "arm_id": str(g.get("ctgov_group_code", gid) or gid),
            "label": title, "role": _group_role(title),
            "N": cmap.get(gid), "total_events": None, "follow_up_max": None,
            "km_points": [], "nar_points": [], "median": None,
        }
        rawKM, nar = [], []
        for _, mrow in gmeas.iterrows():
            ptype = str(mrow.get("param_type", "") or "")
            mtitle = str(mrow.get("title", "") or "")
            cls = mrow.get("classification", None)
            units = str(mrow.get("units", "") or default_unit or "")
            blob = f"{ptype} {mtitle} {units}".lower()
            t_months = parse_timepoint_to_months(cls, default_unit) if cls is not None else None
            if t_months is None:
                t_months = parse_timepoint_to_months(mtitle, default_unit)
            if "at risk" in blob:
                n = parse_count(mrow.get("param_value"), context=mtitle)
                if t_months is not None and n is not None:
                    nar.append({"t": round(t_months, 4), "n": n})
            elif "median" in blob and t_months is None:
                try:
                    mv = float(mrow.get("param_value"))
                    arm["median"] = {"value": mv * _TIME_UNIT_TO_MONTHS.get(units.lower(), 1.0),
                                     "ci_low": None, "ci_high": None}
                except (TypeError, ValueError):
                    pass
            elif ("number of events" in blob or "participants with event" in blob):
                ev = parse_count(mrow.get("param_value"), context=mtitle)
                if ev is not None:
                    arm["total_events"] = ev
            elif _SURV_RE.search(blob) and t_months is not None:
                try:
                    v = float(mrow.get("param_value"))
                    if not math.isnan(v):
                        rawKM.append((t_months, v))
                except (TypeError, ValueError):
                    pass
        arm["km_points"] = orient_to_survival(rawKM)
        nar.sort(key=lambda p: p["t"]); arm["nar_points"] = nar
        if nar:
            arm["follow_up_max"] = max(p["t"] for p in nar)
        elif arm["km_points"]:
            arm["follow_up_max"] = max(p["t"] for p in arm["km_points"])
        arms.append(arm)
    return arms


def resolve_hr_direction(hr: Optional[dict], arms: list[dict]) -> Optional[dict]:
    """Attach favors_arm_id from arm roles WITHOUT ever inverting the HR.
    If exactly one experimental arm exists and HR<1 conventionally favours it, point
    favors_arm_id at the experimental arm; if ambiguous, leave None (engine flags it)."""
    if hr is None:
        return hr
    exp = [a for a in arms if a["role"] == "experimental"]
    if len(exp) == 1:
        hr = dict(hr)
        hr["favors_arm_id"] = exp[0]["arm_id"] if hr["value"] < 1 else \
            (next((a["arm_id"] for a in arms if a["role"] == "comparator"), None))
    return hr


# ----------------------------------------------------------------- I/O shell

def harvest_trial(nct_id: str, location=None, outcome_id=None) -> dict:
    """Pull the survival anchors for one trial from AACT into the trial JSON contract.
    Requires a resolvable AACT snapshot (aact-kit). Raises RuntimeError if none."""
    from aact_kit import load_table  # imported lazily so the pure functions need no snapshot

    where = {"nct_id": nct_id}
    outcomes = load_table("outcomes", location=location, where=where)
    if outcomes.empty:
        raise ValueError(f"No outcomes for {nct_id} (no posted results?)")
    measurements = load_table("outcome_measurements", location=location, where=where)
    analyses = load_table("outcome_analyses", location=location, where=where)
    groups = load_table("result_groups", location=location, where=where)
    if "result_type" in groups:
        groups = groups[groups["result_type"].astype(str).str.lower().str.contains("outcome")]
    try:
        counts = load_table("outcome_counts", location=location, where=where)  # N per arm per outcome
    except Exception:
        counts = pd.DataFrame()

    # pick the time-to-event outcome with the most parseable KM timepoints
    tte = _pick_tte_outcome(outcomes, measurements, analyses) if outcome_id is None else outcome_id
    arms = assemble_arms(outcomes, measurements, groups, counts, tte)
    hr = parse_hazard_ratio(analyses[analyses["outcome_id"] == tte]) if "outcome_id" in analyses else parse_hazard_ratio(analyses)
    hr = resolve_hr_direction(hr, arms)
    return {
        "nct_id": nct_id,
        "source_url": f"https://clinicaltrials.gov/study/{nct_id}",
        "outcome_id": int(tte) if tte is not None else None,
        "time_unit": "months",
        "arms": arms,
        "hr": hr,
    }


def _pick_tte_outcome(outcomes: pd.DataFrame, measurements: pd.DataFrame, analyses: pd.DataFrame):
    """Choose the outcome with the most parseable KM timepoints (best reconstructable curve).
    Falls back to a hazard-ratio outcome, then to the first outcome."""
    if measurements is not None and not measurements.empty:
        blob = (measurements["param_type"].fillna("") + " " + measurements["title"].fillna("")
                + " " + measurements.get("units", pd.Series([""] * len(measurements))).fillna("")).str.lower()
        surv = measurements[blob.str.contains(_SURV_RE)]
        best, best_n = None, 0
        for oid, sub in surv.groupby("outcome_id"):
            tps = set()
            for cls, title in zip(sub["classification"].fillna(""), sub["title"].fillna("")):
                t = parse_timepoint_to_months(cls) or parse_timepoint_to_months(title)
                if t is not None:
                    tps.add(round(t, 3))
            if len(tps) > best_n:
                best, best_n = oid, len(tps)
        if best is not None and best_n >= 3:
            return best
    if "outcome_id" in analyses and not analyses.empty:
        hr_o = set(analyses[analyses["param_type"].fillna("").str.contains(_HR_PARAM)]["outcome_id"])
        cand_hr = outcomes[outcomes["id"].isin(hr_o)]
        if not cand_hr.empty:
            return cand_hr.iloc[0]["id"]
    return outcomes.iloc[0]["id"]
