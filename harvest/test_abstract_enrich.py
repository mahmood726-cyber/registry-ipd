"""Unit tests for the unified abstract-enrichment entry point (harvest/abstract_enrich.py)."""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import abstract_enrich as A  # noqa: E402


def _trial(with_hr=False):
    t = {"nct_id": "NCTX", "arms": [
        {"arm_id": "exp", "label": "Drug", "N": 205, "total_events": None,
         "km_points": [{"t": 1, "S": 0.9}]},
        {"arm_id": "ctl", "label": "Placebo", "N": 97, "total_events": None,
         "km_points": [{"t": 1, "S": 0.8}]},
    ]}
    if with_hr:
        t["hr"] = {"value": 0.50, "ci_low": 0.35, "ci_high": 0.70, "method": "AACT"}
    return t


# single, CI-backed, unambiguous HR + per-arm counts: both levers fire, HR promoted to trial.hr
SINGLE = ("Death occurred in 107 of 205 in the drug group versus 77 of 97 in placebo. "
          "Treatment reduced the risk of death (hazard ratio 0.48, 95% CI 0.35-0.67). "
          "Median overall survival was 11.0 months versus 3.9 months.")


def test_events_and_confident_hr_and_median():
    t = _trial()
    s = A.enrich_from_abstract(t, SINGLE)
    assert s["events"]["patched"] == 2
    assert t["arms"][0]["total_events"] == 107 and t["arms"][1]["total_events"] == 77
    # one CI-backed HR -> promoted to trial.hr
    assert s["hr"]["set_as"] == "trial.hr" and s["hr"]["confident"] is True
    assert t["hr"]["value"] == 0.48 and t["hr"]["source"] == "pubmed_abstract"
    # median attached as a cross-check
    assert s["median"]["attached"] is True and t["median_abstract"]["medians"] == [11.0, 3.9]


def test_hr_not_promoted_when_aact_already_has_one():
    t = _trial(with_hr=True)
    s = A.enrich_from_abstract(t, SINGLE)
    assert t["hr"]["method"] == "AACT" and t["hr"]["value"] == 0.50   # AACT HR untouched
    assert s["hr"]["set_as"] == "trial.hr_abstract" and t["hr_abstract"]["value"] == 0.48


def test_ambiguous_multi_hr_not_promoted():
    # two endpoint HRs -> ambiguous -> attached as cross-check, NOT promoted to trial.hr
    ab = ("Death occurred in 107 of 205 versus 77 of 97. "
          "Progression-free survival favoured treatment (hazard ratio 0.48, 95% CI 0.35-0.67); "
          "overall survival also improved (hazard ratio 0.64, 95% CI 0.40-1.05).")
    t = _trial()
    s = A.enrich_from_abstract(t, ab)
    assert s["hr"]["confident"] is False
    assert "hr" not in t and t["hr_abstract"]["value"] == 0.48     # first HR kept as flagged cross-check


def test_noop_on_empty_abstract():
    t = _trial()
    s = A.enrich_from_abstract(t, "")
    assert s["events"]["patched"] == 0 and "hr" not in t and "median_abstract" not in t


def test_adverse_only_abstract_sets_no_events_but_keeps_hr():
    # RADIANT-4 shape: AE "X of N" counts only (rejected) + a clean HR -> HR promoted, no events
    ab = ("Grade 3 or 4 adverse events included stomatitis (18 of 202 vs 0 of 98). "
          "Everolimus reduced the risk of progression or death (hazard ratio 0.48, 95% CI 0.35-0.67).")
    t = _trial()
    s = A.enrich_from_abstract(t, ab)
    assert s["events"]["patched"] == 0
    assert t["arms"][0]["total_events"] is None
    assert s["hr"]["set_as"] == "trial.hr" and t["hr"]["value"] == 0.48
