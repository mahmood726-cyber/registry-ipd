"""Unit tests for the per-arm event-count extractor (harvest/abstract_events.py).

The event count is the QP's censoring lever; a WRONG count corrupts the reconstruction, so these tests
weight precision (no false positives from enrolment/negation/response fractions) over recall.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import abstract_events as E  # noqa: E402


def test_two_arm_x_of_n_deaths():
    ab = ("Death occurred in 107 of 205 patients in the everolimus group versus "
          "77 of 97 patients in the placebo group.")
    r = E.extract_events(ab)
    assert r["events"] == [107, 77] and r["ns"] == [205, 97] and r["n_fractions"] == 2


def test_progression_or_death_fraction():
    ab = ("Disease progression or death occurred in 51 of 234 patients in arm A "
          "and 73 of 232 in arm B.")
    r = E.extract_events(ab)
    assert r["events"] == [51, 73] and r["ns"] == [234, 232]


def test_slash_form():
    r = E.extract_events("There were deaths in 40/120 versus 60/118.")
    assert r["events"] == [40, 60] and r["ns"] == [120, 118]


def test_negation_skipped():
    # "no deaths" -> the 0-of-50 fraction must not be read as an event count
    r = E.extract_events("No deaths occurred in 0 of 50 patients; the study was stopped early.")
    assert r is None


def test_enrolment_fraction_not_events():
    # "randomly assigned" enrolment fraction must be ignored even though it is X of N
    r = E.extract_events("We randomly assigned 205 of 410 screened patients to everolimus.")
    assert r is None


def test_response_fraction_not_events():
    r = E.extract_events("An objective response was seen in 51 of 205 patients versus 30 of 97.")
    assert r is None


def test_count_cannot_exceed_n():
    # a malformed "of" pair where count > N is not an event fraction
    r = E.extract_events("Across 205 of 97 analyses the death rate varied.")
    assert r is None


def test_endpoint_scopes_to_os():
    ab = ("Death occurred in 107 of 205 and 77 of 97 (overall survival); "
          "recurrence occurred in 60 of 205 and 70 of 97 (disease-free survival).")
    os_r = E.extract_events(ab, endpoint="OS")
    dfs_r = E.extract_events(ab, endpoint="DFS")
    assert os_r["events"] == [107, 77]
    assert dfs_r["events"] == [60, 70]


def test_percent_inside_does_not_become_count():
    # "107 of 205 (52%)" -> count is 107, N is 205; the 52 must not be picked up
    r = E.extract_events("Death occurred in 107 of 205 (52%) versus 77 of 97 (79%).")
    assert r["events"] == [107, 77] and r["ns"] == [205, 97]


def test_adverse_events_not_counted():
    # the dominant real-world false positive: "serious adverse events" fractions are toxicity, not efficacy
    ab = ("Serious adverse events occurred in 83 of 218 patients in the pertuzumab group "
          "versus 36 of 110 in the control group.")
    assert E.extract_events(ab) is None


def test_grade34_toxicity_not_counted():
    assert E.extract_events("Grade 3-4 toxicity was seen in 40 of 120 versus 25 of 118.") is None


def test_mortality_rate_x_of_n_is_counted():
    # the true-positive shape seen in the real cache (PMID 38319812)
    ab = "The mortality rate was 20.2% (19 of 94) for sabizabulin versus 45.1% (23 of 51) for placebo."
    r = E.extract_events(ab)
    assert r["events"] == [19, 23] and r["ns"] == [94, 51]


def test_drug_name_slash_not_a_fraction():
    # "CDK4/6 inhibitor" must not parse "4/6" as an event count (real cache false positive, PMID 38861871)
    ab = "The combination of a CDK4/6 inhibitor with ET significantly improved PFS and overall survival."
    assert E.extract_events(ab) is None


def test_none_when_no_event_fraction():
    assert E.extract_events("The trial enrolled adults with advanced disease across 12 centres.") is None


def test_match_to_arms_by_n_regardless_of_order():
    ab = "Death occurred in 77 of 97 in placebo and 107 of 205 in everolimus."
    r = E.extract_events(ab)
    # arms given in the opposite order (everolimus N=205 first) -> matched by N, not appearance
    mapping = E.match_to_arms(r, arm_ns=[205, 97])
    assert mapping == {0: 107, 1: 77}


def test_match_to_arms_returns_empty_when_n_mismatch():
    r = E.extract_events("Death occurred in 107 of 205 versus 77 of 97.")
    assert E.match_to_arms(r, arm_ns=[300, 400]) == {}


def _trial():
    return {"nct_id": "NCTX", "arms": [
        {"arm_id": "exp", "label": "Everolimus", "N": 205, "total_events": None,
         "km_points": [{"t": 1, "S": 0.9}]},
        {"arm_id": "ctl", "label": "Placebo", "N": 97, "total_events": None,
         "km_points": [{"t": 1, "S": 0.8}]},
    ]}


def test_enrich_fills_total_events_with_provenance():
    ab = "Death occurred in 107 of 205 in everolimus versus 77 of 97 in placebo."
    t = _trial()
    res = E.enrich_trial_events(t, ab)
    assert res["patched"] == 2
    assert t["arms"][0]["total_events"] == 107 and t["arms"][0]["events_source"] == "pubmed_abstract"
    assert t["arms"][1]["total_events"] == 77 and t["arms"][1]["events_source"] == "pubmed_abstract"


def test_enrich_never_overwrites_aact_events():
    ab = "Death occurred in 107 of 205 versus 77 of 97."
    t = _trial()
    t["arms"][0]["total_events"] = 99          # pretend AACT participant-flow already supplied this
    res = E.enrich_trial_events(t, ab)
    assert t["arms"][0]["total_events"] == 99 and "events_source" not in t["arms"][0]   # untouched
    assert t["arms"][1]["total_events"] == 77   # the missing arm still gets filled
    assert res["patched"] == 1


def test_enrich_noop_when_no_abstract_events():
    t = _trial()
    res = E.enrich_trial_events(t, "The study enrolled adults across 12 centres.")
    assert res["patched"] == 0 and t["arms"][0]["total_events"] is None
