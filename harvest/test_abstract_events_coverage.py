"""Unit test for the cohort coverage analyzer (harvest/abstract_events_coverage.py::analyze)."""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import abstract_events_coverage as C  # noqa: E402


def test_analyze_tallies_levers_and_marginal_gain():
    rows = [
        # events that N-match + fill an AACT gap; registry posts an HR (no HR marginal gain)
        {"nct": "NCT_A", "pmid": "1", "curve_endpoint": None, "registry_HR": 0.5},
        # HR only, and the registry posts NO HR -> HR marginal gain
        {"nct": "NCT_B", "pmid": "2", "curve_endpoint": None, "registry_HR": None},
        {"nct": "NCT_C", "pmid": "3", "curve_endpoint": None, "registry_HR": None},   # no usable lever
        {"nct": "NCT_D", "pmid": "9", "curve_endpoint": None, "registry_HR": None},   # no abstract -> ignored
    ]
    abstracts = {
        "1": "Death occurred in 107 of 205 in drug versus 77 of 97 in placebo.",
        "2": "Treatment reduced the risk of death (hazard ratio 0.62, 95% CI 0.50-0.77).",
        "3": "The trial enrolled adults with advanced disease across 12 centres.",
    }
    trials = {
        "NCT_A": {"arms": [{"label": "d", "N": 205, "total_events": None},
                           {"label": "p", "N": 97, "total_events": None}]},
        "NCT_B": {"arms": [{"label": "d", "N": 300, "total_events": None},
                           {"label": "p", "N": 300, "total_events": None}]},
        "NCT_C": {"arms": [{"label": "d", "N": 50, "total_events": None}]},
    }
    out = C.analyze(rows, abstracts, trials)
    s = out["summary"]
    assert s["with_abstract"] == 3                 # NCT_D dropped (no abstract)
    assert s["events_available"] == 1 and s["hr_available"] == 1
    assert s["any_lever"] == 2                      # A (events) + B (HR); C has none
    assert s["n_matched_both_arms"] == 1
    assert s["marginal_gain"] == 1                  # A: abstract fills the AACT event gap, N-matched
    assert out["gains"][0]["nct"] == "NCT_A" and out["gains"][0]["events"] == [107, 77]
    # B: registry posts no HR but the abstract supplies one -> HR marginal gain (A's registry has one)
    assert s["registry_hr_present"] == 1 and s["hr_marginal_gain"] == 1


def test_analyze_no_marginal_gain_when_aact_already_has_events():
    rows = [{"nct": "NCT_A", "pmid": "1", "curve_endpoint": None}]
    abstracts = {"1": "Death occurred in 107 of 205 versus 77 of 97."}
    trials = {"NCT_A": {"arms": [{"label": "d", "N": 205, "total_events": 107},
                                 {"label": "p", "N": 97, "total_events": 77}]}}
    out = C.analyze(rows, abstracts, trials)
    assert out["summary"]["aact_already_has_events"] == 1
    assert out["summary"]["marginal_gain"] == 0    # AACT already had the counts -> no marginal gain
