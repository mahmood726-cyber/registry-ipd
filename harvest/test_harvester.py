"""Unit tests for the pure AACT parsing logic (no live snapshot required).
Run:  python -m pytest harvest/test_harvester.py -q
"""
import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(__file__))
import harvester as H  # noqa: E402


# --------------------------------------------------------------- time parsing
@pytest.mark.parametrize("text,expected", [
    ("Month 12", 12.0),
    ("6 Months", 6.0),
    ("Year 1", 12.0),
    ("Week 24", pytest.approx(24 * 7 / 30.4375, rel=1e-6)),
    ("Day 28", pytest.approx(28 / 30.4375, rel=1e-6)),
    ("garbage", None),
    (None, None),
])
def test_parse_timepoint(text, expected):
    assert H.parse_timepoint_to_months(text) == expected


def test_parse_timepoint_bare_number_with_default_unit():
    assert H.parse_timepoint_to_months("18", default_unit="months") == 18.0
    assert H.parse_timepoint_to_months("18", default_unit=None) is None


# --------------------------------------------------------------- KM normalisation
def test_orient_does_not_flip_flat_low_survival():
    # genuine high-risk survival curve (flat, low) must NOT be flipped to incidence
    flat_low = [(6, 0.45), (12, 0.46), (24, 0.45)]
    out = H.orient_to_survival(flat_low)
    assert out[0]["S"] == pytest.approx(0.45, abs=1e-6)  # kept as survival, not 0.55
    # single point never flips
    assert H.orient_to_survival([(12, 0.4)])[0]["S"] == pytest.approx(0.4)


def test_orient_flips_clear_cumulative_incidence():
    # clearly rising from ~0 => incidence => S = 1 - value
    inc = [(6, 0.05), (12, 0.20), (24, 0.40)]
    out = H.orient_to_survival(inc)
    assert [round(p["S"], 2) for p in out] == [0.95, 0.80, 0.60]


def test_orient_percent_units():
    # percent declared via units flag => divide by 100 regardless of magnitude
    out = H.orient_to_survival([(6, 95.0), (12, 80.0)], is_percent=True)
    assert out[0]["S"] == pytest.approx(0.95) and out[1]["S"] == pytest.approx(0.80)


def test_km_value_percent_vs_proportion():
    assert H.km_value_to_survival(82, "Kaplan-Meier Estimate", "Percentage") == pytest.approx(0.82)
    assert H.km_value_to_survival(0.82, "Kaplan-Meier Estimate", "Proportion") == pytest.approx(0.82)
    assert H.km_value_to_survival(95.0, "", "") == pytest.approx(0.95)  # >1.5 => percent
    assert H.km_value_to_survival(150, "", "") is None                  # impossible survival
    assert H.km_value_to_survival(None) is None


# --------------------------------------------------------------- count + negation guard
def test_parse_count_negation_guard():
    assert H.parse_count("1,807", context="Randomized") == 1807
    assert H.parse_count("1807", context="Not Randomized") is None   # negation guard
    assert H.parse_count("12", context="Number of Events") == 12
    assert H.parse_count(None) is None


# --------------------------------------------------------------- HR parsing (never invert)
def test_parse_hazard_ratio_basic():
    analyses = pd.DataFrame([{
        "param_type": "Hazard Ratio (HR)", "param_value": "0.74",
        "ci_lower_limit": "0.60", "ci_upper_limit": "0.91", "ci_percent": "95",
        "p_value": "0.004", "method": "Cox Proportional Hazards", "ci_n_sides": "2",
    }])
    hr = H.parse_hazard_ratio(analyses)
    assert hr["value"] == pytest.approx(0.74)
    assert hr["ci_low"] == pytest.approx(0.60)
    assert hr["ci_high"] == pytest.approx(0.91)
    assert hr["one_sided"] is False


def test_parse_hazard_ratio_cox_method_fallback():
    analyses = pd.DataFrame([{
        "param_type": "Log Hazard Ratio", "param_value": "0.88",
        "ci_lower_limit": "0.7", "ci_upper_limit": "1.1", "ci_percent": "95",
        "p_value": "0.3", "method": "Cox", "ci_n_sides": "1",
    }])
    hr = H.parse_hazard_ratio(analyses)
    assert hr is not None and hr["value"] == pytest.approx(0.88) and hr["one_sided"] is True


def test_parse_hazard_ratio_none_when_absent():
    analyses = pd.DataFrame([{"param_type": "Odds Ratio", "param_value": "1.2", "method": "Logistic"}])
    assert H.parse_hazard_ratio(analyses) is None


def _analyses_curve_vs_sibling():
    # curve outcome (10) carries NO HR; sibling survival outcome (20) carries the HR;
    # a non-survival sibling (30) carries an unrelated odds ratio that must NOT be picked.
    return pd.DataFrame([
        {"outcome_id": 10, "param_type": "Survival Probability", "param_value": "0.62", "method": ""},
        {"outcome_id": 20, "param_type": "Hazard Ratio (HR)", "param_value": "0.68",
         "ci_lower_limit": "0.50", "ci_upper_limit": "0.92", "method": "Cox"},
        {"outcome_id": 30, "param_type": "Odds Ratio", "param_value": "1.4", "method": "Logistic"},
    ])


def test_select_trial_hr_falls_back_to_survival_sibling():
    an = _analyses_curve_vs_sibling()
    hr, from_sibling = H.select_trial_hr(an, tte=10, survival_outcome_ids={10, 20})
    assert hr is not None and hr["value"] == pytest.approx(0.68)
    assert from_sibling is True            # flagged: HR came from a sibling outcome


def test_select_trial_hr_prefers_curve_outcome_when_present():
    an = _analyses_curve_vs_sibling().copy()
    an.loc[len(an)] = {"outcome_id": 10, "param_type": "Hazard Ratio (HR)", "param_value": "0.80",
                       "ci_lower_limit": "0.6", "ci_upper_limit": "1.05", "method": "Cox"}
    hr, from_sibling = H.select_trial_hr(an, tte=10, survival_outcome_ids={10, 20})
    assert hr["value"] == pytest.approx(0.80) and from_sibling is False


def test_select_trial_hr_sibling_must_be_survival_typed():
    # only a NON-survival sibling (30) has the HR -> not in survival_outcome_ids -> no fallback
    an = pd.DataFrame([
        {"outcome_id": 10, "param_type": "Survival Probability", "param_value": "0.62", "method": ""},
        {"outcome_id": 30, "param_type": "Hazard Ratio (HR)", "param_value": "0.9", "method": "Cox"},
    ])
    hr, from_sibling = H.select_trial_hr(an, tte=10, survival_outcome_ids={10})
    assert hr is None and from_sibling is False


# --------------------------------------------------------------- arm assembly (Tier A shape)
def _synthetic_tables():
    outcomes = pd.DataFrame([{
        "id": 10, "nct_id": "NCTX", "outcome_type": "Primary",
        "title": "Overall Survival", "param_type": "Kaplan-Meier Estimate", "units": "Months",
        "time_frame": "up to 36 months",
    }])
    groups = pd.DataFrame([
        {"id": 1, "ctgov_group_code": "OG000", "title": "Drug A", "result_type": "Outcome"},
        {"id": 2, "ctgov_group_code": "OG001", "title": "Placebo", "result_type": "Outcome"},
    ])
    counts = pd.DataFrame([
        {"result_group_id": 1, "count": "200"},
        {"result_group_id": 2, "count": "200"},
    ])
    rows = []
    # KM estimates + number-at-risk + events for each group
    for gid, S in [(1, [0.9, 0.75, 0.6, 0.5]), (2, [0.85, 0.65, 0.45, 0.35])]:
        for t, s in zip([6, 12, 18, 24], S):
            rows.append({"outcome_id": 10, "result_group_id": gid, "classification": f"Month {t}",
                         "title": "Overall Survival", "param_type": "Kaplan-Meier Estimate",
                         "param_value": s * 100, "units": "Percentage"})
        for t, n in [(6, 180), (12, 150), (24, 90)]:
            rows.append({"outcome_id": 10, "result_group_id": gid, "classification": f"Month {t}",
                         "title": "Number at Risk", "param_type": "Number at Risk",
                         "param_value": n, "units": "Participants"})
        rows.append({"outcome_id": 10, "result_group_id": gid, "classification": None,
                     "title": "Number of Events", "param_type": "Number of events",
                     "param_value": 110 if gid == 1 else 130, "units": "Participants"})
    measurements = pd.DataFrame(rows)
    return outcomes, measurements, groups, counts


def test_assemble_arms_tier_a():
    outcomes, measurements, groups, counts = _synthetic_tables()
    arms = H.assemble_arms(outcomes, measurements, groups, counts, outcome_id=10)
    assert len(arms) == 2
    drug = next(a for a in arms if a["label"] == "Drug A")
    plac = next(a for a in arms if a["label"] == "Placebo")
    assert drug["role"] == "experimental" and plac["role"] == "comparator"
    assert drug["N"] == 200 and drug["total_events"] == 110
    assert len(drug["km_points"]) == 4 and len(drug["nar_points"]) == 3
    assert drug["km_points"][0] == {"t": 6.0, "S": 0.9}
    assert drug["nar_points"][0] == {"t": 6.0, "n": 180}
    assert drug["follow_up_max"] == 24.0


def test_resolve_hr_direction_never_inverts():
    outcomes, measurements, groups, counts = _synthetic_tables()
    arms = H.assemble_arms(outcomes, measurements, groups, counts, outcome_id=10)
    hr = {"value": 0.74, "ci_low": 0.6, "ci_high": 0.91, "favors_arm_id": None}
    out = H.resolve_hr_direction(hr, arms)
    # value MUST be unchanged (never inverted); direction points to experimental arm
    assert out["value"] == 0.74
    assert out["favors_arm_id"] == next(a["arm_id"] for a in arms if a["role"] == "experimental")


def test_assembled_trial_classifies_tier_a_in_engine(tmp_path):
    """Cross-check: the harvested JSON shape is consumable by the JS engine contract.
    We assert the tiering fields are present and well-formed (engine tested separately)."""
    outcomes, measurements, groups, counts = _synthetic_tables()
    arms = H.assemble_arms(outcomes, measurements, groups, counts, outcome_id=10)
    for a in arms:
        assert len(a["km_points"]) >= 3 and len(a["nar_points"]) >= 2 and a["total_events"] is not None
