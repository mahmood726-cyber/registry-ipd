"""Unit tests for the published-HR-from-abstract extractor (harvest/abstract_hr.py)."""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import abstract_hr as A  # noqa: E402


def test_basic_hr_ci_bracket_form():
    hr = A.extract_hr("reduced risk by 37% (hazard ratio [HR] 0.63, 95% CI 0.48-0.83; p=0.003) versus")
    assert hr["value"] == 0.63 and hr["ci_low"] == 0.48 and hr["ci_high"] == 0.83 and hr["has_ci"]


def test_confidence_interval_to_form():
    hr = A.extract_hr("hazard ratio [HR] 0.77, 95% confidence interval [CI] 0.54 to 1.09; studies = 2")
    assert hr["value"] == 0.77 and hr["ci_low"] == 0.54 and hr["ci_high"] == 1.09


def test_equals_paren_form():
    hr = A.extract_hr("Val; hazard ratio=0.54 (95% confidence interval: 0.33-0.89; log-rank p=0.014)")
    assert hr["value"] == 0.54 and hr["ci_low"] == 0.33 and hr["ci_high"] == 0.89


def test_adjusted_no_ci_label():
    hr = A.extract_hr("compared with placebo (adjusted HR 0.88, 0.75-1.03; p=0.123). In patients")
    assert hr["value"] == 0.88 and hr["ci_low"] == 0.75 and hr["ci_high"] == 1.03


def test_first_hr_is_chosen_and_count_flagged():
    ab = ("primary endpoint HR 0.70 (95% CI 0.55-0.89); secondary HR 0.85 (95% CI 0.60-1.20); "
          "subgroup HR 0.60 (95% CI 0.40-0.90)")
    hr = A.extract_hr(ab)
    assert hr["value"] == 0.70 and hr["n_hr_candidates"] == 3   # first chosen, ambiguity flagged


def test_none_when_no_hr_reported():
    assert A.extract_hr("Median overall survival was 12.3 vs 9.8 months (log-rank p=0.04).") is None


def test_html_entities_and_unicode_dash():
    hr = A.extract_hr("placebo (hazard ratio 0.92; 95% CI 0.69&#x2013;1.22; P&#x2009;=&#x2009;0.276)")
    assert hr["value"] == 0.92 and hr["ci_low"] == 0.69 and hr["ci_high"] == 1.22


def test_implausible_value_rejected():
    # a stray "HR 45.0" (e.g. heart rate 45 bpm) must not be taken as a hazard ratio
    assert A.extract_hr("resting HR 45.0 beats per minute at baseline") is None


def test_covariate_prognostic_hr_skipped():
    # MRD-as-predictor HR is a covariate effect, not the treatment-arm comparison -> must be skipped
    assert A.extract_hr("MRD was an independent adverse predictor (hazard ratio, 3.1; 95% CI, 1.36-7.07).") is None


def test_covariate_skipped_but_treatment_hr_kept():
    ab = ("Treatment reduced events (HR 0.72, 95% CI 0.58-0.90). Age was a predictor "
          "(hazard ratio 1.04 per year, 95% CI 1.01-1.07).")
    hr = A.extract_hr(ab)
    assert hr["value"] == 0.72 and hr["ci_low"] == 0.58        # treatment HR kept, covariate skipped


def test_parenthetical_hr_label_and_bracket_ci():
    # "(HR) = 0.75" parenthetical label + "[95% CI, X to Y]" bracket CI must both parse
    hr = A.extract_hr("overall hazard ratio (HR) = 0.75, 95% CI 0.58-0.96 across arms")
    assert hr["value"] == 0.75 and hr["ci_low"] == 0.58 and hr["ci_high"] == 0.96


def test_bracket_ci_form():
    hr = A.extract_hr("NIVO+IPI versus NIVO monotherapy (hazard ratio, 0.78 [95% CI, 0.67 to 0.91])")
    assert hr["value"] == 0.78 and hr["ci_low"] == 0.67 and hr["ci_high"] == 0.91


def test_endpoint_aware_prefers_matching_endpoint_hr():
    ab = ("Overall survival favored treatment (HR 0.70, 95% CI 0.55-0.89). Progression-free survival "
          "also improved (HR 0.55, 95% CI 0.42-0.72).")
    os_hr = A.extract_hr(ab, endpoint="OS")
    pfs_hr = A.extract_hr(ab, endpoint="PFS")
    assert os_hr["value"] == 0.70 and os_hr["endpoint_matched"] is True
    assert pfs_hr["value"] == 0.55 and pfs_hr["endpoint_matched"] is True


def test_endpoint_aware_flags_when_no_matching_endpoint():
    ab = "Progression-free survival improved (HR 0.55, 95% CI 0.42-0.72)."
    hr = A.extract_hr(ab, endpoint="OS")            # no OS HR present
    assert hr["value"] == 0.55 and hr["endpoint_matched"] is False   # fallback, flagged


def test_endpoint_labeled_hrs_counted_as_ambiguous():
    # 3 endpoint-labeled HRs; only the clean one parses to a value, but ambiguity must be flagged (>2)
    ab = ("(hazard ratio for death, 0.61; 95% CI, 0.49-0.77) and (hazard ratio for disease progression "
          "or death, 0.44; 95% CI, 0.36-0.54) and overall (hazard ratio, 0.54; 95% CI, 0.43 to 0.67).")
    hr = A.extract_hr(ab)
    assert hr["n_hr_candidates"] >= 3        # flagged multi-HR -> excluded from high-confidence tier
