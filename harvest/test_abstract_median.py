"""Unit tests for the published-median-survival extractor (harvest/abstract_median.py)."""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import abstract_median as M  # noqa: E402


def test_two_arm_medians_with_ci_stripped():
    ab = ("Median overall survival was 8.4 months (95% CI 6.8-10.0) in the tivantinib group and "
          "9.1 months (7.3-10.4) in the placebo group.")
    r = M.extract_medians(ab)
    assert r["medians"] == [8.4, 9.1] and r["n_numbers"] == 2     # CI bounds not counted as medians


def test_middle_dot_decimal_normalised():
    r = M.extract_medians("median overall survival was 8·4 months versus 9·1 months")
    assert r["medians"] == [8.4, 9.1]


def test_mojibake_replacement_char_decimal():
    r = M.extract_medians("median OS was 11�2 months and 10�0 months")
    assert r["medians"] == [11.2, 10.0]


def test_difference_by_is_skipped():
    # "improved ... by 10.7 months" is a DIFFERENCE, not an arm median -> not extracted
    r = M.extract_medians("treatment improved median progression-free survival by 10.7 months overall")
    assert r is None or r["medians"] == []


def test_weeks_converted_to_months():
    r = M.extract_medians("median PFS was 26 weeks versus 13 weeks")
    assert r["medians"][0] == round(26 / 4.345, 2) and r["medians"][1] == round(13 / 4.345, 2)


def test_not_reached_flagged():
    r = M.extract_medians("median overall survival was not reached versus 18.2 months in the control arm")
    assert r["not_reached"] is True


def test_none_when_no_median_survival():
    assert M.extract_medians("The objective response rate was 42% versus 31% (p=0.01).") is None


def test_endpoint_filter_matches_requested_endpoint():
    ab = ("median overall survival was 18.3 months versus 21.1 months; "
          "median progression-free survival was 4.3 months versus 5.2 months.")
    os_r = M.extract_medians(ab, endpoint="OS")
    pfs_r = M.extract_medians(ab, endpoint="PFS")
    assert os_r["medians"] == [18.3, 21.1]      # OS window read for OS
    assert pfs_r["medians"] == [4.3, 5.2]       # PFS window read for PFS (not the OS numbers)


def test_endpoint_filter_returns_none_when_absent():
    ab = "median progression-free survival was 4.3 months versus 5.2 months."
    assert M.extract_medians(ab, endpoint="OS") is None   # no OS median -> nothing (no PFS fallback)
