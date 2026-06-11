"""Phase 3c regression: de-bias + identification interval AND consistency on the NMA (ADNMAPooler).

Asserts the committed Monte-Carlo result: per-edge reconstruction bias makes the NAIVE pool flag
spurious inconsistency on a consistent network at a higher rate than the all-IPD gold standard, and the
Phase-2c object (de-bias the identifiable offset + carry the identification half-width as inflated
variance) returns the spurious-inconsistency rate to the gold-standard baseline while keeping the network
contrasts calibrated. Reads the committed JSON (no advanced-nma-pooling dependency at test time); skips if
the artifact is absent.
"""
import json
import os

HERE = os.path.dirname(__file__)
RESULTS = os.path.join(HERE, "..", "validate", "phase3c_results.json")


def test_naive_inflates_spurious_inconsistency_honest_restores_baseline():
    if not os.path.exists(RESULTS):
        import pytest
        pytest.skip("phase3c_results.json not generated")
    s = json.load(open(RESULTS, encoding="utf-8"))
    gold = s["true"]["spurious_inconsistency_rate"]
    naive = s["naive"]["spurious_inconsistency_rate"]
    honest = s["honest"]["spurious_inconsistency_rate"]

    # the homogeneous + consistent control keeps the gold-standard false-positive rate near alpha
    alpha = s["config"]["alpha"]
    assert gold <= max(3 * alpha, 0.15), f"gold spurious-inconsistency {gold} should be near alpha {alpha}"

    # ignoring reconstruction bias manufactures spurious inconsistency above the gold baseline
    assert naive > gold, f"naive {naive} should exceed gold {gold}"

    # the Phase-2c de-bias + identification interval brings the flag rate back toward the gold baseline
    assert honest < naive, f"honest {honest} should fall below naive {naive}"
    assert abs(honest - gold) <= abs(naive - gold), \
        f"honest {honest} should be at least as close to gold {gold} as naive {naive}"

    # and the inflated variance shrinks the inconsistency Q-statistic relative to naive
    assert s["honest"]["mean_q_inconsistency"] < s["naive"]["mean_q_inconsistency"]

    # the engine ran every rep and the network contrasts stay calibrated for all three variants
    for m in ("true", "naive", "honest"):
        assert s[m]["n_fits"] >= 0.9 * s["config"]["reps"]
        for cov in (s[m]["coverage_B"], s[m]["coverage_C"]):
            assert 0.85 <= cov <= 0.99
