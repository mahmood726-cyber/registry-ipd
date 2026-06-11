"""Phase 3b regression: weighted_likelihood (spec-collapse-atlas) is the honest within-trial aggregator.

Asserts the committed result of pooling each trial's M reconstruction imputations: weighted_likelihood
keeps honest coverage of the true effect where naive inverse-variance pooling collapses the variance
(false precision). Reads the committed results JSON (no spec-collapse / realipd dependency at test time);
skips cleanly if the artifact is absent.
"""
import json
import os

HERE = os.path.dirname(__file__)
RESULTS = os.path.join(HERE, "..", "validate", "phase3b_weighted_likelihood_results.json")


def test_weighted_likelihood_beats_naive_iv_collapse():
    if not os.path.exists(RESULTS):
        import pytest
        pytest.skip("phase3b_weighted_likelihood_results.json not generated")
    s = json.load(open(RESULTS, encoding="utf-8"))
    # weighted_likelihood keeps honest coverage; naive IV-pooling collapses and (almost) never covers
    assert s["weighted_likelihood_covers_true"] > s["naive_ivre_covers_true"]
    assert s["naive_ivre_covers_true"] <= 1
    # the honest interval is many-fold wider than the collapsed IV interval (the cardinal-sin factor)
    assert s["mean_collapse_ratio"] > 5
    # every cohort's IV interval is narrower than its weighted-likelihood interval
    for c in s["cohorts"]:
        assert c["iv"]["width_log"] < c["wl"]["width_log"]
