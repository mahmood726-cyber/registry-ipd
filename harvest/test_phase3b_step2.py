"""Phase 3b step 2 regression: granularity-mixed survival NMA (advanced-nma-pooling ADNMAPooler).

Asserts the committed Monte-Carlo result: ignoring the reconstructed studies' reconstruction variance
inflates the network heterogeneity tau; propagating it recovers the all-IPD value. Reads the committed
results JSON (no advanced-nma-pooling dependency at test time); skips if the artifact is absent.
"""
import json
import os

HERE = os.path.dirname(__file__)
RESULTS = os.path.join(HERE, "..", "validate", "phase3b_step2_results.json")


def test_network_tau_inflated_by_naive_recovered_by_honest():
    if not os.path.exists(RESULTS):
        import pytest
        pytest.skip("phase3b_step2_results.json not generated")
    s = json.load(open(RESULTS, encoding="utf-8"))
    true_t, naive_t, honest_t = s["true"]["mean_tau"], s["naive"]["mean_tau"], s["honest"]["mean_tau"]
    # ignoring reconstruction variance inflates network heterogeneity
    assert naive_t > honest_t, f"naive tau {naive_t} should exceed honest {honest_t}"
    # honest (propagated) recovers the all-IPD tau more closely than naive
    assert abs(honest_t - true_t) < abs(naive_t - true_t), \
        f"honest {honest_t} should be closer to gold {true_t} than naive {naive_t}"
    # the engine actually ran every rep, and network contrast coverage stays reasonable for all
    for m in ("true", "naive", "honest"):
        assert s[m]["n_fits"] >= 0.9 * s["config"]["reps"]
        for cov in (s[m]["coverage_B"], s[m]["coverage_C"]):
            assert 0.85 <= cov <= 0.99
