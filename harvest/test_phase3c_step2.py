"""Phase 3c step 2 regression: the reconstructed curve unlocks a non-PH survival NMA (SurvivalNPHPooler).

Asserts the committed Monte-Carlo result: a single-log-HR PH pool cannot represent the time-varying
(non-PH) effect -- it misses the early->late gap entirely -- whereas feeding the reconstructed curve's
per-interval events/at-risk into SurvivalNPHPooler recovers the interval-specific effects and the gap. The
new failure mode is the least-identified late interval; the Phase-2c LOO de-bias recovers the late contrast
and restores its coverage. Reads the committed JSON (no advanced-nma-pooling dependency at test time);
skips if the artifact is absent.
"""
import json
import math
import os

HERE = os.path.dirname(__file__)
RESULTS = os.path.join(HERE, "..", "validate", "phase3c_step2_results.json")


def test_curve_unlocks_nonph_and_debias_recovers_late_interval():
    if not os.path.exists(RESULTS):
        import pytest
        pytest.skip("phase3c_step2_results.json not generated")
    s = json.load(open(RESULTS, encoding="utf-8"))
    ph, naive, honest = s["ph"], s["nph_naive"], s["nph_honest"]
    true_gap = s["config"]["true_gap_B_early_to_late"]

    # the non-PH signal is real and substantial
    assert abs(true_gap) > 0.5

    # a single pooled log-HR (PH) misses the early->late gap almost entirely; the non-PH pool recovers it
    assert abs(ph["gap_B_bias"]) > 0.5, "PH pool should miss the whole gap"
    assert abs(naive["gap_B_bias"]) < abs(ph["gap_B_bias"]), "non-PH recovers the gap PH cannot see"
    assert abs(honest["gap_B_bias"]) < abs(ph["gap_B_bias"])

    # the new failure mode is the late interval: de-bias reduces the late-contrast bias below naive
    assert abs(honest["late_B_bias"]) < abs(naive["late_B_bias"]), "Sec 4d de-bias recovers the late interval"

    # ... and restores late-interval coverage toward nominal, above naive
    assert honest["late_B_coverage"] >= naive["late_B_coverage"]
    assert honest["late_B_coverage"] >= 0.90

    # both non-PH variants recover the early interval well (reconstruction bias is late-only)
    for m in (naive, honest):
        assert abs(m["early_B_bias"]) < 0.10

    # the engine ran (nearly) every rep for all three methods
    for m in (ph, naive, honest):
        assert m["n_fits"] >= 0.9 * s["config"]["reps"]
