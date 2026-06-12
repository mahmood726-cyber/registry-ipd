"""Phase 3c step 3 regression: ML-NMR with an effect modifier, time-to-event via RMST (MLNMRPooler).

Asserts the committed Monte-Carlo result: reconstructed-curve trials enter ML-NMR as AD, and ignoring their
reconstruction variance r^2 over-weights them so both the effect-modifier interaction g_B and the
population-adjusted RMST contrast become over-confident (coverage below nominal); propagating r^2 restores
calibrated coverage toward the all-IPD gold. Reads the committed JSON (no advanced-nma-pooling dependency at
test time); skips if the artifact is absent.
"""
import json
import os

HERE = os.path.dirname(__file__)
RESULTS = os.path.join(HERE, "..", "validate", "phase3c_step3_results.json")


def test_ignoring_recon_variance_overconfident_propagating_restores_calibration():
    if not os.path.exists(RESULTS):
        import pytest
        pytest.skip("phase3c_step3_results.json not generated")
    s = json.load(open(RESULTS, encoding="utf-8"))
    gold, naive, honest = s["true"], s["naive"], s["honest"]

    # the all-IPD gold is well-calibrated for both estimands
    assert 0.90 <= gold["g_B_coverage"] <= 0.99
    assert 0.90 <= gold["contrast_coverage"] <= 0.99

    # ignoring r^2 makes the contrast over-confident (coverage materially below nominal and below gold)
    assert naive["contrast_coverage"] < gold["contrast_coverage"] - 0.05
    assert naive["contrast_coverage"] < 0.90

    # ignoring r^2 likewise hurts the effect-modifier interaction coverage
    assert naive["g_B_coverage"] < gold["g_B_coverage"] - 0.05

    # propagating r^2 restores calibration toward the gold level for both estimands
    assert honest["contrast_coverage"] > naive["contrast_coverage"]
    assert abs(honest["contrast_coverage"] - gold["contrast_coverage"]) < abs(naive["contrast_coverage"] - gold["contrast_coverage"])
    assert honest["g_B_coverage"] > naive["g_B_coverage"]
    assert honest["g_B_coverage"] >= 0.88

    # the engine ran every rep for all three variants
    for m in (gold, naive, honest):
        assert m["n_fits"] >= 0.9 * s["config"]["reps"]
