"""Phase 3c step 5 regression: native survival-likelihood ML-NMR (piecewise-exponential Poisson).

(1) cross-validates the IRLS Poisson solver against an independent scipy.optimize fit of the same
log-likelihood; (2) asserts the committed Monte-Carlo result: the native survival ML-NMR recovers the
effect-modifier interaction from mixed IPD + reconstructed-curve survival data (unbiased, ~95% coverage),
and a modifier-ignorant pool is badly biased for an off-average target population.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, "..", "validate"))
RESULTS = os.path.join(HERE, "..", "validate", "phase3c_step5_results.json")


def _toy_rows(seed=3):
    rng = np.random.default_rng(seed)
    rows = []
    for j in range(8):
        x = 0.1 + 0.6 * rng.random()
        t2 = "B" if j % 2 == 0 else "C"
        for s in ("e", "l"):
            base = -1.0 if s == "e" else -0.7
            for trt in ("A", t2):
                logh = base + (0.5 + 1.4 * x if trt == "B" else (0.3 if trt == "C" else 0.0))
                pt = 100.0
                rows.append({"study": f"S{j}", "interval": s, "treatment": trt, "x": x,
                             "events": int(rng.poisson(np.exp(logh) * pt)), "person_time": pt})
    return rows


def test_irls_matches_independent_scipy_optimizer():
    try:
        from scipy.optimize import minimize
    except Exception:
        import pytest
        pytest.skip("scipy not available for the cross-validation")
    from survival_mlnmr import build_design, fit_poisson_irls, neg_loglik
    rows = _toy_rows()
    X, y, offset, _, _ = build_design(rows, "A", include_interaction=True)
    b_irls, _, _, converged = fit_poisson_irls(X, y, offset)
    assert converged, "IRLS should converge on the toy dataset"
    res = minimize(lambda b: neg_loglik(b, X, y, offset), np.zeros(X.shape[1]),
                   method="BFGS", options={"maxiter": 5000, "gtol": 1e-9})
    assert float(np.max(np.abs(b_irls - res.x))) < 1e-3, "IRLS must match the independent scipy MLE"


def test_solver_recovers_population_adjusted_contrast():
    from survival_mlnmr import SurvivalMLNMR
    fit = SurvivalMLNMR("A").fit(_toy_rows(seed=11), include_interaction=True)
    # the population-adjusted B-vs-A log-HR at the covariate centroid is recovered within a loose tolerance
    eff, se = fit.treatment_logHR_at("B", 0.4)
    assert abs(eff - (0.5 + 1.4 * 0.4)) < 0.5, f"contrast {eff} should be near truth {0.5 + 1.4 * 0.4}"
    assert se > 0


def test_committed_result_recovers_modifier_and_beats_ignorant_pool():
    if not os.path.exists(RESULTS):
        import pytest
        pytest.skip("phase3c_step5_results.json not generated")
    s = json.load(open(RESULTS, encoding="utf-8"))
    aware, ignorant = s["modifier_aware"], s["modifier_ignorant"]

    # CLAIM 1: the effect modifier is recovered unbiased with ~nominal coverage
    assert abs(aware["gamma_B_bias"]) < 0.25, f"gamma_B should be ~unbiased, got {aware['gamma_B_bias']}"
    assert 0.88 <= aware["gamma_B_coverage"] <= 0.99

    # CLAIM 2: in the off-average target population the modifier-ignorant pool is materially more biased
    assert abs(ignorant["logHR_high_bias"]) > abs(aware["logHR_high_bias"]) + 0.1
    assert aware["logHR_high_coverage"] > ignorant["logHR_high_coverage"]
    assert abs(aware["logHR_high_bias"]) < 0.15, "modifier-aware prediction is ~unbiased at the target population"

    assert s["n_fits"] >= 0.9 * s["config"]["reps"]
