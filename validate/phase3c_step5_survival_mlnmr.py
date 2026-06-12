#!/usr/bin/env python
"""PHASE 3c (step 5): a NATIVE survival-likelihood ML-NMR with an effect modifier.

Phase 3c step 3 (Sec 4k) ran ML-NMR on a survival question via the RMST-as-Gaussian route. This step uses a
*native survival likelihood* -- a piecewise-exponential (Poisson) ML-NMR (validate/survival_mlnmr.py, the one
synthesis engine the portfolio did not already have) -- combining IPD survival trials (per-patient interval
rows) and reconstructed-curve trials (per-arm per-interval aggregate events/person-time, exactly what a
reconstruction yields). A treatment-by-covariate interaction lets the A-vs-B log-HR vary with an effect
modifier x; C is unmodified.

Two claims:
  1. the native survival ML-NMR RECOVERS the effect-modifier interaction gamma_B (unbiased, ~95% coverage)
     from MIXED IPD + reconstructed-curve survival data -- so a curve-only survival trial joins natively;
  2. a modifier-IGNORANT piecewise-exponential pool (no interaction) returns a single pooled log-HR that is
     materially BIASED when applied to a population whose modifier differs from the corpus average -- the
     value the effect-modifier model adds.

The Poisson IRLS solver is cross-validated against an independent scipy.optimize fit of the same
log-likelihood in harvest/test_phase3c_step5.py. Run from repo root:
  python validate/phase3c_step5_survival_mlnmr.py  ->  validate/phase3c_step5_results.json
"""
import io
import json
import math
import os
import sys

HERE = os.path.dirname(__file__)
sys.path.insert(0, HERE)
from survival_mlnmr import SurvivalMLNMR  # noqa: E402

# true piecewise-exponential model (A vs B), 3 intervals; B has an effect modifier on the log-HR.
INTERVALS = [("i1", 0.0, 1.0), ("i2", 1.0, 2.0), ("i3", 2.0, 4.0)]
BASE_LOGH = {"i1": -1.4, "i2": -1.0, "i3": -0.7}     # reference (A) interval baselines
D_B = 0.5                                            # B main effect (log-HR at x=0)
GAMMA_B = 1.4                                        # effect modifier: log-HR_B(x) = D_B + GAMMA_B * x
XBAR = 0.4                                           # corpus-average modifier
X_HIGH = 0.7                                         # a target population with a higher modifier


def _simulate_patient_rows(rng, study, treatment, x, n, is_B):
    """Simulate n patients' piecewise-exponential survival and expand to (interval) Poisson rows.
    Returns a list of {study, interval, treatment, x, events, person_time} at PATIENT granularity."""
    import numpy as np
    rows = []
    for p in range(n):
        for (iv, t0, t1) in INTERVALS:
            width = t1 - t0
            logh = BASE_LOGH[iv] + (D_B + GAMMA_B * x if is_B else 0.0)
            lam = math.exp(logh)
            # time to event within this interval (exponential), capped at interval width
            te = rng.exponential(1.0 / lam) if lam > 0 else width + 1.0
            if te < width:
                rows.append({"study": study, "interval": iv, "treatment": treatment, "x": x,
                             "events": 1, "person_time": te})
                break                                # patient has the event; not at risk later
            else:
                rows.append({"study": study, "interval": iv, "treatment": treatment, "x": x,
                             "events": 0, "person_time": width})
    return rows


def _aggregate(rows):
    """Collapse patient rows to per-(study,interval,treatment) aggregate (the reconstructed-AD view)."""
    agg = {}
    for r in rows:
        key = (r["study"], r["interval"], r["treatment"])
        if key not in agg:
            agg[key] = {"study": r["study"], "interval": r["interval"], "treatment": r["treatment"],
                        "x": r["x"], "events": 0, "person_time": 0.0}
        agg[key]["events"] += r["events"]
        agg[key]["person_time"] += r["person_time"]
    return list(agg.values())


def run(reps=200, k=12, n_arm=120, seed=20260612):
    import numpy as np
    rng = np.random.default_rng(seed)
    acc = {"gamma_bias": [], "gamma_cov": 0,
           "aware_high_bias": [], "aware_high_cov": 0,
           "ignorant_high_bias": [], "ignorant_high_cov": 0, "n": 0}
    true_logHR_high = D_B + GAMMA_B * X_HIGH
    for _ in range(reps):
        rows = []
        for j in range(k):
            study = f"S{j}"
            x = 0.1 + 0.6 * rng.random()              # modifier spread across the corpus
            reconstructed = (j % 2 == 1)
            patient = (_simulate_patient_rows(rng, study, "A", x, n_arm, is_B=False)
                       + _simulate_patient_rows(rng, study, "B", x, n_arm, is_B=True))
            rows += _aggregate(patient) if reconstructed else patient
        try:
            aware = SurvivalMLNMR("A").fit(rows, include_interaction=True)
            ignorant = SurvivalMLNMR("A").fit(rows, include_interaction=False)
            g, gse = aware.gamma("B")
            ah, ahse = aware.treatment_logHR_at("B", X_HIGH)
            ih, ihse = ignorant.treatment_logHR_at("B", X_HIGH)
        except Exception:
            continue
        if not all(math.isfinite(v) for v in (g, gse, ah, ahse, ih, ihse)):
            continue
        acc["gamma_bias"].append(g - GAMMA_B)
        acc["gamma_cov"] += int(abs(g - GAMMA_B) <= 1.96 * gse)
        acc["aware_high_bias"].append(ah - true_logHR_high)
        acc["aware_high_cov"] += int(abs(ah - true_logHR_high) <= 1.96 * ahse)
        acc["ignorant_high_bias"].append(ih - true_logHR_high)
        acc["ignorant_high_cov"] += int(abs(ih - true_logHR_high) <= 1.96 * ihse)
        acc["n"] += 1

    n = max(acc["n"], 1)
    mean = lambda a: round(float(np.mean(a)), 4) if a else float("nan")
    return {
        "config": {"reps": reps, "k": k, "n_arm": n_arm,
                   "network": "A vs B piecewise-exponential, 6 IPD (patient rows) + 6 reconstructed-curve (aggregate)",
                   "true_gamma_B": GAMMA_B, "true_d_B": D_B, "xbar": XBAR, "x_high": X_HIGH,
                   "true_logHR_B_at_high": round(true_logHR_high, 4),
                   "engine": "validate/survival_mlnmr.py (piecewise-exponential Poisson ML-NMR, native survival likelihood)"},
        "modifier_aware": {
            "gamma_B_bias": mean(acc["gamma_bias"]), "gamma_B_coverage": round(acc["gamma_cov"] / n, 3),
            "logHR_high_bias": mean(acc["aware_high_bias"]), "logHR_high_coverage": round(acc["aware_high_cov"] / n, 3)},
        "modifier_ignorant": {
            "logHR_high_bias": mean(acc["ignorant_high_bias"]), "logHR_high_coverage": round(acc["ignorant_high_cov"] / n, 3)},
        "n_fits": acc["n"],
    }


def main():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    out = run()
    json.dump(out, open(os.path.join(HERE, "phase3c_step5_results.json"), "w"), indent=2)
    c, aw, ig = out["config"], out["modifier_aware"], out["modifier_ignorant"]
    print("=== Phase 3c step 5: native survival-likelihood ML-NMR with an effect modifier ===\n")
    print(f"  engine: {c['engine']}")
    print(f"  network {c['network']}")
    print(f"  true gamma_B {c['true_gamma_B']} (log-HR_B(x) = {c['true_d_B']} + {c['true_gamma_B']}*x); "
          f"target population x_high={c['x_high']} (true log-HR_B {c['true_logHR_B_at_high']}), {c['reps']} reps\n")
    print("  CLAIM 1 -- the native survival ML-NMR recovers the effect modifier from mixed IPD + reconstructed data:")
    print(f"    gamma_B bias {aw['gamma_B_bias']}   coverage {aw['gamma_B_coverage']}")
    print("  CLAIM 2 -- log-HR_B in the high-modifier target population (true "
          f"{c['true_logHR_B_at_high']}):")
    print(f"    modifier-AWARE   bias {aw['logHR_high_bias']:<9} coverage {aw['logHR_high_coverage']}")
    print(f"    modifier-IGNORANT bias {ig['logHR_high_bias']:<9} coverage {ig['logHR_high_coverage']}")
    print(f"  ({out['n_fits']} fits)\n")
    print("  Reading: a piecewise-exponential Poisson ML-NMR -- a NATIVE survival likelihood, not the RMST-as-")
    print("  Gaussian route of step 3 -- recovers the treatment-by-covariate effect modifier (unbiased, ~95%")
    print("  coverage) from a mix of IPD survival trials and reconstructed-curve trials that contribute only")
    print("  per-interval events/at-risk. A modifier-IGNORANT pool reports a single log-HR that is badly biased")
    print("  for a population whose modifier differs from the corpus average -- the value the effect modifier")
    print("  adds. registry-IPD supplies the reconstructed survival input the engine runs on.")
    print("\n  wrote validate/phase3c_step5_results.json")
    return out


if __name__ == "__main__":
    main()
