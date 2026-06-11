#!/usr/bin/env python
"""PHASE 3b (step 2): a granularity-mixed survival NETWORK meta-analysis on the lab's own ADNMAPooler.

The capstone of SYNTHESIS-VISION.md. A 3-treatment survival network (A/B/C, log-HR scale) is pooled with
advanced-nma-pooling's ADNMAPooler (contrast-based random-effects NMA). Some studies are full IPD
(identified point); some are reconstructed-from-curve, carrying the reconstruction variance r^2 the earlier
phases established. We compare pooling the reconstructed studies with their reconstruction variance IGNORED
(naive -- treat pseudo-IPD as exact) vs PROPAGATED (honest), against an all-IPD truth, over a seeded
Monte-Carlo. The claim, now at the network level: ignoring r^2 mis-reads reconstruction noise as between-
study heterogeneity (tau inflated); propagating it recovers tau and keeps the network contrasts calibrated.

Reuses C:\\Projects\\advanced-nma-pooling (ADNMAPooler) -- the engine already exists; registry-IPD supplies
the partially-identified AD trials it was never given. Run from repo root:
  python validate/phase3b_step2_survival_nma.py  ->  validate/phase3b_step2_results.json
"""
import io
import json
import math
import os
import sys

ADNMA = r"C:\Projects\advanced-nma-pooling\src"
HERE = os.path.dirname(__file__)


def _load_engine():
    sys.path.insert(0, ADNMA)
    from nma_pool.data.builder import DatasetBuilder
    from nma_pool.models.core_ad import ADNMAPooler
    from nma_pool.models.spec import ModelSpec
    return DatasetBuilder, ADNMAPooler, ModelSpec


# true network (log-HR vs A): B better, C best; consistent (d_CB = d_C - d_B)
D = {"B": math.log(0.7), "C": math.log(0.5)}
COMPARISONS = [("A", "B"), ("A", "C"), ("B", "C")]          # the three edges of the triangle
TRUE_CONTRAST = {("A", "B"): D["B"], ("A", "C"): D["C"], ("B", "C"): D["C"] - D["B"]}


def _payload(studies):
    """studies: list of (sid, t1, t2, y, se). Encode log-HR y as a continuous contrast: arms (t1:0, t2:y),
    each arm se = se/sqrt(2) so the within-study contrast variance is se^2."""
    S, A, O = [], [], []
    for sid, t1, t2, y, se in studies:
        S.append({"study_id": sid, "design": "rct", "year": 2020, "source_id": sid, "rob_domain_summary": "low"})
        a1, a2 = sid + "_1", sid + "_2"
        A.append({"study_id": sid, "arm_id": a1, "treatment_id": t1, "n": 100})
        A.append({"study_id": sid, "arm_id": a2, "treatment_id": t2, "n": 100})
        sa = se / math.sqrt(2)
        O.append({"study_id": sid, "arm_id": a1, "outcome_id": "loghr", "measure_type": "continuous", "value": 0.0, "se": sa})
        O.append({"study_id": sid, "arm_id": a2, "outcome_id": "loghr", "measure_type": "continuous", "value": y, "se": sa})
    return {"studies": S, "arms": A, "outcomes_ad": O}


def run(reps=300, k=12, tau_true=0.1, s_mean=0.15, r_recon=0.25, seed=20260611):
    import numpy as np
    DatasetBuilder, ADNMAPooler, ModelSpec = _load_engine()
    rng = np.random.default_rng(seed)
    spec = ModelSpec(outcome_id="loghr", measure_type="continuous", reference_treatment="A", random_effects=True)

    def fit(studies):
        ds = DatasetBuilder().from_payload(_payload(studies))
        f = ADNMAPooler().fit(ds, spec)
        return f

    acc = {m: {"tau": 0.0, "covB": 0, "covC": 0, "n": 0} for m in ("true", "naive", "honest")}
    for _ in range(reps):
        base = []  # per study: (sid, t1, t2, theta, s, ipd?, sampling_noise, recon_noise)
        for j in range(k):
            t1, t2 = COMPARISONS[j % 3]
            theta = TRUE_CONTRAST[(t1, t2)] + rng.normal(0, math.sqrt(tau_true))
            s = s_mean * (0.7 + 0.6 * rng.random())
            base.append((f"S{j}", t1, t2, theta, s, (j % 2 == 0),
                         rng.normal(0, s), rng.normal(0, r_recon)))
        builds = {
            "true":   [(sid, t1, t2, th + sn, s) for (sid, t1, t2, th, s, ipd, sn, rn) in base],
            "naive":  [(sid, t1, t2, th + sn + (0 if ipd else rn), s) for (sid, t1, t2, th, s, ipd, sn, rn) in base],
            "honest": [(sid, t1, t2, th + sn + (0 if ipd else rn), s if ipd else math.sqrt(s * s + r_recon * r_recon))
                       for (sid, t1, t2, th, s, ipd, sn, rn) in base],
        }
        for m, studies in builds.items():
            try:
                f = fit(studies)
            except Exception:
                continue
            eB, seB = f.contrast("B", "A")
            eC, seC = f.contrast("C", "A")
            acc[m]["tau"] += f.tau
            acc[m]["covB"] += int(abs(eB - D["B"]) <= 1.96 * seB)
            acc[m]["covC"] += int(abs(eC - D["C"]) <= 1.96 * seC)
            acc[m]["n"] += 1
    out = {"config": {"reps": reps, "k": k, "tau_true_var": tau_true,
                      "tau_true_sd": round(math.sqrt(tau_true), 3), "r_recon": r_recon,
                      "network": "A/B/C log-HR, 6 IPD + 6 reconstructed-curve"}}
    for m in ("true", "naive", "honest"):
        n = max(acc[m]["n"], 1)
        out[m] = {"mean_tau": round(acc[m]["tau"] / n, 4),
                  "coverage_B": round(acc[m]["covB"] / n, 3),
                  "coverage_C": round(acc[m]["covC"] / n, 3), "n_fits": acc[m]["n"]}
    return out


def main():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    out = run()
    json.dump(out, open(os.path.join(HERE, "phase3b_step2_results.json"), "w"), indent=2)
    print("=== Phase 3b step 2: granularity-mixed survival NMA (advanced-nma-pooling ADNMAPooler) ===\n")
    sd = out['config']['tau_true_sd']
    print(f"  network A/B/C, {out['config']['k']} studies (half IPD, half reconstructed-curve), "
          f"true heterogeneity SD {sd}, r_recon {out['config']['r_recon']}, {out['config']['reps']} reps\n")
    print(f"                         mean tau-SD (true {sd})   coverage B(0.95)   coverage C(0.95)")
    for m, lab in [("true", "all-IPD (gold)"), ("naive", "NAIVE (ignore recon var)"), ("honest", "HONEST (propagate r^2)")]:
        d = out[m]
        print(f"  {lab:<26} {d['mean_tau']:<22} {d['coverage_B']:<18} {d['coverage_C']}")
    print("\n  Reading: at the NETWORK level, ignoring reconstruction variance inflates tau (reconstruction")
    print("  noise read as between-study heterogeneity); propagating r^2 recovers it. The reconstructed AD")
    print("  trials join the NMA honestly -- the engine is the lab's ADNMAPooler, the new input is the")
    print("  partially-identified survival trial registry-IPD supplies.")
    print("\n  wrote validate/phase3b_step2_results.json")
    return out


if __name__ == "__main__":
    main()
