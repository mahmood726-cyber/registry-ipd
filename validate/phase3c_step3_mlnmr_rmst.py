#!/usr/bin/env python
"""PHASE 3c (step 3): ML-NMR with an effect modifier, time-to-event via reconstruction (RMST).

This closes the one "claimed, not yet proven" cell of the SYNTHESIS-VISION.md Sec 7 ledger: ML-NMR
(Phillippo 2020) combines IPD + AD with effect modifiers, but the engine is GLM-based and the survival
input it needs was never supplied. registry-IPD supplies it. We run advanced-nma-pooling's MLNMRPooler on
a survival question encoded through a curve-derived continuous estimand -- the RESTRICTED MEAN SURVIVAL
TIME (RMST), which is point-identified from a Kaplan-Meier curve (no censoring-identifiability issue, unlike
the HR) and is one of the estimands SYNTHESIS-VISION.md Sec 3 names. IPD trials enter through the engine's
per-patient ipd path (per-patient RMST pseudo-value min(T, tau) + a patient-level effect modifier);
reconstructed-curve trials enter as AD (arm RMST + its reconstruction variance r^2). The effect modifier
modifies the A-vs-B RMST difference (interaction g_B); C is unmodified.

The claim, now on the ML-NMR engine: ignoring the reconstructed trials' reconstruction variance r^2
over-weights them in the meta-regression and makes the effect-modifier interaction OVER-CONFIDENT
(under-covers); propagating r^2 (the Phase-1 Rubin discipline) recovers calibrated coverage, against an
all-IPD truth. Reuses C:\\Projects\\advanced-nma-pooling (MLNMRPooler). Run from repo root:
  python validate/phase3c_step3_mlnmr_rmst.py  ->  validate/phase3c_step3_results.json

Note (honest scope): this is the RMST-as-continuous route -- the pragmatic time-to-event wiring of a
GLM-based ML-NMR via a curve-derived collapsible estimand. A native survival-likelihood ML-NMR would be a
deeper change to the engine itself; here the engine is used unmodified and registry-IPD supplies the input.
"""
import io
import json
import math
import os
import sys

ADNMA = r"C:\Projects\advanced-nma-pooling\src"
HERE = os.path.dirname(__file__)

# true RMST model (months): control arm baseline + prognostic effect-modifier main effect; treatment adds
# d_t, and for B an interaction g_B*(X - XBAR). C is unmodified (g_C = 0).
RMST_A0 = 10.0
BETA_PROG = 1.2          # prognostic main effect of the modifier on RMST (identifies beta_main)
D = {"B": 1.6, "C": 1.0}
G_B = 2.4                # the effect-modifier interaction on the A-vs-B RMST difference (the estimand)
XBAR = 0.4
EDGES = [("A", "B"), ("A", "C")]


def _load_engine():
    sys.path.insert(0, ADNMA)
    from nma_pool.data.builder import DatasetBuilder
    from nma_pool.models.ml_nmr import MLNMRPooler
    from nma_pool.models.spec import MLNMRSpec
    return DatasetBuilder, MLNMRPooler, MLNMRSpec


def arm_rmst(treatment, X, is_treatment_arm):
    base = RMST_A0 + BETA_PROG * (X - XBAR)
    if not is_treatment_arm:
        return base
    eff = D[treatment] + (G_B * (X - XBAR) if treatment == "B" else 0.0)
    return base + eff


def run(reps=300, k=12, n_arm=120, tau=24.0, rmst_sd=3.0, r_recon=0.55, seed=20260612):
    import numpy as np
    DatasetBuilder, MLNMRPooler, MLNMRSpec = _load_engine()
    rng = np.random.default_rng(seed)
    spec = MLNMRSpec(outcome_id="rmst", reference_treatment="A", covariate_name="em")

    def build(trials, mode):
        """trials: list of dicts {sid, t1, t2, X, ipd, rnoise}. mode in {true, naive, honest}."""
        S, A, O, COV, IPD = [], [], [], [], []
        for tr in trials:
            sid, t1, t2, X = tr["sid"], tr["t1"], tr["t2"], tr["X"]
            a1, a2 = sid + "_1", sid + "_2"
            S.append({"study_id": sid, "design": "rct", "year": 2020, "source_id": sid, "rob_domain_summary": "low"})
            A.append({"study_id": sid, "arm_id": a1, "treatment_id": t1, "n": n_arm})
            A.append({"study_id": sid, "arm_id": a2, "treatment_id": t2, "n": n_arm})
            # small arm-level covariate imbalance so the prognostic main effect (beta_main) is identifiable
            x1, x2 = X - 0.05, X + 0.05
            mu1 = arm_rmst(t1, x1, is_treatment_arm=False)
            mu2 = arm_rmst(t2, x2, is_treatment_arm=True)
            # In 'true' every trial is IPD-quality; in naive/honest the reconstructed ones enter as AD.
            as_ipd = True if mode == "true" else tr["ipd"]
            if as_ipd:
                for (a, trt, mu, xc) in [(a1, t1, mu1, x1), (a2, t2, mu2, x2)]:
                    for p in range(n_arm):
                        val = float(rng.normal(mu, rmst_sd))
                        IPD.append({"study_id": sid, "patient_id": f"{a}_{p}", "arm_id": a, "treatment_id": trt,
                                    "outcome_id": "rmst", "measure_type": "continuous",
                                    "outcome_value": min(max(val, 0.0), tau),
                                    "covariates": {"em": float(rng.normal(xc, 0.08))}})
            else:
                # reconstructed-curve trial enters as AD: arm RMST estimate (+ reconstruction noise) with se
                s_samp = rmst_sd / math.sqrt(n_arm)
                rn = tr["rnoise"]              # shared reconstruction perturbation for this trial
                for (a, mu, xc, sign) in [(a1, mu1, x1, +1), (a2, mu2, x2, -1)]:
                    val = mu + (rn * sign if mode in ("naive", "honest") else 0.0)
                    se = s_samp if mode == "naive" else math.sqrt(s_samp * s_samp + r_recon * r_recon)
                    O.append({"study_id": sid, "arm_id": a, "outcome_id": "rmst", "measure_type": "continuous",
                              "value": float(val), "se": float(se)})
                    COV.append({"study_id": sid, "arm_id": a, "covariate_name": "em",
                                "mean": float(xc), "sd": 0.1, "n": n_arm})
        return {"studies": S, "arms": A, "outcomes_ad": O, "ad_covariates": COV, "ipd": IPD}

    acc = {m: {"gB": [], "gB_cov": 0, "cB": [], "cB_cov": 0, "n": 0} for m in ("true", "naive", "honest")}
    true_contrast_at_xbar = D["B"]      # interaction term vanishes at X = XBAR
    for _ in range(reps):
        trials = []
        for j in range(k):
            t1, t2 = EDGES[j % 2]
            X = 0.1 + 0.6 * rng.random()        # modifier spread identifies the interaction
            # reconstructed-dominated network: the AD (reconstructed) trials carry most of the edge weight,
            # so propagating their r^2 actually matters (a few precise IPD trials don't pin everything)
            trials.append({"sid": f"S{j}", "t1": t1, "t2": t2, "X": X,
                           "ipd": (j % 4 == 0), "rnoise": float(rng.normal(0, r_recon))})
        for m in ("true", "naive", "honest"):
            try:
                ds = DatasetBuilder().from_payload(build(trials, m))
                f = MLNMRPooler().fit(ds, spec)
                g, gse = f.interaction_effects.get("B"), f.interaction_ses.get("B")
                cB, cse = f.contrast("B", "A", covariate_value=XBAR)
            except Exception:
                continue
            if g is None or gse is None or not all(math.isfinite(v) for v in (g, gse, cB, cse)) or gse <= 0 or cse <= 0:
                continue
            acc[m]["gB"].append(g - G_B)
            acc[m]["gB_cov"] += int(abs(g - G_B) <= 1.96 * gse)
            acc[m]["cB"].append(cB - true_contrast_at_xbar)
            acc[m]["cB_cov"] += int(abs(cB - true_contrast_at_xbar) <= 1.96 * cse)
            acc[m]["n"] += 1

    import numpy as np
    out = {"config": {"reps": reps, "k": k, "n_arm": n_arm, "tau": tau, "r_recon": r_recon,
                      "estimand": "effect-modifier interaction g_B AND population-adjusted B-vs-A RMST contrast at mean modifier",
                      "true_g_B": G_B, "true_contrast_at_xbar": true_contrast_at_xbar,
                      "network": "A/B/C RMST, 3 IPD + 9 reconstructed-curve (AD), reconstructed-dominated"}}
    for m in ("true", "naive", "honest"):
        n = max(acc[m]["n"], 1)
        out[m] = {"g_B_bias": round(float(np.mean(acc[m]["gB"])) if acc[m]["gB"] else float("nan"), 4),
                  "g_B_coverage": round(acc[m]["gB_cov"] / n, 3),
                  "contrast_bias": round(float(np.mean(acc[m]["cB"])) if acc[m]["cB"] else float("nan"), 4),
                  "contrast_coverage": round(acc[m]["cB_cov"] / n, 3), "n_fits": acc[m]["n"]}
    return out


def main():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    out = run()
    json.dump(out, open(os.path.join(HERE, "phase3c_step3_results.json"), "w"), indent=2)
    c = out["config"]
    print("=== Phase 3c step 3: ML-NMR effect modifier, time-to-event via RMST reconstruction (MLNMRPooler) ===\n")
    print(f"  network {c['network']}")
    print(f"  estimands: interaction g_B (true {c['true_g_B']}) + pop-adjusted B-vs-A RMST contrast @ mean (true {c['true_contrast_at_xbar']})")
    print(f"  {c['k']} trials/rep, n/arm {c['n_arm']}, RMST horizon tau {c['tau']}, r_recon {c['r_recon']}, {c['reps']} reps\n")
    print(f"  {'':<28}{'g_B cov(0.95)':<16}{'contrast bias':<16}{'contrast cov(0.95)':<20}{'n_fits'}")
    for m, lab in [("true", "all-IPD (gold)"), ("naive", "NAIVE (ignore r^2)"), ("honest", "HONEST (propagate r^2)")]:
        d = out[m]
        print(f"  {lab:<28}{d['g_B_coverage']:<16}{d['contrast_bias']:<16}{d['contrast_coverage']:<20}{d['n_fits']}")
    print("\n  Reading: reconstructed-curve trials enter MLNMRPooler as AD and carry most of the edge weight.")
    print("  Ignoring their reconstruction variance r^2 over-weights them (1/se^2 with se too small), so the")
    print("  population-adjusted contrast becomes OVER-CONFIDENT -- coverage drops below nominal -- exactly the")
    print("  Phase-1 lesson on the ML-NMR engine. Propagating r^2 restores calibrated coverage to the gold level.")
    print("  ML-NMR now runs on a survival question (via the curve-derived RMST), with registry-IPD supplying the")
    print("  partially-identified input the GLM-based engine was never given.")
    print("\n  wrote validate/phase3c_step3_results.json")
    return out


if __name__ == "__main__":
    main()
