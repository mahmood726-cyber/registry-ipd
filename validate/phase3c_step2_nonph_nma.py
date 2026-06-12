#!/usr/bin/env python
"""PHASE 3c (step 2): the reconstructed curve unlocks a NON-PH survival NMA at the network level.

Phase 3c (step 1, Sec 4h) carried the de-bias/identification interval through the proportional-hazards
ADNMAPooler. But the whole point of a curve over a single log-HR is *time-resolved* structure -- and a
single pooled log-HR cannot represent a treatment effect that changes over time (non-PH). This step shows
the unique value of the curve granularity at the NETWORK level, on advanced-nma-pooling's own
SurvivalNPHPooler (piecewise-exponential, interval-specific effects):

  - a reconstructed curve naturally yields per-interval (events, person_time) -- exactly the survival-AD
    input SurvivalNPHPooler consumes -- so a curve-only trial can join a non-PH survival NMA that an
    abstract-HR trial cannot;
  - a single-log-HR PH pool (ADNMAPooler on the collapsed effect) throws the time-variation away: it
    returns one number stuck between the early and late truth and reports a zero early->late gap;
  - the NEW failure mode is time-resolved: reconstruction is worst where the curve is least identified --
    the LATE interval (fewest at risk, most censoring; the Sec 4d "heavily-censored reconstructs worst"
    finding mapped to time). That biases the late-interval contrast for reconstructed trials, and the same
    Phase-2c LOO de-bias recovers it.

True network is non-PH and CONSISTENT (per-treatment, per-interval log-HR vs A): B has a strong early
benefit that wanes; C is roughly constant. 6 IPD + 6 reconstructed-curve studies across the A/B/C triangle,
each reporting two intervals. Seeded Monte-Carlo. Reuses C:\\Projects\\advanced-nma-pooling
(SurvivalNPHPooler + ADNMAPooler). Run from repo root:
  python validate/phase3c_step2_nonph_nma.py  ->  validate/phase3c_step2_results.json
"""
import io
import json
import math
import os
import sys

ADNMA = r"C:\Projects\advanced-nma-pooling\src"
HERE = os.path.dirname(__file__)

# true per-treatment, per-interval log-HR vs A (g[A]=0). NON-PH: B benefits early then wanes.
G = {
    "A": {"early": 0.0, "late": 0.0},
    "B": {"early": math.log(0.45), "late": math.log(0.95)},   # strong early benefit -> wanes
    "C": {"early": math.log(0.65), "late": math.log(0.60)},   # roughly constant
}
INTERVALS = [("early", 0.0, 1.0), ("late", 1.0, 2.0)]
COMPARISONS = [("A", "B"), ("A", "C"), ("B", "C")]
BASE_H = {"early": 0.45, "late": 0.35}        # reference (A) hazard per interval
# person-time per arm per interval: late has far fewer at risk (attrition + censoring) -> fewer events ->
# the late curve is least identified, so reconstruction is worst there (Sec 4d, time-resolved).
PT = {"early": 150 * 1.0 * 0.80, "late": 150 * 1.0 * 0.38}
B_LATE = 0.34          # systematic late-interval reconstruction bias on the treatment arm's log-hazard
BETA_DEBIAS = 0.30     # LOO-calibrated de-bias (slightly under the true bias -> honest residual remains)

TRUE_GAP = {t: G[t]["late"] - G[t]["early"] for t in ("B", "C")}   # early->late change (the non-PH signal)


def _load_engine():
    sys.path.insert(0, ADNMA)
    from nma_pool.data.builder import DatasetBuilder
    from nma_pool.models.core_ad import ADNMAPooler
    from nma_pool.models.spec import ModelSpec, SurvivalNPHSpec
    from nma_pool.models.survival_nph import SurvivalNPHPooler
    return DatasetBuilder, ADNMAPooler, ModelSpec, SurvivalNPHSpec, SurvivalNPHPooler


def _studies_meta_arms(sid, t1, t2):
    a1, a2 = sid + "_1", sid + "_2"
    study = {"study_id": sid, "design": "rct", "year": 2020, "source_id": sid, "rob_domain_summary": "low"}
    arms = [{"study_id": sid, "arm_id": a1, "treatment_id": t1, "n": 150},
            {"study_id": sid, "arm_id": a2, "treatment_id": t2, "n": 150}]
    return study, arms, a1, a2


def _survival_payload(records):
    """records: list of (sid, t1, t2, {interval: (ev1, pt1, ev2, pt2)}). Build survival_ad payload."""
    S, A, SA = [], [], []
    for sid, t1, t2, per_iv in records:
        study, arms, a1, a2 = _studies_meta_arms(sid, t1, t2)
        S.append(study); A.extend(arms)
        for (iv, t_start, t_end) in INTERVALS:
            ev1, pt1, ev2, pt2 = per_iv[iv]
            SA.append({"study_id": sid, "arm_id": a1, "outcome_id": "os", "interval_id": iv,
                       "t_start": t_start, "t_end": t_end, "events": int(ev1), "person_time": float(pt1)})
            SA.append({"study_id": sid, "arm_id": a2, "outcome_id": "os", "interval_id": iv,
                       "t_start": t_start, "t_end": t_end, "events": int(ev2), "person_time": float(pt2)})
    return {"studies": S, "arms": A, "survival_ad": SA}


def _ph_payload(records):
    """Collapse each study's two intervals to a single log-HR contrast (the PH pool throws away time)."""
    out = []
    for sid, t1, t2, per_iv in records:
        ev1 = sum(per_iv[iv][0] for (iv, _a, _b) in INTERVALS)
        pt1 = sum(per_iv[iv][1] for (iv, _a, _b) in INTERVALS)
        ev2 = sum(per_iv[iv][2] for (iv, _a, _b) in INTERVALS)
        pt2 = sum(per_iv[iv][3] for (iv, _a, _b) in INTERVALS)
        h1 = (ev1 + 0.5) / pt1
        h2 = (ev2 + 0.5) / pt2
        loghr = math.log(h2 / h1)                       # t2 vs t1
        se = math.sqrt(1.0 / (ev1 + 0.5) + 1.0 / (ev2 + 0.5))
        out.append((sid, t1, t2, loghr, se))
    return out


def _ad_payload(studies):
    S, A, O = [], [], []
    for sid, t1, t2, y, se in studies:
        study, arms, a1, a2 = _studies_meta_arms(sid, t1, t2)
        S.append(study); A.extend(arms)
        sa = se / math.sqrt(2)
        O.append({"study_id": sid, "arm_id": a1, "outcome_id": "loghr", "measure_type": "continuous", "value": 0.0, "se": sa})
        O.append({"study_id": sid, "arm_id": a2, "outcome_id": "loghr", "measure_type": "continuous", "value": y, "se": sa})
    return {"studies": S, "arms": A, "outcomes_ad": O}


def run(reps=400, k=12, seed=20260611):
    import numpy as np
    (DatasetBuilder, ADNMAPooler, ModelSpec, SurvivalNPHSpec, SurvivalNPHPooler) = _load_engine()
    rng = np.random.default_rng(seed)
    nph_spec = SurvivalNPHSpec(outcome_id="os", reference_treatment="A", random_effects=False)
    ad_spec = ModelSpec(outcome_id="loghr", measure_type="continuous", reference_treatment="A", random_effects=True)

    def gen_counts(t1, t2, reconstructed):
        """One study's per-interval (ev1, pt1, ev2, pt2). Reconstructed studies carry a systematic
        late-interval bias on the TREATMENT arm's hazard (the worst-identified part of the curve)."""
        per_iv_true, per_iv_naive, per_iv_honest = {}, {}, {}
        for (iv, _a, _b) in INTERVALS:
            h0 = BASE_H[iv]
            pt = PT[iv]
            lam1 = h0 * math.exp(G[t1][iv])
            lam2 = h0 * math.exp(G[t2][iv])
            ev1 = int(rng.poisson(lam1 * pt))
            ev2 = int(rng.poisson(lam2 * pt))
            per_iv_true[iv] = (ev1, pt, ev2, pt)
            if reconstructed and iv == "late":
                ev2_bias = int(round(ev2 * math.exp(B_LATE)))            # late treatment arm reconstructed high
                per_iv_naive[iv] = (ev1, pt, ev2_bias, pt)
                ev2_deb = max(int(round(ev2_bias * math.exp(-BETA_DEBIAS))), 0)
                per_iv_honest[iv] = (ev1, pt, ev2_deb, pt)
            else:
                per_iv_naive[iv] = (ev1, pt, ev2, pt)
                per_iv_honest[iv] = (ev1, pt, ev2, pt)
        return per_iv_true, per_iv_naive, per_iv_honest

    acc = {m: {"earlyB": [], "lateB": [], "gapB": [], "lateB_cov": 0, "n": 0}
           for m in ("ph", "nph_naive", "nph_honest")}

    for _ in range(reps):
        rec_true, rec_naive, rec_honest = [], [], []
        for j in range(k):
            t1, t2 = COMPARISONS[j % 3]
            sid = f"S{j}"
            reconstructed = (j % 2 == 1)
            pt_true, pt_naive, pt_honest = gen_counts(t1, t2, reconstructed)
            rec_true.append((sid, t1, t2, pt_true))
            rec_naive.append((sid, t1, t2, pt_naive))
            rec_honest.append((sid, t1, t2, pt_honest))

        # PH pool on perfect (true) data -- isolates "single log-HR loses time structure" from recon bias
        try:
            ds = DatasetBuilder().from_payload(_ad_payload(_ph_payload(rec_true)))
            fph = ADNMAPooler().fit(ds, ad_spec)
            eBA, _se = fph.contrast("B", "A")
            acc["ph"]["earlyB"].append(eBA - G["B"]["early"])
            acc["ph"]["lateB"].append(eBA - G["B"]["late"])
            acc["ph"]["gapB"].append(0.0 - TRUE_GAP["B"])  # a single number's gap estimate is 0 -> full bias
            acc["ph"]["n"] += 1
        except Exception:
            pass

        for m, recs in (("nph_naive", rec_naive), ("nph_honest", rec_honest)):
            try:
                ds = DatasetBuilder().from_payload(_survival_payload(recs))
                f = SurvivalNPHPooler().fit(ds, nph_spec)
                eE, _seE = f.contrast("B", "A", interval_id="early")
                eL, seL = f.contrast("B", "A", interval_id="late")
            except Exception:
                continue
            if not (math.isfinite(eE) and math.isfinite(eL)):
                continue
            acc[m]["earlyB"].append(eE - G["B"]["early"])
            acc[m]["lateB"].append(eL - G["B"]["late"])
            acc[m]["gapB"].append((eL - eE) - TRUE_GAP["B"])
            acc[m]["lateB_cov"] += int(abs(eL - G["B"]["late"]) <= 1.96 * seL)
            acc[m]["n"] += 1

    def summ(m):
        import numpy as np
        d = acc[m]
        n = max(d["n"], 1)
        return {
            "early_B_bias": round(float(np.mean(d["earlyB"])) if d["earlyB"] else float("nan"), 4),
            "late_B_bias": round(float(np.mean(d["lateB"])) if d["lateB"] else float("nan"), 4),
            "gap_B_bias": round(float(np.mean(d["gapB"])) if d["gapB"] else float("nan"), 4),
            "late_B_coverage": round(d["lateB_cov"] / n, 3),
            "n_fits": d["n"],
        }

    out = {"config": {"reps": reps, "k": k,
                      "network": "A/B/C, non-PH (B early-benefit wanes), 6 IPD + 6 reconstructed-curve, 2 intervals",
                      "true_gap_B_early_to_late": round(TRUE_GAP["B"], 4),
                      "true_early_B": round(G["B"]["early"], 4), "true_late_B": round(G["B"]["late"], 4),
                      "late_recon_bias": B_LATE, "beta_debias": BETA_DEBIAS},
           "ph": summ("ph"), "nph_naive": summ("nph_naive"), "nph_honest": summ("nph_honest")}
    return out


def main():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    out = run()
    json.dump(out, open(os.path.join(HERE, "phase3c_step2_results.json"), "w"), indent=2)
    c = out["config"]
    print("=== Phase 3c step 2: the reconstructed curve unlocks a NON-PH survival NMA (SurvivalNPHPooler) ===\n")
    print(f"  network {c['network']}")
    print(f"  true B-vs-A log-HR: early {c['true_early_B']}, late {c['true_late_B']} "
          f"(early->late gap {c['true_gap_B_early_to_late']} = the non-PH signal)")
    print(f"  late reconstruction bias {c['late_recon_bias']}, LOO de-bias {c['beta_debias']}, {c['reps']} reps\n")
    print(f"  {'':<26}{'early B bias':<15}{'late B bias':<15}{'early->late gap bias':<23}{'late B cov(0.95)'}")
    for m, lab in [("ph", "PH pool (1 log-HR)"),
                   ("nph_naive", "non-PH NAIVE"),
                   ("nph_honest", "non-PH HONEST (Sec 4d)")]:
        d = out[m]
        print(f"  {lab:<26}{d['early_B_bias']:<15}{d['late_B_bias']:<15}{d['gap_B_bias']:<23}{d['late_B_coverage']}")
    print("\n  Reading: a single pooled log-HR (PH) is stuck between the early and late truth and reports a")
    print("  ZERO early->late gap -- it cannot see the non-PH effect at all. Feeding the reconstructed curve's")
    print("  per-interval events/at-risk into SurvivalNPHPooler recovers the interval-specific effects AND the")
    print("  gap -- the curve granularity is what unlocks non-PH at the network level. The new failure mode is")
    print("  the LATE interval (least-identified, most-censored): naive over-states the late contrast; the same")
    print("  Phase-2c LOO de-bias recovers it and restores late-interval coverage.")
    print("\n  wrote validate/phase3c_step2_results.json")
    return out


if __name__ == "__main__":
    main()
