# Truth-Recovery Validation — registry-ipd

**Engine:** `src/engine.js` (1206-line UMD module). Re-exported VERBATIM via
`truth-recovery/engine.mjs` (createRequire; zero logic copied or edited). API used:
`reconstruct`, `reconstructEnsemble`, `classifyTier`, and internals `_.coxLogHR`,
`_.kmFromIPD`, `_.medianFromKM`. This is a genuine registry-native pseudo-IPD
reconstructor: tier classification, Guyot inverse-KM, anchor-exact, and a Titman
QP censoring-informed method, with a Cox partial-likelihood (Newton-Raphson) HR.

**Question:** Does the reconstructed pseudo-IPD recover the truth — median, number
of events, and the log hazard ratio + its CI coverage — when fed a registry summary
(KM points + numbers-at-risk + total_events) generated from a KNOWN survival model?

## Method (known-truth DGP)

`dgp-survival.mjs` (standalone, seeded mulberry32) simulates true individual
event/censor times for two arms from Weibull survival (k=1 exponential by default)
with a KNOWN, proportional hazard ratio, then emits exactly the registry summary
the engine consumes: km_points [{t,S}] at 3-month report times, nar_points [{t,n}],
N and total_events per arm. Truth: control median 12 mo, HR 0.65 (true log-HR
-0.431), n=600/arm, independent + administrative censoring at 36 mo.

`harness.mjs` runs the repo's OWN `reconstruct` -> pseudo-IPD -> `coxLogHR`, and
`reconstructEnsemble` for the HR credible interval. Reference truth uses the SAME
Cox estimator on the TRUE individual data. Coverage = fraction of reconstructed
95% HR-CIs (40 seeded replicates) containing the true theoretical log-HR.

## Results (measured, 40-replicate study)

| Quantity | Reconstructed | Truth | Error |
|---|---|---|---|
| Tier / method | A / Titman QP | (rich registry) | correct routing |
| Events, Exp | rel err 0.16% | (registry count) | near-exact |
| Events, Ctrl | rel err 0.12% | (registry count) | near-exact |
| Median (both arms) | abs err 0.15-0.18 mo | sample KM median | excellent |
| **mean HR** | **0.657** | **0.657** (true IPD) / 0.650 (theory) | **log-HR bias 0.0005 vs true** |
| HR-CI coverage | **0.775** | nominal 0.95 | under-coverage |

Single replicate (seed 42): events 390/390 (Exp), 446/445 (Ctrl); median 18.83 vs
18.69 (Exp), 12.95 vs 12.85 (Ctrl); HR 0.691 rec vs 0.698 true vs 0.650 theory.

### Branch + robustness sweep
- **QP (event counts known):** log-HR bias vs true 0.0016 — essentially exact.
- **Curve-only (total_events stripped):** falls back to anchor-exact; HR 0.672 vs
  0.653 true — mild attenuation toward null (bias 0.030), exactly as the engine's
  own code comments predict ("censoring is invisible to the anchors").
- **Weibull k=1.5 (non-exponential, increasing hazard):** rec HR 0.6455 vs true
  0.6450 — robust to PH-preserving non-exponential shapes.

## Verdict

GENUINE METHODS ENGINE — VALIDATED, STRONG RECOVERY.

- Tier/method routing: correct. PASS
- Event-count recovery (QP pins posted count): near-exact (<1%). PASS
- Median recovery: excellent (<0.2 mo). PASS
- **HR point estimate: ~unbiased (log-HR bias 0.0005 vs true IPD).** PASS
- HR-CI coverage: 0.775 (honest negative — interval too narrow). PARTIAL
- Curve-only attenuation toward null when no event count posted. EXPECTED/DOCUMENTED

This is a materially stronger reconstruction than a naive endpoint-placement
digitizer: the censoring-informed Titman QP, combined with a proper Cox partial
likelihood, recovers the HR almost exactly when the registry posts event counts —
consistent with the advanced-stats requirement to verify events/median/HR against
truth rather than claiming exactness.

## Key finding (honest negative)

The HR POINT estimate is excellent, but the ensemble 95% HR-CI covers the true
log-HR only ~78% of the time (nominal 95%) — the credible interval is too narrow.
The ensemble (`reconstructEnsemble` with pinEvents) samples a +/-10% band around the
posted event count plus +/-0.005 KM-anchor rounding and Guyot-vs-anchor model
uncertainty; this under-propagates the true sampling variability of the HR. The
interval understates uncertainty even though the central estimate is trustworthy.

## Recommendation

1. **HR point estimates from rich (Tier-A, event-count) registry trials are
   trustworthy** for pseudo-IPD meta-analysis; near-exact recovery confirmed.
2. **Widen / recalibrate the ensemble HR-CI** so coverage approaches nominal 0.95.
   Candidates: add Cox-model sampling variance (the dominant missing source) to the
   Monte-Carlo band rather than only event-count/anchor jitter, or bootstrap the
   reconstructed IPD. Test A6 is written to fail loudly when coverage is fixed.
3. **Curve-only (no event count) HR carries ~3% null-attenuation** — flag this in
   exported bundles so downstream users do not treat it as equivalent to the QP HR.

## DGP bug found + fixed during validation (additive, in truth-recovery/ only)

The harness initially hung. Root cause was IN MY DGP (not the engine): `buildArm`
destructured `k` from a cfg object that carried `shape` (not `k`), so `k=undefined`
-> `Math.pow(x, 1/undefined)=NaN` -> all simulated times NaN -> the KM step loop
`while(times[idx[j]]===t)` never advanced (`NaN!==NaN`) -> infinite loop. Fixed by
destructuring `{ shape: k }` and adding a fail-closed `Number.isFinite` guard +
anti-stall `if(j===i)j++` in `kmSteps` (a lessons.md-style infinite-loop trap on
NaN-poisoned sorted arrays). Engine source untouched.

## Files
- engine.mjs — verbatim re-export of src/engine.js
- dgp-survival.mjs — seeded known-truth 2-arm Weibull/exponential survival DGP
- harness.mjs — wires repo reconstruct/ensemble/coxLogHR; measures recovery + coverage
- test.mjs — 7 assertions (all pass; A6/A7 lock the honest negatives)
- REPORT.md — this file
