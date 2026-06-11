#!/usr/bin/env node
/*
 * PHASE 3: granularity-mixed evidence synthesis with partial-identification-aware intervals.
 * =====================================================================================================
 *
 * The capstone of SYNTHESIS-VISION.md: one pooled survival contrast over a corpus whose trials sit at
 * DIFFERENT data granularities, each contributing exactly what it identifies. For the pooled-HR estimand:
 *
 *   - IPD trial          -> an IDENTIFIED POINT  (logHR, s^2), delta = 0
 *   - HR-only trial      -> an IDENTIFIED POINT  (posted logHR, s^2), delta = 0   (no time-resolved info,
 *                           but for a pooled HR it contributes a clean point)
 *   - reconstructed CURVE -> a PARTIALLY-IDENTIFIED contribution: a de-biased point (Phase 2c) plus a
 *                           residual-bias identification half-width delta_i, and an inflated within-trial
 *                           variance s_rec^2 + r^2 (Phase 1/2b)
 *
 * The pool is a REML random-effects meta-analysis (HKSJ, t_{k-1} PI). Because the CURVE trials carry
 * identification intervals (not points), the OUTPUT is an identified SET: the pooled HR and tau^2 range
 * over all residual-bias configurations of the curve trials within their boxes (greedy extremes); the IPD
 * and HR-only trials are fixed points. All curve calibration (de-bias offset + delta) is LEAVE-ONE-OUT
 * against the true-IPD gold standard, so it is out-of-sample.
 *
 * Two demonstrations:
 *   (1) a realistic mixed corpus -> pooled HR + identified set; check it brackets the all-IPD truth.
 *   (2) the EVIDENCE-COMPLETENESS CURVE: sweep the fraction of trials that are curve-only (rest IPD) and
 *       show the identified-set width grow smoothly from a point (all IPD) to its widest (all curve) -- a
 *       quantitative "granularity -> synthesis precision" trade-off, the new artifact.
 *
 * Run from repo root (needs realipd/cbio_*.csv):
 *   node validate/phase3_granularity_mixed.js  ->  validate/phase3_granularity_mixed_results.json
 */
const fs = require('fs');
const path = require('path');
const { metaRE } = require('./phase2_real_pooling.js');
const { buildCorpus, mean, variance } = require('./phase2c_bias_offset.js');

// Per-trial calibrated contributions (LOO). Each trial gets BOTH an identified-point form (IPD/HR-only:
// true logHR + s_true^2) and a partially-identified curve form (de-biased recon point, inflated variance,
// residual-bias half-width delta).
function calibrate(rows) {
  return rows.map((r, j) => {
    const eO = rows.filter((_, i) => i !== j).map(o => o.e);
    return {
      ds: r.ds,
      point_y: r.yTrue, point_v: r.sTrue2,                       // IPD / HR-only contribution
      curve_y: r.yRec - mean(eO),                                // de-biased reconstruction point
      curve_v: r.sRec2 + r.r2ens,                                // within-trial variance (sampling + recon)
      curve_delta: 1.64 * Math.sqrt(variance(eO)),               // residual-bias identification half-width
    };
  });
}

// Pool a corpus given a granularity assignment (per trial: 'IPD' | 'HRONLY' | 'CURVE'). Returns the central
// estimate and the identified SET over the curve trials' residual-bias boxes.
function pool(cal, assign) {
  const y = [], v = [], delta = [];
  cal.forEach((c, i) => {
    if (assign[i] === 'CURVE') { y.push(c.curve_y); v.push(c.curve_v); delta.push(c.curve_delta); }
    else { y.push(c.point_y); v.push(c.point_v); delta.push(0); }     // IPD or HR-only = identified point
  });
  const central = metaRE(y, v);
  const mu = central.pooled_logHR;
  const shift = (dir) => metaRE(y.map((yi, i) => yi + dir * Math.sign(yi - mu) * delta[i]), v);
  const apart = shift(+1), together = metaRE(
    y.map((yi, i) => yi - Math.sign(yi - mu) * Math.min(delta[i], Math.abs(yi - mu))), v);
  const hiMu = metaRE(y.map((yi, i) => yi + delta[i]), v).pooled_logHR;
  const loMu = metaRE(y.map((yi, i) => yi - delta[i]), v).pooled_logHR;
  return {
    pooled_HR: central.pooled_HR, tau2: +central.tau2.toFixed(3),
    HR_set: [+Math.exp(loMu).toFixed(2), +Math.exp(hiMu).toFixed(2)],
    tau2_set: [+together.tau2.toFixed(3), +apart.tau2.toFixed(3)],
    HR_set_width: +(Math.exp(hiMu) - Math.exp(loMu)).toFixed(3),
  };
}

function run() {
  const rows = buildCorpus();
  const cal = calibrate(rows);
  const k = cal.length;
  const allIPD = pool(cal, cal.map(() => 'IPD'));

  // (1) realistic mixed corpus: alternate IPD / CURVE / HR-only by index (deterministic, no RNG)
  const mixAssign = cal.map((_, i) => ['IPD', 'CURVE', 'HRONLY'][i % 3]);
  const mixed = pool(cal, mixAssign);
  const counts = mixAssign.reduce((m, g) => (m[g] = (m[g] || 0) + 1, m), {});

  // (2) evidence-completeness curve: first f*k trials are CURVE, the rest IPD
  const curve = [];
  for (const f of [0, 0.25, 0.5, 0.75, 1]) {
    const nCurve = Math.round(f * k);
    const assign = cal.map((_, i) => (i < nCurve ? 'CURVE' : 'IPD'));
    const p = pool(cal, assign);
    curve.push({ curve_fraction: f, n_curve: nCurve, pooled_HR: p.pooled_HR,
      HR_set: p.HR_set, HR_set_width: p.HR_set_width, tau2_set: p.tau2_set });
  }
  return {
    k,
    true_all_ipd: { pooled_HR: allIPD.pooled_HR, tau2: allIPD.tau2 },
    mixed_corpus: { composition: counts, ...mixed,
      brackets_truth: allIPD.pooled_HR >= mixed.HR_set[0] && allIPD.pooled_HR <= mixed.HR_set[1] },
    completeness_curve: curve,
  };
}

if (require.main === module) {
  const out = run();
  fs.writeFileSync(path.join(__dirname, 'phase3_granularity_mixed_results.json'), JSON.stringify(out, null, 2));
  console.log(`=== Phase 3: granularity-mixed synthesis (k=${out.k} trials) ===\n`);
  console.log(`  all-IPD reference     : HR ${out.true_all_ipd.pooled_HR.toFixed(2)}  tau^2 ${out.true_all_ipd.tau2}\n`);
  const m = out.mixed_corpus;
  console.log(`  MIXED corpus ${JSON.stringify(m.composition)}:`);
  console.log(`    pooled HR ${m.pooled_HR.toFixed(2)}  identified set [${m.HR_set[0]}, ${m.HR_set[1]}]  (brackets all-IPD HR? ${m.brackets_truth})`);
  console.log(`    tau^2 set [${m.tau2_set[0]}, ${m.tau2_set[1]}]\n`);
  console.log('  EVIDENCE-COMPLETENESS CURVE (more curve-only trials -> wider identified set):');
  console.log('    curve frac  pooled HR   HR identified set        set width');
  for (const c of out.completeness_curve) {
    console.log(`      ${c.curve_fraction.toFixed(2)}      ${c.pooled_HR.toFixed(2)}       [${c.HR_set[0]}, ${c.HR_set[1]}]`.padEnd(54) + `${c.HR_set_width}`);
  }
  console.log('\n  wrote validate/phase3_granularity_mixed_results.json');
}
module.exports = { run, pool, calibrate };
