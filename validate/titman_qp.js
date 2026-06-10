#!/usr/bin/env node
/*
 * TITMAN-2026-style QUADRATIC-PROGRAM reconstruction (Stat Med 2026;45(6-7):e70474, PMID 41775249).
 *
 * Titman formulates KM-IPD reconstruction as a QP with linear constraints, shining "when numbers at
 * risk and marked censoring times are available". We adapt it to the AACT regime (no number-at-risk,
 * but a posted total-event count) with a key linearisation:
 *
 *   On the cumulative-hazard scale the posted curve fixes the per-interval discrete hazards
 *      h_j = 1 - S_j/S_{j-1}   (known from the anchors),
 *   so events in interval j are d_j = h_j * n_j and the at-risk recursion
 *      n_{j+1} = n_j (1 - h_j) - c_j
 *   is LINEAR in the unknown censoring counts c_j. The total-event count is a LINEAR constraint
 *      E = sum_j h_j n_j,
 *   and the one remaining degree of freedom (how censoring is distributed) is resolved by the convex
 *   QP   min ½‖c‖²  s.t.  E(c) = E_reported,  c ≥ 0,  n_j ≥ 0.
 *   Because ∂E/∂c_j ≡ A_j is constant, the minimum-norm non-negative solution is closed-form:
 *      c_j = max(0, λ A_j),   λ = (E_reported − E0) / Σ A_j²,   E0 = N(1 − S_K).
 *
 * This is a faithful, exactly-solvable QP for the structured-registry setting. We benchmark it against
 * the engine's current censoring-informed (anchor-exact) reconstruction on the true-IPD gold standard.
 *
 * Run from repo root: node validate/titman_qp.js   ->   realipd/titman_qp_results.json
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');
const _ = RIPD._;
const GS = require('./goldstandard.js');

const coxHR = (a, b) => _.coxLogHR(a.map(r => ({ time: r.time, status: r.status, x: 1 }))
  .concat(b.map(r => ({ time: r.time, status: r.status, x: 0 })))).hr;

// Titman-style QP reconstruction of one arm. Needs km_points (t,S), N, total_events.
function reconstructArmQP(arm) {
  const pts = (arm.km_points || []).slice().sort((a, b) => a.t - b.t);
  if (!pts.length) return [];
  if (pts[0].t > 0) pts.unshift({ t: 0, S: 1 });
  const K = pts.length - 1;
  const N = arm.N;
  const E = arm.total_events;
  const hh = [null];
  for (let j = 1; j <= K; j++) { let hj = 1 - pts[j].S / pts[j - 1].S; hh[j] = Math.min(0.999, Math.max(0, hj)); }
  // curve-only at-risk + event ceiling E0 = N(1 - S_K)
  const n0 = [null, N];
  for (let j = 1; j <= K; j++) n0[j + 1] = n0[j] * (1 - hh[j]);
  let E0 = 0; for (let j = 1; j <= K; j++) E0 += hh[j] * n0[j];
  // sensitivities A_j = ∂E/∂c_j = -Σ_{m>j} h_m Π_{l=j+1}^{m-1}(1-h_l)  (constant, ≤ 0)
  const A = [null];
  for (let j = 1; j <= K; j++) {
    let a = 0;
    for (let m = j + 1; m <= K; m++) { let prod = 1; for (let l = j + 1; l <= m - 1; l++) prod *= (1 - hh[l]); a += hh[m] * prod; }
    A[j] = -a;
  }
  const sumA2 = A.slice(1).reduce((s, a) => s + a * a, 0) || 1;
  const Etarget = (E != null ? Math.min(E, E0) : E0);
  const lambda = (Etarget - E0) / sumA2;            // ≤ 0
  const c = [null];
  for (let j = 1; j <= K; j++) c[j] = Math.max(0, lambda * A[j]);
  // realise at-risk with c, guard ≥ 0
  const n = [null, N];
  for (let j = 1; j <= K; j++) n[j + 1] = Math.max(0, n[j] * (1 - hh[j]) - c[j]);
  // build pseudo-IPD: d_j events spread across interval j, c_j censorings at the interval end
  const ipd = [];
  for (let j = 1; j <= K; j++) {
    const dj = Math.round(hh[j] * n[j]);
    const t0 = pts[j - 1].t, t1 = pts[j].t;
    for (let i = 0; i < dj; i++) ipd.push({ time: t0 + (t1 - t0) * (i + 0.5) / Math.max(1, dj), status: 1 });
    const cj = Math.round(c[j]);
    for (let i = 0; i < cj; i++) ipd.push({ time: t1, status: 0 });
  }
  const rem = Math.max(0, N - ipd.length);
  for (let i = 0; i < rem; i++) ipd.push({ time: pts[K].t, status: 0 });
  return ipd;
}

function coarse(ipd, K) {
  const km = _.kmFromIPD(ipd), tmax = 0.95 * Math.max(...ipd.map(r => r.time));
  const pts = [{ t: 0, S: 1 }];
  for (let i = 1; i <= K; i++) { const t = tmax * i / K; pts.push({ t: +t.toFixed(2), S: +_.evalKM(km, t).toFixed(4) }); }
  return { km_points: pts, nar_points: [], N: ipd.length, total_events: ipd.filter(r => r.status === 1).length, follow_up_max: +tmax.toFixed(2) };
}

const rows = [];
for (const cfg of GS.CONFIGS) {
  let arms; try { arms = GS.loadArms(cfg); } catch (e) { continue; }
  const { expT, ctlT } = arms;
  if (expT.length < 20 || ctlT.length < 20) continue;
  const trueHR = coxHR(expT, ctlT);
  const aE = coarse(expT, 8), aC = coarse(ctlT, 8);
  // QP reconstruction
  const hrQP = coxHR(reconstructArmQP(Object.assign({ role: 'exp' }, aE)), reconstructArmQP(Object.assign({ role: 'ctl' }, aC)));
  // previous engine censoring-informed method (anchor-exact) for comparison — the QP is now the
  // engine default, so we force anchor-exact here to show what the QP replaced.
  const trial = { nct_id: cfg.ds, time_unit: 'd',
    arms: [Object.assign({ arm_id: 'exp', role: 'experimental' }, aE), Object.assign({ arm_id: 'ctl', role: 'comparator' }, aC)] };
  const rci = RIPD.reconstruct(trial, { method: 'anchor-exact' });
  const hrCI = rci.arms ? coxHR(rci.arms[0].ipd, rci.arms[1].ipd) : NaN;
  const fe = (h) => isFinite(h) ? +Math.exp(Math.abs(Math.log(h) - Math.log(trueHR))).toFixed(3) : null;
  rows.push({ ds: cfg.ds, tcga: cfg.ds.startsWith('cbio_'), n_exp: expT.length, n_ctl: ctlT.length,
    true_HR: +trueHR.toFixed(3), qp_HR: +hrQP.toFixed(3), censinf_HR: +hrCI.toFixed(3),
    qp_fold: fe(hrQP), censinf_fold: fe(hrCI) });
}

function agg(set) {
  const med = (k) => { const xs = set.map(r => r[k]).filter(x => x != null).sort((a, b) => a - b); return +xs[xs.length >> 1].toFixed(3); };
  const w20 = (k) => set.filter(r => r[k] != null && r[k] < 1.2).length + '/' + set.length;
  return { n: set.length, qp: { median_fold: med('qp_fold'), within20: w20('qp_fold') }, censinf: { median_fold: med('censinf_fold'), within20: w20('censinf_fold') } };
}
const out = { summary: { all: agg(rows), tcga: agg(rows.filter(r => r.tcga)), big: agg(rows.filter(r => r.n_exp >= 100 && r.n_ctl >= 100)),
  note: 'Titman-2026-style QP reconstruction vs the engine censoring-informed (anchor-exact), both using the registry total-event count, scored as HR fold-error vs true IPD.' }, per_dataset: rows };
fs.writeFileSync(path.join(GS.dir, 'titman_qp_results.json'), JSON.stringify(out, null, 2));
console.log(JSON.stringify(out.summary, null, 2));
console.log('\nTCGA (true | QP fold | cens-informed fold):');
for (const r of rows.filter(r => r.tcga)) console.log('  ' + r.ds.replace('cbio_', '').padEnd(10) + r.true_HR.toFixed(2).padStart(6) + ' | ' + r.qp_fold.toFixed(2) + ' | ' + r.censinf_fold.toFixed(2));
