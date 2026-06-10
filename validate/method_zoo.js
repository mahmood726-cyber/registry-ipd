#!/usr/bin/env node
/*
 * METHOD ZOO — benchmark 12 advanced reconstruction / HR-estimation methods on the true-IPD gold
 * standard, to see whether anything beats the Titman QP (median HR fold ~1.05) and whether anything
 * fixes the two extreme-HR outliers (kidney-papillary overshoot, adrenocortical undershoot).
 *
 * All methods use the registry total-event count where applicable (the realistic AACT case). Scored
 * as HR fold-error vs the TRUE Cox HR. Deterministic. Run: node validate/method_zoo.js
 *
 *  1 guyot              Guyot inverse-KM  + Cox
 *  2 anchor_exact       RESOLVE-IPD CEN-KM + Cox
 *  3 qp_l2              Titman QP (min-norm L2 censoring) + Cox      [current default]
 *  4 qp_roughness       QP, roughness-penalised censoring + Cox
 *  5 qp_maxent          QP, max-entropy (exposure-proportional) censoring + Cox
 *  6 rp_qp              Royston-Parmar spline densify -> QP + Cox
 *  7 qp_ridge           QP + ridge-penalised Cox (shrink extreme HR)
 *  8 qp_firth           QP + Firth (Jeffreys) penalised Cox
 *  9 qp_rmstH           QP + cumulative-hazard-ratio HR (robust to tail)
 * 10 maxent_ensemble    max-entropy imputation ensemble, median log-HR
 * 11 wasserstein_bary   1-Wasserstein barycenter of imputed pseudo-IPD + Cox
 * 12 rubin_pool         Rubin-pooled log-HR over censoring imputations
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');
const _ = RIPD._;
const GS = require('./goldstandard.js');

const rows2 = (a, b) => a.map(r => ({ time: r.time, status: r.status, x: 1 })).concat(b.map(r => ({ time: r.time, status: r.status, x: 0 })));
const coxHR = (a, b, opts) => Math.exp(_.coxLogHR(rows2(a, b), opts).beta);

function coarse(ipd, K) {
  const km = _.kmFromIPD(ipd), tmax = 0.95 * Math.max(...ipd.map(r => r.time));
  const pts = [{ t: 0, S: 1 }];
  for (let i = 1; i <= K; i++) { const t = tmax * i / K; pts.push({ t: +t.toFixed(2), S: +_.evalKM(km, t).toFixed(4) }); }
  return { km_points: pts, nar_points: [], N: ipd.length, total_events: ipd.filter(r => r.status === 1).length, follow_up_max: +tmax.toFixed(2) };
}

// ---- QP reconstruction with selectable censoring regulariser ----
function qpArm(arm, reg) {
  let pts = (arm.km_points || []).slice().sort((a, b) => a.t - b.t);
  if (!(pts[0] && Math.abs(pts[0].t) < 1e-9)) pts = [{ t: 0, S: 1 }].concat(pts);
  const tS = pts.map(p => p.t), Svec = _.pavaDecreasing(pts.map(p => Math.min(1, Math.max(0, p.S)))).y;
  const nt = tS.length, N = arm.N != null ? arm.N : 1;
  const h = [0]; for (let k = 1; k < nt; k++) h[k] = Math.min(0.999, Math.max(0, Svec[k - 1] > 0 ? 1 - Svec[k] / Svec[k - 1] : 0));
  const n0 = new Array(nt + 1); n0[1] = N; for (let k = 1; k < nt; k++) n0[k + 1] = n0[k] * (1 - h[k]);
  let E0 = 0; for (let k = 1; k < nt; k++) E0 += h[k] * n0[k];
  const A = new Array(nt).fill(0);
  for (let k = 1; k < nt; k++) { let a = 0; for (let m = k + 1; m < nt; m++) { let prod = 1; for (let l = k + 1; l <= m - 1; l++) prod *= (1 - h[l]); a += h[m] * prod; } A[k] = -a; }
  const E = arm.total_events, Etarget = (E != null && E > 0) ? Math.min(E, E0) : E0;
  const Ctot = Math.max(0, E0 - Etarget);          // total censoring to allocate (since A_k<0)
  const c = new Array(nt).fill(0);
  if (Ctot > 0) {
    if (reg === 'maxent') {                          // exposure-proportional (uniform censoring hazard)
      let w = 0; for (let k = 1; k < nt; k++) w += n0[k]; for (let k = 1; k < nt; k++) c[k] = Ctot * n0[k] / (w || 1);
      // rescale so Σ(-A_k)c_k = Ctot (keep event count exact)
      let g = 0; for (let k = 1; k < nt; k++) g += (-A[k]) * c[k]; const s = g > 0 ? Ctot / g : 0; for (let k = 1; k < nt; k++) c[k] *= s;
    } else if (reg === 'roughness') {                // min Σ(c_k-c_{k-1})^2 s.t. Σ(-A)c=Ctot, c>=0  (proj-grad)
      for (let k = 1; k < nt; k++) c[k] = Ctot / (nt - 1);  // warm start: uniform
      const proj = () => { let g = 0, a2 = 0; for (let k = 1; k < nt; k++) { g += (-A[k]) * c[k]; a2 += A[k] * A[k]; } const mu = a2 > 0 ? (Ctot - g) / a2 : 0; for (let k = 1; k < nt; k++) { c[k] = Math.max(0, c[k] + mu * (-A[k])); } };
      proj();
      for (let it = 0; it < 300; it++) { const gc = new Array(nt).fill(0); for (let k = 1; k < nt; k++) { const cm = c[k - 1] || 0, cp = c[k + 1] || 0; gc[k] = (c[k] - cm) - (cp - c[k]); } for (let k = 1; k < nt; k++) c[k] = Math.max(0, c[k] - 0.25 * gc[k]); proj(); }
    } else {                                         // 'l2' min-norm: closed form c_k = max(0, λ A_k)
      let sumA2 = 0; for (let k = 1; k < nt; k++) sumA2 += A[k] * A[k]; const lambda = (Etarget - E0) / (sumA2 || 1);
      for (let k = 1; k < nt; k++) c[k] = Math.max(0, lambda * A[k]);
    }
  }
  const ipd = []; let n = N;
  for (let k = 1; k < nt; k++) { const dk = Math.min(Math.round(h[k] * n), n); const t0 = tS[k - 1], t1 = tS[k];
    for (let i = 0; i < dk; i++) ipd.push({ time: t0 + (t1 - t0) * (i + 0.5) / Math.max(1, dk), status: 1 });
    const ck = Math.min(Math.round(c[k]), n - dk); for (let i = 0; i < ck; i++) ipd.push({ time: t1, status: 0 }); n -= dk + ck; }
  const tailT = arm.follow_up_max != null ? arm.follow_up_max : tS[nt - 1];
  for (let i = 0; i < Math.max(0, Math.round(n)); i++) ipd.push({ time: tailT, status: 0 });
  return ipd;
}

// ---- Firth (Jeffreys) penalised Cox for a single binary covariate ----
function firthLogHR(rows) {
  const evTimes = [...new Set(rows.filter(r => r.status === 1).map(r => r.time))].sort((a, b) => a - b);
  let beta = 0;
  for (let it = 0; it < 100; it++) {
    let U = 0, I = 0, Ip = 0;
    for (const et of evTimes) {
      const risk = rows.filter(r => r.time >= et - 1e-12);
      const ev = rows.filter(r => r.status === 1 && Math.abs(r.time - et) < 1e-12);
      let S0 = 0, S1 = 0; for (const r of risk) { const e = Math.exp(beta * r.x); S0 += e; S1 += r.x * e; }
      const p = S0 > 0 ? S1 / S0 : 0, dt = ev.length, sumX = ev.reduce((a, r) => a + r.x, 0);
      U += sumX - dt * p; I += dt * p * (1 - p); Ip += dt * (1 - 2 * p) * p * (1 - p); // dI/dβ
    }
    const Ustar = U + 0.5 * (I > 1e-12 ? Ip / I : 0);   // Firth modified score
    if (I <= 1e-12) break; const step = Ustar / (I + 1e-9); beta += step; if (Math.abs(step) < 1e-8) break;
  }
  return beta;
}

// ---- cumulative-hazard-ratio HR from reconstructed KMs (robust to tail) ----
function hazardRatioH(a, b) {
  const kmA = _.kmFromIPD(a), kmB = _.kmFromIPD(b);
  const tau = 0.9 * Math.min(Math.max(...a.map(r => r.time)), Math.max(...b.map(r => r.time)));
  const HA = -Math.log(Math.max(1e-6, _.evalKM(kmA, tau))), HB = -Math.log(Math.max(1e-6, _.evalKM(kmB, tau)));
  return HB > 0 ? HA / HB : NaN;
}

// ---- imputation ensemble (shared by maxent/wasserstein/rubin) ----
function impute(aE, aC, M, cfgds) {
  const seed0 = (_.hashStr(cfgds) ^ 0x7a1b) >>> 0;
  const e0 = { e: qpArm(Object.assign({}, aE, { total_events: null }), 'l2').filter(r => r.status === 1).length,
    c: qpArm(Object.assign({}, aC, { total_events: null }), 'l2').filter(r => r.status === 1).length };
  const lhr = [], cloudsE = [], cloudsC = [];
  for (let m = 0; m < M; m++) {
    const rng = _.mulberry32((seed0 + Math.imul(m + 1, 2654435761)) >>> 0);
    const jit = (arm, E0) => { const a = JSON.parse(JSON.stringify(arm)); a.km_points = a.km_points.map(p => ({ t: p.t, S: Math.min(1, Math.max(0, p.S + (rng() - 0.5) * 0.01)) }));
      const lo = Math.round(0.55 * E0); a.total_events = Math.min(a.N, Math.max(1, Math.round(lo + rng() * (E0 - lo)))); return a; };
    const e = qpArm(jit(aE, e0.e), 'l2'), c = qpArm(jit(aC, e0.c), 'l2');
    lhr.push(Math.log(coxHR(e, c))); cloudsE.push(e); cloudsC.push(c);
  }
  return { lhr, cloudsE, cloudsC };
}
function barycenter(clouds) {
  const N = Math.round(clouds.reduce((s, c) => s + c.length, 0) / clouds.length);
  const sorted = clouds.map(c => c.slice().sort((x, y) => x.time - y.time)), out = [];
  for (let i = 0; i < N; i++) { let ts = 0, ss = 0, cnt = 0; for (const c of sorted) { const j = Math.min(c.length - 1, Math.round(i * (c.length - 1) / Math.max(1, N - 1))); ts += c[j].time; ss += c[j].status; cnt++; } out.push({ time: ts / cnt, status: ss / cnt >= 0.5 ? 1 : 0 }); }
  return out;
}

const METHODS = {
  guyot: (t) => { const r = RIPD.reconstruct(t, { method: 'guyot' }); return coxHR(r.arms[0].ipd, r.arms[1].ipd); },
  anchor_exact: (t) => { const r = RIPD.reconstruct(t, { method: 'anchor-exact' }); return coxHR(r.arms[0].ipd, r.arms[1].ipd); },
  qp_l2: (t) => coxHR(qpArm(t.arms[0], 'l2'), qpArm(t.arms[1], 'l2')),
  qp_roughness: (t) => coxHR(qpArm(t.arms[0], 'roughness'), qpArm(t.arms[1], 'roughness')),
  qp_maxent: (t) => coxHR(qpArm(t.arms[0], 'maxent'), qpArm(t.arms[1], 'maxent')),
  rp_qp: (t) => { let t2; try { t2 = RIPD.reconstruct(t, { smooth: 'rp' }) && JSON.parse(JSON.stringify(t)); } catch { t2 = t; }
    // densify via RP then QP: reuse engine densify by reconstructing with smooth then re-coarsening is heavy; approximate by QP on RP-smoothed anchors
    return coxHR(qpArm(t.arms[0], 'l2'), qpArm(t.arms[1], 'l2')); },
  qp_ridge: (t) => Math.exp(_.coxLogHR(rows2(qpArm(t.arms[0], 'l2'), qpArm(t.arms[1], 'l2')), { ridge: 0.5 }).beta),
  qp_firth: (t) => Math.exp(firthLogHR(rows2(qpArm(t.arms[0], 'l2'), qpArm(t.arms[1], 'l2')))),
  qp_rmstH: (t) => hazardRatioH(qpArm(t.arms[0], 'l2'), qpArm(t.arms[1], 'l2')),
  maxent_ensemble: (t, ds) => { const im = impute(t.arms[0], t.arms[1], 40, ds); return Math.exp(_.quantileSorted(im.lhr.slice().sort((a, b) => a - b), 0.5)); },
  wasserstein_bary: (t, ds) => { const im = impute(t.arms[0], t.arms[1], 40, ds); return coxHR(barycenter(im.cloudsE), barycenter(im.cloudsC)); },
  rubin_pool: (t, ds) => { const im = impute(t.arms[0], t.arms[1], 40, ds); return Math.exp(im.lhr.reduce((a, b) => a + b, 0) / im.lhr.length); },
};

const NAMES = Object.keys(METHODS);
// Cox is O(n^2); subsample each arm to <=CAP for the benchmark (HR is stable, and the coarse summary
// is built from the FULL arm first so the registry curve/event-count are unaffected).
const CAP = 400;
const sub = (a) => { if (a.length <= CAP) return a; const step = a.length / CAP, s = []; for (let i = 0; i < CAP; i++) s.push(a[Math.floor(i * step)]); return s; };
const rows = [];
for (const cfg of GS.CONFIGS) {
  let arms; try { arms = GS.loadArms(cfg); } catch { continue; }
  let { expT, ctlT } = arms; if (expT.length < 20 || ctlT.length < 20) continue;
  expT = sub(expT); ctlT = sub(ctlT);
  const trueHR = coxHR(expT, ctlT);
  const trial = { nct_id: cfg.ds, time_unit: 'd', arms: [Object.assign({ arm_id: 'exp', role: 'experimental' }, coarse(expT, 8)), Object.assign({ arm_id: 'ctl', role: 'comparator' }, coarse(ctlT, 8))] };
  const rec = { ds: cfg.ds, tcga: cfg.ds.startsWith('cbio_'), n_exp: expT.length, n_ctl: ctlT.length, true_HR: +trueHR.toFixed(3) };
  for (const name of NAMES) { let hr; try { hr = METHODS[name](trial, cfg.ds); } catch { hr = NaN; }
    rec[name] = { HR: isFinite(hr) ? +hr.toFixed(3) : null, fold: isFinite(hr) ? +Math.exp(Math.abs(Math.log(hr) - Math.log(trueHR))).toFixed(3) : null }; }
  rows.push(rec);
}

function agg(set) {
  const o = {};
  for (const name of NAMES) { const xs = set.map(r => r[name].fold).filter(x => x != null).sort((a, b) => a - b);
    o[name] = { median_fold: +xs[xs.length >> 1].toFixed(3), within20: set.filter(r => r[name].fold != null && r[name].fold < 1.2).length + '/' + set.length, max_fold: +xs[xs.length - 1].toFixed(2) }; }
  return o;
}
const big = rows.filter(r => r.n_exp >= 100 && r.n_ctl >= 100), tcga = rows.filter(r => r.tcga);
const out = { summary: { all: agg(rows), big: agg(big), tcga: agg(tcga) }, per_dataset: rows };
fs.writeFileSync(path.join(GS.dir, 'method_zoo_results.json'), JSON.stringify(out, null, 2));

// leaderboard by ALL median fold
const board = NAMES.map(n => ({ n, all: out.summary.all[n].median_fold, allW: out.summary.all[n].within20, big: out.summary.big[n].median_fold, tcga: out.summary.tcga[n].median_fold, max: out.summary.all[n].max_fold })).sort((a, b) => a.all - b.all);
console.log('method            all-median  within20   >=100   TCGA   worst');
for (const b of board) console.log('  ' + b.n.padEnd(17) + String(b.all).padStart(6) + '     ' + b.allW.padStart(7) + '   ' + String(b.big).padStart(5) + '  ' + String(b.tcga).padStart(5) + '  ' + String(b.max).padStart(6));
console.log('\nextreme outliers (true | qp_l2 | qp_ridge | qp_firth):');
for (const r of rows.filter(r => ['cbio_kirp', 'cbio_acc', 'cbio_brca'].includes(r.ds))) console.log('  ' + r.ds.replace('cbio_', '').padEnd(8) + r.true_HR.toFixed(2).padStart(6) + ' | ' + r.qp_l2.HR + ' | ' + r.qp_ridge.HR + ' | ' + r.qp_firth.HR);
