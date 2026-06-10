#!/usr/bin/env node
/*
 * ADVANCED ESTIMATORS for the curve-only (no registry event count) regime.
 *
 * The registry KM curve + N do NOT identify the censoring level, and the Cox HR depends on the
 * at-risk sets, which depend on censoring. So the curve-only point estimate is genuinely ambiguous.
 * The current default censors to the tail, which UNDERESTIMATES large HRs under heavy censoring
 * (the TCGA finding). This benchmarks three principled point estimators on the gold standard:
 *
 *   A. censor-to-tail        : the current curve-only reconstruct() (baseline).
 *   B. max-entropy ensemble  : median log-HR over M imputations of the under-identified censoring
 *                              level (the engine's ensemble summary) — a model-averaging estimate.
 *   C. Wasserstein barycenter: the 1-Wasserstein barycenter of the M reconstructed IPD point-clouds
 *                              per arm. In 1-D the optimal coupling is rank-matching, so the
 *                              barycenter is the rank-averaged event times WITH rank-averaged
 *                              censoring status — uniquely, this preserves the at-risk structure the
 *                              HR depends on, unlike averaging the KM curves. Cox on the two
 *                              barycentric arms. (Cutting-edge / optimal-transport point estimate.)
 *
 * Deterministic (seeded). Run from repo root: node validate/advanced_estimators.js
 * Writes realipd/advanced_estimators_results.json.
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');
const _ = RIPD._;
const GS = require('./goldstandard.js');

const M = 160;          // imputations
const K = 8;            // posted timepoints

function coarse(ipd) {
  const km = _.kmFromIPD(ipd), tmax = 0.95 * Math.max(...ipd.map(r => r.time));
  const pts = [{ t: 0, S: 1 }];
  for (let i = 1; i <= K; i++) { const t = tmax * i / K; pts.push({ t: +t.toFixed(2), S: +_.evalKM(km, t).toFixed(4) }); }
  return { km_points: pts, nar_points: [], N: ipd.length, total_events: ipd.filter(r => r.status === 1).length, follow_up_max: +tmax.toFixed(2) };
}
const coxHR = (a, b) => _.coxLogHR(a.map(r => ({ time: r.time, status: r.status, x: 1 }))
  .concat(b.map(r => ({ time: r.time, status: r.status, x: 0 })))).hr;

// Wasserstein-1 barycenter of M labelled point-clouds for ONE arm. Each cloud is an IPD (time,status)
// of (near-)equal size N. Optimal 1-D coupling = sort each cloud, match by rank, average per rank.
function barycenterArm(clouds) {
  const N = Math.round(clouds.reduce((s, c) => s + c.length, 0) / clouds.length);
  const sorted = clouds.map(c => c.slice().sort((x, y) => x.time - y.time));
  const out = [];
  for (let i = 0; i < N; i++) {
    let tsum = 0, ssum = 0, cnt = 0;
    for (const c of sorted) {
      // map rank i in [0,N) to this cloud's index
      const j = Math.min(c.length - 1, Math.round(i * (c.length - 1) / Math.max(1, N - 1)));
      tsum += c[j].time; ssum += c[j].status; cnt++;
    }
    out.push({ time: tsum / cnt, status: ssum / cnt >= 0.5 ? 1 : 0 });
  }
  return out;
}

function runDataset(cfg) {
  let arms; try { arms = GS.loadArms(cfg); } catch (e) { return null; }
  const { expT, ctlT } = arms;
  if (expT.length < 20 || ctlT.length < 20) return null;
  const trueHR = coxHR(expT, ctlT);
  const trial = { nct_id: cfg.ds, time_unit: 'd',
    arms: [Object.assign({ arm_id: 'exp', role: 'experimental' }, coarse(expT)),
           Object.assign({ arm_id: 'ctl', role: 'comparator' }, coarse(ctlT)) ] };
  // curve-only event ceiling E0 per arm (no intermediate censoring)
  const e0 = {};
  { const cr = RIPD.reconstruct(trial, { ignoreTotalEvents: true });
    if (cr.arms) cr.arms.forEach(a => { e0[a.arm_id] = a.ipd.filter(x => x.status === 1).length; }); }

  // A. censor-to-tail baseline (curve-only: strip events)
  const t0 = JSON.parse(JSON.stringify(trial)); t0.arms.forEach(a => { a.total_events = null; });
  const rA = RIPD.reconstruct(t0, {});
  const hrA = rA.arms ? coxHR(rA.arms[0].ipd, rA.arms[1].ipd) : NaN;

  // B + C: M imputations of the censoring level (max-entropy band [0.55 E0, E0])
  const seed0 = (_.hashStr(cfg.ds) ^ 0x5151abcd) >>> 0;
  const lhr = [];
  const clouds = { exp: [], ctl: [] };
  for (let m = 0; m < M; m++) {
    const rng = _.mulberry32((seed0 + Math.imul(m + 1, 2654435761)) >>> 0);
    const t2 = JSON.parse(JSON.stringify(trial));
    for (const a of t2.arms) {
      a.km_points = a.km_points.map(p => ({ t: p.t, S: Math.min(1, Math.max(0, p.S + (rng() - 0.5) * 0.01)) }));
      const E0 = e0[a.arm_id];
      if (E0 != null && E0 > 0) { const lo = Math.round(0.55 * E0); a.total_events = Math.min(a.N, Math.max(1, Math.round(lo + rng() * (E0 - lo)))); }
    }
    const method = rng() < 0.5 ? 'guyot' : 'anchor-exact';
    let r; try { r = RIPD.reconstruct(t2, { method }); } catch { continue; }
    if (!r.arms) continue;
    const e = r.arms.find(a => a.role === 'experimental'), c = r.arms.find(a => a.role === 'comparator');
    lhr.push(Math.log(coxHR(e.ipd, c.ipd)));
    clouds.exp.push(e.ipd); clouds.ctl.push(c.ipd);
  }
  const hrB = lhr.length ? Math.exp(_.quantileSorted(lhr.slice().sort((a, b) => a - b), 0.5)) : NaN;
  const hrC = (clouds.exp.length && clouds.ctl.length)
    ? coxHR(barycenterArm(clouds.exp), barycenterArm(clouds.ctl)) : NaN;

  const fe = (h) => isFinite(h) ? Math.exp(Math.abs(Math.log(h) - Math.log(trueHR))) : null;
  return { ds: cfg.ds, n_exp: expT.length, n_ctl: ctlT.length, true_HR: +trueHR.toFixed(3),
    tcga: cfg.ds.startsWith('cbio_'),
    censor_to_tail: { HR: +hrA.toFixed(3), fold: +fe(hrA).toFixed(3) },
    maxent_ensemble: { HR: +hrB.toFixed(3), fold: +fe(hrB).toFixed(3) },
    wasserstein_barycenter: { HR: +hrC.toFixed(3), fold: +fe(hrC).toFixed(3) } };
}

const rows = GS.CONFIGS.map(runDataset).filter(Boolean);
function agg(set) {
  const med = (key) => { const xs = set.map(r => r[key].fold).filter(x => x != null).sort((a, b) => a - b); return +xs[xs.length >> 1].toFixed(3); };
  const within20 = (key) => set.filter(r => r[key].fold != null && r[key].fold < 1.2).length + '/' + set.length;
  return {
    n: set.length,
    censor_to_tail: { median_fold: med('censor_to_tail'), within20: within20('censor_to_tail') },
    maxent_ensemble: { median_fold: med('maxent_ensemble'), within20: within20('maxent_ensemble') },
    wasserstein_barycenter: { median_fold: med('wasserstein_barycenter'), within20: within20('wasserstein_barycenter') },
  };
}
const out = {
  summary: {
    all: agg(rows),
    tcga_heavily_censored: agg(rows.filter(r => r.tcga)),
    note: 'Curve-only regime (no registry event count). censor_to_tail = current default; '
      + 'maxent_ensemble = median log-HR over censoring-level imputations; wasserstein_barycenter = '
      + 'Cox HR on the 1-Wasserstein barycenter (rank-averaged, censoring-preserving) of the imputed IPD.',
  },
  per_dataset: rows,
};
fs.writeFileSync(path.join(GS.dir, 'advanced_estimators_results.json'), JSON.stringify(out, null, 2));
console.log(JSON.stringify(out.summary, null, 2));
console.log('\nTCGA per-dataset (true | censor-tail | maxent | wasserstein):');
for (const r of rows.filter(r => r.tcga)) console.log('  ' + r.ds.replace('cbio_', '').padEnd(10)
  + r.true_HR.toFixed(2).padStart(6) + ' | ' + r.censor_to_tail.fold.toFixed(2) + ' | '
  + r.maxent_ensemble.fold.toFixed(2) + ' | ' + r.wasserstein_barycenter.fold.toFixed(2));
