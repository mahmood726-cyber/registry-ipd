#!/usr/bin/env node
/*
 * NAR-FUSION experiment — registry-exact curve + figure-OCR'd numbers-at-risk.
 *
 * registry-IPD's binding accuracy limit is the ABSENT number-at-risk (AACT has 0 NAR rows; censoring is
 * unidentified from the curve alone, so curve-only underestimates large HRs ~1.5-fold). The sibling
 * `kmcurve` project OCRs the numbers-at-risk table off a published figure. This tests the synergy idea
 * (KMCURVE-SYNERGY.md, idea 2): for a trial that is in both AACT and a publication, FUSE the
 * registry-exact curve anchors with the figure's NAR — and recover the HR WITHOUT a registry event
 * count.
 *
 * On the true-IPD gold standard we simulate the fusion honestly: exact KM anchors (registry) + the NAR
 * at a SPARSE set of times (as a published risk table prints, ~4 columns) with OCR-style rounding noise
 * — but NO total-event count. The engine's anchor-exact reconstructor consumes nar_points natively, so
 * this is the fusion. Compared against (a) curve-only (no NAR, no events) and (b) the Titman QP (curve +
 * registry event count). If fusion ≈ QP, the figure's NAR substitutes for the missing registry event
 * count — best-of-both registry+figure reconstruction.
 *
 * Deterministic. Run: node validate/nar_fusion.js  ->  realipd/nar_fusion_results.json
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');
const _ = RIPD._;
const GS = require('./goldstandard.js');

const coxHR = (a, b) => Math.exp(_.coxLogHR(a.map(r => ({ time: r.time, status: r.status, x: 1 }))
  .concat(b.map(r => ({ time: r.time, status: r.status, x: 0 }))).slice()).beta);
const atRisk = (ipd, t) => ipd.filter(r => r.time >= t - 1e-9).length;

// coarse registry-style curve + (optionally) a sparse OCR'd NAR table at NEVERY-th anchor
function buildArm(ipd, K, narEvery, rng) {
  const km = _.kmFromIPD(ipd), tmax = 0.95 * Math.max(...ipd.map(r => r.time));
  const pts = [{ t: 0, S: 1 }], nar = [{ t: 0, n: ipd.length }];
  for (let i = 1; i <= K; i++) {
    const t = +(tmax * i / K).toFixed(2);
    pts.push({ t, S: +_.evalKM(km, t).toFixed(4) });
    if (narEvery && i % narEvery === 0) {
      let n = atRisk(ipd, t);
      if (rng) n = Math.max(0, Math.round(n * (1 + (rng() - 0.5) * 0.06))); // ±3% OCR-style noise
      nar.push({ t, n });
    }
  }
  return { km_points: pts, nar_points: nar, N: ipd.length, total_events: ipd.filter(r => r.status === 1).length, follow_up_max: +tmax.toFixed(2) };
}

function mkTrial(ds, e, c) {
  return { nct_id: ds, time_unit: 'd', arms: [Object.assign({ arm_id: 'exp', role: 'experimental' }, e), Object.assign({ arm_id: 'ctl', role: 'comparator' }, c)] };
}

// NAR-aware QP-style reconstruction: take the at-risk path n(t) from the figure NAR (interpolated to
// every anchor), compute events d_k = h_k·n_k from the curve, and SPREAD events within each interval
// (the QP's key fix). Uses the figure NAR INSTEAD of a registry event count.
function narAwareReconstruct(arm) {
  let pts = arm.km_points.slice().sort((a, b) => a.t - b.t);
  if (!(pts[0] && Math.abs(pts[0].t) < 1e-9)) pts = [{ t: 0, S: 1 }].concat(pts);
  const tS = pts.map(p => p.t), Sv = _.pavaDecreasing(pts.map(p => Math.min(1, Math.max(0, p.S)))).y, nt = tS.length, N = arm.N;
  const nar = (arm.nar_points || []).slice().filter(p => isFinite(p.t) && isFinite(p.n)).sort((a, b) => a.t - b.t);
  if (!nar.length || nar[0].t > 1e-9) nar.unshift({ t: 0, n: N });
  const interpN = (t) => { // piecewise-linear at-risk from the (sparse) NAR table, monotone non-increasing
    if (t <= nar[0].t) return nar[0].n;
    for (let i = 1; i < nar.length; i++) if (t <= nar[i].t) { const a = nar[i - 1], b = nar[i]; return a.n + (b.n - a.n) * (t - a.t) / Math.max(1e-9, b.t - a.t); }
    return nar[nar.length - 1].n;
  };
  const n = tS.map(interpN);
  const h = [0]; for (let k = 1; k < nt; k++) h[k] = Math.min(0.999, Math.max(0, Sv[k - 1] > 0 ? 1 - Sv[k] / Sv[k - 1] : 0));
  const ipd = [];
  for (let k = 1; k < nt; k++) {
    const dk = Math.min(Math.round(h[k] * n[k - 1]), Math.round(n[k - 1]));
    const t0 = tS[k - 1], t1 = tS[k];
    for (let i = 0; i < dk; i++) ipd.push({ time: t0 + (t1 - t0) * (i + 0.5) / Math.max(1, dk), status: 1 });
    const ck = Math.max(0, Math.round(n[k - 1] - n[k] - dk));   // censoring = at-risk drop minus events
    for (let i = 0; i < ck; i++) ipd.push({ time: t1, status: 0 });
  }
  const rem = Math.max(0, Math.round(N - ipd.length)), tailT = arm.follow_up_max != null ? arm.follow_up_max : tS[nt - 1];
  for (let i = 0; i < rem; i++) ipd.push({ time: tailT, status: 0 });
  return ipd;
}
function hrOf(trial, opts) { const r = RIPD.reconstruct(trial, opts); return r.arms ? coxHR(r.arms[0].ipd, r.arms[1].ipd) : NaN; }

const rows = [];
for (const cfg of GS.CONFIGS) {
  let arms; try { arms = GS.loadArms(cfg); } catch { continue; }
  const { expT, ctlT } = arms; if (expT.length < 20 || ctlT.length < 20) continue;
  const trueHR = coxHR(expT, ctlT);
  const rng = _.mulberry32(_.hashStr(cfg.ds) ^ 0x4e41524f);
  // curve-only: no NAR, no events
  const cuE = buildArm(expT, 8, 0, null), cuC = buildArm(ctlT, 8, 0, null);
  const tCu = mkTrial(cfg.ds, cuE, cuC); tCu.arms.forEach(a => { a.total_events = null; a.nar_points = []; });
  // fusion: curve + sparse OCR'd NAR (every 2nd anchor → 4 NAR columns), NO event count
  const fuE = buildArm(expT, 8, 2, rng), fuC = buildArm(ctlT, 8, 2, rng);
  const tFu = mkTrial(cfg.ds, fuE, fuC); tFu.arms.forEach(a => { a.total_events = null; }); // keep NAR, drop events
  // QP: curve + registry event count (current best)
  const qpE = buildArm(expT, 8, 0, null), qpC = buildArm(ctlT, 8, 0, null);
  const tQp = mkTrial(cfg.ds, qpE, qpC); tQp.arms.forEach(a => { a.nar_points = []; });

  const fe = (h) => isFinite(h) ? +Math.exp(Math.abs(Math.log(h) - Math.log(trueHR))).toFixed(3) : null;
  const curveOnly = hrOf(tCu, {});                       // no events, no NAR → guyot/anchor-exact best-of
  const fusion = hrOf(tFu, { method: 'anchor-exact' });  // curve + sparse OCR'd NAR via anchor-exact
  const fusionQP = coxHR(narAwareReconstruct(tFu.arms[0]), narAwareReconstruct(tFu.arms[1])); // curve + sparse NAR via NAR-aware QP
  const qp = hrOf(tQp, {});                               // curve + events → QP default
  rows.push({ ds: cfg.ds, tcga: cfg.ds.startsWith('cbio'), n_exp: expT.length, n_ctl: ctlT.length, true_HR: +trueHR.toFixed(3),
    curve_only: { HR: +curveOnly.toFixed(3), fold: fe(curveOnly) },
    fusion_nar_anchorexact: { HR: +fusion.toFixed(3), fold: fe(fusion) },
    fusion_nar_qp: { HR: +fusionQP.toFixed(3), fold: fe(fusionQP) },
    qp_curve_plus_events: { HR: +qp.toFixed(3), fold: fe(qp) } });
}

function agg(set, key) {
  const xs = set.map(r => r[key].fold).filter(x => x != null).sort((a, b) => a - b);
  return { median_fold: +xs[xs.length >> 1].toFixed(3), within20: set.filter(r => r[key].fold != null && r[key].fold < 1.2).length + '/' + set.length };
}
const KEYS = ['curve_only', 'fusion_nar_anchorexact', 'fusion_nar_qp', 'qp_curve_plus_events'];
const aggAll = (set) => { const o = {}; KEYS.forEach(k => o[k] = agg(set, k)); return o; };
const big = rows.filter(r => r.n_exp >= 100 && r.n_ctl >= 100);
const heavy = rows.filter(r => r.tcga);   // the heavily-censored TCGA cohorts where curve-only fails
const out = { summary: {
  n: rows.length, all: aggAll(rows), ge100: aggAll(big), heavily_censored_tcga: aggAll(heavy),
  note: 'Fusion = registry-exact curve + numbers-at-risk (figure risk table), NO registry event count. '
    + 'sparse = ~4 OCR-noised NAR columns (realistic risk table); dense = NAR at every anchor, clean '
    + '(upper bound). QP = curve + registry event count. If fusion approaches QP, the figure NAR '
    + 'substitutes for the missing registry event count.',
}, per_dataset: rows };
fs.writeFileSync(path.join(GS.dir, 'nar_fusion_results.json'), JSON.stringify(out, null, 2));
console.log('NAR fusion (' + rows.length + ' datasets) — median HR fold-error (within-20%):');
console.log('                   curve-only      NAR+anchorexact  NAR+QP(no evts)  QP(+events)');
for (const [lab, s] of [['all', out.summary.all], ['>=100/arm', out.summary.ge100], ['heavily-cens TCGA', out.summary.heavily_censored_tcga]])
  console.log('  ' + lab.padEnd(18) + KEYS.map(k => (String(s[k].median_fold) + ' (' + s[k].within20 + ')').padStart(15)).join(' '));
