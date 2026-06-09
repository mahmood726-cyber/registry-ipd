#!/usr/bin/env node
/*
 * Robust-estimand validation: median survival and RMST. These are CURVE-DERIVED quantities, so
 * (unlike the HR) they should be recovered far more faithfully — this test quantifies that and is
 * not circular against the registry HR.
 *
 *  1. MEDIAN accuracy  — reconstructed arm median vs the registry-reported median (external ground
 *     truth, only for arms where the trial posts a median measure).
 *  2. RMST fidelity     — reconstructed-IPD RMST vs the RMST of the registry anchor curve itself,
 *     to a common horizon tau. Measures whether the pseudo-IPD preserves the curve area. Available
 *     for every Tier-A arm.
 *
 * Reconstructs CURVE-ONLY (the default method) via ignoreTotalEvents.
 * Usage: node validate/validate_rmst.js [cohort_dir]
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');
const _ = RIPD._;

const dir = process.argv[2] && !process.argv[2].startsWith('-') ? process.argv[2] : path.join(__dirname, '..', 'cohort');
function median(a) { if (!a.length) return null; const s = a.slice().sort((x, y) => x - y); const m = s.length >> 1; return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; }
function quant(a, p) { if (!a.length) return null; const s = a.slice().sort((x, y) => x - y); const i = p * (s.length - 1), lo = Math.floor(i), hi = Math.ceil(i); return lo === hi ? s[lo] : s[lo] + (s[hi] - s[lo]) * (i - lo); }

// RMST of the registry anchor curve (right-continuous step on the posted (t,S) points) to tau
function anchorRMST(km_points, tau) {
  const p = km_points.slice().sort((a, b) => a.t - b.t);
  let area = 0, prevT = 0, prevS = 1;
  for (const q of p) {
    const t = Math.min(q.t, tau);
    if (t > prevT) area += prevS * (t - prevT);
    prevT = q.t; prevS = q.S;
    if (q.t >= tau) return area;
  }
  if (tau > prevT) area += prevS * (tau - prevT);
  return area;
}

const files = fs.readdirSync(dir).filter(f => f.endsWith('.json') && !f.startsWith('manifest') && !f.startsWith('validation') && !f.startsWith('registry'));
// external registry medians per arm (ctgov_group_code), harvested separately
let regMed = {};
try { regMed = JSON.parse(fs.readFileSync(path.join(dir, 'registry_medians.json'), 'utf8')); } catch { /* optional */ }
const medErr = [], rmstErr = [], medFidErr = [];
let arms = 0, medN = 0, medFidN = 0;

// median of the registry ANCHOR curve itself (0.5-crossing of the posted step) — same endpoint,
// uncontaminated; tests whether the pseudo-IPD reproduces the curve's own median.
function anchorMedian(km_points) {
  const p = km_points.slice().sort((a, b) => a.t - b.t);
  for (const q of p) if (q.S <= 0.5 + 1e-9) return q.t;
  return null;
}

for (const f of files) {
  let t; try { t = JSON.parse(fs.readFileSync(path.join(dir, f), 'utf8')); } catch { continue; }
  let r; try { r = RIPD.reconstruct(t, { ignoreTotalEvents: true }); } catch { continue; }
  if (r.tier !== 'A' || !r.arms) continue;
  // common horizon = min over arms of last anchor time
  const taus = t.arms.filter(a => a.km_points && a.km_points.length).map(a => a.km_points[a.km_points.length - 1].t);
  if (!taus.length) continue;
  const tau = Math.min(...taus);
  for (const a of t.arms) {
    const rec = r.arms.find(x => x.arm_id === a.arm_id);
    if (!rec || !a.km_points || !a.km_points.length) continue;
    arms++;
    const km = _.kmFromIPD(rec.ipd);
    // RMST fidelity
    const rRec = _.rmst(km, tau), rAnc = anchorRMST(a.km_points, tau);
    if (rAnc > 0) rmstErr.push(Math.abs(rRec - rAnc) / rAnc);
    // median-from-curve fidelity (uncontaminated): recon median vs anchor-curve median
    const aMed = anchorMedian(a.km_points);
    if (aMed != null && aMed > 0) {
      const mr = _.medianFromKM(km);
      if (mr != null) { medFidErr.push(Math.abs(mr - aMed) / aMed); medFidN++; }
    }
    // median accuracy vs registry-reported median (from the trial's own median measure, external)
    const regM = (regMed[t.nct_id] && regMed[t.nct_id][a.arm_id]) != null ? regMed[t.nct_id][a.arm_id]
      : (a.median && a.median.value > 0 ? a.median.value : null);
    if (regM != null && regM > 0 && regM <= tau) {   // only where the median is within follow-up
      const mRec = _.medianFromKM(km);
      if (mRec != null) { medErr.push(Math.abs(mRec - regM) / regM); medN++; }
    }
  }
}

const report = {
  cohort_dir: dir, arms_evaluated: arms,
  RMST_fidelity: {
    n: rmstErr.length,
    median_pct_err: rmstErr.length ? +(100 * median(rmstErr)).toFixed(2) : null,
    p90_pct_err: rmstErr.length ? +(100 * quant(rmstErr, 0.9)).toFixed(2) : null,
    pct_within_5pct: rmstErr.length ? Math.round(100 * rmstErr.filter(e => e <= 0.05).length / rmstErr.length) : null,
  },
  median_from_curve_fidelity: {
    n: medFidN,
    median_pct_err: medFidErr.length ? +(100 * median(medFidErr)).toFixed(2) : null,
    p90_pct_err: medFidErr.length ? +(100 * quant(medFidErr, 0.9)).toFixed(2) : null,
    pct_within_10pct: medFidErr.length ? Math.round(100 * medFidErr.filter(e => e <= 0.10).length / medFidErr.length) : null,
    note: 'recon median vs the registry curve\'s own 0.5-crossing (same endpoint, uncontaminated)',
  },
  median_accuracy_vs_registry: {
    n_arms_with_registry_median: medN,
    median_pct_err: medErr.length ? +(100 * median(medErr)).toFixed(2) : null,
    p90_pct_err: medErr.length ? +(100 * quant(medErr, 0.9)).toFixed(2) : null,
    pct_within_10pct: medErr.length ? Math.round(100 * medErr.filter(e => e <= 0.10).length / medErr.length) : null,
    note: 'vs registry-reported median — only trustworthy after survival-endpoint filtering (see harvest_medians)',
  },
  interpretation: 'RMST/median are curve-derived, so errors should be far below the HR errors '
    + '(median ~8-12%). RMST fidelity ~0 confirms the pseudo-IPD preserves the registry curve area.',
};
fs.writeFileSync(path.join(dir, 'validation_rmst_report.json'), JSON.stringify(report, null, 2));
console.log(JSON.stringify(report, null, 2));
