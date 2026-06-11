#!/usr/bin/env node
/*
 * PHASE 3b (step 1): export per-trial reconstruction IMPUTATIONS as specs for spec-collapse's
 * weighted_likelihood — the principled within-trial aggregator.
 * =====================================================================================================
 *
 * The M reconstructions of ONE curve (sampling the under-identified censoring level) are a "multiverse of
 * one dataset": correlated specs, NOT M independent studies. Aggregating them by inverse-variance pooling
 * collapses the variance by ~M (the cardinal sin advanced-stats.md / spec-collapse name). The honest
 * aggregator is spec-collapse-atlas's weighted_likelihood (Gaussian/t mixture, never narrower than one
 * draw). This exports, per cohort, the M imputation specs {theta=logHR, var=se^2, k} plus the held-out true
 * logHR, so the Python side (phase3b_weighted_likelihood.py) can pool them honestly and check coverage.
 *
 * Run from repo root (needs realipd/cbio_*.csv):
 *   node validate/phase3b_export_imputations.js  ->  validate/phase3b_imputations.json
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');
const GS = require('./goldstandard.js');
const { coxBetaSE, coarse } = require('./phase2_real_pooling.js');

function mulberry32(a) {
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

const M = 120;
// a few representative cohorts: heavily-censored large-effect + a cleaner one
const PICK = ['cbio_kirc', 'cbio_coadread', 'cbio_luad', 'cbio_hnsc'];

function imputationsFor(cfg) {
  const arms = GS.loadArms(cfg);
  const { expT, ctlT } = arms;
  const tTrue = coxBetaSE(expT, ctlT);
  const trial = { nct_id: cfg.ds, time_unit: 'd', arms: [
    Object.assign({ arm_id: 'exp', role: 'experimental' }, coarse(expT, 8)),
    Object.assign({ arm_id: 'ctl', role: 'comparator' }, coarse(ctlT, 8))] };
  // curve-only event ceiling E0 per arm (no intermediate censoring)
  const e0 = {};
  const cr = RIPD.reconstruct(trial, { ignoreTotalEvents: true });
  cr.arms.forEach(a => { e0[a.arm_id] = a.ipd.filter(x => x.status === 1).length; });

  const specs = [];
  let seed = 777;
  for (let m = 0; m < M; m++) {
    seed = (seed + 0x6D2B79F5) >>> 0;
    const rng = mulberry32(seed);
    const t2 = JSON.parse(JSON.stringify(trial));
    for (const a of t2.arms) {
      a.km_points = a.km_points.map(p => ({ t: p.t, S: Math.min(1, Math.max(0, p.S + (rng() - 0.5) * 0.01)) }));
      const E0 = e0[a.arm_id];
      if (E0 > 0) {                                     // sample the censoring level across the full band
        const lo = Math.round(0.55 * E0);
        a.total_events = Math.min(a.N, Math.max(1, Math.round(lo + rng() * (E0 - lo))));
      }
    }
    const method = rng() < 0.5 ? 'guyot' : 'anchor-exact';
    let r; try { r = RIPD.reconstruct(t2, { method }); } catch { continue; }
    if (!r.arms) continue;
    const rc = coxBetaSE(r.arms.find(a => a.role === 'experimental').ipd,
      r.arms.find(a => a.role === 'comparator').ipd);
    if (rc.se == null || !isFinite(rc.logHR)) continue;
    const events = r.arms.reduce((s, a) => s + a.ipd.filter(x => x.status === 1).length, 0);
    specs.push({ theta: +rc.logHR.toFixed(5), var: +(rc.se * rc.se).toFixed(6), k: events,
      significant: (rc.logHR - 1.96 * rc.se > 0 || rc.logHR + 1.96 * rc.se < 0) });
  }
  return { ds: cfg.ds.replace('cbio_', ''), true_logHR: +tTrue.logHR.toFixed(5),
    true_HR: +Math.exp(tTrue.logHR).toFixed(3), true_se: +tTrue.se.toFixed(4),
    n_imputations: specs.length, specs };
}

const out = {};
for (const ds of PICK) {
  const cfg = GS.CONFIGS.find(c => c.ds === ds);
  if (!cfg) continue;
  try { out[ds] = imputationsFor(cfg); } catch (e) { out[ds] = { error: String(e) }; }
}
fs.writeFileSync(path.join(__dirname, 'phase3b_imputations.json'), JSON.stringify(out, null, 2));
const done = Object.keys(out).filter(k => out[k].specs);
console.log(`wrote validate/phase3b_imputations.json — ${done.length} cohorts, M=${M} imputations each`);
for (const k of done) console.log(`  ${out[k].ds}: true HR ${out[k].true_HR}, ${out[k].n_imputations} specs`);
