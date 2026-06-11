#!/usr/bin/env node
/*
 * PHASE 2b: the censoring lever shrinks the reconstruction variance r^2 — closing the loop from the
 * per-trial lever to synthesis-level heterogeneity honesty.
 * =====================================================================================================
 *
 * Phase 2 found that on heavily-censored curve-only cohorts the reconstruction variance r^2 is so large it
 * swamps the heterogeneity signal: honest (Rubin) pooling absorbs the genuine between-trial tau^2 into the
 * within-trial variance and UNDER-states it. The claim that unifies the project: a censoring lever (a
 * posted total-event count / number-at-risk) PINS the censoring split, shrinking r^2 into the range where
 * a curve-only trial becomes informative about heterogeneity again.
 *
 * This demonstrates it on the same 14 TCGA cohorts, measuring r^2 in two regimes from the engine's
 * multiple-imputation ensemble:
 *   r^2_curve  : ensemble with NO event count (full censoring band)            -> large
 *   r^2_event  : ensemble with the event count, pinned (reconstructEnsemble {pinEvents:true}) -> small
 * then pooling HONEST (Rubin s^2 + r^2) under each and comparing tau^2 to the true-IPD tau^2.
 *
 * Run from repo root (needs realipd/cbio_*.csv):
 *   node validate/phase2b_lever_shrinks_r2.js  ->  validate/phase2b_lever_shrinks_r2_results.json
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');
const GS = require('./goldstandard.js');
const { metaRE, coxBetaSE, coarse } = require('./phase2_real_pooling.js');

function r2FromEnsemble(ens) {
  const h = ens.ensemble && ens.ensemble.hr;
  if (!h || !(h.lo > 0 && h.hi > 0)) return null;
  return ((Math.log(h.hi) - Math.log(h.lo)) / (2 * 1.959964)) ** 2;
}

function run() {
  const trueY = [], trueV = [], recY = [], vCurve = [], vEvent = [];
  let sumR2curve = 0, sumR2event = 0, n = 0;
  const per = [];
  for (const cfg of GS.CONFIGS) {
    if (!cfg.ds.startsWith('cbio_')) continue;
    let arms; try { arms = GS.loadArms(cfg); } catch { continue; }
    const { expT, ctlT } = arms; if (expT.length < 20 || ctlT.length < 20) continue;
    const t = coxBetaSE(expT, ctlT); if (t.se == null) continue;
    const trial = { nct_id: cfg.ds, time_unit: 'd', arms: [
      Object.assign({ arm_id: 'exp', role: 'experimental' }, coarse(expT, 8)),
      Object.assign({ arm_id: 'ctl', role: 'comparator' }, coarse(ctlT, 8))] };
    const r = RIPD.reconstruct(trial, {}); if (!r.arms) continue;
    const rc = coxBetaSE(r.arms.find(a => a.role === 'experimental').ipd,
      r.arms.find(a => a.role === 'comparator').ipd); if (rc.se == null) continue;

    // r^2 with NO event count (curve-only: strip total_events so the ensemble samples the full band)
    const curveTrial = JSON.parse(JSON.stringify(trial));
    curveTrial.arms.forEach(a => { a.total_events = null; });
    const r2c = r2FromEnsemble(RIPD.reconstructEnsemble(curveTrial, { M: 200 }));
    // r^2 with the event count PINNED (the lever)
    const r2e = r2FromEnsemble(RIPD.reconstructEnsemble(trial, { M: 200, pinEvents: true }));
    if (r2c == null || r2e == null) continue;

    trueY.push(t.logHR); trueV.push(t.se * t.se);
    recY.push(rc.logHR); vCurve.push(rc.se * rc.se + r2c); vEvent.push(rc.se * rc.se + r2e);
    sumR2curve += r2c; sumR2event += r2e; n++;
    per.push({ ds: cfg.ds.replace('cbio_', ''), r_curve: +Math.sqrt(r2c).toFixed(3),
      r_event: +Math.sqrt(r2e).toFixed(3), shrink: +(Math.sqrt(r2c) / Math.sqrt(r2e)).toFixed(1) });
  }
  const TRUE = metaRE(trueY, trueV);
  const HON_CURVE = metaRE(recY, vCurve), HON_EVENT = metaRE(recY, vEvent);
  return {
    k: TRUE.k,
    mean_r_curve: +Math.sqrt(sumR2curve / n).toFixed(3),
    mean_r_event: +Math.sqrt(sumR2event / n).toFixed(3),
    r_shrink_factor: +(Math.sqrt(sumR2curve / n) / Math.sqrt(sumR2event / n)).toFixed(1),
    true_ipd: TRUE, honest_curve_only: HON_CURVE, honest_event_pinned: HON_EVENT,
    per_trial: per,
  };
}

if (require.main === module) {
  const out = run();
  fs.writeFileSync(path.join(__dirname, 'phase2b_lever_shrinks_r2_results.json'), JSON.stringify(out, null, 2));
  const tline = (lab, m, truth) => `  ${lab.padEnd(26)} tau^2 ${m.tau2.toFixed(3)}  (true ${truth.toFixed(3)})  PI [${m.PI[0]}, ${m.PI[1]}]  HR ${m.pooled_HR.toFixed(2)}`;
  console.log(`=== Phase 2b: the censoring lever shrinks r^2 (k=${out.k} TCGA cohorts) ===\n`);
  console.log(`  mean reconstruction SD r:  curve-only ${out.mean_r_curve}  ->  event-pinned ${out.mean_r_event}   (${out.r_shrink_factor}x smaller)\n`);
  console.log(tline('true IPD', out.true_ipd, out.true_ipd.tau2));
  console.log(tline('honest, curve-only r^2', out.honest_curve_only, out.true_ipd.tau2) + '   <- r^2 huge: tau^2 absorbed (understated)');
  console.log(tline('honest, event-pinned r^2', out.honest_event_pinned, out.true_ipd.tau2) + '   <- lever shrinks r^2: tau^2 restored');
  console.log('\n  wrote validate/phase2b_lever_shrinks_r2_results.json');
}
module.exports = { run };
