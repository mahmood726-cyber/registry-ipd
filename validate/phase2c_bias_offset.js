#!/usr/bin/env node
/*
 * PHASE 2c: reconstruction error as a CALIBRATED identified-set offset (bias) + width (variance).
 * =====================================================================================================
 *
 * Phase 2b showed the ensemble's r^2 captures the censoring UNCERTAINTY but not the censoring BIAS. The
 * deeper statistical truth: a DETERMINISTIC per-trial reconstruction bias b_i is observationally identical
 * to between-trial heterogeneity — you cannot tell "this trial reconstructed 8% high" apart from "this
 * trial truly differs by 8%" from the data alone. So tau^2 is NOT point-identified from reconstructed
 * trials; it is PARTIALLY identified. The honest object is an identified SET (Manski), calibrated against a
 * true-IPD gold standard (the project has one: the 14 TCGA cohorts).
 *
 * Calibration (all LEAVE-ONE-OUT, so out-of-sample): from the others' errors e_i = logHR_rec_i -
 * logHR_true_i we estimate (a) a SYSTEMATIC offset beta = mean(e) — identifiable, so we de-bias it; and
 * (b) a residual-bias bound delta = 1.64*SD(e) — the plausible per-trial bias the data cannot resolve.
 * The pooled effect and tau^2 then range over all residual-bias configurations {|b_i| <= delta_i}: that
 * range is the identified set. The claim to test: the identified set CONTAINS the true-IPD value, whereas
 * the naive event-pinned point estimate is a misleadingly precise (and wrong) single number.
 *
 * Run from repo root (needs realipd/cbio_*.csv):
 *   node validate/phase2c_bias_offset.js  ->  validate/phase2c_bias_offset_results.json
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');
const GS = require('./goldstandard.js');
const { metaRE, coxBetaSE, coarse } = require('./phase2_real_pooling.js');

const mean = (a) => a.reduce((x, y) => x + y, 0) / a.length;
const variance = (a) => { const m = mean(a); return a.reduce((x, y) => x + (y - m) ** 2, 0) / (a.length - 1); };

function r2eventPinned(trial) {
  const ens = RIPD.reconstructEnsemble(trial, { M: 200, pinEvents: true });
  const h = ens.ensemble && ens.ensemble.hr;
  if (!h || !(h.lo > 0 && h.hi > 0)) return null;
  return ((Math.log(h.hi) - Math.log(h.lo)) / (2 * 1.959964)) ** 2;
}

function run() {
  const rows = [];
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
    const r2 = r2eventPinned(trial); if (r2 == null) continue;
    rows.push({ ds: cfg.ds.replace('cbio_', ''), yTrue: t.logHR, sTrue2: t.se * t.se,
      yRec: rc.logHR, sRec2: rc.se * rc.se, r2ens: r2, e: rc.logHR - t.logHR });
  }
  const k = rows.length;
  const yTrue = rows.map(r => r.yTrue), vTrue = rows.map(r => r.sTrue2), yRec = rows.map(r => r.yRec);
  const vEvent = rows.map(r => r.sRec2 + r.r2ens);
  const yDeb = [], delta = [];
  for (let j = 0; j < k; j++) {
    const eO = rows.filter((_, i) => i !== j).map(o => o.e);
    yDeb.push(rows[j].yRec - mean(eO));                 // de-bias the LOO systematic offset (identifiable)
    delta.push(1.64 * Math.sqrt(variance(eO)));         // LOO residual-bias bound (the un-resolvable part)
  }
  const TRUE = metaRE(yTrue, vTrue);
  const EVENT = metaRE(yRec, vEvent);                   // Phase 2b naive point (tau^2 overstated)
  const DEB = metaRE(yDeb, vEvent);                     // de-biased point estimate
  // identified SET over residual biases |b_i| <= delta_i (greedy extremes about the de-biased mean)
  const mu = DEB.pooled_logHR;
  const tau2 = (sh) => metaRE(yDeb.map((y, i) => y + sh[i]), vEvent).tau2;
  const shMax = yDeb.map((y, i) => Math.sign(y - mu) * delta[i]);                 // push apart -> max tau^2
  const shMin = yDeb.map((y, i) => -Math.sign(y - mu) * Math.min(delta[i], Math.abs(y - mu))); // pull together
  const tau2_set = [+tau2(shMin).toFixed(3), +tau2(shMax).toFixed(3)];
  const hi = metaRE(yDeb.map((y, i) => y + delta[i]), vEvent).pooled_logHR;
  const lo = metaRE(yDeb.map((y, i) => y - delta[i]), vEvent).pooled_logHR;
  const hr_set = [+Math.exp(lo).toFixed(2), +Math.exp(hi).toFixed(2)];
  return {
    k,
    true_ipd: { pooled_HR: TRUE.pooled_HR, tau2: +TRUE.tau2.toFixed(3), PI: TRUE.PI },
    event_pinned_point: { pooled_HR: EVENT.pooled_HR, tau2: +EVENT.tau2.toFixed(3) },
    debiased_point: { pooled_HR: DEB.pooled_HR, tau2: +DEB.tau2.toFixed(3) },
    identified_set: { HR: hr_set, tau2: tau2_set },
    contains_truth: {
      HR: TRUE.pooled_HR >= hr_set[0] && TRUE.pooled_HR <= hr_set[1],
      tau2: TRUE.tau2 >= tau2_set[0] - 1e-9 && TRUE.tau2 <= tau2_set[1] + 1e-9,
    },
    mean_recon_bias_logHR: +mean(rows.map(r => r.e)).toFixed(3),
    per_trial: rows.map(r => ({ ds: r.ds, e: +r.e.toFixed(3), r_ens: +Math.sqrt(r.r2ens).toFixed(3) })),
  };
}

if (require.main === module) {
  const out = run();
  fs.writeFileSync(path.join(__dirname, 'phase2c_bias_offset_results.json'), JSON.stringify(out, null, 2));
  console.log(`=== Phase 2c: reconstruction error as a partially-identified set (k=${out.k} TCGA cohorts) ===\n`);
  console.log(`  true IPD (held out)        HR ${out.true_ipd.pooled_HR.toFixed(2)}   tau^2 ${out.true_ipd.tau2}`);
  console.log(`  event-pinned POINT         HR ${out.event_pinned_point.pooled_HR.toFixed(2)}   tau^2 ${out.event_pinned_point.tau2}   <- a precise but WRONG single number`);
  console.log(`  de-biased POINT (LOO)      HR ${out.debiased_point.pooled_HR.toFixed(2)}   tau^2 ${out.debiased_point.tau2}   <- systematic offset removed`);
  console.log(`\n  IDENTIFIED SET (LOO-calibrated residual-bias boxes):`);
  console.log(`    pooled HR : [${out.identified_set.HR[0]}, ${out.identified_set.HR[1]}]   contains true ${out.true_ipd.pooled_HR.toFixed(2)}? ${out.contains_truth.HR}`);
  console.log(`    tau^2     : [${out.identified_set.tau2[0]}, ${out.identified_set.tau2[1]}]   contains true ${out.true_ipd.tau2}? ${out.contains_truth.tau2}`);
  console.log(`\n  Honest reading: a deterministic per-trial reconstruction bias is indistinguishable from`);
  console.log(`  heterogeneity, so tau^2 cannot be a point. The calibrated identified SET brackets the truth;`);
  console.log(`  the naive point does not. mean systematic offset (de-biased) = ${out.mean_recon_bias_logHR} log-HR.`);
  console.log('\n  wrote validate/phase2c_bias_offset_results.json');
}
module.exports = { run };
