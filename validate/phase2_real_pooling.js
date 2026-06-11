#!/usr/bin/env node
/*
 * PHASE 2: the honest-pooling linchpin on REAL reconstructions (not a simulation).
 * ================================================================================
 *
 * honest_pooling_sim.js proved, in Monte-Carlo, that pooling reconstructed effects on sampling variance
 * alone mis-reads reconstruction noise as between-trial heterogeneity (tau^2 inflated, PI too wide), and
 * that propagating the reconstruction variance via Rubin's rules recovers it. This runs the SAME three-way
 * comparison on the 14 TCGA stage cohorts where we hold the true IPD, using the engine's own
 * multiple-imputation ensemble (`reconstructEnsemble`) to measure each trial's REAL reconstruction
 * (censoring-level) variance r_i^2 = ((ln hi - ln lo)/(2*1.96))^2 from its HR credible interval.
 *
 * Per cohort i:
 *   true IPD            : (logHR_i, s_i^2)             from Cox on the true patient data
 *   reconstructed point : (logHR_i^rec, s_i^rec 2)     from Cox on the central reconstruction
 *   reconstruction var  : r_i^2                        from the M-imputation HR credible interval
 * Pool three ways (REML tau^2, HKSJ, PI on t_{k-1}):
 *   TRUE  = metaRE(true logHR, s^2)
 *   NAIVE = metaRE(recon logHR, s_rec^2)               <- treats pseudo-IPD as exact
 *   HONEST= metaRE(recon logHR, s_rec^2 + r^2)         <- Rubin total variance
 *
 * Helpers (coxBetaSE/metaRE/coarse) mirror ipd_meta_fidelity.js; kept local so that headline-pinned
 * script is untouched. Run from repo root (needs realipd/cbio_*.csv):
 *   node validate/phase2_real_pooling.js  ->  validate/phase2_real_pooling_results.json
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');
const _ = RIPD._;
const GS = require('./goldstandard.js');

function coxBetaSE(a, b) {
  const rows = a.map(r => ({ time: r.time, status: r.status, x: 1 }))
    .concat(b.map(r => ({ time: r.time, status: r.status, x: 0 })));
  const beta = _.coxLogHR(rows).beta;
  const evTimes = [...new Set(rows.filter(r => r.status === 1).map(r => r.time))];
  let I = 0;
  for (const et of evTimes) {
    const risk = rows.filter(r => r.time >= et - 1e-12);
    const ev = rows.filter(r => r.status === 1 && Math.abs(r.time - et) < 1e-12);
    let S0 = 0, S1 = 0;
    for (const r of risk) { const e = Math.exp(beta * r.x); S0 += e; S1 += r.x * e; }
    const p = S0 > 0 ? S1 / S0 : 0;
    I += ev.length * (p - p * p);
  }
  return { logHR: beta, se: I > 1e-9 ? 1 / Math.sqrt(I) : null };
}

function qt975(df) { const z = 1.959964; return z + (z ** 3 + z) / (4 * df) + (5 * z ** 5 + 16 * z ** 3 + 3 * z) / (96 * df * df); }

function metaRE(y, v) {
  const k = y.length;
  let tau2 = 0;
  for (let it = 0; it < 200; it++) {
    const w = v.map(vi => 1 / (vi + tau2)), sw = w.reduce((a, b) => a + b, 0);
    const mu = y.reduce((s, yi, i) => s + w[i] * yi, 0) / sw;
    let num = 0, den = 0;
    for (let i = 0; i < k; i++) { num += w[i] * w[i] * ((y[i] - mu) ** 2 - v[i]); den += w[i] * w[i]; }
    num += 1 / sw;
    const t2 = Math.max(0, num / den);
    if (Math.abs(t2 - tau2) < 1e-10) { tau2 = t2; break; }
    tau2 = t2;
  }
  const w = v.map(vi => 1 / (vi + tau2)), sw = w.reduce((a, b) => a + b, 0);
  const mu = y.reduce((s, yi, i) => s + w[i] * yi, 0) / sw;
  let q = 0; for (let i = 0; i < k; i++) q += w[i] * (y[i] - mu) ** 2;
  const seHK = Math.sqrt(Math.max(1, q / (k - 1)) / sw);
  const I2 = Math.max(0, (q - (k - 1)) / q) * 100;
  const tcrit = qt975(k - 1);
  const ci = [mu - tcrit * seHK, mu + tcrit * seHK];
  const pi = [mu - tcrit * Math.sqrt(tau2 + seHK * seHK), mu + tcrit * Math.sqrt(tau2 + seHK * seHK)];
  return { k, pooled_logHR: mu, pooled_HR: +Math.exp(mu).toFixed(3), tau2, I2: +I2.toFixed(1),
    seHK, CI: ci.map(x => +Math.exp(x).toFixed(3)), PI: pi.map(x => +Math.exp(x).toFixed(3)),
    PI_width_log: +(pi[1] - pi[0]).toFixed(3) };
}

function coarse(ipd, K) {
  const km = _.kmFromIPD(ipd), tmax = 0.95 * Math.max(...ipd.map(r => r.time)), pts = [{ t: 0, S: 1 }];
  for (let i = 1; i <= K; i++) { const t = tmax * i / K; pts.push({ t: +t.toFixed(2), S: +_.evalKM(km, t).toFixed(4) }); }
  return { km_points: pts, nar_points: [], N: ipd.length, total_events: ipd.filter(r => r.status === 1).length, follow_up_max: +tmax.toFixed(2) };
}

function run() {
  const trueY = [], trueV = [], recY = [], naiveV = [], honestV = [], per = [];
  for (const cfg of GS.CONFIGS) {
    if (!cfg.ds.startsWith('cbio_')) continue;                  // the 14 TCGA stage cohorts
    let arms; try { arms = GS.loadArms(cfg); } catch { continue; }
    const { expT, ctlT } = arms; if (expT.length < 20 || ctlT.length < 20) continue;
    const t = coxBetaSE(expT, ctlT); if (t.se == null) continue;
    const trial = { nct_id: cfg.ds, time_unit: 'd', arms: [
      Object.assign({ arm_id: 'exp', role: 'experimental' }, coarse(expT, 8)),
      Object.assign({ arm_id: 'ctl', role: 'comparator' }, coarse(ctlT, 8))] };
    const r = RIPD.reconstruct(trial, {}); if (!r.arms) continue;
    const rc = coxBetaSE(r.arms.find(a => a.role === 'experimental').ipd,
      r.arms.find(a => a.role === 'comparator').ipd); if (rc.se == null) continue;
    // REAL reconstruction variance from the engine's multiple-imputation ensemble
    const ens = RIPD.reconstructEnsemble(trial, { M: 200 });
    const h = ens.ensemble && ens.ensemble.hr;
    if (!h || !(h.lo > 0 && h.hi > 0)) continue;
    const r2 = ((Math.log(h.hi) - Math.log(h.lo)) / (2 * 1.959964)) ** 2;
    trueY.push(t.logHR); trueV.push(t.se * t.se);
    recY.push(rc.logHR); naiveV.push(rc.se * rc.se); honestV.push(rc.se * rc.se + r2);
    per.push({ ds: cfg.ds.replace('cbio_', ''), true_HR: +Math.exp(t.logHR).toFixed(2),
      recon_HR: +Math.exp(rc.logHR).toFixed(2), s_rec: +rc.se.toFixed(3), r_recon: +Math.sqrt(r2).toFixed(3),
      var_inflation: +((rc.se * rc.se + r2) / (rc.se * rc.se)).toFixed(2) });
  }
  const TRUE = metaRE(trueY, trueV), NAIVE = metaRE(recY, naiveV), HONEST = metaRE(recY, honestV);
  return { k: TRUE.k, true_ipd: TRUE, naive: NAIVE, honest: HONEST, per_trial: per };
}

if (require.main === module) {
  const out = run();
  fs.writeFileSync(path.join(__dirname, 'phase2_real_pooling_results.json'), JSON.stringify(out, null, 2));
  const fmt = (m) => `HR ${m.pooled_HR.toFixed(2)}  CI [${m.CI[0]}, ${m.CI[1]}]  PI [${m.PI[0]}, ${m.PI[1]}]  tau2 ${m.tau2.toFixed(3)}  I2 ${m.I2}%`;
  console.log(`=== Phase 2: honest pooling on ${out.k} REAL TCGA reconstructions (late vs early stage OS) ===\n`);
  console.log('  true IPD                  ', fmt(out.true_ipd));
  console.log('  NAIVE (ignore recon var)  ', fmt(out.naive));
  console.log('  HONEST (Rubin s^2 + r^2)  ', fmt(out.honest));
  const meanInfl = (out.per_trial.reduce((a, p) => a + p.var_inflation, 0) / out.per_trial.length).toFixed(2);
  console.log(`\n  mean per-trial variance inflation from reconstruction (honest/naive): ${meanInfl}x`);
  console.log(`  tau^2:  true ${out.true_ipd.tau2.toFixed(3)} | naive ${out.naive.tau2.toFixed(3)} | honest ${out.honest.tau2.toFixed(3)}`);
  console.log(`  PI width (log): true ${out.true_ipd.PI_width_log} | naive ${out.naive.PI_width_log} | honest ${out.honest.PI_width_log}`);
  console.log('\n  wrote validate/phase2_real_pooling_results.json');
}
module.exports = { run, metaRE, coxBetaSE, coarse };
