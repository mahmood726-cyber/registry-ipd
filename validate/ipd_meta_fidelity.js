#!/usr/bin/env node
/*
 * IPD META-ANALYSIS FIDELITY — does reconstruction error propagate into the POOLED estimate?
 *
 * The method's purpose is to enable IPD meta-analysis. We have only ever validated per-trial HR
 * recovery; this asks the use-case question: if you pool reconstructed pseudo-IPD across trials, do you
 * get the same meta-analytic answer as pooling the TRUE IPD? Reconstruction error could (a) wash out in
 * pooling, or (b) bias/inflate the pooled estimate and heterogeneity.
 *
 * Test set: the 14 TCGA stage cohorts (a real meta-analytic question — the pooled prognostic effect of
 * advanced vs early stage across cancers, with genuine between-cancer heterogeneity). For each we take
 * the true Cox (logHR, SE) and the QP-reconstructed Cox (logHR, SE), then run a proper random-effects
 * meta-analysis on BOTH and compare pooled HR, τ², 95% CI and prediction interval.
 *
 * Methods (per advanced-stats rules): REML τ² (k≥10), HKSJ variance with floor max(1, q), CI/PI on the
 * Student t_{k-1} distribution, pooling on the log-HR scale and back-transforming.
 *
 * Deterministic. Run: node validate/ipd_meta_fidelity.js  ->  realipd/ipd_meta_fidelity_results.json
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');
const _ = RIPD._;
const GS = require('./goldstandard.js');

// Cox logHR + SE (1/sqrt information) for a single binary covariate (exp=1, ctl=0)
function coxBetaSE(a, b) {
  const rows = a.map(r => ({ time: r.time, status: r.status, x: 1 })).concat(b.map(r => ({ time: r.time, status: r.status, x: 0 })));
  const beta = _.coxLogHR(rows.slice()).beta;
  const evTimes = [...new Set(rows.filter(r => r.status === 1).map(r => r.time))].sort((x, y) => x - y);
  let I = 0;
  for (const et of evTimes) {
    const risk = rows.filter(r => r.time >= et - 1e-12), ev = rows.filter(r => r.status === 1 && Math.abs(r.time - et) < 1e-12);
    let S0 = 0, S1 = 0; for (const r of risk) { const e = Math.exp(beta * r.x); S0 += e; S1 += r.x * e; }
    const p = S0 > 0 ? S1 / S0 : 0; I += ev.length * p * (1 - p);
  }
  const se = I > 1e-9 ? 1 / Math.sqrt(I) : null;
  return { logHR: beta, se };
}

// Student-t 0.975 quantile via Cornish–Fisher (3 terms) — accurate to <0.5% for df≥5
function qt975(df) {
  const z = 1.959963985; const g1 = (z * z * z + z) / 4, g2 = (5 * z ** 5 + 16 * z ** 3 + 3 * z) / 96, g3 = (3 * z ** 7 + 19 * z ** 5 + 17 * z ** 3 - 15 * z) / 384;
  return z + g1 / df + g2 / (df * df) + g3 / (df ** 3);
}

// random-effects meta-analysis: REML τ², HKSJ CI (floored), prediction interval (t_{k-1})
function metaRE(y, v) {
  const k = y.length;
  // REML iteration for τ²
  let tau2 = 0;
  for (let it = 0; it < 200; it++) {
    const w = v.map(vi => 1 / (vi + tau2)), sw = w.reduce((a, b) => a + b, 0);
    const mu = y.reduce((a, yi, i) => a + w[i] * yi, 0) / sw;
    let num = 0, den = 0; for (let i = 0; i < k; i++) { num += w[i] * w[i] * ((y[i] - mu) ** 2 - v[i]); den += w[i] * w[i]; }
    const t2 = Math.max(0, num / den + 1 / sw);
    if (Math.abs(t2 - tau2) < 1e-9) { tau2 = t2; break; } tau2 = t2;
  }
  const w = v.map(vi => 1 / (vi + tau2)), sw = w.reduce((a, b) => a + b, 0);
  const mu = y.reduce((a, yi, i) => a + w[i] * yi, 0) / sw;
  // Cochran Q and I²
  const wFE = v.map(vi => 1 / vi), swFE = wFE.reduce((a, b) => a + b, 0), muFE = y.reduce((a, yi, i) => a + wFE[i] * yi, 0) / swFE;
  const Q = y.reduce((a, yi, i) => a + wFE[i] * (yi - muFE) ** 2, 0); const I2 = Math.max(0, (Q - (k - 1)) / Q) * 100;
  // HKSJ variance with floor (prevents narrowing below standard RE when q<1)
  const q = Math.max(1, y.reduce((a, yi, i) => a + w[i] * (yi - mu) ** 2, 0) / (k - 1));
  const seHK = Math.sqrt(q / sw), tcrit = qt975(k - 1);
  const ci = [mu - tcrit * seHK, mu + tcrit * seHK];
  const pi = [mu - tcrit * Math.sqrt(tau2 + seHK * seHK), mu + tcrit * Math.sqrt(tau2 + seHK * seHK)];
  return { k, pooled_logHR: mu, tau2, I2: +I2.toFixed(1), Q: +Q.toFixed(2), seHKSJ: seHK,
    pooled_HR: +Math.exp(mu).toFixed(3), CI: ci.map(x => +Math.exp(x).toFixed(3)), PI: pi.map(x => +Math.exp(x).toFixed(3)) };
}

function coarse(ipd, K) {
  const km = _.kmFromIPD(ipd), tmax = 0.95 * Math.max(...ipd.map(r => r.time)), pts = [{ t: 0, S: 1 }];
  for (let i = 1; i <= K; i++) { const t = tmax * i / K; pts.push({ t: +t.toFixed(2), S: +_.evalKM(km, t).toFixed(4) }); }
  return { km_points: pts, nar_points: [], N: ipd.length, total_events: ipd.filter(r => r.status === 1).length, follow_up_max: +tmax.toFixed(2) };
}

const trueY = [], trueV = [], recY = [], recV = [], per = [];
for (const cfg of GS.CONFIGS) {
  if (!cfg.ds.startsWith('cbio_')) continue;                 // the 14 TCGA stage cohorts
  let arms; try { arms = GS.loadArms(cfg); } catch { continue; }
  const { expT, ctlT } = arms; if (expT.length < 20 || ctlT.length < 20) continue;
  const t = coxBetaSE(expT, ctlT);
  const trial = { nct_id: cfg.ds, time_unit: 'd', arms: [Object.assign({ arm_id: 'exp', role: 'experimental' }, coarse(expT, 8)), Object.assign({ arm_id: 'ctl', role: 'comparator' }, coarse(ctlT, 8))] };
  const r = RIPD.reconstruct(trial, {}); if (!r.arms || t.se == null) continue;
  const rc = coxBetaSE(r.arms.find(a => a.role === 'experimental').ipd, r.arms.find(a => a.role === 'comparator').ipd);
  if (rc.se == null) continue;
  trueY.push(t.logHR); trueV.push(t.se * t.se); recY.push(rc.logHR); recV.push(rc.se * rc.se);
  per.push({ ds: cfg.ds.replace('cbio_', ''), true_HR: +Math.exp(t.logHR).toFixed(2), true_se: +t.se.toFixed(3), recon_HR: +Math.exp(rc.logHR).toFixed(2), recon_se: +rc.se.toFixed(3) });
}

const trueMA = metaRE(trueY, trueV), reconMA = metaRE(recY, recV);
const out = {
  summary: {
    k: trueMA.k, comparison: 'late vs early stage, overall survival, across TCGA cancers',
    true_IPD: { pooled_HR: trueMA.pooled_HR, CI_HKSJ: trueMA.CI, PI: trueMA.PI, tau2: +trueMA.tau2.toFixed(3), I2: trueMA.I2 },
    reconstructed_QP: { pooled_HR: reconMA.pooled_HR, CI_HKSJ: reconMA.CI, PI: reconMA.PI, tau2: +reconMA.tau2.toFixed(3), I2: reconMA.I2 },
    pooled_HR_fold_diff: +Math.exp(Math.abs(trueMA.pooled_logHR - reconMA.pooled_logHR)).toFixed(3),
    tau2_abs_diff: +Math.abs(trueMA.tau2 - reconMA.tau2).toFixed(3),
    methods: 'REML τ², HKSJ CI (floored), prediction interval on t_{k-1}, log-HR pooling. k=' + trueMA.k + '.',
  },
  per_trial: per,
};
fs.writeFileSync(path.join(GS.dir, 'ipd_meta_fidelity_results.json'), JSON.stringify(out, null, 2));
console.log('IPD meta-analysis fidelity (k=' + trueMA.k + ' TCGA cohorts, late vs early stage OS):');
console.log('                  pooled HR     95% CI (HKSJ)          95% PI               τ²      I²');
const fmt = (m) => `${m.pooled_HR.toFixed(2).padStart(5)}   [${m.CI[0]}, ${m.CI[1]}]`.padEnd(28) + `[${m.PI[0]}, ${m.PI[1]}]`.padEnd(20) + `  ${m.tau2.toFixed(3)}  ${m.I2}%`;
console.log('  true IPD       ', fmt(trueMA));
console.log('  reconstructed  ', fmt(reconMA));
console.log('\n  pooled-HR fold difference:', out.summary.pooled_HR_fold_diff, '| τ² abs diff:', out.summary.tau2_abs_diff);
