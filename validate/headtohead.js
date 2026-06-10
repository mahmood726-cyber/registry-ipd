#!/usr/bin/env node
/*
 * HEAD-TO-HEAD on REAL IPD: registry-structured anchors vs figure-digitization, both vs truth.
 *
 * The paper's central value claim is that registry-native reconstruction carries *zero figure-
 * digitization error*. This quantifies that on real data instead of asserting it. For each open-IPD
 * dataset we take the SAME true Kaplan-Meier curve and build two posted-data scenarios:
 *
 *   REGISTRY  : exact KM survival at the posted anchor timepoints (what ct.gov/AACT exposes) — no
 *               pixel error.
 *   DIGITIZED : the same curve as a *plotted figure* read by a digitizer — coordinates sampled along
 *               the curve with realistic pixel noise on both axes (Gaussian; survival sigma=1 pp,
 *               time sigma=0.5% of t_max), then monotone-enforced. This is how Guyot/IPDfromKM/etc.
 *               ingest a curve.
 *
 * Both reconstruct through the SAME engine and are scored against the true IPD (Cox HR, median, RMST).
 * Two comparisons:
 *   (A) EQUAL DENSITY (K points each): isolates the *pure digitization-error* cost — same information,
 *       one exact and one pixel-noised.
 *   (B) REALISTIC (registry few-exact K=8 vs digitized many-noisy K=25): the real-world trade-off
 *       (registries post few timepoints; a digitizer clicks many points off the figure).
 *
 * Deterministic: noise seeded per dataset via the engine's mulberry32(hashStr(ds)). No new tolerances
 * are asserted here — it is a reporting script. Run: node validate/headtohead.js
 * Writes realipd/headtohead_results.json.
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');
const _ = RIPD._;
const GS = require('./goldstandard.js');

const DIG_S_SIGMA = 0.01;     // default survival-axis digitization noise (1 percentage point)
const DIG_T_FRAC = 0.005;     // default time-axis noise as fraction of t_max (0.5%)
const K_REG = 8;              // realistic registry anchor count
const K_DIG = 25;             // realistic digitizer click count

function gauss(rng) { // Box-Muller from a uniform PRNG
  let u = 0, v = 0; while (u === 0) u = rng(); while (v === 0) v = rng();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}

// exact KM survival sampled at K evenly-spaced timepoints up to 0.95*tmax
function anchors(km, ipd, K) {
  const tmax = 0.95 * Math.max(...ipd.map(r => r.time));
  const pts = [{ t: 0, S: 1 }];
  for (let i = 1; i <= K; i++) { const t = tmax * i / K; pts.push({ t: +t.toFixed(3), S: +_.evalKM(km, t).toFixed(4) }); }
  return { pts, tmax };
}

// same timepoints, but survival & time perturbed by digitization pixel noise, then monotone-fixed
function digitize(km, ipd, K, rng, sSigma, tFrac) {
  sSigma = sSigma == null ? DIG_S_SIGMA : sSigma;
  tFrac = tFrac == null ? DIG_T_FRAC : tFrac;
  const tmax = 0.95 * Math.max(...ipd.map(r => r.time));
  const tnoise = tFrac * tmax;
  const raw = [{ t: 0, S: 1 }];
  for (let i = 1; i <= K; i++) {
    const t0 = tmax * i / K;
    const t = Math.max(raw[raw.length - 1].t + 1e-6, t0 + gauss(rng) * tnoise);
    const S = Math.min(1, Math.max(0, _.evalKM(km, t0) + gauss(rng) * sSigma));
    raw.push({ t: +t.toFixed(3), S });
  }
  // enforce non-increasing survival (a real digitizer / the reconstruction PAVA-fixes this)
  const Svec = _.pavaDecreasing(raw.map(p => p.S)).y;
  return { pts: raw.map((p, i) => ({ t: p.t, S: +Svec[i].toFixed(4) })), tmax };
}

function armFrom(summary, ipd, tmax) {
  return { km_points: summary, nar_points: [], N: ipd.length,
    total_events: ipd.filter(r => r.status === 1).length, follow_up_max: +tmax.toFixed(2) };
}

function buildTrial(ds, eA, cA) {
  return { nct_id: ds, time_unit: 'days',
    arms: [Object.assign({ arm_id: 'exp', label: 'Experimental', role: 'experimental' }, eA),
           Object.assign({ arm_id: 'ctl', label: 'Comparator', role: 'comparator' }, cA)] };
}

function scoreRecon(trial, truth, tau) {
  const r = RIPD.reconstruct(trial, {});
  if (!r.arms) return null;
  const e = r.arms.find(a => a.role === 'experimental'), c = r.arms.find(a => a.role === 'comparator');
  const hr = GS.coxHR(e.ipd, c.ipd).hr;
  const medE = _.medianFromKM(_.kmFromIPD(e.ipd), { interpolate: true });
  const rmstd = _.rmst(_.kmFromIPD(e.ipd), tau) - _.rmst(_.kmFromIPD(c.ipd), tau);
  return {
    HR: +hr.toFixed(3),
    logHR_err: +Math.abs(Math.log(hr) - Math.log(truth.HR)).toFixed(4),
    median_pcterr: truth.medE ? +(100 * Math.abs(medE - truth.medE) / truth.medE).toFixed(1) : null,
    RMSTdiff_err: +Math.abs(rmstd - truth.RMSTd).toFixed(2),
  };
}

function runOne(cfg, opts) {
  opts = opts || {};
  const sSigma = opts.sSigma == null ? DIG_S_SIGMA : opts.sSigma;
  const tFrac = opts.tFrac == null ? DIG_T_FRAC : opts.tFrac;
  const { expT, ctlT } = GS.loadArms(cfg);
  if (expT.length < 20 || ctlT.length < 20) return null;
  const kmE = _.kmFromIPD(expT), kmC = _.kmFromIPD(ctlT);
  const tau = 0.9 * Math.min(Math.max(...expT.map(r => r.time)), Math.max(...ctlT.map(r => r.time)));
  const truth = { HR: GS.coxHR(expT, ctlT).hr, medE: _.medianFromKM(kmE, { interpolate: true }),
    RMSTd: _.rmst(kmE, tau) - _.rmst(kmC, tau) };
  const rng = _.mulberry32(_.hashStr(cfg.ds));

  // (A) equal density K_REG: exact vs noisy
  const regA_e = anchors(kmE, expT, K_REG), regA_c = anchors(kmC, ctlT, K_REG);
  const digA_e = digitize(kmE, expT, K_REG, rng, sSigma, tFrac), digA_c = digitize(kmC, ctlT, K_REG, rng, sSigma, tFrac);
  const regEqual = scoreRecon(buildTrial('REG=' + cfg.ds, armFrom(regA_e.pts, expT, regA_e.tmax), armFrom(regA_c.pts, ctlT, regA_c.tmax)), truth, tau);
  const digEqual = scoreRecon(buildTrial('DIG=' + cfg.ds, armFrom(digA_e.pts, expT, digA_e.tmax), armFrom(digA_c.pts, ctlT, digA_c.tmax)), truth, tau);

  // (B) realistic: registry few-exact (K_REG) vs digitized many-noisy (K_DIG)
  const digB_e = digitize(kmE, expT, K_DIG, rng, sSigma, tFrac), digB_c = digitize(kmC, ctlT, K_DIG, rng, sSigma, tFrac);
  const digReal = scoreRecon(buildTrial('DIGR=' + cfg.ds, armFrom(digB_e.pts, expT, digB_e.tmax), armFrom(digB_c.pts, ctlT, digB_c.tmax)), truth, tau);

  return { ds: cfg.ds, label: cfg.label, true_HR: +truth.HR.toFixed(3),
    equal_density_K: K_REG, registry_exact: regEqual, digitized_noisy: digEqual,
    realistic: { registry_K: K_REG, digitized_K: K_DIG, registry_exact: regEqual, digitized_noisy: digReal } };
}

function aggregate(rows) {
  const ok = rows.filter(Boolean);
  const mean = (xs) => xs.length ? xs.reduce((a, b) => a + b, 0) / xs.length : null;
  const med = (xs) => { if (!xs.length) return null; const s = xs.slice().sort((a, b) => a - b); const m = Math.floor(s.length / 2); return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; };
  const regE = ok.map(r => r.registry_exact.logHR_err);
  const digE = ok.map(r => r.digitized_noisy.logHR_err);
  const digR = ok.map(r => r.realistic.digitized_noisy.logHR_err);
  const regWinsEqual = ok.filter(r => r.registry_exact.logHR_err <= r.digitized_noisy.logHR_err).length;
  const regWinsReal = ok.filter(r => r.registry_exact.logHR_err <= r.realistic.digitized_noisy.logHR_err).length;
  return {
    n_datasets: ok.length,
    equal_density: {
      registry_exact_mean_logHR_err: +mean(regE).toFixed(4), registry_exact_median_logHR_err: +med(regE).toFixed(4),
      digitized_noisy_mean_logHR_err: +mean(digE).toFixed(4), digitized_noisy_median_logHR_err: +med(digE).toFixed(4),
      registry_le_digitized: `${regWinsEqual}/${ok.length}`,
      interpretation: 'Same K points, one exact one pixel-noised: isolates the pure digitization-error cost.',
    },
    realistic: {
      registry_K8_exact_median_logHR_err: +med(regE).toFixed(4),
      digitized_K25_noisy_median_logHR_err: +med(digR).toFixed(4),
      registry_le_digitized: `${regWinsReal}/${ok.length}`,
      interpretation: 'Registry posts few exact anchors; a digitizer reads many noisy points. Real-world trade-off.',
    },
    noise_model: { survival_sigma_pp: DIG_S_SIGMA * 100, time_sigma_frac: DIG_T_FRAC,
      note: 'Conservative: ~1 percentage-point survival error and 0.5% time error per digitized point, '
        + 'within the range reported for KM-figure digitization. Deterministic (seeded per dataset).' },
  };
}

if (require.main === module) {
  const rows = GS.CONFIGS.map(runOne);
  const out = { summary: aggregate(rows), per_dataset: rows.filter(Boolean) };
  fs.writeFileSync(path.join(GS.dir, 'headtohead_results.json'), JSON.stringify(out, null, 2));
  console.log(JSON.stringify(out.summary, null, 2));
  console.log(`\nwrote realipd/headtohead_results.json (${out.per_dataset.length} datasets)`);
} else {
  module.exports = { runOne, aggregate };
}
