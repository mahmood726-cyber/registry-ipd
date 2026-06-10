#!/usr/bin/env node
/*
 * UNCERTAINTY CALIBRATION across nominal levels — is the credible interval honest at every width, or
 * only at 95%?
 *
 * `goldstandard_uncertainty.js` reports marginal 95% coverage (true HR inside the 95% credible
 * interval). That is necessary but weak: an interval can hit 95% coverage while being badly
 * miscalibrated elsewhere (e.g. far too wide at 50%). This script reconstructs the multiple-imputation
 * log-HR distribution per dataset (replicating the engine ensemble's censoring-level sampling) and, for
 * EACH nominal level p ∈ {50,80,90,95}, checks whether the TRUE log-HR falls in the central p% credible
 * interval. Empirical coverage ≈ nominal at every level ⇒ well-calibrated; consistently above ⇒
 * conservative (intervals too wide); below ⇒ anti-conservative (overconfident).
 *
 * Adequate-N (≥100/arm) gold-standard datasets. Deterministic. Run: node validate/uncertainty_calibration.js
 * Writes realipd/uncertainty_calibration_results.json.
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');
const _ = RIPD._;
const GS = require('./goldstandard.js');

const M = 400;                              // imputations per dataset
const LEVELS = [50, 80, 90, 95];
const coxHR = (a, b) => Math.exp(_.coxLogHR(a.map(r => ({ time: r.time, status: r.status, x: 1 }))
  .concat(b.map(r => ({ time: r.time, status: r.status, x: 0 }))).slice()).beta);

function coarse(ipd, K) {
  const km = _.kmFromIPD(ipd), tmax = 0.95 * Math.max(...ipd.map(r => r.time)), pts = [{ t: 0, S: 1 }];
  for (let i = 1; i <= K; i++) { const t = tmax * i / K; pts.push({ t: +t.toFixed(2), S: +_.evalKM(km, t).toFixed(4) }); }
  return { km_points: pts, nar_points: [], N: ipd.length, total_events: ipd.filter(r => r.status === 1).length, follow_up_max: +tmax.toFixed(2) };
}

// replicate the engine ensemble's per-arm censoring-level sampling and collect the log-HR distribution
function lhrSamples(trial) {
  const seed0 = (_.hashStr(trial.nct_id || '') ^ 0x1234abcd) >>> 0;
  const e0 = {};
  { const cr = RIPD.reconstruct(trial, { ignoreTotalEvents: true }); if (cr.arms) cr.arms.forEach(a => { e0[a.arm_id] = a.ipd.filter(x => x.status === 1).length; }); }
  const ei = trial.arms.findIndex(a => a.role === 'experimental'), ci = trial.arms.findIndex(a => a.role === 'comparator');
  const out = [];
  for (let m = 0; m < M; m++) {
    const rng = _.mulberry32((seed0 + Math.imul(m + 1, 2654435761)) >>> 0);
    const t2 = JSON.parse(JSON.stringify(trial));
    for (const a of t2.arms) {
      a.km_points = a.km_points.map(p => ({ t: p.t, S: Math.min(1, Math.max(0, p.S + (rng() - 0.5) * 0.01)) }));
      const E0 = e0[a.arm_id];
      if (E0 != null && E0 > 0) { const reg = a.total_events, lo = Math.round(0.55 * E0); let e = Math.round(lo + rng() * (E0 - lo));
        if (reg != null) e = rng() < 0.5 ? Math.max(0, Math.round(reg * (1 + (rng() - 0.5) * 0.30))) : e;
        a.total_events = Math.min(a.N != null ? a.N : e, Math.max(1, e)); }
    }
    const method = rng() < 0.5 ? 'guyot' : 'anchor-exact';
    let r; try { r = RIPD.reconstruct(t2, { method }); } catch { continue; }
    if (!r.arms) continue;
    out.push(Math.log(coxHR(r.arms[ei].ipd, r.arms[ci].ipd)));
  }
  return out.sort((a, b) => a - b);
}

const perLevel = {}; LEVELS.forEach(p => perLevel[p] = { covered: 0, n: 0 });
const rows = [];
for (const cfg of GS.CONFIGS) {
  let arms; try { arms = GS.loadArms(cfg); } catch { continue; }
  const { expT, ctlT } = arms; if (expT.length < 100 || ctlT.length < 100) continue;
  const trueLHR = Math.log(coxHR(expT, ctlT));
  const trial = { nct_id: cfg.ds, time_unit: 'd', arms: [Object.assign({ arm_id: 'exp', role: 'experimental' }, coarse(expT, 8)), Object.assign({ arm_id: 'ctl', role: 'comparator' }, coarse(ctlT, 8))] };
  const s = lhrSamples(trial); if (s.length < 50) continue;
  const rec = { ds: cfg.ds, true_HR: +Math.exp(trueLHR).toFixed(3), levels: {} };
  for (const p of LEVELS) {
    const lo = _.quantileSorted(s, (1 - p / 100) / 2), hi = _.quantileSorted(s, 1 - (1 - p / 100) / 2);
    const inside = trueLHR >= lo && trueLHR <= hi;
    perLevel[p].covered += inside ? 1 : 0; perLevel[p].n++;
    rec.levels[p] = inside;
  }
  rows.push(rec);
}

const calibration = LEVELS.map(p => ({ nominal: p, empirical: +(100 * perLevel[p].covered / perLevel[p].n).toFixed(1), covered: perLevel[p].covered, n: perLevel[p].n }));
const out = {
  summary: {
    n_datasets: rows.length, M,
    calibration,
    interpretation: 'Empirical ≈ nominal at every level ⇒ well-calibrated. Consistently above nominal ⇒ '
      + 'conservative (intervals wider than needed). Below ⇒ overconfident. This is the multi-level check '
      + 'behind the marginal 95% coverage reported in goldstandard_uncertainty.js.',
  },
  per_dataset: rows,
};
fs.writeFileSync(path.join(GS.dir, 'uncertainty_calibration_results.json'), JSON.stringify(out, null, 2));
console.log('Uncertainty calibration across nominal levels (' + rows.length + ' datasets, M=' + M + '):');
console.log('  nominal  empirical  (covered/n)');
for (const c of calibration) console.log('   ' + (c.nominal + '%').padStart(5) + '     ' + (c.empirical + '%').padStart(6) + '    ' + c.covered + '/' + c.n);
const mae = +(calibration.reduce((a, c) => a + Math.abs(c.empirical - c.nominal), 0) / calibration.length).toFixed(1);
console.log('  mean |empirical − nominal| =', mae, 'pp', mae <= 8 ? '(well-calibrated)' : '(miscalibrated)');
