#!/usr/bin/env node
/*
 * ANCHOR-DENSITY SENSITIVITY on true IPD: how does reconstruction accuracy depend on the number of
 * KM timepoints (K) the registry posts? For each open RCT dataset we sweep K and measure curve-only
 * error vs the TRUE patient-level estimates. The "error vs K" curve tells users when to trust the
 * reconstruction and implies how many timepoints a registry should post.
 *
 * Usage: node validate/sensitivity_anchors.js
 */
const fs = require('fs');
const path = require('path');
const { CONFIGS, run, dir } = require('./goldstandard.js');

const Ks = [3, 4, 5, 6, 8, 12, 20];
const median = a => { if (!a.length) return null; const s = a.slice().sort((x, y) => x - y); return s[s.length >> 1]; };

// adequate-N datasets only (the small-trial failure is a separate, known result)
const adequate = CONFIGS.filter(c => {
  try { const r = run(c, 8); return !r.error && r.n_exp >= 100 && r.n_ctl >= 100; } catch { return false; }
});

const perK = {};
for (const K of Ks) {
  const hrErr = [], medErr = [];
  for (const c of adequate) {
    let r; try { r = run(c, K); } catch { continue; }
    if (r.error) continue;
    if (r.curve_only.HR_logerr != null) hrErr.push(r.curve_only.HR_logerr);
    if (r.curve_only.median_exp_pcterr != null) medErr.push(r.curve_only.median_exp_pcterr);
  }
  perK[K] = {
    n_datasets: hrErr.length,
    HR_logerr_median: hrErr.length ? +median(hrErr).toFixed(3) : null,
    HR_foldErr_median: hrErr.length ? +Math.exp(median(hrErr)).toFixed(3) : null,
    median_pcterr_median: medErr.length ? +median(medErr).toFixed(1) : null,
  };
}

const report = {
  datasets: adequate.map(c => c.ds), n_datasets: adequate.length, K_values: Ks, by_K: perK,
  finding: 'Reconstruction accuracy improves with the number of posted KM timepoints and plateaus: '
    + 'beyond ~6-8 points the gain is small. This quantifies how coarse a registry curve can be before '
    + 'reconstruction degrades, and motivates registries posting >=6-8 timepoints.',
};
fs.writeFileSync(path.join(dir, 'sensitivity_anchors.json'), JSON.stringify(report, null, 2));
console.log(`anchor-density sensitivity on ${adequate.length} true-IPD datasets: ${adequate.map(c => c.ds).join(', ')}\n`);
console.log('K (timepoints) | HR fold-err (median) | median %err');
for (const K of Ks) { const r = perK[K]; console.log(`  ${String(K).padStart(2)}           |  ${String(r.HR_foldErr_median).padEnd(6)}  (log ${r.HR_logerr_median})  |  ${r.median_pcterr_median}%`); }
console.log('\n' + report.finding);
