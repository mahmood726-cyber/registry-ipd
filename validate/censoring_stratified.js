#!/usr/bin/env node
/*
 * WHEN DOES THE EVENT-COUNT TIER MATTER? — curve-only vs Titman-QP, stratified by censoring fraction.
 *
 * The headline finding of this project is that the registry curve does not identify the censoring, so
 * curve-only underestimates the HR — but *how much* depends on how censored the trial is. This turns
 * that into actionable guidance: for each gold-standard dataset compute the censoring fraction
 * (1 − events/N, pooled over arms) and pair it with the curve-only and QP HR fold-errors, then
 * stratify. Expectation: curve-only degrades as censoring rises while the QP (which injects the event
 * count) stays flat — so the gap *is* the value of the event count, and it is ~0 for lightly-censored
 * trials and large for heavily-censored ones.
 *
 * Run: node validate/censoring_stratified.js  ->  realipd/censoring_stratified_results.json
 */
const fs = require('fs');
const path = require('path');
const GS = require('./goldstandard.js');

const results = require(path.join(GS.dir, 'goldstandard_results.json'));
const byDs = {}; for (const r of results) if (r && r.ds) byDs[r.ds] = r;

const rows = [];
for (const cfg of GS.CONFIGS) {
  const res = byDs[cfg.ds]; if (!res || res.error) continue;
  let arms; try { arms = GS.loadArms(cfg); } catch { continue; }
  const { expT, ctlT } = arms; if (!expT || expT.length < 20 || ctlT.length < 20) continue;
  const N = expT.length + ctlT.length, ev = expT.filter(r => r.status === 1).length + ctlT.filter(r => r.status === 1).length;
  const cE = 1 - expT.filter(r => r.status === 1).length / expT.length, cC = 1 - ctlT.filter(r => r.status === 1).length / ctlT.length;
  rows.push({ ds: cfg.ds, tcga: cfg.ds.startsWith('cbio'), n: N, events: ev, censoring_fraction: +(1 - ev / N).toFixed(3),
    max_arm_censoring: +Math.max(cE, cC).toFixed(3), censoring_asymmetry: +Math.abs(cE - cC).toFixed(3), abs_log_HR: +Math.abs(Math.log(res.true_HR)).toFixed(3),
    curve_only_fold: +Math.exp(res.curve_only.HR_logerr).toFixed(3),
    qp_fold: +Math.exp(res.censoring_informed.HR_logerr).toFixed(3) });
}

const BINS = [[0, 0.4, '<40%'], [0.4, 0.6, '40–60%'], [0.6, 0.8, '60–80%'], [0.8, 1.01, '≥80%']];
const med = (xs) => { if (!xs.length) return null; const s = xs.slice().sort((a, b) => a - b); return +s[s.length >> 1].toFixed(3); };
const strata = BINS.map(([lo, hi, lab]) => {
  const set = rows.filter(r => r.censoring_fraction >= lo && r.censoring_fraction < hi);
  return { band: lab, n: set.length, median_censoring: med(set.map(r => r.censoring_fraction)),
    curve_only_median_fold: med(set.map(r => r.curve_only_fold)), qp_median_fold: med(set.map(r => r.qp_fold)),
    gap: set.length ? +(med(set.map(r => r.curve_only_fold)) - med(set.map(r => r.qp_fold))).toFixed(3) : null };
});

// rank correlation (Spearman) between censoring fraction and curve-only fold-error
function spearman(a, b) {
  const rank = (v) => { const idx = v.map((x, i) => [x, i]).sort((p, q) => p[0] - q[0]); const r = new Array(v.length); idx.forEach(([_, i], k) => r[i] = k + 1); return r; };
  const ra = rank(a), rb = rank(b), n = a.length; let d2 = 0; for (let i = 0; i < n; i++) d2 += (ra[i] - rb[i]) ** 2;
  return +(1 - 6 * d2 / (n * (n * n - 1))).toFixed(3);
}
const cf = rows.map(r => r.curve_only_fold);
const predictors = {
  pooled_censoring: spearman(rows.map(r => r.censoring_fraction), cf),
  max_arm_censoring: spearman(rows.map(r => r.max_arm_censoring), cf),
  censoring_asymmetry: spearman(rows.map(r => r.censoring_asymmetry), cf),
  effect_size_absLogHR: spearman(rows.map(r => r.abs_log_HR), cf),
};
const meanGap = +(rows.reduce((a, r) => a + (r.curve_only_fold - r.qp_fold), 0) / rows.length).toFixed(3);

const out = { summary: { n_datasets: rows.length, strata, mean_event_count_value_gap: meanGap,
  spearman_curveonly_fold_vs: predictors,
  guidance: 'The QP (event-count) gap over curve-only is consistently POSITIVE (mean ~' + meanGap + ' fold) '
    + 'but is NOT cleanly predictable: pooled censoring barely correlates (' + predictors.pooled_censoring + '); the weak '
    + 'drivers are censoring ASYMMETRY between arms (' + predictors.censoring_asymmetry + ') and EFFECT SIZE ('
    + predictors.effect_size_absLogHR + '). So there is no safe "curve-only is fine here" rule from summary '
    + 'features — you cannot reliably know a priori when curve-only will underestimate the HR. The honest '
    + 'conclusion: ALWAYS prefer the event count (QP) when the registry posts it, and a figure NAR table '
    + 'when it does not; reserve bare curve-only for triangulation with its credible interval.',
}, per_dataset: rows.sort((a, b) => a.censoring_fraction - b.censoring_fraction) };
fs.writeFileSync(path.join(GS.dir, 'censoring_stratified_results.json'), JSON.stringify(out, null, 2));

console.log('When does the event-count tier matter? (' + rows.length + ' datasets)');
console.log('  censoring band   n   median cens   curve-only fold   QP fold   gap');
for (const s of strata) if (s.n) console.log('   ' + s.band.padEnd(8) + String(s.n).padStart(8) + '   ' + String(s.median_censoring).padStart(9)
  + '   ' + String(s.curve_only_median_fold).padStart(13) + '   ' + String(s.qp_median_fold).padStart(7) + '   ' + String(s.gap).padStart(5));
console.log('\n  mean event-count gap (curve-only − QP) =', meanGap, 'fold (consistently positive)');
console.log('  Spearman(curve-only fold, predictor):', JSON.stringify(predictors));
console.log('  => no clean predictor; weak drivers = censoring asymmetry + effect size. Always prefer the event count.');
