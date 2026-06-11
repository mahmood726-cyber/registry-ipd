// Phase 2: honest pooling on the 14 REAL TCGA reconstructions. Needs realipd/cbio_*.csv (gitignored);
// skips cleanly when the data isn't present so a fresh clone's test run does not fail.
const { test } = require('node:test');
const assert = require('node:assert');
const { run } = require('../validate/phase2_real_pooling.js');

test('honest pooling on real reconstructions: HR recovered, naive overstates heterogeneity', (t) => {
  const r = run();
  if (r.k < 5) { t.skip(`realipd/cbio_*.csv not present (k=${r.k})`); return; }
  // the honest (Rubin) pooled HR is closer to the true-IPD pooled HR than the naive pool
  const dTrueHonest = Math.abs(r.honest.pooled_HR - r.true_ipd.pooled_HR);
  const dTrueNaive = Math.abs(r.naive.pooled_HR - r.true_ipd.pooled_HR);
  assert.ok(dTrueHonest <= dTrueNaive + 1e-9,
    `honest HR ${r.honest.pooled_HR} should be >= as close to true ${r.true_ipd.pooled_HR} as naive ${r.naive.pooled_HR}`);
  // ignoring reconstruction variance overstates between-trial heterogeneity
  assert.ok(r.naive.tau2 > r.true_ipd.tau2,
    `naive tau2 ${r.naive.tau2} should exceed true ${r.true_ipd.tau2}`);
  // reconstruction materially adds per-trial variance on these heavily-censored cohorts
  const meanInfl = r.per_trial.reduce((a, p) => a + p.var_inflation, 0) / r.per_trial.length;
  assert.ok(meanInfl > 1.5, `mean variance inflation ${meanInfl} should exceed 1.5x`);
});
