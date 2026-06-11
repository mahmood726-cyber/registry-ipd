// Phase 2c: reconstruction bias is partially-identified heterogeneity. The LOO-calibrated identified set
// must bracket the true-IPD pooled HR and tau^2; de-biasing the systematic offset recovers the point HR.
// Needs realipd/cbio_*.csv (gitignored); skips cleanly when absent.
const { test } = require('node:test');
const assert = require('node:assert');
const { run } = require('../validate/phase2c_bias_offset.js');

test('the calibrated identified set brackets the true HR and tau^2; de-bias recovers the point HR', (t) => {
  const r = run();
  if (r.k < 5) { t.skip(`realipd/cbio_*.csv not present (k=${r.k})`); return; }
  // the identified set contains the held-out true-IPD values (partial identification done right)
  assert.ok(r.contains_truth.HR, `HR set [${r.identified_set.HR}] should contain true ${r.true_ipd.pooled_HR}`);
  assert.ok(r.contains_truth.tau2, `tau^2 set [${r.identified_set.tau2}] should contain true ${r.true_ipd.tau2}`);
  // de-biasing the LOO systematic offset moves the pooled HR closer to truth than the naive point
  const dDeb = Math.abs(r.debiased_point.pooled_HR - r.true_ipd.pooled_HR);
  const dNaive = Math.abs(r.event_pinned_point.pooled_HR - r.true_ipd.pooled_HR);
  assert.ok(dDeb <= dNaive + 1e-9,
    `de-biased HR ${r.debiased_point.pooled_HR} should be >= as close to true ${r.true_ipd.pooled_HR} as naive ${r.event_pinned_point.pooled_HR}`);
  // the set is a genuine interval (a point would be false precision)
  assert.ok(r.identified_set.tau2[1] > r.identified_set.tau2[0],
    'tau^2 identified set should have positive width');
});
