// The synthesis linchpin: pooling reconstructed effects must propagate reconstruction variance, or the
// heterogeneity (tau^2) and prediction interval are distorted. Monte-Carlo with relaxed (stochastic)
// tolerances per the project's testing rules; the effects are large (~2x), so the assertions are robust.
const { test } = require('node:test');
const assert = require('node:assert');
const { run } = require('../validate/honest_pooling_sim.js');

test('ignoring reconstruction variance inflates tau^2; Rubin propagation recovers it', () => {
  const r = run({ reps: 4000, seed: 12345 });
  const trueT = r.true_ipd.mean_tau2, naive = r.naive.mean_tau2, honest = r.honest.mean_tau2;
  // NAIVE mis-reads reconstruction noise as between-trial heterogeneity -> tau^2 markedly inflated
  assert.ok(naive > 1.5 * trueT, `naive tau^2 ${naive} should be >1.5x true ${trueT}`);
  // HONEST (Rubin total variance) recovers the true-IPD tau^2 within ~25%
  assert.ok(Math.abs(honest - trueT) / trueT < 0.25, `honest tau^2 ${honest} should be near true ${trueT}`);
  assert.ok(honest < naive * 0.7, `honest tau^2 ${honest} should be well below naive ${naive}`);
});

test('naive pooling inflates the prediction interval; honest pooling matches true-IPD width', () => {
  const r = run({ reps: 4000, seed: 777 });
  const trueW = r.true_ipd.mean_PI_width, naiveW = r.naive.mean_PI_width, honestW = r.honest.mean_PI_width;
  assert.ok(naiveW > 1.2 * trueW, `naive PI width ${naiveW} should be >1.2x true ${trueW}`);
  assert.ok(Math.abs(honestW - trueW) / trueW < 0.15, `honest PI width ${honestW} should match true ${trueW}`);
});

test('pooled-mean coverage stays near nominal for all three (REML tau^2 self-compensates)', () => {
  const r = run({ reps: 4000, seed: 2024 });
  for (const m of [r.true_ipd, r.naive, r.honest]) {
    assert.ok(m.coverage_mu > 0.92 && m.coverage_mu < 0.98,
      `pooled-mean coverage ${m.coverage_mu} should be near nominal 0.95`);
  }
});
