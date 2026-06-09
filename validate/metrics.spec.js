/* Validation-metric tests. Run: node --test validate/metrics.spec.js */
const { test } = require('node:test');
const assert = require('node:assert');
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');
const M = require('./metrics.js');

const fx = (n) => JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'test', 'fixtures', n), 'utf8'));

test('AACT-only fidelity: dense-anchor Tier A drives anchor sup-error ~0', () => {
  const t = fx('fixture_exp_known.json');
  const r = RIPD.reconstruct(t);
  const f = M.fidelity(t, r);
  assert.ok(f.anchor_sup_error <= 1e-3, `anchor sup-error ${f.anchor_sup_error}`);
  assert.ok(f.logHR_err != null && f.logHR_err < 0.3, `logHR err ${f.logHR_err}`);
  assert.ok(f.wasserstein_to_anchors < 0.02, `W1 ${f.wasserstein_to_anchors}`);
});

test('head-to-head: AACT-only beats simulated digitization on anchor fidelity', () => {
  const t = fx('fixture_exp_known.json');
  // simulate digitization: perturb the KM anchors by pixel-like noise (deterministic)
  const rng = RIPD._.mulberry32(99);
  const dig = JSON.parse(JSON.stringify(t));
  for (const a of dig.arms) {
    a.km_points = a.km_points.map(p => ({ t: p.t, S: Math.min(1, Math.max(0, p.S + (rng() - 0.5) * 0.06)) }));
  }
  const h = M.headToHead(t, { digitizedTrial: dig });
  assert.ok(h.aact_only.anchor_sup_error <= 1e-3, `aact sup ${h.aact_only.anchor_sup_error}`);
  assert.ok(h.digitization.anchor_sup_error > h.aact_only.anchor_sup_error,
    `digitization sup-error ${h.digitization.anchor_sup_error} should exceed AACT-only ${h.aact_only.anchor_sup_error}`);
  assert.ok(h.digitization.wasserstein_to_anchors > h.aact_only.wasserstein_to_anchors,
    'digitization W1 should exceed AACT-only W1');
});

test('aggregate computes cohort RMSE and tier counts', () => {
  const trials = ['fixture_exp_known.json'].map(fx);
  const fids = trials.map(t => M.fidelity(t, RIPD.reconstruct(t)));
  const agg = M.aggregate(fids);
  assert.strictEqual(agg.n, 1);
  assert.ok(agg.logHR_RMSE != null);
  assert.ok(agg.anchor_sup_error_max <= 1e-3);
  assert.deepStrictEqual(agg.tier_counts, { A: 1 });
});
