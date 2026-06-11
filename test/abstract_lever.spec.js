// Regression: the censoring lever (per-arm total_events -> QP) recovers RADIANT-4's posted HR.
// Grounds the count->QP->truth half of the abstract-event-count lever on a real trial (NCT01524783).
const { test } = require('node:test');
const assert = require('node:assert');
const { run } = require('../validate/abstract_lever_realtrial.js');

test('RADIANT-4: the event-count lever pulls the HR from outside the posted CI to inside it', () => {
  const r = run();
  assert.strictEqual(r.posted_hr, 0.48);
  // curve-only (no event count) sits OUTSIDE the posted 95% CI — the identifiability trap
  assert.strictEqual(r.curve_only.inside_ci, false);
  // censoring-informed (107/77 fed to the QP) lands INSIDE the posted CI
  assert.strictEqual(r.censoring_informed.inside_ci, true);
  // and is strictly closer to truth than curve-only
  assert.ok(r.censoring_informed.fold < r.curve_only.fold,
    `informed fold ${r.censoring_informed.fold} should beat curve-only ${r.curve_only.fold}`);
});

test('RADIANT-4: calibrating to the abstract HR reproduces it (the high-coverage lever)', () => {
  const r = run();
  const c = r.abstract_hr_calibrated;
  assert.strictEqual(c.abstract_hr, 0.48);
  // calibration imposes the abstract HR almost exactly while preserving the anchors -> inside the CI
  assert.ok(c.inside_ci, `calibrated HR ${c.hr} should be inside the posted CI`);
  assert.ok(c.fold <= 1.05, `calibrated fold ${c.fold} should be ~1 (HR imposed)`);
  // it solves a plausible experimental event count near the true 107
  assert.ok(c.exp_total_events > 80 && c.exp_total_events < 140,
    `solved exp events ${c.exp_total_events} should be near the true 107`);
});
