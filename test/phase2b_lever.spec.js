// Phase 2b: a posted event count (pinned in the ensemble) shrinks the reconstruction variance r^2.
// Needs realipd/cbio_*.csv (gitignored); skips cleanly when absent.
const { test } = require('node:test');
const assert = require('node:assert');
const { run } = require('../validate/phase2b_lever_shrinks_r2.js');

test('the censoring lever shrinks reconstruction variance r^2 (the mechanism)', (t) => {
  const r = run();
  if (r.k < 5) { t.skip(`realipd/cbio_*.csv not present (k=${r.k})`); return; }
  // event-pinned reconstruction SD is materially smaller than curve-only
  assert.ok(r.mean_r_event < r.mean_r_curve,
    `event-pinned r ${r.mean_r_event} should be < curve-only r ${r.mean_r_curve}`);
  assert.ok(r.r_shrink_factor > 2,
    `lever should shrink reconstruction SD >2x (got ${r.r_shrink_factor}x)`);
  // per-trial: the shrink holds for the large majority of cohorts, not just on average
  const shrunk = r.per_trial.filter(p => p.r_event < p.r_curve).length;
  assert.ok(shrunk >= 0.8 * r.per_trial.length,
    `event count should shrink r for >=80% of trials (${shrunk}/${r.per_trial.length})`);
});
