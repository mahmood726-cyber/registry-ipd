// Phase 3: granularity-mixed synthesis. The mixed-corpus identified set must bracket the all-IPD truth,
// and the evidence-completeness curve (more curve-only trials -> wider set) must be monotone, starting at a
// point. Needs realipd/cbio_*.csv (gitignored); skips cleanly when absent.
const { test } = require('node:test');
const assert = require('node:assert');
const { run } = require('../validate/phase3_granularity_mixed.js');

test('mixed-granularity pool brackets truth; completeness curve grows monotonically from a point', (t) => {
  const r = run();
  if (r.k < 5) { t.skip(`realipd/cbio_*.csv not present (k=${r.k})`); return; }
  // a corpus mixing IPD + reconstructed-curve + HR-only trials produces a set that contains the all-IPD HR
  assert.ok(r.mixed_corpus.brackets_truth,
    `mixed set [${r.mixed_corpus.HR_set}] should bracket all-IPD HR ${r.true_all_ipd.pooled_HR}`);
  const c = r.completeness_curve;
  // all-IPD (fraction 0) is fully identified: zero-width set (a point)
  assert.strictEqual(c[0].HR_set_width, 0, 'all-IPD corpus should give a zero-width (point) set');
  // identified-set width is monotone non-decreasing in the curve-only fraction
  for (let i = 1; i < c.length; i++) {
    assert.ok(c[i].HR_set_width >= c[i - 1].HR_set_width - 1e-9,
      `set width should not shrink as curve fraction rises (${c[i - 1].HR_set_width} -> ${c[i].HR_set_width})`);
  }
  // an all-curve corpus carries a materially wide set (the honest cost of low granularity)
  assert.ok(c[c.length - 1].HR_set_width > 1,
    `all-curve set width ${c[c.length - 1].HR_set_width} should be materially > 1`);
});
