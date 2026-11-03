// P63 validation-lane regression tests (local branch f1000/p63-registryipd).
// 1. Kaplan-Meier median on an exact S = 0.5 plateau: survival::survfit returns the midpoint between the time the
//    curve reaches 0.5 and the next time it drops (and the plateau start when it never drops again). The engine's
//    step median returned the plateau start, so 89 of 1,772 arm medians on the committed registry pseudo-IPD
//    disagreed with survfit.
// 2. The offline web page dist/registry-ipd.html must carry the same engine as src/engine.js (it did not: the
//    engine fixes of 19d0426 were never rebuilt into it, so the offline page still lacked the Cox divergence guard).
'use strict';
const test = require('node:test');
const assert = require('node:assert');
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');
const _ = RIPD._;

const rows = (spec) => spec.flatMap(([t, d, c]) => [...Array(d).fill({ time: t, status: 1 }), ...Array(c).fill({ time: t, status: 0 })]);

test('step median on an exact 0.5 plateau is the midpoint to the next drop (survfit rule)', () => {
  // NCT00248040 comparison 39, experimental arm of the committed pseudo-IPD: 13 events at 1, then 8 events + 5 censored at 12
  const km = _.kmFromIPD(rows([[1, 13, 0], [12, 8, 5]]));
  assert.strictEqual(_.medianFromKM(km), 6.5);
  // two events at 2 and 4 of n = 2: S(2) = 0.5, S(4) = 0  ->  3
  assert.strictEqual(_.medianFromKM(_.kmFromIPD(rows([[2, 1, 0], [4, 1, 0]]))), 3);
  // censoring between the plateau start and the next event does not move the end of the plateau
  assert.strictEqual(_.medianFromKM(_.kmFromIPD(rows([[1, 3, 0], [3, 0, 1], [5, 0, 1], [10, 1, 0]]))), 5.5);
});

test('step median on a 0.5 plateau that never drops again is the plateau start (survfit rule)', () => {
  assert.strictEqual(_.medianFromKM(_.kmFromIPD(rows([[1, 1, 0], [2, 1, 0], [3, 0, 1], [4, 0, 1]]))), 2);
});

test('step median without a plateau is unchanged: first time with S < 0.5', () => {
  assert.strictEqual(_.medianFromKM(_.kmFromIPD(rows([[1, 1, 0], [2, 1, 0], [3, 1, 0]]))), 2);
  assert.strictEqual(_.medianFromKM(_.kmFromIPD(rows([[1, 1, 0], [2, 0, 3]]))), null);
});

test('interpolated median is unchanged by the plateau rule', () => {
  const km = _.kmFromIPD(rows([[1, 13, 0], [12, 8, 5]]));
  assert.strictEqual(_.medianFromKM(km, { interpolate: true }), 1);
});

test('offline web page dist/registry-ipd.html inlines exactly src/engine.js', () => {
  const root = path.join(__dirname, '..');
  const html = fs.readFileSync(path.join(root, 'dist', 'registry-ipd.html'), 'utf8');
  const engine = fs.readFileSync(path.join(root, 'src', 'engine.js'), 'utf8').replace(/<\/script>/g, '<\\/script>');
  assert.ok(html.includes('<script>\n' + engine + '\n</script>'), 'dist/registry-ipd.html is stale: run node build.js');
});
