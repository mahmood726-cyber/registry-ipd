// R path (index.html inline script, inlined again into dist/registry-ipd.html by build.js): the generated R
// script carries the page's own pseudo-IPD and estimator values at full precision; a Tier C trial (no
// pseudo-IPD) gives NOT CHECKED instead of a script.
'use strict';
const test = require('node:test');
const assert = require('node:assert');
const fs = require('fs');
const path = require('path');
const vm = require('vm');
const root = path.join(__dirname, '..');
function loadPage(rel) {
  const html = fs.readFileSync(path.join(root, rel), 'utf8');
  const ctx = { console, document: { getElementById: () => null, createElement: () => ({}) }, getComputedStyle: () => ({ getPropertyValue: () => '' }) };
  ctx.window = ctx; ctx.self = ctx; vm.createContext(ctx);
  for (const m of html.matchAll(/<script([^>]*)>([\s\S]*?)<\/script>/gi)) {
    const src = (m[1].match(/src="([^"]+)"/) || [])[1];
    try { vm.runInContext(src ? fs.readFileSync(path.join(root, src), 'utf8') : m[2], ctx); } catch (e) { /* UI script needs a DOM */ }
  }
  return ctx;
}
for (const rel of ['index.html', 'dist/registry-ipd.html']) {
  test(`${rel}: R script embeds the page's Cox estimate and pseudo-IPD exactly`, () => {
    const c = loadPage(rel);
    const trial = c.RIPD_EXAMPLES.realOnc.trial, res = c.RIPD.reconstruct(trial, { bootstrap: 300 });
    const b = c.RIPD_RPATH.buildScript(trial, res);
    assert.ok(b.ok);
    const hr = res.audit.checks.C4_hr.recon_hr;
    assert.ok(b.script.includes('recon HR = exp(beta) = ' + String(hr)));
    assert.ok(b.script.includes('time = c(' + res.arms[0].ipd.map((r) => String(r.time)).join(', ') + ')'));
    assert.match(b.script, /MATCH %d\/%d/); assert.match(b.script, /MISMATCH/);
    assert.match(b.script, /quadprog not available/);          // named when missing, never faked
  });
  test(`${rel}: Tier C gives NOT CHECKED, not a script`, () => {
    const c = loadPage(rel);
    const trial = c.RIPD_EXAMPLES.tierC.trial;
    const b = c.RIPD_RPATH.buildScript(trial, c.RIPD.reconstruct(trial, {}));
    assert.strictEqual(b.ok, false); assert.match(b.reason, /Tier C/);
  });
}
