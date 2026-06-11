#!/usr/bin/env node
/*
 * REPRODUCIBILITY HARNESS — pin every paper-cited validation headline and re-derive the row-based ones,
 * so a future engine/extractor change cannot silently drift the numbers the manuscript reports.
 *
 * Two guards per result JSON:
 *   (1) GOLDEN: the committed summary still contains the exact value the paper/README/GALLERY cite
 *       (fractions exact, floats within a small tolerance).
 *   (2) CONSISTENCY: for files that ship per-row data, re-derive the headline FROM the rows and assert it
 *       matches the committed summary -- catching a summary-computation bug even if the golden value was
 *       updated to match it. (No external data needed: rows are committed.)
 *
 * Exits non-zero on any mismatch. Run from repo root:  node validate/verify_headlines.js
 *   (also wired as part of `npm run verify`).
 */
const fs = require('fs');
const path = require('path');
const RIPD = path.join(__dirname, '..', 'realipd');
const load = (f) => JSON.parse(fs.readFileSync(path.join(RIPD, f), 'utf8'));
const get = (o, p) => p.split('.').reduce((a, k) => (a == null ? a : a[k]), o);

let pass = 0, fail = 0; const fails = [];
function check(label, actual, expected, tol) {
  let ok;
  if (typeof expected === 'number' && typeof actual === 'number') ok = Math.abs(actual - expected) <= (tol == null ? 1e-9 : tol);
  else ok = String(actual) === String(expected);
  if (ok) { pass++; } else { fail++; fails.push(`${label}: expected ${expected}, got ${actual}`); }
}
const medJS = (xs) => { const s = xs.slice().sort((a, b) => a - b); return s.length ? +s[s.length >> 1].toFixed(3) : null; };
const frac = (rows, pred, base) => { const b = rows.filter(base); return `${b.filter(pred).length}/${b.length}`; };

// ---------------------------------------------------------------- (1) GOLDEN headline values
const GOLDEN = {
  'census_full_aact.json': { 'validation_grade_curve_and_hr_broad': 112, 'validation_grade_curve_and_hr_strict': 77 },
  'gallery_expanded.json': { 'summary.scored_against_registry_HR': 57, 'summary.all_scored.in_registry_CI': '47/57', 'summary.by_source.curve.in_registry_CI': '24/25' },
  'cohort_pubmed_validation.json': { 'summary.published_HR_high_confidence': 20, 'summary.of_which_no_registry_HR': 12, 'summary.HR.recon_within_published_95CI': '17/20' },
  'pubmed_validation.json': { 'summary.high_confidence_extractions': 16, 'summary.high_confidence.recon_within_published_95CI': '13/16' },
  'pubmed_median_validation.json': { 'summary.trials_validated': 5, 'summary.median_arm_fold': 1.071 },
  'cohort_uncertainty_validation.json': { 'summary.published_HR_in_reconstructed_95CI': '17/20' },
  'registry_median_validation.json': { 'summary.trials_curve_consistent': 11, 'summary.median_arm_fold': 1.133, 'summary.registry_curve_median_inconsistent': 2 },
  'tierb_scale.json': { 'reconstructed_tierB': 1144 },
};
for (const file in GOLDEN) {
  let d; try { d = load(file); } catch { fail++; fails.push(`${file}: MISSING (run its generator)`); continue; }
  for (const p in GOLDEN[file]) check(`${file}:${p}`, get(d, p), GOLDEN[file][p], typeof GOLDEN[file][p] === 'number' ? 0.005 : null);
}

// ---------------------------------------------------------------- (2) CONSISTENCY: re-derive from rows
function consistency(file, fn) {
  let d; try { d = load(file); } catch { fail++; fails.push(`${file}: MISSING`); return; }
  try { fn(d); } catch (e) { fail++; fails.push(`${file}: consistency threw ${e.message}`); }
}
consistency('gallery_expanded.json', (d) => {
  const curve = d.rows.filter(r => r.source === 'curve');
  check('gallery_expanded curve.median_fold (rederived)', medJS(curve.map(r => r.fold)), d.summary.by_source.curve.median_fold, 1e-6);
  check('gallery_expanded all in_CI (rederived)', frac(d.rows, r => r.in_registry_CI, r => r.in_registry_CI !== null), d.summary.all_scored.in_registry_CI);
});
consistency('cohort_pubmed_validation.json', (d) => {
  const hc = d.rows.filter(r => r.hr_high_confidence);
  check('cohort_pubmed hc count (rederived)', hc.length, d.summary.published_HR_high_confidence);
  const ci = hc.filter(r => 'recon_in_published_HR_CI' in r);
  check('cohort_pubmed recon_in_CI (rederived)', `${ci.filter(r => r.recon_in_published_HR_CI).length}/${ci.length}`, d.summary.HR.recon_within_published_95CI);
});
consistency('cohort_uncertainty_validation.json', (d) => {
  check('uncertainty coverage (rederived)', `${d.rows.filter(r => r.published_in_recon_CI).length}/${d.rows.length}`, d.summary.published_HR_in_reconstructed_95CI);
});
consistency('registry_median_validation.json', (d) => {
  const folds = d.rows.flatMap(r => r.arm_folds);
  check('registry_median fold (rederived)', medJS(folds), d.summary.median_arm_fold, 1e-6);
});
consistency('pubmed_median_validation.json', (d) => {
  const folds = d.rows.flatMap(r => r.arm_folds);
  check('pubmed_median fold (rederived)', medJS(folds), d.summary.median_arm_fold, 1e-6);
});

// ---------------------------------------------------------------- report
console.log(`headline verification: ${pass} passed, ${fail} failed`);
if (fail) { console.error('\nFAILURES:'); for (const f of fails) console.error('  ✗ ' + f); process.exit(1); }
console.log('all pinned headlines + row-derived summaries consistent.');
