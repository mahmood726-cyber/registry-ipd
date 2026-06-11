#!/usr/bin/env node
/*
 * UNCERTAINTY VALIDATION ON REAL DATA — does the reconstruction's own 95% credible interval cover the
 * INDEPENDENT published HR, on real (coarse, NAR-less) registry curves?
 *
 * The gold-standard uncertainty check (goldstandard_uncertainty.js) showed the multiple-imputation HR
 * interval covers the TRUE HR 14/14 — but on clean OPEN IPD. This asks the harder, real-world question:
 * take each cohort trial with a high-confidence PUBLISHED HR (cohort_pubmed_validation.json), run
 * reconstructEnsemble on the actual registry curve, and check whether the published HR falls inside the
 * reconstructed 95% credible interval. Independent truth (PubMed), real registry input, the engine's own
 * UQ — a direct test of whether the uncertainty is calibrated where it matters.
 *
 * Run from repo root: node validate/cohort_uncertainty_validation.js
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');

const COHORT = path.join(__dirname, '..', 'cohort');
const RIPDdir = path.join(__dirname, '..', 'realipd');
if (!fs.existsSync(COHORT)) { console.error('cohort/ missing (gitignored; re-harvest)'); process.exit(2); }
const cpv = JSON.parse(fs.readFileSync(path.join(RIPDdir, 'cohort_pubmed_validation.json'), 'utf8'));
const rowsIn = cpv.rows.filter(r => r.hr_high_confidence && r.published_HR != null);

const med = (xs) => { if (!xs.length) return null; const s = xs.slice().sort((a, b) => a - b); return +s[s.length >> 1].toFixed(3); };
const out = [];
let covered = 0, n = 0;
const widthRatios = [];
for (const r of rowsIn) {
  const fp = path.join(COHORT, `${r.nct}.json`);
  if (!fs.existsSync(fp)) continue;
  let t; try { t = JSON.parse(fs.readFileSync(fp, 'utf8')); } catch { continue; }
  let e; try { e = RIPD.reconstructEnsemble(t, { M: 200 }); } catch { continue; }
  if (!e.ensemble || !e.ensemble.hr) continue;
  const { lo, hi, est } = e.ensemble.hr;
  if (!(lo > 0 && hi > lo)) continue;
  n++;
  const inCI = r.published_HR >= lo && r.published_HR <= hi;
  if (inCI) covered++;
  widthRatios.push(hi / lo);
  out.push({ nct: r.nct, condition: r.condition, published_HR: r.published_HR,
    recon_HR_est: est, recon_95CI: [lo, hi], published_in_recon_CI: inCI,
    registry_HR: r.registry_HR == null ? null : r.registry_HR });
}

let goldRef = 'goldstandard_uncertainty.js (clean open IPD, true HR)';
try {
  const gu = JSON.parse(fs.readFileSync(path.join(RIPDdir, 'goldstandard_uncertainty.json'), 'utf8'));
  goldRef = `clean open IPD: interval covers the TRUE HR ${gu.true_HR_in_95pct_credible_interval}, median CI width fold ${gu.median_CI_width_fold} (goldstandard_uncertainty.js)`;
} catch { /* keep default */ }
const summary = {
  trials: n,
  published_HR_in_reconstructed_95CI: `${covered}/${n}`,
  coverage_pct: n ? +(100 * covered / n).toFixed(0) : null,
  median_CI_width_ratio_hi_over_lo: med(widthRatios),
  gold_standard_reference: goldRef,
  note: 'Reconstructed multiple-imputation 95% credible interval (reconstructEnsemble, M=200) vs the '
    + 'INDEPENDENT published HR, on real coarse registry curves. Coverage here folds in reconstruction '
    + 'uncertainty AND the published-vs-registry effect divergence (~1 in 5 trials, see GALLERY.md), so '
    + 'it is a conservative real-world floor, not a clean coverage probability. A wide-ish median CI '
    + 'width ratio reflects the genuine censoring-level under-identification the interval is meant to span.',
};
fs.writeFileSync(path.join(RIPDdir, 'cohort_uncertainty_validation.json'), JSON.stringify({ summary, rows: out }, null, 2));
console.log(JSON.stringify(summary, null, 2));
console.log(`\nPublished HR inside the reconstructed 95% CI: ${covered}/${n} on real registry curves`);
console.log('wrote realipd/cohort_uncertainty_validation.json');
