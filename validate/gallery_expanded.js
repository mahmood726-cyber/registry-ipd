#!/usr/bin/env node
/*
 * EXPANDED PRODUCTION VALIDATION — score the reconstructed HR against a held-out registry HR across the
 * FULL validation-grade population, not just the 30 trials whose HR the curve-scoped harvester captured.
 *
 * The census finds 112 trials posting both a curve and an HR; the original gallery validated 30 because
 * the harvester dropped HRs posted in a sibling outcome (now fixed: harvester.select_trial_hr). This
 * consumes realipd/validation_hr_backfill.json (the sibling-aware HRs recovered in one snapshot pass by
 * harvest/backfill_validation_hr.py) and re-scores every validation-grade trial that reconstructs as a
 * 2-arm Tier-A curve, reporting median fold-error by HR source:
 *   - curve  : HR from the curve's own outcome (the original 30 + any newly curve-sourced)
 *   - sibling: HR recovered from a survival sibling outcome (flagged: endpoint may be OS-vs-PFS)
 *
 * Honest: the registry HR is a coarse held-out truth and its arm-reference is assumed to match the
 * reconstruction's exp/ctl convention (same caveat as the original gallery); sibling-sourced rows carry
 * the extra endpoint caveat and are reported separately so neither inflates the headline.
 *
 * Run from repo root: node validate/gallery_expanded.js  ->  realipd/gallery_expanded.json + console.
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');

const COHORT = path.join(__dirname, '..', 'cohort');
const RIPDdir = path.join(__dirname, '..', 'realipd');
const load = (f) => JSON.parse(fs.readFileSync(path.join(RIPDdir, f), 'utf8'));
if (!fs.existsSync(COHORT)) { console.error('cohort/ missing (gitignored; re-harvest)'); process.exit(2); }

const cen = load('census_full_aact.json');
const vg = cen.validation_grade_ncts_broad || [];
let bf; try { bf = load('validation_hr_backfill.json'); } catch { bf = { hr: {} }; }
const backfillHR = bf.hr || {};

const coxHR = (a, b) => Math.exp(RIPD._.coxLogHR(a.map(r => ({ time: r.time, status: r.status, x: 1 }))
  .concat(b.map(r => ({ time: r.time, status: r.status, x: 0 }))).slice()).beta);
const fold = (recon, reg) => +Math.exp(Math.abs(Math.log(recon) - Math.log(reg))).toFixed(3);
const med = (xs) => { if (!xs.length) return null; const s = xs.slice().sort((a, b) => a - b); return +s[s.length >> 1].toFixed(3); };

const rows = [];
let n2arm = 0, nNoHR = 0, nNot2arm = 0;
for (const nct of vg) {
  const fp = path.join(COHORT, `${nct}.json`);
  if (!fs.existsSync(fp)) continue;
  let t; try { t = JSON.parse(fs.readFileSync(fp, 'utf8')); } catch { continue; }
  let r; try { r = RIPD.reconstruct(t); } catch { continue; }
  if (r.tier !== 'A' || !r.arms || r.arms.length !== 2) { nNot2arm++; continue; }
  n2arm++;
  // registry HR (+ its posted CI): prefer the cohort-native HR (curve outcome, direction-resolved);
  // else the backfill. The CI lets us ask the stronger question — does the reconstructed HR fall
  // within the registry's OWN reported uncertainty — instead of treating the point HR as exact truth.
  let regHR = null, ciLow = null, ciHigh = null, source = 'curve';
  if (t.hr && t.hr.value != null) { regHR = t.hr.value; ciLow = t.hr.ci_low; ciHigh = t.hr.ci_high; }
  else if (backfillHR[nct] && backfillHR[nct].value != null) {
    const b = backfillHR[nct]; regHR = b.value; ciLow = b.ci_low; ciHigh = b.ci_high;
    source = b.from_sibling_outcome ? 'sibling' : 'curve';
  }
  if (regHR == null) { nNoHR++; continue; }
  const exp = r.arms.find(a => a.role === 'experimental') || r.arms[1];
  const ctl = r.arms.find(a => a.role === 'comparator') || r.arms[0];
  const recon = +coxHR(exp.ipd, ctl.ipd).toFixed(3);
  const ev = exp.ipd.filter(x => x.status === 1).length + ctl.ipd.filter(x => x.status === 1).length;
  const hasCI = Number.isFinite(ciLow) && Number.isFinite(ciHigh) && ciLow > 0 && ciHigh > ciLow;
  const inCI = hasCI ? (recon >= ciLow && recon <= ciHigh) : null;
  rows.push({ nct, condition: (t.condition || '').split(';')[0].trim().slice(0, 40),
    badge: r.audit && r.audit.badge, source, registry_HR: regHR,
    registry_CI: hasCI ? [ciLow, ciHigh] : null, recon_HR: recon,
    fold: fold(recon, regHR), in_registry_CI: inCI, events: ev });
}

const curve = rows.filter(r => r.source === 'curve');
const sibling = rows.filter(r => r.source === 'sibling');
const within20 = (set) => `${set.filter(r => r.fold < 1.2).length}/${set.length}`;
// CI coverage: of rows with a posted registry CI, how many contain the reconstructed HR.
const ciCov = (set) => { const w = set.filter(r => r.in_registry_CI !== null); return `${w.filter(r => r.in_registry_CI).length}/${w.length}`; };
const summary = {
  validation_grade_population: vg.length,
  reconstructed_2arm_tierA: n2arm,
  scored_against_registry_HR: rows.length,
  original_gallery_n: 30,
  by_source: {
    curve: { n: curve.length, median_fold: med(curve.map(r => r.fold)), within_1_2: within20(curve), in_registry_CI: ciCov(curve) },
    sibling: { n: sibling.length, median_fold: med(sibling.map(r => r.fold)), within_1_2: within20(sibling), in_registry_CI: ciCov(sibling) },
  },
  all_scored: { n: rows.length, median_fold: med(rows.map(r => r.fold)), within_1_2: within20(rows), in_registry_CI: ciCov(rows) },
  no_hr_recovered: nNoHR, not_2arm_tierA: nNot2arm,
  note: 'fold = reconstructed Cox HR vs registry HR (coarse held-out truth; arm-reference assumed to '
    + 'match the exp/ctl convention, same as the original gallery). "sibling" rows took the HR from a '
    + 'SAME-ENDPOINT survival sibling outcome (harvester.select_trial_hr is endpoint-aware: an OS curve '
    + 'is never scored against a PFS HR -- explicit-mismatch siblings are dropped). Reported separately '
    + 'from curve-sourced only because the HR sits in a different outcome record, not a different endpoint.',
};
fs.writeFileSync(path.join(RIPDdir, 'gallery_expanded.json'), JSON.stringify({ summary, rows }, null, 2));

console.log(JSON.stringify(summary, null, 2));
console.log(`\nExpanded validation-grade scoring: ${rows.length} trials (was 30) — curve ${curve.length} (median fold ${summary.by_source.curve.median_fold}), sibling ${sibling.length} (median fold ${summary.by_source.sibling.median_fold})`);
console.log(`reconstructed HR within the registry's posted 95% CI: ${summary.all_scored.in_registry_CI}`);
console.log('wrote realipd/gallery_expanded.json');
