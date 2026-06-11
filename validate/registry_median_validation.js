#!/usr/bin/env node
/*
 * REGISTRY-SIDE MEDIAN CROSS-CHECK — reconstructed per-arm median vs the AACT-POSTED median, at scale.
 *
 * Complements the small-n independent PubMed median check (5 trials) with a large-n registry-side one:
 * the registry posts a median survival for many trials (in a sibling survival outcome -- recovered
 * endpoint-matched by harvest/registry_medians.py). The reconstruction does NOT use the posted median
 * (Tier A reconstructs from the KM timepoints), so comparing the reconstructed median to the
 * separately-posted median is a genuine recovery check. Matched by sorted magnitude (arm labelling
 * irrelevant). NOT independent of the registry -- read alongside the PubMed median check, not instead.
 *
 * Run from repo root: node validate/registry_median_validation.js
 */
const fs = require('fs');
const path = require('path');

const RIPDdir = path.join(__dirname, '..', 'realipd');
const recon = JSON.parse(fs.readFileSync(path.join(RIPDdir, 'cohort_recon.json'), 'utf8'));
const posted = JSON.parse(fs.readFileSync(path.join(RIPDdir, 'registry_medians.json'), 'utf8'));
const reconBy = {}; for (const r of recon) reconBy[r.nct] = r;

const COHORT = path.join(__dirname, '..', 'cohort');
const fold = (a, b) => +Math.exp(Math.abs(Math.log(a) - Math.log(b))).toFixed(3);
const med = (xs) => { if (!xs.length) return null; const s = xs.slice().sort((a, b) => a - b); return +s[s.length >> 1].toFixed(3); };
// last observed KM timepoint per trial — the reconstructed median (= the curve's S=0.5 crossing, by
// construction) is OBSERVED within the data when it is below this. A registry-posted median far beyond
// the observed follow-up, while the curve has already dropped under 0.5, contradicts the curve itself.
const lastTimepoint = (nct) => {
  try {
    const t = JSON.parse(fs.readFileSync(path.join(COHORT, `${nct}.json`), 'utf8'));
    return Math.max(...(t.arms || []).flatMap(a => (a.km_points || []).map(p => p.t)));
  } catch { return null; }
};

const rows = [], flagged = [];
const armFolds = [];
for (const nct in posted) {
  const p = posted[nct], r = reconBy[nct];
  if (!r || !r.recon_medians || r.recon_medians.length !== 2) continue;
  if (!(p.medians && p.medians.length === 2)) continue;
  const ps = p.medians.slice().sort((a, b) => a - b);
  const rs = r.recon_medians.slice().sort((a, b) => a - b);
  if (!(ps[0] > 0 && rs[0] > 0)) continue;
  const f = [fold(rs[0], ps[0]), fold(rs[1], ps[1])];
  const tmax = lastTimepoint(nct);
  // curve-consistency gate: the posted median must be plausible given the observed curve. If the posted
  // median exceeds ~1.5x the last observed timepoint while the reconstructed (curve) median sits well
  // inside it, the registry's posted CURVE and posted MEDIAN disagree -> not a reconstruction error.
  const inconsistent = tmax != null && (ps[1] > 1.5 * tmax) && (rs[1] < tmax);
  const rec = { nct, endpoint: p.endpoint, posted_medians: ps, reconstructed_medians: rs, arm_folds: f, last_curve_timepoint: tmax };
  if (inconsistent) { flagged.push(rec); continue; }
  armFolds.push(...f);
  rows.push(rec);
}

const summary = {
  trials_curve_consistent: rows.length,
  arm_medians: armFolds.length,
  median_arm_fold: med(armFolds),
  within_10pct: `${armFolds.filter(f => f < 1.1).length}/${armFolds.length}`,
  within_20pct: `${armFolds.filter(f => f < 1.2).length}/${armFolds.length}`,
  registry_curve_median_inconsistent: flagged.length,
  inconsistent_ncts: flagged.map(r => r.nct),
  note: 'Reconstructed per-arm median vs the AACT-POSTED same-endpoint median (harvest/registry_medians.py), '
    + 'matched by sorted magnitude. Registry-provenance (NOT independent), but the reconstruction does not '
    + 'use the posted median (Tier A reconstructs from the KM timepoints), so this is a genuine recovery '
    + 'check at larger n than the independent PubMed median validation. A curve-consistency gate excludes '
    + 'trials whose registry-posted median grossly exceeds the observed curve follow-up while the curve is '
    + 'already below 0.5 -- there the registry’s OWN curve and median disagree (a registry data-quality '
    + 'finding the cross-check surfaces; the reconstruction faithfully follows the posted curve). Read '
    + 'with the PubMed median check: this for breadth, PubMed for independence.',
};
fs.writeFileSync(path.join(RIPDdir, 'registry_median_validation.json'), JSON.stringify({ summary, rows, flagged }, null, 2));
console.log(JSON.stringify(summary, null, 2));
console.log(`\nRegistry-posted median recovery: ${armFolds.length} arm-medians across ${rows.length} curve-consistent trials, median fold ${summary.median_arm_fold} (within 10%: ${summary.within_10pct}); ${flagged.length} flagged registry curve/median inconsistencies`);
console.log('wrote realipd/registry_median_validation.json');
