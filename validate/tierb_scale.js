#!/usr/bin/env node
/*
 * TIER-B SCALE — reconstruct the harvested Tier-B population (median + HR + N, NO curve) and quantify
 * how far it extends the reconstructable registry beyond the curve-based Tier-A cohort.
 *
 * The Tier-A curve cohort is SATURATED (595 harvested >= 514 broad-census curves). The expansion lives
 * in Tier B: trials that post a survival median + a hazard ratio but no KM curve. harvest/harvest_tierb.py
 * builds tierb_cohort/<nct>.json from the AACT snapshot; this reconstructs each via the engine's
 * parametric (exponential) Tier-B path and reports the count + a fidelity-honesty note.
 *
 * Tier B is intentionally LOW-fidelity: the median pins one point of an assumed-exponential curve, the
 * HR is IMPOSED (not recovered), and RMST carries ~7% error with worse tails on non-exponential survival
 * (see VALIDATION.md). It is a coverage/triage tier, not an HR-recovery tier. We report it as such.
 *
 * Run from repo root: node validate/tierb_scale.js  ->  realipd/tierb_scale.json
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');

const DIR = path.join(__dirname, '..', 'tierb_cohort');
if (!fs.existsSync(DIR)) { console.error('tierb_cohort/ missing — run python harvest/harvest_tierb.py'); process.exit(2); }
const files = fs.readdirSync(DIR).filter(f => f.endsWith('.json') && f.startsWith('NCT'));
if (!files.length) { console.error('tierb_cohort/ empty — run python harvest/harvest_tierb.py'); process.exit(2); }

let nB = 0, nC = 0, nErr = 0, nSkip = 0, badge = { silver: 0, bronze: 0, none: 0 };
const hrs = [], medRatios = [];
let idx = 0;
for (const f of files) {
  let t; try { t = JSON.parse(fs.readFileSync(path.join(DIR, f), 'utf8')); } catch { nErr++; continue; }
  // Guard against pathological inputs that make the exponential cutoff solver diverge: sane HR and
  // median spread, finite/positive medians and N. (Tier B couples median <-> HR on the exponential;
  // an extreme HR or near-equal/degenerate medians can prevent convergence.)
  const hv = t.hr && t.hr.value;
  const me = t.arms && t.arms[0] && t.arms[0].median && t.arms[0].median.value;
  const mc = t.arms && t.arms[1] && t.arms[1].median && t.arms[1].median.value;
  const ratio = (me > 0 && mc > 0) ? me / mc : 0;
  if (!(hv > 0.05 && hv < 20) || !(me > 0 && mc > 0) || ratio > 50 || ratio < 0.02) { nSkip++; continue; }
  if (process.env.TB_TRACE && (++idx % 100 === 0)) process.stderr.write(`#${idx} ${t.nct_id}\n`);
  // bootstrap:0 -> skip the 1000-draw uncertainty envelope (O(n^2) Cox each). We only need tier
  // classification + badge for the coverage count, not per-trial CIs; the central draw is built either way.
  let r; try { r = RIPD.reconstruct(t, { bootstrap: 0 }); } catch { nErr++; continue; }
  if (r.tier === 'B') {
    nB++;
    badge[r.audit.badge] = (badge[r.audit.badge] || 0) + 1;
    if (t.hr && t.hr.value != null) hrs.push(t.hr.value);
    const me = t.arms[0].median.value, mc = t.arms[1].median.value;
    if (me > 0 && mc > 0) medRatios.push(me / mc);
  } else if (r.tier === 'C') { nC++; } else { nB++; }  // (a Tier-A here would be unexpected)
}
const sorted = hrs.slice().sort((a, b) => a - b);
const medHR = sorted.length ? sorted[sorted.length >> 1] : null;

const out = {
  generated: 'Tier-B scale: parametric reconstruction of median+HR trials with no posted KM curve',
  harvested_tierb_trials: files.length,
  reconstructed_tierB: nB,
  rejected_tierC: nC,
  skipped_pathological: nSkip,
  errors: nErr,
  badge_distribution: badge,
  median_registry_HR: medHR != null ? +medHR.toFixed(3) : null,
  badge_interpretation: 'Most Tier-B trials self-audit to badge "none" because the exponential PH ' +
    'coupling cannot honor BOTH the posted experimental median AND the posted HR at once — the ' +
    'reconstructed exp median is forced to cMed/HR, so it disagrees with the independently-posted exp ' +
    'median whenever the trial is non-exponential (the common case). silver/bronze flag the minority ' +
    'where they happen to be consistent. This badge split is itself the honest signal: Tier B is ' +
    'coverage, not fidelity.',
  honest_note: 'Tier B is a COVERAGE tier, not an HR-recovery tier: the curve is assumed exponential, ' +
    'the HR is imposed (not recovered), RMST carries ~7% error and worse on non-exponential survival. ' +
    'Use Tier-B pseudo-IPD for triage/scoping only; the curve-based Tier-A cohort remains the ' +
    'validated-fidelity deliverable.',
  reproduce: 'python harvest/harvest_tierb.py && node validate/tierb_scale.js',
};
fs.writeFileSync(path.join(__dirname, '..', 'realipd', 'tierb_scale.json'), JSON.stringify(out, null, 2));
console.log(`Tier-B reconstructed: ${nB}/${files.length} (rejected C: ${nC}, errors: ${nErr})`);
console.log(`  badges: ${JSON.stringify(badge)}`);
console.log(`  median registry HR: ${medHR != null ? medHR.toFixed(3) : 'n/a'}`);
console.log('  wrote realipd/tierb_scale.json');
