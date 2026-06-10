#!/usr/bin/env node
/*
 * SCALE RUN — reconstruct EVERY pairwise comparison across the full reconstructable AACT cohort.
 *
 * The production gallery (gallery.js) characterised the 2-arm trials. This scales to the whole cohort,
 * including the multi-arm trials (3–15 arms), where each pair of arms is a separate reconstructable
 * comparison. For every trial with ≥2 curve-bearing arms we form the arm pairs (capped), build a 2-arm
 * sub-trial, and reconstruct it with the shipped engine (Titman-QP default). Output: a registry-wide
 * INDEX (metadata only — NCT, condition, badge, method, HR; no patient rows) + scale summary + SCALE.md.
 *
 * This turns the validated method into a registry-scale resource: how many comparisons can be
 * reconstructed from ClinicalTrials.gov today, with what audit quality. Run from repo root:
 *   node validate/scale_run.js   ->   realipd/scale_index.json  +  SCALE.md
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');

const COHORT = path.join(__dirname, '..', 'cohort');
const files = fs.existsSync(COHORT) ? fs.readdirSync(COHORT).filter(f => f.endsWith('.json') && f.startsWith('NCT')) : [];
if (!files.length) { console.error('cohort/ empty (gitignored; re-harvest with harvest/harvest_cohort.py)'); process.exit(2); }

const coxHR = (a, b) => Math.exp(RIPD._.coxLogHR(a.map(r => ({ time: r.time, status: r.status, x: 1 }))
  .concat(b.map(r => ({ time: r.time, status: r.status, x: 0 }))).slice()).beta);
const MAX_PAIRS = 6;   // cap pairwise comparisons per trial (avoid combinatorial blow-up on 15-arm trials)

let nTrials = 0, nReconstructable = 0, nComparisons = 0, nMultiArm = 0, nSingleArm = 0;
const badge = { gold: 0, silver: 0, bronze: 0, none: 0 };
const method = { qp: 0, 'anchor-exact': 0, guyot: 0 };
const byCondition = {};
const index = [];

for (const f of files) {
  let t; try { t = JSON.parse(fs.readFileSync(path.join(COHORT, f), 'utf8')); } catch { continue; }
  nTrials++;
  const arms = (t.arms || []).filter(a => (a.km_points || []).length >= 3 && a.N != null);
  if (arms.length < 2) { if (arms.length === 1) nSingleArm++; continue; }
  if (arms.length > 2) nMultiArm++;
  const cond = (t.condition || 'unspecified').split(';')[0].trim().slice(0, 40);
  const comparisons = [];
  let pairCount = 0;
  for (let i = 0; i < arms.length && pairCount < MAX_PAIRS; i++) {
    for (let j = i + 1; j < arms.length && pairCount < MAX_PAIRS; j++) {
      const trial2 = { nct_id: t.nct_id, time_unit: t.time_unit,
        arms: [Object.assign({}, arms[j], { role: 'experimental', arm_id: 'e' }), Object.assign({}, arms[i], { role: 'comparator', arm_id: 'c' })],
        hr: arms.length === 2 ? t.hr : null };
      let r; try { r = RIPD.reconstruct(trial2); } catch { continue; }
      if (r.tier !== 'A' || !r.arms) continue;
      pairCount++; nComparisons++;
      badge[r.audit.badge] = (badge[r.audit.badge] || 0) + 1;
      method[r.method] = (method[r.method] || 0) + 1;
      const hr = coxHR(r.arms[0].ipd, r.arms[1].ipd);
      comparisons.push({ exp: (arms[j].label || 'arm' + j).slice(0, 28), ctl: (arms[i].label || 'arm' + i).slice(0, 28),
        recon_HR: +hr.toFixed(3), badge: r.audit.badge, method: r.method, exportable: r.exportable });
    }
  }
  if (!comparisons.length) continue;
  nReconstructable++;
  byCondition[cond] = (byCondition[cond] || 0) + 1;
  index.push({ nct: t.nct_id, url: t.source_url, condition: cond, n_arms: (t.arms || []).length,
    registry_HR: t.hr && t.hr.value != null ? t.hr.value : null, comparisons });
}

const topConditions = Object.entries(byCondition).sort((a, b) => b[1] - a[1]).slice(0, 15);
const summary = {
  cohort_trials_scanned: nTrials, single_arm: nSingleArm, multi_arm: nMultiArm,
  reconstructable_trials: nReconstructable, total_pairwise_comparisons: nComparisons,
  badge_distribution: badge, method_distribution: method,
  exportable_comparisons: index.reduce((a, t) => a + t.comparisons.filter(c => c.exportable).length, 0),
  distinct_conditions: Object.keys(byCondition).length, top_conditions: topConditions,
  note: 'Registry-scale reconstruction: every pairwise arm comparison (capped at ' + MAX_PAIRS + '/trial) '
    + 'across the harvested AACT Tier-A cohort, via the shipped engine (Titman-QP default). Index carries '
    + 'metadata only (NCT, condition, badge, method, reconstructed HR) — no patient-level rows.',
};
fs.writeFileSync(path.join(__dirname, '..', 'realipd', 'scale_index.json'), JSON.stringify({ summary, index }, null, 2));

const L = [];
L.push('# Scale — registry-wide reconstruction across ClinicalTrials.gov', '');
L.push('*The validated method run at scale: every pairwise arm comparison across the harvested AACT',
  'Tier-A cohort (Titman-QP default, `validate/scale_run.js`). Metadata index at `realipd/scale_index.json`;',
  'the pseudo-IPD itself is reproducible per trial (no patient rows committed). This quantifies what',
  'ClinicalTrials.gov can yield as reconstructed survival evidence today.*', '');
L.push('## At a glance', '');
L.push('- **' + nTrials + '** cohort trials scanned · **' + nSingleArm + '** single-arm (no comparison) · **' + nMultiArm + '** multi-arm.');
L.push('- **' + nReconstructable + '** trials yield ≥1 reconstructable comparison → **' + nComparisons + '** pairwise pseudo-IPD comparisons.');
L.push('- Self-audit badges: ' + Object.entries(badge).filter(([, v]) => v).map(([k, v]) => `${v} ${k}`).join(' · ') + ' · **' + summary.exportable_comparisons + ' exportable**.');
L.push('- Reconstruction method: ' + Object.entries(method).filter(([, v]) => v).map(([k, v]) => `${v} ${k}`).join(' · ') + ' (QP fires only where a total-event count is posted).');
L.push('- **' + summary.distinct_conditions + '** distinct conditions — survival evidence well beyond oncology.', '');
L.push('## Most-represented conditions', '');
L.push('| condition | reconstructable trials |', '|---|---|');
for (const [c, n] of topConditions) L.push(`| ${c} | ${n} |`);
L.push('', '*The multi-arm trials are where scale compounds: a 4-arm trial is up to 6 pairwise comparisons.',
  'The binding limit remains coverage (only a fraction of AACT posts a curve; see the census) and the',
  'posted event count (only a minority enables the QP; see the production gallery).*', '');
fs.writeFileSync(path.join(__dirname, '..', 'SCALE.md'), L.join('\n'));

console.log(JSON.stringify(summary, null, 2));
console.log('\nwrote realipd/scale_index.json + SCALE.md');
