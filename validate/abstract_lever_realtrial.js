#!/usr/bin/env node
/*
 * REAL-TRIAL grounding of the censoring lever: RADIANT-4 (NCT01524783), validated vs the POSTED HR.
 * ====================================================================================================
 *
 * The abstract event-count lever (harvest/abstract_events.py) supplies the QP its per-arm total_events
 * from in-scope data. This script grounds the OTHER half of the chain — total_events -> QP -> recovered
 * HR — on a real trial where we hold the curve AND the truth:
 *
 *   - curve: this repo's harvested AACT record NCT01524783.json (exact posted KM-estimate anchors + N)
 *   - truth: the registry-/publication-posted Cox HR 0.48 (95% CI 0.35-0.67) for PFS (held out)
 *   - lever: the per-arm event count 107/77 (everolimus/placebo)
 *
 * We reconstruct two ways from the SAME exact anchors and score the Cox HR against 0.48:
 *   curve-only           : reconstruct with total_events stripped  (the identifiability trap)
 *   censoring-informed QP : reconstruct with the 107/77 event counts (the lever)
 *
 * Honest provenance note on the lever's SOURCE for THIS trial: RADIANT-4's Lancet abstract
 * (PMID 26703889) reports only ADVERSE-event "X of N" counts (correctly rejected by abstract_events ->
 * None) and the HR 0.48 (extracted by abstract_hr). So for RADIANT-4 the in-scope event count comes from
 * AACT participant-flow (107/77), and the abstract contributes the HR. The end-to-end abstract->count
 * extraction itself is grounded separately on DAPA-HF (PMID 31535829), whose abstract DOES post per-arm
 * outcome counts (386/2373 vs 502/2371) — see VALIDATION.md "Abstract event-count lever".
 *
 * Run from repo root:  node validate/abstract_lever_realtrial.js
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');
const _ = RIPD._;

const coxHR = (e, c) => _.coxLogHR(
  e.map(r => ({ time: r.time, status: r.status, x: 1 }))
    .concat(c.map(r => ({ time: r.time, status: r.status, x: 0 })))).hr;

function loadTrial(p) {
  const t = JSON.parse(fs.readFileSync(p, 'utf8'));
  const seen = new Set(), arms = [];
  for (const a of t.arms) {
    if (seen.has(a.arm_id) || !(a.km_points && a.N)) continue;
    seen.add(a.arm_id);
    const fu = Math.max(...a.km_points.map(p => p.t));
    arms.push({ arm_id: a.arm_id, label: a.label, role: a.role, N: a.N,
                total_events: a.total_events, follow_up_max: fu,
                km_points: a.km_points, nar_points: a.nar_points || [] });
  }
  return { nct_id: t.nct_id, time_unit: t.time_unit || 'months', arms, hr: t.hr };
}

function reconHR(trial, useEvents) {
  const t2 = JSON.parse(JSON.stringify(trial));
  if (!useEvents) t2.arms.forEach(a => { a.total_events = null; });
  const r = RIPD.reconstruct(t2, {});
  const e = r.arms.find(a => a.role === 'experimental');
  const c = r.arms.find(a => a.role === 'comparator');
  return coxHR(e.ipd, c.ipd);
}

function fold(h, truth) { return +(Math.max(h, truth) / Math.min(h, truth)).toFixed(3); }
function inCI(h, lo, hi) { return lo != null && hi != null && h >= lo && h <= hi; }

function run() {
  const trial = loadTrial(path.join(__dirname, '..', 'NCT01524783.json'));
  const posted = trial.hr.value, lo = trial.hr.ci_low, hi = trial.hr.ci_high;
  const co = reconHR(trial, false);
  const inf = reconHR(trial, true);
  const out = {
    trial: trial.nct_id,
    arms: trial.arms.map(a => ({ label: a.label, N: a.N, total_events: a.total_events })),
    posted_hr: posted, posted_ci: [lo, hi],
    curve_only: { hr: +co.toFixed(3), fold: fold(co, posted), inside_ci: inCI(co, lo, hi) },
    censoring_informed: { hr: +inf.toFixed(3), fold: fold(inf, posted), inside_ci: inCI(inf, lo, hi) },
  };
  const outPath = path.join(__dirname, 'abstract_lever_realtrial.json');
  fs.mkdirSync(path.dirname(outPath), { recursive: true });
  fs.writeFileSync(outPath, JSON.stringify(out, null, 2));
  console.log(`\n=== RADIANT-4 (${out.trial}) — censoring lever vs posted HR ===`);
  out.arms.forEach(a => console.log(`  ${a.label}: N=${a.N}, events=${a.total_events}`));
  console.log(`\n  posted Cox HR (truth)          : ${posted}  (95% CI ${lo}-${hi})`);
  console.log(`  curve-only (no event count)    : HR ${out.curve_only.hr}  fold ${out.curve_only.fold}  inside CI? ${out.curve_only.inside_ci}`);
  console.log(`  censoring-informed (107/77 QP) : HR ${out.censoring_informed.hr}  fold ${out.censoring_informed.fold}  inside CI? ${out.censoring_informed.inside_ci}`);
  console.log(`\n  wrote ${outPath}`);
  return out;
}

if (require.main === module) run();
module.exports = { run, loadTrial, reconHR };
