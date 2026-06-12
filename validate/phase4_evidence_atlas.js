/*
 * PHASE 4 — the EVIDENCE-COMPLETENESS ATLAS (SYNTHESIS-VISION.md Sec 5).
 *
 * The capstone artefact of the granularity-manifold reframing: the "information map" of a body of
 * evidence for a target estimand (here the hazard ratio), computed BEFORE any pooling. For every
 * harvested trial we ask not "did it publish?" but "what does what-it-published *identify* about the
 * HR?", and place it on the manifold:
 *
 *   point_nar    curve + numbers-at-risk  -> HR point-identified (Guyot)         delta ~ 0.02
 *   point_qp     curve + total events     -> HR point-identified (Titman-QP)     delta ~ 0.042  (Phase 2b)
 *   point_hr     a posted HR (+ CI)       -> logHR point-identified directly     delta_recon = 0
 *   partial_curve  curve only, no lever   -> HR PARTIALLY identified (a region)  delta ~ 0.219 x anchorFactor(K)
 *   none         < 2 arms / no curve & no HR -> not contributable to an HR synthesis
 *
 * Each contributable trial gets an information score s^2/(s^2 + delta_recon^2) in (0,1]: the fraction of
 * its full-IPD HR information that survives the reconstruction/identification uncertainty (delta_recon is
 * the EXTRA width beyond sampling). The atlas then reports what fraction of the evidence -- by trial count
 * AND information-weighted (inverse-variance) -- is point- vs partially-identified, plus the fraction that
 * is time-resolvable (a curve is present, so a non-PH / Sec 4i analysis is even possible).
 *
 * All reconstruction half-widths are the MEASURED values from earlier phases (Phase 2b SDs, the operating
 * curve e(K) in PAPER.md), not invented. Reads the same harvested cohort the coverage census uses.
 *
 * Run from repo root:  node validate/phase4_evidence_atlas.js [cohortDir]   (default: cohort)
 * Writes: realipd/evidence_atlas.json
 */
const fs = require('fs');
const path = require('path');

const COHORT = process.argv[2] || path.join(__dirname, '..', 'cohort');
const OUT = path.join(__dirname, '..', 'realipd', 'evidence_atlas.json');

// --- measured reconstruction half-widths (log-HR SD), with provenance ---
const DELTA_NAR = 0.02;      // dense-anchor Guyot near-exact (validation HR median ~3%)
const DELTA_QP = 0.042;      // Phase 2b event-pinned reconstruction SD (validate/phase2b_lever_shrinks_r2.js)
const DELTA_CURVE = 0.219;   // Phase 2b curve-only reconstruction SD (the under-identified regime)
// operating curve e(K) from PAPER.md: K>=6 reliable (~1.10), K=5 borderline (~1.15), K3-4 weak (~1.3-1.4)
function anchorFactor(kBind) {
  if (kBind >= 6) return 1.0;
  if (kBind === 5) return 1.15;
  if (kBind >= 3) return 1.4;
  return 2.0;                // K<3: barely reconstructable, region is widest
}

function bindingK(arms) {
  const a = arms.filter(x => x && x.N != null && (x.km_points || []).length > 0);
  if (!a.length) return 0;
  return Math.min(...a.map(x => x.km_points.length));
}
function lastS(arm) {
  const km = arm.km_points || [];
  if (!km.length) return null;
  let best = km[0];
  for (const p of km) if (p.t != null && best.t != null && p.t >= best.t) best = p;
  return best && best.S != null ? best.S : null;
}
// arm event count: posted total_events, else curve-implied N*(1-S_last) (>=1), else null
function armEvents(arm) {
  if (arm.total_events != null && arm.total_events >= 0) return Math.max(arm.total_events, 0.5);
  const S = lastS(arm);
  if (arm.N != null && S != null && S >= 0 && S <= 1) return Math.max(arm.N * (1 - S), 0.5);
  return null;
}
// sampling variance of logHR ~ 1/E1 + 1/E2 (arm-based); null if either arm's events are unknown
function samplingVar(arms) {
  const ev = arms.map(armEvents).filter(e => e != null);
  if (ev.length < 2) return null;
  ev.sort((a, b) => a - b);
  const [e1, e2] = [ev[0], ev[ev.length - 1]];   // weakest + strongest arm of the contrast
  return 1 / e1 + 1 / e2;
}

function hrSE(trial) {     // posted-HR sampling SE from its own reported CI (log scale)
  const h = trial.hr;
  if (!h || h.value == null || h.ci_low == null || h.ci_high == null) return null;
  if (h.ci_low <= 0 || h.ci_high <= 0 || h.ci_high <= h.ci_low) return null;
  return (Math.log(h.ci_high) - Math.log(h.ci_low)) / (2 * 1.96);
}

function classify(trial) {
  const arms = (trial.arms || []);
  const twoArm = arms.length >= 2;
  const anyNar = arms.some(a => (a.nar_points || []).length > 0);
  const allEvents = twoArm && arms.every(a => a.total_events != null);
  const hasHR = !!(trial.hr && trial.hr.value != null);
  const allKm = twoArm && arms.every(a => (a.km_points || []).length > 0);
  const kBind = bindingK(arms);
  const sVar = samplingVar(arms);

  let tier, deltaRecon, sUse;
  if (twoArm && anyNar) { tier = 'point_nar'; deltaRecon = DELTA_NAR; sUse = sVar; }
  else if (allEvents) { tier = 'point_qp'; deltaRecon = DELTA_QP; sUse = sVar; }
  else if (hasHR) {
    tier = 'point_hr'; deltaRecon = 0.0;
    const se = hrSE(trial); sUse = se != null ? se * se : (sVar != null ? sVar : 0.04 * 0.04);
  } else if (allKm) {
    tier = 'partial_curve'; deltaRecon = DELTA_CURVE * anchorFactor(kBind); sUse = sVar;
  } else {
    tier = 'none'; deltaRecon = null; sUse = null;
  }

  const contributable = tier !== 'none' && sUse != null;
  const infoScore = contributable ? sUse / (sUse + deltaRecon * deltaRecon) : null;
  const weight = contributable ? 1 / (sUse + deltaRecon * deltaRecon) : null;   // inverse-variance HR weight
  const timeResolvable = allKm;     // a curve exists -> non-PH / Sec 4i analysis is possible

  return {
    nct: trial.nct_id || null, condition: trial.condition || null,
    nArms: arms.length, kBind, anyNar, allEvents, hasHR, allKm,
    tier, deltaRecon: deltaRecon == null ? null : round(deltaRecon, 4),
    samplingVar: sUse == null ? null : round(sUse, 5),
    infoScore: infoScore == null ? null : round(infoScore, 4),
    weight: weight == null ? null : round(weight, 3),
    timeResolvable, contributable,
  };
}

function round(x, n) { const f = 10 ** n; return Math.round(x * f) / f; }

function summarize(rows) {
  const POINT = new Set(['point_nar', 'point_qp', 'point_hr']);
  const contributable = rows.filter(r => r.contributable);
  const tierCounts = {};
  for (const t of ['point_nar', 'point_qp', 'point_hr', 'partial_curve', 'none']) tierCounts[t] = rows.filter(r => r.tier === t).length;
  const nContrib = contributable.length || 1;
  const pointTrials = contributable.filter(r => POINT.has(r.tier)).length;
  const partialTrials = contributable.filter(r => r.tier === 'partial_curve').length;
  const Wtot = contributable.reduce((s, r) => s + r.weight, 0) || 1;
  const Wpoint = contributable.filter(r => POINT.has(r.tier)).reduce((s, r) => s + r.weight, 0);
  const Wpartial = contributable.filter(r => r.tier === 'partial_curve').reduce((s, r) => s + r.weight, 0);
  const scores = contributable.map(r => r.infoScore).sort((a, b) => a - b);
  const median = scores.length ? scores[Math.floor((scores.length - 1) / 2)] : null;
  const timeResolvable = rows.filter(r => r.timeResolvable).length;
  return {
    n_trials: rows.length,
    n_contributable_to_HR: contributable.length,
    tier_counts: tierCounts,
    point_identified_trials: pointTrials,
    partial_identified_trials: partialTrials,
    point_identified_trial_pct: round((100 * pointTrials) / nContrib, 1),
    partial_identified_trial_pct: round((100 * partialTrials) / nContrib, 1),
    point_identified_info_weighted_pct: round((100 * Wpoint) / Wtot, 1),
    partial_identified_info_weighted_pct: round((100 * Wpartial) / Wtot, 1),
    median_info_score: median == null ? null : round(median, 3),
    time_resolvable_trials: timeResolvable,
    time_resolvable_pct: round((100 * timeResolvable) / (rows.length || 1), 1),
  };
}

function main() {
  if (!fs.existsSync(COHORT)) {
    console.error(`ERROR: cohort dir not found: ${COHORT} (gitignored; re-harvest or pass a dir).`);
    process.exit(2);
  }
  const files = fs.readdirSync(COHORT).filter(f => f.endsWith('.json'));
  const rows = [];
  let parseErrors = 0, skippedNonTrial = 0;
  for (const f of files) {
    let trial;
    try { trial = JSON.parse(fs.readFileSync(path.join(COHORT, f), 'utf8')); }
    catch (e) { parseErrors++; continue; }
    if (!trial || typeof trial !== 'object' || Array.isArray(trial) || !('arms' in trial)) { skippedNonTrial++; continue; }
    rows.push(classify(trial));
  }

  const overall = summarize(rows);
  // per-condition "one review question" view (only conditions with >= 4 contributable trials)
  const byCond = {};
  for (const r of rows) { if (!r.condition) continue; (byCond[r.condition] = byCond[r.condition] || []).push(r); }
  const perCondition = Object.entries(byCond)
    .map(([cond, rs]) => ({ condition: cond, ...summarize(rs) }))
    .filter(c => c.n_contributable_to_HR >= 4)
    .sort((a, b) => b.n_contributable_to_HR - a.n_contributable_to_HR);

  const out = {
    summary: {
      generated_from: path.relative(path.join(__dirname, '..'), COHORT).replace(/\\/g, '/'),
      files_scanned: files.length, parse_errors: parseErrors, skipped_non_trial: skippedNonTrial,
      estimand: 'hazard ratio (log-HR)',
      reconstruction_halfwidths_logHR: { point_nar: DELTA_NAR, point_qp: DELTA_QP, curve_only_base: DELTA_CURVE,
        anchor_factor: 'K>=6:1.0, K=5:1.15, K3-4:1.4, K<3:2.0 (operating curve e(K), PAPER.md)' },
      ...overall,
      note: 'Information map for the HR estimand BEFORE pooling. delta_recon = measured reconstruction SD by '
        + 'granularity (Phase 2b: curve-only 0.219, event-pinned 0.042; Guyot dense ~0.02; posted-HR 0). '
        + 'infoScore = s^2/(s^2+delta_recon^2) = fraction of full-IPD HR information surviving identification. '
        + 'Binding finding: structured numbers-at-risk are essentially absent in AACT (point_nar ~ 0), so HR '
        + 'point-identification rests on a posted HR or a posted event count; the curve-only majority is only '
        + 'PARTIALLY identified for HR -- the exact population reconstruction-with-UQ (Phases 1-3c) contributes '
        + 'honestly. Denominator caveat: registry->publication linkage is ~63.6% (lessons.md), so the true '
        + 'evidence base is larger than any registry census.',
    },
    per_condition: perCondition,
    per_trial: rows,
  };
  fs.mkdirSync(path.dirname(OUT), { recursive: true });
  fs.writeFileSync(OUT, JSON.stringify(out, null, 2));

  const s = overall;
  console.log('=== Phase 4: evidence-completeness atlas (HR estimand, before pooling) ===\n');
  console.log(`  corpus: ${s.n_trials} trials (${out.summary.files_scanned} files), ${s.n_contributable_to_HR} contributable to an HR contrast\n`);
  console.log('  identification tier (granularity manifold):');
  for (const [t, lab] of [['point_nar', 'curve + numbers-at-risk (Guyot)'], ['point_qp', 'curve + total events (Titman-QP)'],
    ['point_hr', 'posted HR + CI'], ['partial_curve', 'curve only -> PARTIALLY identified'], ['none', 'not contributable']]) {
    console.log(`    ${t.padEnd(15)} ${String(s.tier_counts[t]).padStart(4)}   ${lab}`);
  }
  console.log('');
  console.log(`  HR point-identified:   ${s.point_identified_trials} trials (${s.point_identified_trial_pct}%),  ${s.point_identified_info_weighted_pct}% of the poolable HR information`);
  console.log(`  HR partially-identified: ${s.partial_identified_trials} trials (${s.partial_identified_trial_pct}%),  ${s.partial_identified_info_weighted_pct}% of the poolable HR information`);
  console.log(`  median per-trial information score: ${s.median_info_score}   |   time-resolvable (curve present): ${s.time_resolvable_trials} (${s.time_resolvable_pct}%)`);
  console.log('\n  Reading: the map is computed BEFORE any pooling -- it says, per trial, what the posted statistics');
  console.log('  identify about the HR. Numbers-at-risk are essentially absent, so most trials are only PARTIALLY');
  console.log('  identified for HR (a region, not a point) -- the population the reconstruction-with-UQ pipeline');
  console.log('  (Phases 1-3c) lets join a synthesis honestly instead of being silently dropped to AD.');
  console.log(`\n  per-condition breakdown for ${out.per_condition.length} questions (>=4 contributable trials):`);
  for (const c of out.per_condition.slice(0, 8)) {
    console.log(`    ${c.condition.slice(0, 26).padEnd(28)} k=${String(c.n_contributable_to_HR).padStart(2)}  point ${String(c.point_identified_trial_pct).padStart(5)}%  partial ${String(c.partial_identified_trial_pct).padStart(5)}%  med-info ${c.median_info_score}`);
  }
  console.log(`\n  wrote ${path.relative(process.cwd(), OUT).replace(/\\/g, '/')}`);
  return out;
}

if (require.main === module) main();
module.exports = { classify, summarize };
