/*
 * COVERAGE CENSUS over the harvested AACT cohort (cohort/*.json, gitignored).
 *
 * Answers the binding-limitation question with a precise number instead of "hundreds of trials":
 *   - Of trials that posted *some* survival result, how many are actually reconstructable?
 *   - What is the anchor-density (KM-timepoint) distribution, and how many clear the >=5-6
 *     timepoint reliability standard established by the operating curve in PAPER.md?
 *   - What does the engine's own self-audit say (gold/silver/bronze/none, exportable fraction)?
 *
 * Uses the SAME engine logic the tool ships (classifyTier + reconstruct + selfAudit) — no
 * reinvented classifier. Deterministic (reconstruct() with ensemble disabled).
 *
 * Run from repo root:  node validate/census_cohort.js [cohortDir]   (default: cohort)
 * Writes: realipd/coverage_census.json
 */
const fs = require('fs');
const path = require('path');
const engine = require('../src/engine.js');

const COHORT = process.argv[2] || path.join(__dirname, '..', 'cohort');
const OUT = path.join(__dirname, '..', 'realipd', 'coverage_census.json');

if (!fs.existsSync(COHORT)) {
  console.error(`ERROR: cohort dir not found: ${COHORT}`);
  console.error('cohort/ is gitignored; re-harvest with harvest/harvest_cohort.py or pass a dir.');
  process.exit(2);
}

const files = fs.readdirSync(COHORT).filter(f => f.endsWith('.json'));
if (!files.length) { console.error(`ERROR: no .json trials in ${COHORT}`); process.exit(2); }

// binding anchor density = fewest KM timepoints among arms that carry data (weakest arm limits
// a two-arm reconstruction). Only count arms with a population N.
function bindingK(trial) {
  const arms = (trial.arms || []).filter(a => a && a.N != null && (a.km_points || []).length > 0);
  if (arms.length < 1) return 0;
  return Math.min(...arms.map(a => a.km_points.length));
}
function maxK(trial) {
  const arms = (trial.arms || []);
  const ks = arms.map(a => (a.km_points || []).length);
  return ks.length ? Math.max(...ks) : 0;
}

const rows = [];
let parseErrors = 0;
for (const f of files) {
  let trial;
  try { trial = JSON.parse(fs.readFileSync(path.join(COHORT, f), 'utf8')); }
  catch (e) { parseErrors++; continue; }
  const tier = engine.classifyTier(trial);
  const nArms = (trial.arms || []).length;
  const kBind = bindingK(trial);
  const kMax = maxK(trial);
  const hasHR = !!(trial.hr && trial.hr.value != null);
  const anyNar = (trial.arms || []).some(a => (a.nar_points || []).length > 0);

  let badge = 'none', exportable = false, hardFail = null, recErr = null;
  try {
    const r = engine.reconstruct(trial, { ensemble: false });
    badge = (r.audit && r.audit.badge) || 'none';
    exportable = !!r.exportable;
    hardFail = r.audit && r.audit.hardFail || false;
  } catch (e) { recErr = String(e && e.message || e).slice(0, 120); }

  rows.push({ nct: trial.nct_id || f.replace('.json', ''), nArms, tier, kBind, kMax,
    hasHR, anyNar, badge, exportable, hardFail, recErr });
}

const N = rows.length;
const by = (pred) => rows.filter(pred).length;
const tierCounts = { A: by(r => r.tier === 'A'), B: by(r => r.tier === 'B'), C: by(r => r.tier === 'C') };
const badgeCounts = { gold: by(r => r.badge === 'gold'), silver: by(r => r.badge === 'silver'),
  bronze: by(r => r.badge === 'bronze'), none: by(r => r.badge === 'none') };

// anchor-density histogram on the binding K, Tier A only (the reconstructable population)
const tierA = rows.filter(r => r.tier === 'A');
const kHist = {};
for (const r of tierA) { const k = r.kBind; const bucket = k >= 12 ? '12+' : String(k); kHist[bucket] = (kHist[bucket] || 0) + 1; }

// reliability bands from the operating curve in PAPER.md (HR fold-error e(K))
//   K>=6 reliable (e~1.10), K 5 borderline (e~1.15), K 3-4 weak (e~1.2-1.4), K<3 not Tier A
const reliable = by(r => r.tier === 'A' && r.kBind >= 6);
const borderline = by(r => r.tier === 'A' && r.kBind === 5);
const weak = by(r => r.tier === 'A' && r.kBind >= 3 && r.kBind <= 4);
const structuredNAR = by(r => r.anyNar);

const pct = (n) => Math.round((1000 * n) / N) / 10;
const summary = {
  generated_from: path.relative(path.join(__dirname, '..'), COHORT).replace(/\\/g, '/'),
  universe_trials: N,
  parse_errors: parseErrors,
  tier_counts: tierCounts,
  tier_pct: { A: pct(tierCounts.A), B: pct(tierCounts.B), C: pct(tierCounts.C) },
  badge_counts: badgeCounts,
  exportable: by(r => r.exportable),
  exportable_pct: pct(by(r => r.exportable)),
  structured_number_at_risk_trials: structuredNAR,
  reconstructable_tierA: tierCounts.A,
  reliable_K_ge_6: reliable,
  reliable_K_ge_6_pct: pct(reliable),
  borderline_K_eq_5: borderline,
  weak_K_3_to_4: weak,
  anchor_density_hist_bindingK_tierA: kHist,
  reconstruction_errors: by(r => r.recErr != null),
  note: 'Binding K = fewest KM timepoints among arms with N (weakest arm limits a 2-arm '
    + 'reconstruction). Reliability bands from the operating curve e(K) in PAPER.md: K>=6 reliable, '
    + 'K=5 borderline, K=3-4 weak. Structured NAR count corroborates the AACT zero-NAR finding. '
    + 'badge/exportable from the shipped engine self-audit (ensemble disabled for determinism).',
};

fs.mkdirSync(path.dirname(OUT), { recursive: true });
fs.writeFileSync(OUT, JSON.stringify({ summary, per_trial: rows }, null, 2));

console.log(JSON.stringify(summary, null, 2));
console.log(`\nwrote ${path.relative(process.cwd(), OUT).replace(/\\/g, '/')} (${N} trials)`);
