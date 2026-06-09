#!/usr/bin/env node
/*
 * Cohort validation runner. Reconstructs every harvested trial JSON in a directory,
 * scores fidelity to the registry anchors, and (optionally) runs the AACT-only vs
 * digitization head-to-head when a matching digitized JSON is supplied.
 *
 * Usage:
 *   node validate/validate.js <cohort_dir> [--digitized <dir>] [-o report.json] [--bootstrap N]
 *
 * cohort_dir : directory of trial JSON files (from harvest_trial.py)
 * --digitized: directory with same-named JSONs whose km_points are digitized-from-figure
 *
 * Reports honest, scoped numbers — no universal-superiority wording. If a trial is Tier C it is
 * reported as "not reconstructable", never silently dropped.
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');
const M = require('./metrics.js');

function parseArgs(argv) {
  const a = { dir: null, digitized: null, out: 'validation_report.json', bootstrap: 300 };
  const rest = [];
  for (let i = 0; i < argv.length; i++) {
    if (argv[i] === '--digitized') a.digitized = argv[++i];
    else if (argv[i] === '-o') a.out = argv[++i];
    else if (argv[i] === '--bootstrap') a.bootstrap = parseInt(argv[++i], 10);
    else rest.push(argv[i]);
  }
  a.dir = rest[0];
  return a;
}

function main() {
  const args = parseArgs(process.argv.slice(2));
  if (!args.dir || !fs.existsSync(args.dir)) {
    console.error('usage: node validate/validate.js <cohort_dir> [--digitized <dir>] [-o report.json]');
    console.error('  (cohort_dir must contain harvested trial JSON files)');
    process.exit(2);
  }
  const files = fs.readdirSync(args.dir).filter(f => f.endsWith('.json'));
  if (!files.length) { console.error('no .json trial files in', args.dir); process.exit(2); }

  const opts = { bootstrap: args.bootstrap };
  const rows = [], fidsAll = [], reconstructable = [], notReconstructable = [];
  const h2h = [];
  for (const f of files) {
    let trial;
    try { trial = JSON.parse(fs.readFileSync(path.join(args.dir, f), 'utf8')); }
    catch (e) { console.error('skip (bad JSON):', f, e.message); continue; }
    const res = RIPD.reconstruct(trial, opts);
    const fid = M.fidelity(trial, res);
    rows.push({ file: f, nct: trial.nct_id, tier: res.tier, badge: res.audit.badge,
      logHR_err: fid.logHR_err, anchor_sup_error: fid.anchor_sup_error,
      wasserstein: fid.wasserstein_to_anchors, exportable: res.exportable });
    if (res.tier === 'C') notReconstructable.push(trial.nct_id);
    else { reconstructable.push(trial.nct_id); fidsAll.push(fid); }

    // head-to-head if a digitized counterpart exists
    if (args.digitized) {
      const dpath = path.join(args.digitized, f);
      if (fs.existsSync(dpath)) {
        const dig = JSON.parse(fs.readFileSync(dpath, 'utf8'));
        const hh = M.headToHead(trial, Object.assign({ digitizedTrial: dig }, opts));
        if (hh.digitization) {
          h2h.push({
            nct: trial.nct_id,
            aact_sup: hh.aact_only.anchor_sup_error, digi_sup: hh.digitization.anchor_sup_error,
            aact_w1: hh.aact_only.wasserstein_to_anchors, digi_w1: hh.digitization.wasserstein_to_anchors,
            aact_logHR_err: hh.aact_only.logHR_err, digi_logHR_err: hh.digitization.logHR_err,
          });
        }
      }
    }
  }

  const report = {
    cohort_dir: args.dir,
    n_trials: files.length,
    n_reconstructable: reconstructable.length,
    n_not_reconstructable: notReconstructable.length,
    not_reconstructable_nct: notReconstructable,
    aggregate: M.aggregate(fidsAll),
    per_trial: rows,
  };
  if (h2h.length) {
    const winsSup = h2h.filter(x => x.aact_sup != null && x.digi_sup != null && x.aact_sup < x.digi_sup).length;
    const winsW1 = h2h.filter(x => x.aact_w1 != null && x.digi_w1 != null && x.aact_w1 < x.digi_w1).length;
    report.head_to_head = {
      n: h2h.length,
      aact_better_anchor_sup: `${winsSup}/${h2h.length}`,
      aact_better_wasserstein: `${winsW1}/${h2h.length}`,
      detail: h2h,
      note: 'AACT-only is expected to win on anchor fidelity (exact registry values); this is the '
        + 'scoped claim, NOT universal superiority on between-anchor shape.',
    };
  }
  fs.writeFileSync(args.out, JSON.stringify(report, null, 2));
  console.log(`trials=${report.n_trials}  reconstructable=${report.n_reconstructable}  `
    + `tierC=${report.n_not_reconstructable}`);
  console.log('aggregate:', JSON.stringify(report.aggregate));
  if (report.head_to_head) console.log('head-to-head:', JSON.stringify({
    n: report.head_to_head.n, aact_better_anchor_sup: report.head_to_head.aact_better_anchor_sup,
    aact_better_wasserstein: report.head_to_head.aact_better_wasserstein }));
  console.log('wrote', args.out);
}

main();
