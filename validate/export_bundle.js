#!/usr/bin/env node
/*
 * EXPORT BUNDLE — materialise the reconstructed pseudo-IPD for every exportable comparison as a
 * downloadable, analysis-ready dataset (the actual deliverable: pseudo-IPD others can pool).
 *
 * For each exportable pairwise comparison in the cohort (Titman-QP default) we emit the reconstructed
 * per-patient rows in TIDY LONG format — the standard IPD-meta-analysis input:
 *     nct, comparison, arm (exp|ctl), time, status (1=event)
 * The full stacked dataset goes to dist/pseudo_ipd/registry_pseudo_ipd.csv (committed, ~6.5MB; also
 * reproducible). A committed MANIFEST (realipd/pseudo_ipd_manifest.json: one row per comparison with
 * NCT/condition/badge/method/HR/N/events/provenance — no patient rows) indexes it, and a small SAMPLE
 * (dist/pseudo_ipd/SAMPLE_*.csv) is written for the first few so the format is inspectable.
 *
 * The pseudo-IPD is RECONSTRUCTED (synthetic), so it is freely shareable; it is not real patient data.
 * Run from repo root: node validate/export_bundle.js [--max N]
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');

const COHORT = path.join(__dirname, '..', 'cohort');
const OUTDIR = path.join(__dirname, '..', 'dist', 'pseudo_ipd');
const files = fs.existsSync(COHORT) ? fs.readdirSync(COHORT).filter(f => f.endsWith('.json') && f.startsWith('NCT')) : [];
if (!files.length) { console.error('cohort/ empty (gitignored; re-harvest with harvest/harvest_cohort.py)'); process.exit(2); }
fs.mkdirSync(OUTDIR, { recursive: true });
const maxIdx = process.argv.indexOf('--max'); const MAX = maxIdx >= 0 ? parseInt(process.argv[maxIdx + 1], 10) : Infinity;
const MAX_PAIRS = 6;

const coxHR = (a, b) => Math.exp(RIPD._.coxLogHR(a.map(r => ({ time: r.time, status: r.status, x: 1 }))
  .concat(b.map(r => ({ time: r.time, status: r.status, x: 0 }))).slice()).beta);
const fnum = (x) => Number.isFinite(x) ? (Math.round(x * 1000) / 1000) : '';

const bigCsv = [path.join(OUTDIR, 'registry_pseudo_ipd.csv')];
const stream = fs.createWriteStream(bigCsv[0]); stream.write('nct,comparison,arm,time,status\n');
const manifest = []; let nComparisons = 0, nRows = 0, nSample = 0;

outer:
for (const f of files) {
  let t; try { t = JSON.parse(fs.readFileSync(path.join(COHORT, f), 'utf8')); } catch { continue; }
  const arms = (t.arms || []).filter(a => (a.km_points || []).length >= 3 && a.N != null);
  if (arms.length < 2) continue;
  const cond = (t.condition || 'unspecified').split(';')[0].trim().slice(0, 50);
  let pc = 0;
  for (let i = 0; i < arms.length && pc < MAX_PAIRS; i++) for (let j = i + 1; j < arms.length && pc < MAX_PAIRS; j++) {
    const trial2 = { nct_id: t.nct_id, time_unit: t.time_unit,
      arms: [Object.assign({}, arms[j], { role: 'experimental', arm_id: 'e' }), Object.assign({}, arms[i], { role: 'comparator', arm_id: 'c' })],
      hr: arms.length === 2 ? t.hr : null };
    let r; try { r = RIPD.reconstruct(trial2); } catch { continue; }
    if (r.tier !== 'A' || !r.arms || !r.exportable) continue;
    pc++; const cmp = nComparisons++;
    const e = r.arms.find(a => a.role === 'experimental'), c = r.arms.find(a => a.role === 'comparator');
    const hr = coxHR(e.ipd, c.ipd);
    const sampleRows = [];
    for (const [arm, ipd] of [['exp', e.ipd], ['ctl', c.ipd]]) for (const row of ipd) {
      const line = `${t.nct_id},${cmp},${arm},${fnum(row.time)},${row.status}\n`;
      stream.write(line); nRows++; if (nSample < 6) sampleRows.push(line);
    }
    if (nSample < 6) { fs.writeFileSync(path.join(OUTDIR, `SAMPLE_${t.nct_id}_cmp${cmp}.csv`), 'nct,comparison,arm,time,status\n' + sampleRows.join('')); nSample++; }
    manifest.push({ comparison: cmp, nct: t.nct_id, url: t.source_url, condition: cond,
      exp_label: (arms[j].label || 'exp').slice(0, 40), ctl_label: (arms[i].label || 'ctl').slice(0, 40),
      n_exp: e.ipd.length, n_ctl: c.ipd.length, events: e.ipd.filter(x => x.status === 1).length + c.ipd.filter(x => x.status === 1).length,
      reconstructed_HR: +hr.toFixed(3), registry_HR: t.hr && t.hr.value != null ? t.hr.value : null,
      badge: r.audit.badge, method: r.method, anchors: Math.min(arms[i].km_points.length, arms[j].km_points.length) });
    if (nComparisons >= MAX) { stream.end(); break outer; }
  }
}
stream.end();

const out = { generated: 'registry-wide reconstructed pseudo-IPD (Titman-QP default)', n_comparisons: nComparisons, n_patient_rows: nRows,
  format: 'tidy long: nct, comparison, arm (exp|ctl), time, status (1=event)',
  full_dataset: 'dist/pseudo_ipd/registry_pseudo_ipd.csv (committed; reproducible via this script)',
  sharing: 'pseudo-IPD is RECONSTRUCTED (synthetic) — freely shareable; not real patient data.',
  comparisons: manifest };
fs.writeFileSync(path.join(__dirname, '..', 'realipd', 'pseudo_ipd_manifest.json'), JSON.stringify(out, null, 2));
console.log(`wrote ${nComparisons} comparisons / ${nRows} pseudo-patient rows`);
console.log(`  full dataset: ${path.relative(process.cwd(), bigCsv[0]).replace(/\\/g, '/')} (committed)`);
console.log(`  manifest: realipd/pseudo_ipd_manifest.json (committed, metadata only)`);
console.log(`  ${nSample} SAMPLE_*.csv written for format inspection`);
