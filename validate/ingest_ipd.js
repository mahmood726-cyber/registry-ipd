#!/usr/bin/env node
/*
 * INGEST EXTERNAL TRUE-IPD (credentialed repositories: Vivli / YODA / Project Data Sphere).
 *
 * Item 3 of the roadmap: extend the gold standard to dozens–hundreds of trials using credentialed
 * patient-level data. The data themselves are DUA-protected and must be exported by the credentialed
 * user — this adapter is the plug-in path so that drop is a one-command operation, NOT a code edit.
 *
 * It reads a manifest (JSON) describing each trial's CSV + column mapping, builds the true IPD, runs
 * the SAME pipeline as the open gold standard (true HR/median/RMST → registry-style coarse summary →
 * curve-only + Titman-QP censoring-informed reconstruction → multiple-imputation uncertainty), and
 * scores each against truth. Output mirrors goldstandard.js so results merge directly.
 *
 * Supports two formats per dataset:
 *   format:"cdisc-adtte"  — CDISC ADaM time-to-event. time=AVAL, ARM=TRTP/TRT01P, and the CENSOR flag
 *                           CNSR is INVERTED (1=censored, 0=event) — handled here. (The #1 ingestion
 *                           bug for Vivli/YODA exports.)
 *   format:"generic"      — time / status (1=event) [/ eventVal] / arm columns, like goldstandard.js.
 *
 * Run:  node validate/ingest_ipd.js <manifest.json> [-o out.json]
 * Data and manifests live under realipd/credentialed/ (gitignored). See CREDENTIALED.md.
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');
const _ = RIPD._;

function parseCSV(file) {
  const txt = fs.readFileSync(file, 'utf8').replace(/^﻿/, '').trim().split(/\r?\n/);
  const head = splitCsvLine(txt[0]);
  return txt.slice(1).map(line => { const c = splitCsvLine(line); const o = {}; head.forEach((h, i) => o[h.trim()] = (c[i] != null ? c[i] : '').trim()); return o; });
}
// minimal RFC-4180-ish splitter (handles quoted commas, common in clinical exports)
function splitCsvLine(line) {
  const out = []; let cur = '', q = false;
  for (let i = 0; i < line.length; i++) { const ch = line[i];
    if (q) { if (ch === '"' && line[i + 1] === '"') { cur += '"'; i++; } else if (ch === '"') q = false; else cur += ch; }
    else { if (ch === '"') q = true; else if (ch === ',') { out.push(cur); cur = ''; } else cur += ch; } }
  out.push(cur); return out.map(s => s.replace(/^"|"$/g, ''));
}
const num = x => { const v = parseFloat(x); return isFinite(v) ? v : null; };

// build a 2-arm true IPD ({time,status:1=event}) from one dataset config
function loadArms(cfg, baseDir) {
  const rows = parseCSV(path.join(baseDir, cfg.csv));
  const cdisc = cfg.format === 'cdisc-adtte';
  const timeCol = cfg.time || (cdisc ? 'AVAL' : 'time');
  const armCol = cfg.arm || (cdisc ? 'TRTP' : 'arm');
  const toEvent = (r) => {
    if (cdisc || cfg.censorFlag) { const cn = num(r[cfg.censorFlag || 'CNSR']); return cn == null ? null : (cn === 0 ? 1 : 0); } // CNSR: 1=censored→0, 0=event→1
    const s = num(r[cfg.status || 'status']); if (s == null) return null;
    return cfg.eventVal != null ? (s === cfg.eventVal ? 1 : 0) : (s >= 1 ? 1 : 0);
  };
  const arm = (which) => rows.filter(r => String(r[armCol]) === String(which))
    .map(r => ({ time: num(r[timeCol]), status: toEvent(r) }))
    .filter(r => r.time != null && r.status != null && r.time > 0);
  const cap = (a) => { const CAP = 2500; if (a.length <= CAP) return a; const step = a.length / CAP, s = []; for (let i = 0; i < CAP; i++) s.push(a[Math.floor(i * step)]); return s; };
  return { expT: cap(arm(cfg.exp)), ctlT: cap(arm(cfg.ctl)) };
}

const rows2 = (a, b) => a.map(r => ({ time: r.time, status: r.status, x: 1 })).concat(b.map(r => ({ time: r.time, status: r.status, x: 0 })));
const coxHR = (a, b) => Math.exp(_.coxLogHR(rows2(a, b)).beta);
function coarse(ipd, K) {
  const km = _.kmFromIPD(ipd), tmax = 0.95 * Math.max(...ipd.map(r => r.time)), pts = [{ t: 0, S: 1 }];
  for (let i = 1; i <= K; i++) { const t = tmax * i / K; pts.push({ t: +t.toFixed(2), S: +_.evalKM(km, t).toFixed(4) }); }
  return { km_points: pts, nar_points: [], N: ipd.length, total_events: ipd.filter(r => r.status === 1).length, follow_up_max: +tmax.toFixed(2) };
}

function runDataset(cfg, baseDir) {
  let expT, ctlT;
  try { ({ expT, ctlT } = loadArms(cfg, baseDir)); } catch (e) { return { id: cfg.id, error: 'load: ' + e.message }; }
  if (expT.length < 20 || ctlT.length < 20) return { id: cfg.id, error: `too few rows (${expT.length}/${ctlT.length} per arm; need ≥20)` };
  const kmE = _.kmFromIPD(expT), kmC = _.kmFromIPD(ctlT);
  const tau = 0.9 * Math.min(Math.max(...expT.map(r => r.time)), Math.max(...ctlT.map(r => r.time)));
  const truth = { HR: coxHR(expT, ctlT), medE: _.medianFromKM(kmE, { interpolate: true }), RMSTd: _.rmst(kmE, tau) - _.rmst(kmC, tau) };
  const trial = { nct_id: cfg.id, time_unit: cfg.time_unit || 'd',
    arms: [Object.assign({ arm_id: 'exp', role: 'experimental' }, coarse(expT, cfg.K || 8)), Object.assign({ arm_id: 'ctl', role: 'comparator' }, coarse(ctlT, cfg.K || 8))] };
  // curve-only (strip events) and censoring-informed (QP default, keeps events)
  const t0 = JSON.parse(JSON.stringify(trial)); t0.arms.forEach(a => { a.total_events = null; });
  const rCu = RIPD.reconstruct(t0, {}), rCi = RIPD.reconstruct(trial, {});
  const fe = (ipds) => { const h = coxHR(ipds[0], ipds[1]); return { HR: +h.toFixed(3), fold: +Math.exp(Math.abs(Math.log(h) - Math.log(truth.HR))).toFixed(3) }; };
  const cu = rCu.arms ? fe([rCu.arms[0].ipd, rCu.arms[1].ipd]) : null;
  const ci = rCi.arms ? fe([rCi.arms[0].ipd, rCi.arms[1].ipd]) : null;
  let cover = null, width = null;
  if (expT.length >= 100 && ctlT.length >= 100) {
    try { const ens = RIPD.reconstructEnsemble(trial, { M: 200 }); const hr = ens.ensemble && ens.ensemble.hr;
      if (hr) { cover = truth.HR >= hr.lo && truth.HR <= hr.hi; width = +(hr.hi / hr.lo).toFixed(2); } } catch {}
  }
  return { id: cfg.id, label: cfg.label || cfg.id, source: cfg.source || null, n_exp: expT.length, n_ctl: ctlT.length,
    true_HR: +truth.HR.toFixed(3), curve_only: cu, censoring_informed_qp: ci, ci_covers_true_HR: cover, ci_width_fold: width };
}

function main(argv) {
  const manifestPath = argv.find(a => !a.startsWith('-'));
  if (!manifestPath) { console.error('usage: node validate/ingest_ipd.js <manifest.json> [-o out.json]'); return 2; }
  if (!fs.existsSync(manifestPath)) { console.error(`manifest not found: ${manifestPath}\nSee CREDENTIALED.md for the format and folder layout.`); return 2; }
  const manifest = JSON.parse(fs.readFileSync(manifestPath, 'utf8'));
  const baseDir = path.dirname(path.resolve(manifestPath));
  const outIdx = argv.indexOf('-o'); const outPath = outIdx >= 0 ? argv[outIdx + 1] : path.join(baseDir, 'ingest_results.json');
  const results = (manifest.datasets || []).map(c => runDataset(c, baseDir));
  const ok = results.filter(r => !r.error);
  const big = ok.filter(r => r.n_exp >= 100 && r.n_ctl >= 100);
  const med = (set, get) => { const xs = set.map(get).filter(x => x != null).sort((a, b) => a - b); return xs.length ? +xs[xs.length >> 1].toFixed(3) : null; };
  const w20 = (set, get) => set.filter(r => get(r) != null && get(r) < 1.2).length + '/' + set.length;
  const covered = ok.filter(r => r.ci_covers_true_HR != null);
  const summary = {
    source: manifest.source || 'mixed', n_datasets: ok.length, n_errored: results.length - ok.length, n_ge100_per_arm: big.length,
    curve_only_median_fold: med(big, r => r.curve_only && r.curve_only.fold), curve_only_within20: w20(big, r => r.curve_only && r.curve_only.fold),
    censoring_informed_qp_median_fold: med(big, r => r.censoring_informed_qp && r.censoring_informed_qp.fold), censoring_informed_qp_within20: w20(big, r => r.censoring_informed_qp && r.censoring_informed_qp.fold),
    uncertainty_coverage: covered.length ? `${covered.filter(r => r.ci_covers_true_HR).length}/${covered.length}` : 'n/a',
    note: 'External credentialed-IPD validation via the same pipeline as the open gold standard. '
      + 'CDISC ADTTE CNSR inversion handled (1=censored). Aggregates over >=100/arm trials.',
  };
  fs.writeFileSync(outPath, JSON.stringify({ summary, per_dataset: results }, null, 2));
  console.log(JSON.stringify(summary, null, 2));
  for (const r of results) console.log('  ' + (r.id || '?').padEnd(16) + (r.error ? 'ERR ' + r.error
    : `N ${r.n_exp}/${r.n_ctl}  true ${r.true_HR}  curve ${r.curve_only && r.curve_only.fold}  QP ${r.censoring_informed_qp && r.censoring_informed_qp.fold}` + (r.ci_covers_true_HR != null ? `  cover ${r.ci_covers_true_HR}` : '')));
  console.log(`\nwrote ${path.relative(process.cwd(), outPath).replace(/\\/g, '/')}`);
  return 0;
}
if (require.main === module) process.exit(main(process.argv.slice(2)));
else module.exports = { loadArms, runDataset };
