#!/usr/bin/env node
/*
 * REAL CANCER-COHORT slice: reconstruction on TCGA / cBioPortal overall-survival IPD.
 *
 * Separate from the main gold standard on purpose. These 8 TCGA cohorts are split by SEX, which in
 * most cancers is a NEAR-NULL survival contrast (true HR ~0.7–1.2). For near-null effects the *fold*
 * (ratio) error used in the main gold standard is the wrong metric — a tiny absolute log-HR wobble
 * looks enormous as a ratio. So this slice scores the **absolute log-HR error** (the honest metric
 * when there is essentially no signal to recover), and it is deliberately NOT folded into the
 * headline fold-error aggregate, which is for adequately-powered non-null contrasts.
 *
 * Purpose: show the method on real, large oncology IPD from an open repository, and characterise its
 * precision in the low-signal regime (where it is least useful and should be least trusted).
 *
 * Data: harvest/fetch_cbioportal.js (open cBioPortal API; CSVs gitignored). Run from repo root:
 *   node validate/goldstandard_cbio.js   ->   realipd/goldstandard_cbio_results.json
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');
const _ = RIPD._;
const GS = require('./goldstandard.js');

const CBIO = [
  { ds: 'cbio_luad', label: 'TCGA lung adenocarcinoma (OS, M/F)' },
  { ds: 'cbio_coadread', label: 'TCGA colorectal (OS, M/F)' },
  { ds: 'cbio_stad', label: 'TCGA stomach (OS, M/F)' },
  { ds: 'cbio_lihc', label: 'TCGA liver HCC (OS, M/F)' },
  { ds: 'cbio_kirc', label: 'TCGA kidney clear-cell (OS, M/F)' },
  { ds: 'cbio_hnsc', label: 'TCGA head & neck (OS, M/F)' },
  { ds: 'cbio_skcm', label: 'TCGA melanoma (OS, M/F)' },
  { ds: 'cbio_blca', label: 'TCGA bladder (OS, M/F)' },
].map(c => Object.assign(c, { time: 'time', status: 'status', arm: 'sex', exp: 'Male', ctl: 'Female' }));

function coarse(ipd, K) {
  const km = _.kmFromIPD(ipd), tmax = 0.95 * Math.max(...ipd.map(r => r.time));
  const pts = [{ t: 0, S: 1 }];
  for (let i = 1; i <= K; i++) { const t = tmax * i / K; pts.push({ t: +t.toFixed(2), S: +_.evalKM(km, t).toFixed(4) }); }
  return { km_points: pts, nar_points: [], N: ipd.length, total_events: ipd.filter(r => r.status === 1).length, follow_up_max: +tmax.toFixed(2) };
}

const rows = [];
for (const cfg of CBIO) {
  let arms; try { arms = GS.loadArms(cfg); } catch (e) { rows.push({ ds: cfg.ds, error: String(e.message).slice(0, 80) }); continue; }
  const { expT, ctlT } = arms;
  if (expT.length < 100 || ctlT.length < 100) { rows.push({ ds: cfg.ds, error: 'too few rows' }); continue; }
  const trueHR = GS.coxHR(expT, ctlT).hr;
  const trial = { nct_id: cfg.ds, time_unit: 'months',
    arms: [Object.assign({ arm_id: 'exp', label: 'Male', role: 'experimental' }, coarse(expT, 8)),
           Object.assign({ arm_id: 'ctl', label: 'Female', role: 'comparator' }, coarse(ctlT, 8))] };
  const r = RIPD.reconstruct(trial, {});
  if (!r.arms) { rows.push({ ds: cfg.ds, error: 'no reconstruction' }); continue; }
  const e = r.arms.find(a => a.role === 'experimental'), c = r.arms.find(a => a.role === 'comparator');
  const reconHR = GS.coxHR(e.ipd, c.ipd).hr;
  rows.push({ ds: cfg.ds, label: cfg.label, n_exp: expT.length, n_ctl: ctlT.length,
    true_HR: +trueHR.toFixed(3), recon_HR: +reconHR.toFixed(3),
    abs_logHR_err: +Math.abs(Math.log(reconHR) - Math.log(trueHR)).toFixed(3),
    direction_agrees: (trueHR - 1) === 0 || ((reconHR - 1) >= 0) === ((trueHR - 1) >= 0) });
}

const ok = rows.filter(r => !r.error);
const errs = ok.map(r => r.abs_logHR_err).sort((a, b) => a - b);
const mean = errs.reduce((a, b) => a + b, 0) / (errs.length || 1);
const summary = {
  n_cohorts: ok.length,
  metric: 'absolute log-HR error (appropriate for near-null effects; NOT fold-error)',
  mean_abs_logHR_err: +mean.toFixed(3),
  median_abs_logHR_err: +(errs[errs.length >> 1] || 0).toFixed(3),
  max_abs_logHR_err: +(errs[errs.length - 1] || 0).toFixed(3),
  direction_agrees: `${ok.filter(r => r.direction_agrees).length}/${ok.length}`,
  true_HR_range: ok.length ? `${Math.min(...ok.map(r => r.true_HR))}–${Math.max(...ok.map(r => r.true_HR))}` : 'n/a',
  finding: 'On 8 real TCGA cohorts split by sex (a near-null contrast, true HR ~0.7-1.2), curve-only '
    + 'reconstruction recovers the log-HR to a mean absolute error of ~' + mean.toFixed(2) + '. The '
    + 'method has limited precision in this low-signal regime and can drift — it should not be trusted '
    + 'to detect small effects. This is a stress test outside the adequately-powered non-null regime '
    + 'the headline targets; reported separately for that reason.',
};
fs.writeFileSync(path.join(GS.dir, 'goldstandard_cbio_results.json'), JSON.stringify({ summary, per_cohort: rows }, null, 2));
console.log(JSON.stringify(summary, null, 2));
console.log('\nper cohort:');
for (const r of rows) console.log('  ' + r.ds.padEnd(15) + (r.error ? 'ERR ' + r.error : `true ${r.true_HR}  recon ${r.recon_HR}  |logHR err| ${r.abs_logHR_err}  dir ${r.direction_agrees ? 'ok' : 'FLIP'}`));
