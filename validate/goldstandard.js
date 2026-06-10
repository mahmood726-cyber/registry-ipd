#!/usr/bin/env node
/*
 * GOLD-STANDARD validation on TRUE patient-level RCT IPD (openly published datasets).
 *
 * For each real dataset we: (1) compute the TRUE estimates (median, Cox HR, RMST) from the full
 * patient-level data; (2) generate the REGISTRY-STYLE coarse summary a sponsor would post to ct.gov
 * — KM survival at ~8 timepoints + N + total events per arm; (3) reconstruct pseudo-IPD from that
 * summary alone; (4) compare reconstruction vs TRUTH. This is real-data validation (real hazard
 * shapes, censoring, non-PH), not synthetic, and not circular: the engine never sees the true IPD.
 *
 * Datasets (open, Rdatasets mirror of R packages). Re-download (public), by source package:
 *   base=https://vincentarelbundock.github.io/Rdatasets/csv
 *   for ds in rotterdam gbsg colon veteran pbc diabetic nwtco myeloid kidtran retinopathy mgus2 \
 *             flchain nafld1 cancer ovarian; do curl -L -o realipd/$ds.csv $base/survival/$ds.csv; done
 *   for ds in alloauto larynx burn pneumon bfeed; do curl -L -o realipd/$ds.csv $base/KMsurv/$ds.csv; done
 *   for ds in prostateSurvival pharmacoSmoking hepatoCellular; do curl -L -o realipd/$ds.csv $base/asaur/$ds.csv; done
 *   # cbio_* (7 TCGA cohorts, late vs early stage): node harvest/fetch_cbioportal.js  (open cBioPortal API)
 *   # remaining (udca2, gehan, tongue, bmt, melanoma, ebmt1/3, aidssi, kidney, etc.) live under
 *   # realipd/ already or in their respective packages; see CONFIGS below for the per-dataset source.
 * Usage: node validate/goldstandard.js [realipd_dir]
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');
const _ = RIPD._;

const dir = process.argv[2] && !process.argv[2].startsWith('-') ? process.argv[2] : path.join(__dirname, '..', 'realipd');

function parseCSV(file) {
  const txt = fs.readFileSync(file, 'utf8').trim().split(/\r?\n/);
  const head = txt[0].split(',');
  return txt.slice(1).map(line => {
    // simple CSV (these datasets have no quoted commas except colon rx which has no comma)
    const cells = line.split(',');
    const o = {}; head.forEach((h, i) => { o[h.replace(/"/g, '')] = (cells[i] || '').replace(/"/g, ''); });
    return o;
  });
}
const num = x => { const v = parseFloat(x); return isFinite(v) ? v : null; };

// config: dataset -> how to extract a 2-arm time-to-event comparison
const CONFIGS = [
  { ds: 'gbsg', label: 'GBSG breast cancer (recurrence-free survival, hormonal Rx)', time: 'rfstime', status: 'status', arm: 'hormon', exp: '1', ctl: '0' },
  { ds: 'veteran', label: 'Veterans lung cancer (overall survival, treatment)', time: 'time', status: 'status', arm: 'trt', exp: '2', ctl: '1' },
  { ds: 'rotterdam', label: 'Rotterdam breast cancer (overall survival, hormonal Rx)', time: 'dtime', status: 'death', arm: 'hormon', exp: '1', ctl: '0' },
  { ds: 'pbc', label: 'PBC primary biliary cirrhosis (OS, D-penicillamine vs placebo)', time: 'time', status: 'status', eventVal: 2, arm: 'trt', exp: '2', ctl: '1' },
  { ds: 'diabetic', label: 'Diabetic retinopathy (time to vision loss, laser vs control)', time: 'time', status: 'status', arm: 'trt', exp: '1', ctl: '0' },
  { ds: 'nwtco', label: 'NWTSG Wilms tumour (relapse, histology group)', time: 'edrel', status: 'rel', arm: 'histol', exp: '2', ctl: '1' },
  { ds: 'myeloid', label: 'Myeloid AML RCT (overall survival, trt A vs B)', time: 'futime', status: 'death', arm: 'trt', exp: 'B', ctl: 'A' },
  { ds: 'kidtran', label: 'Kidney transplant (graft survival, by sex)', time: 'time', status: 'delta', arm: 'gender', exp: '2', ctl: '1' },
  { ds: 'udca2', label: 'UDCA in PBC RCT (progression, UDCA vs placebo)', time: 'futime', status: 'status', arm: 'trt', exp: '1', ctl: '0' },
  { ds: 'retinopathy', label: 'Diabetic retinopathy (laser, time to vision loss)', time: 'futime', status: 'status', arm: 'trt', exp: '1', ctl: '0' },
  { ds: 'mgus2', label: 'MGUS cohort (death, by sex)', time: 'futime', status: 'death', arm: 'sex', exp: 'M', ctl: 'F' },
  { ds: 'flchain', label: 'Free light chain cohort (death, by sex)', time: 'futime', status: 'death', arm: 'sex', exp: 'M', ctl: 'F' },
  { ds: 'nafld1', label: 'NAFLD cohort (death, by sex)', time: 'futime', status: 'status', arm: 'male', exp: '1', ctl: '0' },
  { ds: 'prostateSurvival', label: 'Prostate cancer (cancer death, grade)', time: 'survTime', status: 'status', eventVal: 1, arm: 'grade', exp: 'poor', ctl: 'mode' },
  { ds: 'bmt', label: 'Bone marrow transplant (DFS, risk group)', time: 't2', status: 'd3', arm: 'group', exp: '3', ctl: '1' },
  { ds: 'melanoma', label: 'Melanoma (cancer death, ulceration)', time: 'time', status: 'status', eventVal: 1, arm: 'ulcer', exp: '1', ctl: '0' },
  { ds: 'pharmacoSmoking', label: 'Smoking-cessation RCT (time to relapse)', time: 'ttr', status: 'relapse', arm: 'grp', exp: 'patchOnly', ctl: 'combination' },
  { ds: 'gehan', label: 'Gehan leukemia RCT (remission, 6-MP vs control)', time: 'time', status: 'cens', arm: 'treat', exp: '6-MP', ctl: 'control' },
  { ds: 'tongue', label: 'Tongue cancer (death, ploidy)', time: 'time', status: 'delta', arm: 'type', exp: '2', ctl: '1' },
  { ds: 'cancer', label: 'NCCTG lung cancer (OS, by sex)', time: 'time', status: 'status', eventVal: 2, arm: 'sex', exp: '2', ctl: '1' },
  { ds: 'ebmt3', label: 'EBMT transplant (relapse-free survival, T-cell depletion)', time: 'rfstime', status: 'rfsstat', arm: 'tcd', exp: 'TCD', ctl: 'No TCD' },
  { ds: 'aidssi', label: 'AIDS cohort (time to AIDS, CCR5 genotype)', time: 'time', status: 'status', eventVal: 1, arm: 'ccr5', exp: 'WM', ctl: 'WW' },
  { ds: 'ebmt1', label: 'EBMT transplant (overall survival, risk score)', time: 'srv', status: 'srvstat', arm: 'score', exp: 'High risk', ctl: 'Low risk' },
  { ds: 'bnct', label: 'BNCT brain tumour (OS, treated vs untreated; small)', time: 'time', status: 'death', arm: 'trt', exp: '3', ctl: '1' },
  { ds: 'cgd_fe', label: 'CGD RCT (time to first serious infection, rIFN-g vs placebo; recurrent→first-event)', time: 'time', status: 'status', arm: 'arm', exp: 'rIFN-g', ctl: 'placebo' },
  { ds: 'bladder_fe', label: 'Bladder cancer (time to first recurrence, thiotepa vs placebo; recurrent→first-event)', time: 'time', status: 'status', arm: 'arm', exp: '2', ctl: '1' },
  // --- open packaged datasets added 2026-06-10 (KMsurv / asaur via Rdatasets mirror) ---
  { ds: 'alloauto', label: 'Leukemia transplant (DFS, allogeneic vs autologous)', time: 'time', status: 'delta', arm: 'type', exp: '1', ctl: '2' },
  { ds: 'larynx', label: 'Larynx cancer (death, stage III vs I)', time: 'time', status: 'delta', arm: 'stage', exp: '3', ctl: '1' },
  { ds: 'burn', label: 'Burn RCT (time to staph infection, body cleansing vs routine bathing)', time: 'T3', status: 'D3', arm: 'Z1', exp: '1', ctl: '0' },
  { ds: 'pneumon', label: 'Infant pneumonia (time to hospitalization, mother smoking)', time: 'chldage', status: 'hospital', arm: 'smoke', exp: '1', ctl: '0' },
  { ds: 'bfeed', label: 'Breastfeeding (time to weaning, mother smoking)', time: 'duration', status: 'delta', arm: 'smoke', exp: '1', ctl: '0' },
  { ds: 'hepatoCellular', label: 'Hepatocellular carcinoma (OS, vascular invasion)', time: 'OS', status: 'Death', arm: 'Vascularinvasion', exp: '1', ctl: '0' },
  // --- real cancer-survival IPD from cBioPortal/TCGA, STRONG contrast = late vs early stage ---
  { ds: 'cbio_luad', label: 'TCGA lung adeno (OS, late vs early stage)', time: 'time', status: 'status', arm: 'stage_group', exp: 'late', ctl: 'early' },
  { ds: 'cbio_coadread', label: 'TCGA colorectal (OS, late vs early stage)', time: 'time', status: 'status', arm: 'stage_group', exp: 'late', ctl: 'early' },
  { ds: 'cbio_stad', label: 'TCGA stomach (OS, late vs early stage)', time: 'time', status: 'status', arm: 'stage_group', exp: 'late', ctl: 'early' },
  { ds: 'cbio_kirc', label: 'TCGA kidney clear-cell (OS, late vs early stage)', time: 'time', status: 'status', arm: 'stage_group', exp: 'late', ctl: 'early' },
  { ds: 'cbio_hnsc', label: 'TCGA head & neck (OS, late vs early stage)', time: 'time', status: 'status', arm: 'stage_group', exp: 'late', ctl: 'early' },
  { ds: 'cbio_skcm', label: 'TCGA melanoma (OS, late vs early stage)', time: 'time', status: 'status', arm: 'stage_group', exp: 'late', ctl: 'early' },
  { ds: 'cbio_blca', label: 'TCGA bladder (OS, late vs early stage)', time: 'time', status: 'status', arm: 'stage_group', exp: 'late', ctl: 'early' },
];
const CAP = 2500; // subsample cap per arm (file order; keeps huge cohorts tractable)

function coxHR(expIpd, ctlIpd) {
  return _.coxLogHR(expIpd.map(r => ({ time: r.time, status: r.status, x: 1 }))
    .concat(ctlIpd.map(r => ({ time: r.time, status: r.status, x: 0 }))));
}

// Load the two true patient-level arms for a dataset config (shared by run() and the head-to-head).
function loadArms(cfg) {
  const rows = parseCSV(path.join(dir, cfg.ds + '.csv'));
  const toEvent = (s) => cfg.eventVal != null ? (s === cfg.eventVal ? 1 : 0) : (s >= 1 ? 1 : 0);
  const arm = (which) => rows.filter(r => String(r[cfg.arm]) === which)
    .map(r => ({ time: num(r[cfg.time]), status: num(r[cfg.status]) }))
    .filter(r => r.time != null && r.status != null && r.time > 0)
    .map(r => ({ time: r.time, status: toEvent(r.status) }));
  const subsample = (a) => { if (a.length <= CAP) return a; const step = a.length / CAP, s = []; for (let i = 0; i < CAP; i++) s.push(a[Math.floor(i * step)]); return s; };
  return { expT: subsample(arm(cfg.exp)), ctlT: subsample(arm(cfg.ctl)) };
}

function run(cfg, K) {
  K = K || 8;  // number of registry-style posted KM timepoints
  const { expT, ctlT } = loadArms(cfg);
  if (expT.length < 20 || ctlT.length < 20) return { ds: cfg.ds, error: 'too few rows' };

  // ---- TRUE estimates from full IPD ----
  const kmE = _.kmFromIPD(expT), kmC = _.kmFromIPD(ctlT);
  const tau = 0.9 * Math.min(Math.max(...expT.map(r => r.time)), Math.max(...ctlT.map(r => r.time)));
  const trueMedE = _.medianFromKM(kmE, { interpolate: true }), trueMedC = _.medianFromKM(kmC, { interpolate: true });
  const trueHR = coxHR(expT, ctlT).hr;
  const trueRMSTd = _.rmst(kmE, tau) - _.rmst(kmC, tau);

  // ---- registry-style coarse summary (what ct.gov would post): KM at K timepoints + N + events ----
  function coarse(km, ipd) {
    const tmax = 0.95 * Math.max(...ipd.map(r => r.time));
    const pts = [{ t: 0, S: 1 }];
    for (let i = 1; i <= K; i++) { const t = tmax * i / K; pts.push({ t: +t.toFixed(2), S: +_.evalKM(km, t).toFixed(4) }); }
    return { km_points: pts, nar_points: [], N: ipd.length, total_events: ipd.filter(r => r.status === 1).length, follow_up_max: +tmax.toFixed(2) };
  }
  const aE = coarse(kmE, expT), aC = coarse(kmC, ctlT);
  const trial = {
    nct_id: 'GOLD-' + cfg.ds, time_unit: 'days',
    arms: [
      Object.assign({ arm_id: 'exp', label: 'Experimental', role: 'experimental' }, aE),
      Object.assign({ arm_id: 'ctl', label: 'Comparator', role: 'comparator' }, aC),
    ],
    hr: { value: +trueHR.toFixed(4), favors_arm_id: trueHR < 1 ? 'exp' : 'ctl', method: 'Cox(true)' },
  };

  function evalRecon(opts, useEvents) {
    // curve-only strips total_events; informed keeps them
    const t2 = JSON.parse(JSON.stringify(trial));
    if (!useEvents) t2.arms.forEach(a => { a.total_events = null; });
    const r = RIPD.reconstruct(t2, opts);
    if (!r.arms) return null;
    const e = r.arms.find(a => a.role === 'experimental'), c = r.arms.find(a => a.role === 'comparator');
    const me = _.medianFromKM(_.kmFromIPD(e.ipd), { interpolate: true });
    const mc = _.medianFromKM(_.kmFromIPD(c.ipd), { interpolate: true });
    const hr = coxHR(e.ipd, c.ipd).hr;
    const rmstd = _.rmst(_.kmFromIPD(e.ipd), tau) - _.rmst(_.kmFromIPD(c.ipd), tau);
    return { hr, medE: me, medC: mc, rmstd };
  }
  const co = evalRecon({}, false);       // curve-only
  const inf = evalRecon({}, true);       // censoring-informed (uses registry total events)
  const pe = (a, b) => (a == null || !b) ? null : Math.abs(a - b) / Math.abs(b);
  return {
    ds: cfg.ds, label: cfg.label, n_exp: expT.length, n_ctl: ctlT.length,
    true_HR: +trueHR.toFixed(3), true_median_exp: round1(trueMedE), true_median_ctl: round1(trueMedC), true_RMST_diff: round1(trueRMSTd),
    curve_only: { HR: +co.hr.toFixed(3), HR_logerr: +Math.abs(Math.log(co.hr) - Math.log(trueHR)).toFixed(3), median_exp_pcterr: pct(pe(co.medE, trueMedE)), RMSTdiff: round1(co.rmstd) },
    censoring_informed: { HR: +inf.hr.toFixed(3), HR_logerr: +Math.abs(Math.log(inf.hr) - Math.log(trueHR)).toFixed(3), median_exp_pcterr: pct(pe(inf.medE, trueMedE)), RMSTdiff: round1(inf.rmstd) },
  };
}
const round1 = x => x == null ? null : +x.toFixed(1);
const pct = x => x == null ? null : +(100 * x).toFixed(1);

if (require.main === module) {
  const results = CONFIGS.map(c => run(c, 8));
  fs.writeFileSync(path.join(dir, 'goldstandard_results.json'), JSON.stringify(results, null, 2));
  console.log(JSON.stringify(results, null, 2));
} else {
  module.exports = { CONFIGS, run, dir, loadArms, coxHR };
}
