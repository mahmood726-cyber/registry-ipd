#!/usr/bin/env node
/*
 * GOLD-STANDARD competing-risks validation on TRUE patient-level data (R survival::colon).
 *
 * colon has two rows/patient: etype 1 = time to recurrence, etype 2 = time to death. We build the
 * true cause-labeled IPD for the RECURRENCE endpoint with death as a COMPETING event
 * (cause 1 = recurrence, cause 2 = death-without-recurrence, 0 = censored), per treatment arm.
 *
 * Then: (1) TRUE Aalen–Johansen CIF for recurrence + the (biased) naive 1-KM; (2) the registry-style
 * coarse summary a sponsor would post (cause-specific recurrence curve @8 timepoints + N + cause-1
 * events + competing-death count); (3) reconstruct cause-labeled pseudo-IPD; (4) compare the
 * reconstructed AJ CIF to the TRUE AJ CIF, and confirm the naive 1-KM overestimates on real data.
 *
 * Usage: node validate/goldstandard_cr.js [realipd_dir]
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');
const _ = RIPD._;
const dir = process.argv[2] && !process.argv[2].startsWith('-') ? process.argv[2] : path.join(__dirname, '..', 'realipd');

const txt = fs.readFileSync(path.join(dir, 'colon.csv'), 'utf8').trim().split(/\r?\n/);
const head = txt[0].split(',').map(h => h.replace(/"/g, ''));
const rows = txt.slice(1).map(line => { const c = line.split(','); const o = {}; head.forEach((h, i) => o[h] = (c[i] || '').replace(/"/g, '')); return o; });
const num = x => { const v = parseFloat(x); return isFinite(v) ? v : null; };

// build per-patient cause-labeled record for recurrence (death competing), per arm
function crByArm(rxVal) {
  const byId = {};
  for (const r of rows) {
    if (r.rx !== rxVal) continue;
    (byId[r.id] = byId[r.id] || {})[num(r.etype)] = { time: num(r.time), status: num(r.status) };
  }
  const ipd = [];
  for (const id in byId) {
    const rec = byId[id][1], dth = byId[id][2];
    if (!rec || !dth) continue;
    if (rec.status === 1) ipd.push({ time: rec.time, cause: 1 });            // recurrence first
    else if (dth.status === 1) ipd.push({ time: dth.time, cause: 2 });       // died w/o recurrence (competing)
    else ipd.push({ time: rec.time, cause: 0 });                             // censored
  }
  return ipd;
}

function cifAt(ajSteps, t) { let v = 0; for (const s of ajSteps) { if (s.t <= t + 1e-9) v = s.cif1; else break; } return v; }
function naiveAt(naiveSteps, t) { let v = 0; for (const s of naiveSteps) { if (s.t <= t + 1e-9) v = s.cif1; else break; } return v; }

function run(rxVal, label) {
  const cr = crByArm(rxVal);
  if (cr.length < 50) return { arm: label, error: 'too few' };
  const trueAJ = _.cifAalenJohansen(cr);
  const trueNaive = _.cifNaive1(cr);
  const tmax = 0.9 * Math.max(...cr.map(r => r.time));
  const ts = []; for (let i = 1; i <= 8; i++) ts.push(+(tmax * i / 8).toFixed(0));

  // registry-style coarse summary: cause-specific recurrence "survival" (event=cause1, else censored)
  const causeSpecificIPD = cr.map(r => ({ time: r.time, status: r.cause === 1 ? 1 : 0 }));
  const km = _.kmFromIPD(causeSpecificIPD);
  const km_points = [{ t: 0, S: 1 }].concat(ts.map(t => ({ t, S: +_.evalKM(km, t).toFixed(4) })));
  const N = cr.length;
  const cause1 = cr.filter(r => r.cause === 1).length, cause2 = cr.filter(r => r.cause === 2).length;
  const trial = {
    nct_id: 'GOLDCR-' + rxVal, time_unit: 'days',
    arms: [{ arm_id: 'a', label, role: 'experimental', N, total_events: cause1, competing_events: cause2,
      follow_up_max: +tmax.toFixed(0), km_points, nar_points: [] }],
    hr: null,
  };
  const r = RIPD.reconstructCompetingRisks(trial);
  const reconAJ = r.arms[0].cif;

  // compare reconstructed AJ CIF vs TRUE AJ CIF at the coarse timepoints
  let maxAbs = 0; const tbl = [];
  for (const t of ts) {
    const tv = cifAt(trueAJ, t), rv = cifAt(reconAJ, t), nv = naiveAt(trueNaive, t);
    maxAbs = Math.max(maxAbs, Math.abs(rv - tv));
    tbl.push({ t, true_AJ: +tv.toFixed(3), recon_AJ: +rv.toFixed(3), true_naive_1KM: +nv.toFixed(3) });
  }
  const finalTrueAJ = cifAt(trueAJ, tmax), finalNaive = naiveAt(trueNaive, tmax);
  return {
    arm: label, n: N, recurrences: cause1, competing_deaths: cause2,
    recon_AJ_CIF_max_abs_err_vs_true: +maxAbs.toFixed(3),
    true_final_AJ_CIF: +finalTrueAJ.toFixed(3),
    true_final_naive_1KM_CIF: +finalNaive.toFixed(3),
    naive_overestimates_by_pp: +(100 * (finalNaive - finalTrueAJ)).toFixed(1),
    curve: tbl,
  };
}

const out = [run('Obs', 'Observation'), run('Lev+5FU', 'Levamisole+5FU')];
fs.writeFileSync(path.join(dir, 'goldstandard_cr_results.json'), JSON.stringify(out, null, 2));
out.forEach(o => {
  if (o.error) return console.log(o.arm, o.error);
  console.log(`\n${o.arm} (n=${o.n}, ${o.recurrences} recurrences, ${o.competing_deaths} competing deaths)`);
  console.log(`  reconstructed AJ CIF vs TRUE AJ CIF: max abs err = ${o.recon_AJ_CIF_max_abs_err_vs_true}`);
  console.log(`  TRUE final recurrence CIF: Aalen-Johansen ${o.true_final_AJ_CIF} vs naive 1-KM ${o.true_final_naive_1KM_CIF} (naive overestimates by ${o.naive_overestimates_by_pp} pp)`);
});
