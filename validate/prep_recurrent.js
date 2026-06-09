#!/usr/bin/env node
/*
 * Collapse RECURRENT-EVENT datasets to time-to-FIRST-event (the standard way to bring recurrent
 * survival data into a single-event reconstruction): writes realipd/<ds>_fe.csv with time,status,arm.
 *   cgd      -> time to first serious infection (rIFN-g vs placebo RCT)
 *   bladder  -> time to first tumour recurrence (thiotepa vs placebo)
 * Run: node validate/prep_recurrent.js
 */
const fs = require('fs');
const path = require('path');
const SRC = path.join(__dirname, '..', 'realipd', 'more');
const OUT = path.join(__dirname, '..', 'realipd');
function parse(file) {
  const t = fs.readFileSync(file, 'utf8').trim().split(/\r?\n/);
  const h = t[0].split(',').map(x => x.replace(/"/g, ''));
  return t.slice(1).map(l => { const c = l.split(','); const o = {}; h.forEach((k, i) => o[k] = (c[i] || '').replace(/"/g, '')); return o; });
}
const N = x => { const v = parseFloat(x); return isFinite(v) ? v : null; };

// cgd: per id, first row (by enum) with status==1 -> (tstop,1); else last tstop censored
(function () {
  const rows = parse(path.join(SRC, 'cgd.csv'));
  const byId = {};
  for (const r of rows) (byId[r.id] = byId[r.id] || []).push(r);
  const out = ['time,status,arm'];
  for (const id in byId) {
    const recs = byId[id].sort((a, b) => N(a.enum) - N(b.enum));
    const firstEv = recs.find(r => N(r.status) === 1);
    const time = firstEv ? N(firstEv.tstop) : N(recs[recs.length - 1].tstop);
    const status = firstEv ? 1 : 0;
    out.push(`${time},${status},${recs[0].treat}`);
  }
  fs.writeFileSync(path.join(OUT, 'cgd_fe.csv'), out.join('\n'));
  console.log('cgd_fe.csv:', out.length - 1, 'patients');
})();

// bladder: enum==1 row per patient gives first-recurrence (stop, event); arm = rx
(function () {
  const rows = parse(path.join(SRC, 'bladder.csv')).filter(r => N(r.enum) === 1);
  const out = ['time,status,arm'];
  for (const r of rows) out.push(`${N(r.stop)},${N(r.event)},${r.rx}`);
  fs.writeFileSync(path.join(OUT, 'bladder_fe.csv'), out.join('\n'));
  console.log('bladder_fe.csv:', out.length - 1, 'patients');
})();
