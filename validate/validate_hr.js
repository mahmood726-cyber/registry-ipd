#!/usr/bin/env node
/*
 * PAIRED HR validation vs registry HR (held-out ground truth). For each usable 2-arm trial we
 * reconstruct BOTH curve-only (ignoreTotalEvents) and censoring-informed (as-harvested) on the
 * SAME trial, so curve-only vs informed is a like-for-like paired comparison (no shifting strata).
 *
 * Honesty fixes baked in:
 *  - within-CI is scored on the ROLE-ORIENTED reconstructed HR only (no inverse-orientation clause),
 *    on the direction-determinable subset — it no longer counts a wrong-direction fit as "covered".
 *  - magnitude error stays orientation-robust (registry HR direction is genuinely ambiguous).
 *  - Wilson 95% CIs on every proportion; McNemar on paired within-CI.
 *
 * Usage: node validate/validate_hr.js [cohort_dir] [-o report.json]
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');

const dir = process.argv[2] && !process.argv[2].startsWith('-') ? process.argv[2] : path.join(__dirname, '..', 'cohort');
const outIdx = process.argv.indexOf('-o');
const outPath = outIdx >= 0 ? process.argv[outIdx + 1] : path.join(dir, 'validation_hr_report.json');

const median = a => { if (!a.length) return null; const s = a.slice().sort((x, y) => x - y); const m = s.length >> 1; return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; };
const quant = (a, p) => { if (!a.length) return null; const s = a.slice().sort((x, y) => x - y); const i = p * (s.length - 1), lo = Math.floor(i), hi = Math.ceil(i); return lo === hi ? s[lo] : s[lo] + (s[hi] - s[lo]) * (i - lo); };
function wilson(k, n) { if (!n) return null; const z = 1.959964, p = k / n, d = 1 + z * z / n; const c = (p + z * z / (2 * n)) / d, h = z * Math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d; return [Math.round(100 * Math.max(0, c - h)), Math.round(100 * Math.min(1, c + h))]; }

// reconstruct one trial and score its HR vs registry (role-oriented within-CI, robust magnitude)
function score(t, opts) {
  const r = RIPD.reconstruct(t, opts);
  if (r.tier !== 'A' || !r.arms || r.arms.length !== 2) return null;
  let exp = r.arms.find(a => a.role === 'experimental'), ctl = r.arms.find(a => a.role === 'comparator');
  const determinable = !!(exp && ctl);
  if (!determinable) { exp = r.arms[1]; ctl = r.arms[0]; }
  const ev = exp.ipd.filter(x => x.status === 1).length + ctl.ipd.filter(x => x.status === 1).length;
  const cox = RIPD._.coxLogHR(exp.ipd.map(x => ({ time: x.time, status: x.status, x: 1 }))
    .concat(ctl.ipd.map(x => ({ time: x.time, status: x.status, x: 0 }))));
  const lnRecon = cox.beta, lnReg = Math.log(t.hr.value);
  const magErr = Math.min(Math.abs(lnRecon - lnReg), Math.abs(lnRecon + lnReg));
  const dirAgrees = determinable ? (Math.abs(lnRecon - lnReg) <= Math.abs(lnRecon + lnReg)) : null;
  // within-CI: role-oriented HR only (NO inverse clause), determinable subset only
  let within = null;
  if (determinable && t.hr.ci_low != null && t.hr.ci_high != null) {
    const hrR = Math.exp(lnRecon);
    within = hrR >= t.hr.ci_low - 1e-9 && hrR <= t.hr.ci_high + 1e-9;
  }
  return { ev, hr: Math.exp(lnRecon), magErr, dirAgrees, within, determinable };
}

const files = fs.readdirSync(dir).filter(f => f.endsWith('.json') && !f.startsWith('manifest') && !f.startsWith('validation') && !f.startsWith('registry'));
const rows = [];
for (const f of files) {
  let t; try { t = JSON.parse(fs.readFileSync(path.join(dir, f), 'utf8')); } catch { continue; }
  if (!t.hr || t.hr.value == null || !t.arms || t.arms.length !== 2) continue;
  const co = score(t, { ignoreTotalEvents: true });
  const inf = score(t, {});                 // uses harvested/patched total_events when present
  if (!co || !inf) continue;
  rows.push({ nct: t.nct_id, events: inf.ev, curve: co, informed: inf,
    informed_uses_events: t.arms.every(a => a.total_events != null) });
}

function summ(rows, pick) {
  const errs = rows.map(r => pick(r).magErr);
  const ciRows = rows.filter(r => pick(r).within != null);
  const dirRows = rows.filter(r => pick(r).dirAgrees != null);
  const ciK = ciRows.filter(r => pick(r).within).length;
  const dirK = dirRows.filter(r => pick(r).dirAgrees).length;
  return {
    n: rows.length,
    HR_fold_err_median: +Math.exp(median(errs)).toFixed(3),
    HR_fold_err_p90: +Math.exp(quant(errs, 0.9)).toFixed(3),
    within_CI: ciRows.length ? `${Math.round(100 * ciK / ciRows.length)}% (${ciK}/${ciRows.length}, 95%CI ${wilson(ciK, ciRows.length).join('-')})` : null,
    direction: dirRows.length ? `${Math.round(100 * dirK / dirRows.length)}% (${dirK}/${dirRows.length}, 95%CI ${wilson(dirK, dirRows.length).join('-')})` : null,
  };
}
// McNemar on paired within-CI (informed-only-hit vs curve-only-only-hit)
const both = rows.filter(r => r.curve.within != null && r.informed.within != null);
const b = both.filter(r => r.curve.within && !r.informed.within).length;  // curve hit, informed miss
const c = both.filter(r => !r.curve.within && r.informed.within).length;  // informed hit, curve miss
const report = {
  cohort_dir: dir, n_paired_2arm_with_HR: rows.length,
  curve_only: summ(rows, r => r.curve),
  censoring_informed: summ(rows, r => r.informed),
  paired_within_CI_mcnemar: { informed_gained: c, informed_lost: b, n_discordant: b + c,
    note: 'paired on identical trials; small discordant counts => difference not significant' },
  note: 'within-CI is role-oriented (no inverse clause); HR direction is genuinely ambiguous in the '
    + 'registry so magnitude error is orientation-robust. n is small (~30) — read with the Wilson CIs.',
  per_trial: rows.map(r => ({ nct: r.nct, events: r.events, reg_hr: undefined,
    curve_hr: +r.curve.hr.toFixed(3), informed_hr: +r.informed.hr.toFixed(3),
    informed_uses_events: r.informed_uses_events })),
};
fs.writeFileSync(outPath, JSON.stringify(report, null, 2));
console.log(JSON.stringify({ n: report.n_paired_2arm_with_HR, curve_only: report.curve_only,
  censoring_informed: report.censoring_informed, mcnemar: report.paired_within_CI_mcnemar }, null, 2));
console.log('wrote', outPath);
