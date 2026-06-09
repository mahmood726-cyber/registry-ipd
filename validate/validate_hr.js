#!/usr/bin/env node
/*
 * Empirical "is it good enough?" validation, using the registry-reported hazard ratio as
 * HELD-OUT GROUND TRUTH. For every Tier-A trial that has a curve AND a reported HR:
 *   reconstruct pseudo-IPD from the KM curve alone -> Cox HR on the pseudo-IPD -> compare.
 *
 * HONESTY on direction: ct.gov HR direction (which arm is reference) is frequently ambiguous,
 * so we score MAGNITUDE error orientation-robustly:  min(|lnHRr - lnHRreg|, |lnHRr + lnHRreg|),
 * and report directional agreement separately. We stratify by reconstructed event count, because
 * Cox HR is inherently unstable when events are few — that stratification IS the answer to
 * "good enough for what".
 *
 * Usage: node validate/validate_hr.js [cohort_dir] [-o report.json]
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');

const dir = process.argv[2] && !process.argv[2].startsWith('-') ? process.argv[2]
  : path.join(__dirname, '..', 'cohort');
const outIdx = process.argv.indexOf('-o');
const outPath = outIdx >= 0 ? process.argv[outIdx + 1] : path.join(dir, 'validation_hr_report.json');

function median(a) { if (!a.length) return null; const s = a.slice().sort((x, y) => x - y); const m = s.length >> 1; return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; }
function quant(a, p) { if (!a.length) return null; const s = a.slice().sort((x, y) => x - y); const i = p * (s.length - 1); const lo = Math.floor(i), hi = Math.ceil(i); return lo === hi ? s[lo] : s[lo] + (s[hi] - s[lo]) * (i - lo); }

const files = fs.readdirSync(dir).filter(f => f.endsWith('.json') && !f.startsWith('manifest') && !f.startsWith('validation'));
const rows = [];
let skippedNoHr = 0, skippedArms = 0, skippedRecon = 0;

for (const f of files) {
  let t; try { t = JSON.parse(fs.readFileSync(path.join(dir, f), 'utf8')); } catch { continue; }
  if (!t.hr || t.hr.value == null) { skippedNoHr++; continue; }
  if (!t.arms || t.arms.length !== 2) { skippedArms++; continue; }      // clean 2-arm comparison only
  let res; try { res = RIPD.reconstruct(t); } catch { skippedRecon++; continue; }
  if (res.tier !== 'A' || !res.arms || res.arms.length !== 2) { skippedRecon++; continue; }

  // Orient by ROLE (experimental x=1 vs comparator x=0), the convention registry HR uses, so
  // direction is meaningful. If there is no clear comparator, direction is undeterminable.
  let exp = res.arms.find(a => a.role === 'experimental');
  let ctl = res.arms.find(a => a.role === 'comparator');
  const dirDeterminable = !!(exp && ctl);
  if (!dirDeterminable) { exp = res.arms[1]; ctl = res.arms[0]; }
  const ev = exp.ipd.filter(r => r.status === 1).length + ctl.ipd.filter(r => r.status === 1).length;
  const cox = RIPD._.coxLogHR(exp.ipd.map(r => ({ time: r.time, status: r.status, x: 1 }))
    .concat(ctl.ipd.map(r => ({ time: r.time, status: r.status, x: 0 }))));
  const lnRecon = cox.beta, lnReg = Math.log(t.hr.value);
  // magnitude error orientation-robust; direction scored only when determinable
  const eSame = Math.abs(lnRecon - lnReg), eInv = Math.abs(lnRecon + lnReg);
  const magErr = Math.min(eSame, eInv);
  const directionAgrees = dirDeterminable ? (eSame <= eInv) : null;
  // within registry CI (try both orientations)
  let withinCI = null;
  if (t.hr.ci_low != null && t.hr.ci_high != null) {
    const hrR = Math.exp(lnRecon), inv = 1 / hrR;
    withinCI = (hrR >= t.hr.ci_low - 1e-9 && hrR <= t.hr.ci_high + 1e-9) ||
      (inv >= t.hr.ci_low - 1e-9 && inv <= t.hr.ci_high + 1e-9);
  }
  rows.push({
    nct: t.nct_id, events: ev, recon_hr: +Math.exp(lnRecon).toFixed(4), reg_hr: t.hr.value,
    logHR_mag_err: +magErr.toFixed(4), direction_agrees: directionAgrees, within_ci: withinCI,
    separated: cox.separated, method: res.method, badge: res.audit.badge,
    n_timepoints: Math.max(...t.arms.map(a => a.km_points.length)),
  });
}

function summarize(subset, label) {
  if (!subset.length) return { stratum: label, n: 0 };
  const errs = subset.map(r => r.logHR_mag_err);
  const ci = subset.filter(r => r.within_ci != null);
  return {
    stratum: label, n: subset.length,
    logHR_mag_err_median: +median(errs).toFixed(4),
    logHR_mag_err_p90: +quant(errs, 0.9).toFixed(4),
    HR_fold_err_median: +Math.exp(median(errs)).toFixed(3),  // multiplicative: e.g. 1.15x
    within_registry_CI_pct: ci.length ? Math.round(100 * ci.filter(r => r.within_ci).length / ci.length) : null,
    direction_determinable_n: subset.filter(r => r.direction_agrees !== null).length,
    direction_agree_pct: (() => { const d = subset.filter(r => r.direction_agrees !== null); return d.length ? Math.round(100 * d.filter(r => r.direction_agrees).length / d.length) : null; })(),
  };
}

const lo = rows.filter(r => r.events < 10), mid = rows.filter(r => r.events >= 10 && r.events < 50), hi = rows.filter(r => r.events >= 50);
const report = {
  cohort_dir: dir, n_trials_files: files.length,
  usable_2arm_with_HR: rows.length,
  skipped: { no_registry_hr: skippedNoHr, not_2_arms: skippedArms, reconstruction: skippedRecon },
  overall: summarize(rows, 'all'),
  by_events: [summarize(lo, 'events<10 (Cox unstable)'), summarize(mid, 'events 10-49'), summarize(hi, 'events>=50')],
  interpretation: 'logHR_mag_err = |ln(recon HR) - ln(registry HR)| orientation-robust. '
    + 'HR_fold_err_median ~1.0 = excellent; ~1.1 = within 10%. Cox HR is unstable when events<10, '
    + 'so the events>=50 stratum is the fair test of reconstruction-from-curve quality.',
  per_trial: rows.sort((a, b) => b.events - a.events),
};
fs.writeFileSync(outPath, JSON.stringify(report, null, 2));
console.log(JSON.stringify({ usable: report.usable_2arm_with_HR, overall: report.overall, by_events: report.by_events }, null, 2));
console.log('wrote', outPath);
