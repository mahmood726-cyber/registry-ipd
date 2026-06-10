#!/usr/bin/env node
/*
 * TIER-B VALIDATION — how good is the median+HR parametric reconstruction (no KM curve)?
 *
 * Tier B fires when a trial posts a median + HR + N + events but NO Kaplan-Meier curve. The engine
 * reconstructs each arm as an exponential parametric model with the posted median (and imposes the HR).
 * This tier was never validated against truth (the gold standard tests Tier-A curve reconstruction).
 *
 * Because the HR and median are INPUTS, recovering them is circular — the honest test is the estimand
 * Tier B does NOT receive: **RMST** (restricted mean survival time), which depends on the survival
 * *shape* the exponential assumes. On the true-IPD gold standard we build the Tier-B inputs
 * (median + HR + N + events + follow-up, no curve), reconstruct, and compare the reconstructed
 * RMST-difference to the TRUE RMST-difference. Large errors flag where the exponential (constant-hazard)
 * assumption breaks (early-heavy or plateauing hazards).
 *
 * Deterministic. Run: node validate/tierb_validation.js  ->  realipd/tierb_validation_results.json
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');
const _ = RIPD._;
const GS = require('./goldstandard.js');

const coxFit = (a, b) => _.coxLogHR(a.map(r => ({ time: r.time, status: r.status, x: 1 }))
  .concat(b.map(r => ({ time: r.time, status: r.status, x: 0 }))).slice());

// numerical RMST = ∫_0^tau S(t) dt for a given survival function
function rmstOf(S, tau, steps) { steps = steps || 400; const dt = tau / steps; let a = 0; for (let i = 0; i < steps; i++) { const t = (i + 0.5) * dt; a += S(t) * dt; } return a; }
// closed-form 2-parameter Weibull from a median and S(tmax)≈1-events/N:
//   S(t)=exp(-(t/scale)^shape); S(median)=0.5 and S(tmax)=Stmax ⇒ shape=ln(R/L)/ln(tmax/median),
//   R=-ln(Stmax), L=ln2, scale=median/L^(1/shape). shape≈1 ⇒ exponential.
function weibullFromMedianRate(median, events, N, tmax) {
  const Stmax = Math.min(0.999, Math.max(1e-3, 1 - events / N));
  const A = tmax / median, R = -Math.log(Stmax), L = Math.LN2;
  if (!(A > 1.001) || !(R > 0)) return { shape: 1, scale: median / Math.pow(L, 1) }; // fall back to exponential
  let shape = Math.log(R / L) / Math.log(A); if (!isFinite(shape) || shape <= 0.05) shape = 1; shape = Math.min(8, shape);
  return { shape, scale: median / Math.pow(L, 1 / shape) };
}

const rows = [];
for (const cfg of GS.CONFIGS) {
  let arms; try { arms = GS.loadArms(cfg); } catch { continue; }
  const { expT, ctlT } = arms; if (expT.length < 20 || ctlT.length < 20) continue;
  const kmE = _.kmFromIPD(expT), kmC = _.kmFromIPD(ctlT);
  const trueMedE = _.medianFromKM(kmE, { interpolate: true }), trueMedC = _.medianFromKM(kmC, { interpolate: true });
  if (trueMedE == null || trueMedC == null) continue;   // median undefined (curve never reaches 0.5) → Tier B n/a
  const fit = coxFit(expT, ctlT), trueHR = Math.exp(fit.beta), hrSE = 1 / Math.sqrt(Math.max(1, expT.filter(r => r.status === 1).length + ctlT.filter(r => r.status === 1).length) / 4);
  const tau = 0.9 * Math.min(Math.max(...expT.map(r => r.time)), Math.max(...ctlT.map(r => r.time)));
  const trueRmstE = _.rmst(kmE, tau), trueRmstC = _.rmst(kmC, tau), trueRmstD = trueRmstE - trueRmstC;

  // Tier-B inputs: median + HR + N + events + follow-up, NO km_points
  const arm = (ipd, med, role) => ({ arm_id: role.slice(0, 3), role, N: ipd.length, total_events: ipd.filter(r => r.status === 1).length,
    follow_up_max: +(Math.max(...ipd.map(r => r.time))).toFixed(2), median: { value: +med.toFixed(2), ci_low: +(med * 0.8).toFixed(2), ci_high: +(med * 1.25).toFixed(2) } });
  const trial = { nct_id: 'TB-' + cfg.ds, time_unit: 'd', arms: [arm(expT, trueMedE, 'experimental'), arm(ctlT, trueMedC, 'comparator')],
    hr: { value: +trueHR.toFixed(3), ci_low: +Math.exp(fit.beta - 1.96 * hrSE).toFixed(3), ci_high: +Math.exp(fit.beta + 1.96 * hrSE).toFixed(3), favors_arm_id: trueHR < 1 ? 'exp' : 'com' } };

  const r = RIPD.reconstruct(trial, { bootstrap: 0 });
  if (r.tier !== 'B' || !r.arms) { rows.push({ ds: cfg.ds, error: 'not Tier B (' + r.tier + ')' }); continue; }
  const e = r.arms.find(a => a.role === 'experimental'), c = r.arms.find(a => a.role === 'comparator');
  const recRmstE = _.rmst(_.kmFromIPD(e.ipd), tau), recRmstC = _.rmst(_.kmFromIPD(c.ipd), tau), recRmstD = recRmstE - recRmstC;
  const recMedE = _.medianFromKM(_.kmFromIPD(e.ipd), { interpolate: true });
  const pe = (a, b) => (b && isFinite(a)) ? +(100 * Math.abs(a - b) / Math.abs(b)).toFixed(1) : null;

  // per-arm shape-model comparison (exponential vs Weibull) on RMST, from median + events/N + tmax.
  // The clean shape test: each arm uses its OWN median (no imposed-HR coupling).
  const armShape = (ipd, med) => {
    const ev = ipd.filter(r => r.status === 1).length, N = ipd.length, tmax = Math.max(...ipd.map(r => r.time));
    const lam = Math.LN2 / med; const Sexp = (t) => Math.exp(-lam * t);
    const w = weibullFromMedianRate(med, ev, N, tmax); const Sw = (t) => Math.exp(-Math.pow(t / w.scale, w.shape));
    const trueR = _.rmst(_.kmFromIPD(ipd), tau);
    return { exp: rmstOf(Sexp, tau), weib: rmstOf(Sw, tau), true: trueR, shape: +w.shape.toFixed(2) };
  };
  const aE = armShape(expT, trueMedE), aC = armShape(ctlT, trueMedC);
  const expRmstDiff_err = Math.abs((aE.exp - aC.exp) - trueRmstD), weibRmstDiff_err = Math.abs((aE.weib - aC.weib) - trueRmstD);
  rows.push({ ds: cfg.ds, tcga: cfg.ds.startsWith('cbio'), n_exp: expT.length, n_ctl: ctlT.length,
    true_HR: +trueHR.toFixed(3), true_median_exp: +trueMedE.toFixed(1), true_RMSTdiff: +trueRmstD.toFixed(1),
    recon_RMSTdiff: +recRmstD.toFixed(1), RMSTdiff_abs_err: +Math.abs(recRmstD - trueRmstD).toFixed(1),
    RMST_exp_pcterr: pe(recRmstE, trueRmstE), median_exp_pcterr: pe(recMedE, trueMedE),
    weibull_shape_exp: aE.shape, weibull_shape_ctl: aC.shape,
    exp_RMSTdiff_err: +expRmstDiff_err.toFixed(1), weibull_RMSTdiff_err: +weibRmstDiff_err.toFixed(1),
    exp_arm_RMST_pcterr: pe(aE.exp, aE.true), weibull_arm_RMST_pcterr: pe(aE.weib, aE.true) });
}

const ok = rows.filter(r => !r.error);
const med = (get) => { const xs = ok.map(get).filter(x => x != null).sort((a, b) => a - b); return xs.length ? +xs[xs.length >> 1].toFixed(2) : null; };
const within = (get, thr) => ok.filter(r => get(r) != null && get(r) < thr).length + '/' + ok.length;
const summary = {
  n_datasets: ok.length, n_errored: rows.length - ok.length,
  engine_exponential: { median_RMST_exp_arm_pcterr: med(r => r.RMST_exp_pcterr), within20pct: within(r => r.RMST_exp_pcterr, 20), median_median_exp_pcterr: med(r => r.median_exp_pcterr) },
  shape_model_RMSTdiff_abs_err: { exponential: med(r => r.exp_RMSTdiff_err), weibull: med(r => r.weibull_RMSTdiff_err) },
  shape_model_arm_RMST_pcterr: { exponential: med(r => r.exp_arm_RMST_pcterr), weibull: med(r => r.weibull_arm_RMST_pcterr),
    exp_within20: within(r => r.exp_arm_RMST_pcterr, 20), weibull_within20: within(r => r.weibull_arm_RMST_pcterr, 20) },
  non_exponential_examples: ok.slice().sort((a, b) => (b.exp_arm_RMST_pcterr || 0) - (a.exp_arm_RMST_pcterr || 0)).slice(0, 6)
    .map(r => ({ ds: r.ds, shape: r.weibull_shape_exp, exp_RMST_err: r.exp_arm_RMST_pcterr, weibull_RMST_err: r.weibull_arm_RMST_pcterr })),
  note: 'Tier B = median+HR parametric, NO KM curve. HR/median are inputs (circular); RMST tests the '
    + 'survival-shape assumption. Compared the engine exponential vs a closed-form 2-param Weibull fit '
    + 'from median + events/N (shape≈1 ⇒ exponential; <1 early-heavy, >1 early-light).',
};
fs.writeFileSync(path.join(GS.dir, 'tierb_validation_results.json'), JSON.stringify({ summary, per_dataset: rows }, null, 2));
console.log(JSON.stringify(summary, null, 2));
