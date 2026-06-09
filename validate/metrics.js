/*
 * Validation metrics for the head-to-head: AACT-only reconstruction vs digitization.
 *
 * The HONEST comparison the plan commits to: on trials with rich registry survival data,
 * AACT-only reconstruction should drive ANCHOR error to ~0 (exact registry values) while a
 * digitized curve carries positive pixel error. We quantify fidelity of any reconstruction to
 * the registry-reported HR / median / RMST / KM anchors, so the same yardstick applies to both
 * the AACT-only result and a digitization result.
 *
 * All metrics reuse the audited engine primitives (kmFromIPD, evalKM, rmst, coxLogHR).
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory(require('../src/engine.js'));
  else root.RIPD_METRICS = factory(root.RIPD);
}(typeof self !== 'undefined' ? self : this, function (RIPD) {
  'use strict';
  const _ = RIPD._;

  // Build a right-continuous survival step function from sparse registry KM anchors.
  function anchorStep(km_points) {
    const pts = km_points.slice().sort((a, b) => a.t - b.t);
    return function (t) { let S = 1; for (const p of pts) { if (p.t <= t + 1e-9) S = p.S; else break; } return S; };
  }

  // Cox HR from a reconstruction result (experimental vs comparator).
  function reconHR(result) {
    if (!result.arms || result.expIdx == null || result.ctlIdx == null) {
      // fall back: find by role
      const e = (result.arms || []).find(a => a.role === 'experimental');
      const c = (result.arms || []).find(a => a.role === 'comparator');
      if (!e || !c) return null;
      return _.coxLogHR(e.ipd.map(r => ({ time: r.time, status: r.status, x: 1 }))
        .concat(c.ipd.map(r => ({ time: r.time, status: r.status, x: 0 })))).hr;
    }
    const e = result.arms[result.expIdx], c = result.arms[result.ctlIdx];
    return _.coxLogHR(e.ipd.map(r => ({ time: r.time, status: r.status, x: 1 }))
      .concat(c.ipd.map(r => ({ time: r.time, status: r.status, x: 0 })))).hr;
  }

  // max_k | S_recon(t_k) - S_registry(t_k) | over every arm's registry anchors.
  function anchorSupError(trial, result) {
    if (!result.arms) return null;
    let mx = 0, any = false;
    for (const a of trial.arms) {
      const rec = result.arms.find(x => x.arm_id === a.arm_id);
      if (!rec || !(a.km_points && a.km_points.length)) continue;
      const km = _.kmFromIPD(rec.ipd);
      for (const p of a.km_points) { mx = Math.max(mx, Math.abs(_.evalKM(km, p.t) - p.S)); any = true; }
    }
    return any ? mx : null;
  }

  // 1-Wasserstein (L1 area) between reconstructed KM and the registry anchor step, per arm, averaged.
  function wassersteinToAnchors(trial, result, tau) {
    if (!result.arms) return null;
    const vals = [];
    for (const a of trial.arms) {
      const rec = result.arms.find(x => x.arm_id === a.arm_id);
      if (!rec || !(a.km_points && a.km_points.length)) continue;
      const km = _.kmFromIPD(rec.ipd);
      const aStep = anchorStep(a.km_points);
      const grid = 400, dt = tau / grid; let acc = 0;
      for (let i = 0; i < grid; i++) { const t = (i + 0.5) * dt; acc += Math.abs(_.evalKM(km, t) - aStep(t)) * dt; }
      vals.push(acc);
    }
    return vals.length ? vals.reduce((x, y) => x + y, 0) / vals.length : null;
  }

  function reconMedian(result, role) {
    if (!result.arms) return null;
    const a = result.arms.find(x => x.role === role);
    return a ? _.medianFromKM(_.kmFromIPD(a.ipd)) : null;
  }

  function maxAnchorTime(trial) {
    let t = 0; for (const a of trial.arms) for (const p of (a.km_points || [])) t = Math.max(t, p.t);
    return t || 1;
  }

  // Fidelity of one reconstruction to the registry-reported summaries + anchors.
  function fidelity(trial, result) {
    const tau = maxAnchorTime(trial);
    const out = { tier: result.tier, badge: result.audit && result.audit.badge };
    // HR
    if (trial.hr && trial.hr.value != null) {
      const rh = reconHR(result);
      out.recon_hr = rh;
      out.logHR_err = (rh != null) ? Math.abs(Math.log(rh) - Math.log(trial.hr.value)) : null;
    }
    // median per registry-reported arm
    out.median_abs_err = {};
    for (const a of trial.arms) {
      if (a.median && a.median.value != null) {
        const rm = reconMedian(result, a.role);
        out.median_abs_err[a.arm_id] = (rm != null) ? Math.abs(rm - a.median.value) : null;
      }
    }
    out.anchor_sup_error = anchorSupError(trial, result);
    out.wasserstein_to_anchors = wassersteinToAnchors(trial, result, tau);
    return out;
  }

  // Head-to-head: same trial, AACT-only vs a digitization result (engine run on digitized points).
  // `digitizedTrial` is the same trial with km_points replaced by digitized-from-figure points.
  function headToHead(trial, opts) {
    opts = opts || {};
    const aact = RIPD.reconstruct(trial, opts);
    const res = { nct_id: trial.nct_id, aact_only: fidelity(trial, aact) };
    if (opts.digitizedTrial) {
      const digi = RIPD.reconstruct(opts.digitizedTrial, opts);
      // compare digitized reconstruction against the SAME registry anchors (trial.arms[].km_points)
      res.digitization = fidelity(trial, mergeArms(opts.digitizedTrial, digi));
    }
    return res;
  }
  // helper: attach original registry km_points so anchor metrics use registry truth, not digitized pts
  function mergeArms(digTrial, digResult) { return digResult; }

  // Aggregate a cohort of per-trial fidelity objects.
  function aggregate(fidelities) {
    const logHR = fidelities.map(f => f.logHR_err).filter(x => x != null);
    const sup = fidelities.map(f => f.anchor_sup_error).filter(x => x != null);
    const w = fidelities.map(f => f.wasserstein_to_anchors).filter(x => x != null);
    const rmse = (arr) => arr.length ? Math.sqrt(arr.reduce((a, b) => a + b * b, 0) / arr.length) : null;
    const mean = (arr) => arr.length ? arr.reduce((a, b) => a + b, 0) / arr.length : null;
    return {
      n: fidelities.length,
      logHR_RMSE: rmse(logHR), logHR_mean_abs: mean(logHR),
      anchor_sup_error_mean: mean(sup), anchor_sup_error_max: sup.length ? Math.max(...sup) : null,
      wasserstein_mean: mean(w),
      tier_counts: fidelities.reduce((m, f) => (m[f.tier] = (m[f.tier] || 0) + 1, m), {})
    };
  }

  return { fidelity, headToHead, aggregate, reconHR, anchorSupError, wassersteinToAnchors, reconMedian, anchorStep };
}));
