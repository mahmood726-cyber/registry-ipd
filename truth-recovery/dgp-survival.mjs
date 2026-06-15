// dgp-survival.mjs — STANDALONE seeded known-truth 2-arm survival DGP for registry-ipd.
//
// Simulates true individual event/censor times from a KNOWN parametric model
// (Weibull per arm; exponential when shape k=1) with a KNOWN hazard ratio, then
// emits exactly the REGISTRY SUMMARY the engine consumes:
//   - km_points: [{t, S}]  digitized KM read-off at fixed registry report times
//   - nar_points: [{t, n}] numbers-at-risk table
//   - N, total_events      per arm (registry-reported counts)
//   - role: experimental / comparator
//
// Truth we later check recovery against:
//   - true median per arm   (Weibull: b*(ln2)^(1/k))
//   - true total events      (administrative censoring at tMax)
//   - true log hazard ratio  = ln(hr)  (proportional hazards by construction)
//
// Self-contained mulberry32 PRNG; no engine imports (independent generator).

function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// Weibull(shape k, scale b) inverse-CDF; k=1 => exponential. Proportional hazards
// across arms achieved by scaling the *hazard* by hr (i.e. scale b *= hr^(-1/k)).
function rweibull(rng, k, b) {
  let u = rng(); if (u <= 0) u = 1e-12;
  return b * Math.pow(-Math.log(u), 1 / k);
}

function simArm(rng, n, k, b, cLambda, tMax) {
  const times = [], status = [];
  for (let i = 0; i < n; i++) {
    const tEvent = rweibull(rng, k, b);
    const tCens = cLambda > 0 ? -Math.log(Math.max(1e-12, rng())) / cLambda : Infinity;
    const t = Math.min(tEvent, tCens, tMax);
    const e = (tEvent <= tCens && tEvent <= tMax) ? 1 : 0;
    times.push(t); status.push(e);
  }
  return { times, status };
}

// KM step function from individual data.
function kmSteps(times, status) {
  if (times.some(t => !Number.isFinite(t))) throw new Error('kmSteps: non-finite time in DGP output (check shape/scale params)');
  const idx = times.map((_, i) => i).sort((a, b) => times[a] - times[b]);
  const n0 = times.length;
  const steps = []; let surv = 1, atRisk = n0, i = 0;
  while (i < n0) {
    const t = times[idx[i]];
    let d = 0, tot = 0, j = i;
    while (j < n0 && times[idx[j]] === t) { tot++; if (status[idx[j]] === 1) d++; j++; }
    if (j === i) j++; // defensive: never stall (cannot happen with finite, sorted times)
    if (d > 0) { surv *= (1 - d / atRisk); steps.push({ t, S: surv }); }
    atRisk -= tot; i = j;
  }
  return steps;
}

function evalKM(steps, t) {
  let S = 1; for (const s of steps) { if (s.t <= t + 1e-9) S = s.S; else break; } return S;
}

function round(x, d) { const f = Math.pow(10, d); return Math.round(x * f) / f; }

function buildArm(rng, cfg, armCfg) {
  const { shape: k, cLambda, tMax, reportTimes, narTimes } = cfg;
  const { n, b, role, arm_id, label } = armCfg;
  const { times, status } = simArm(rng, n, k, b, cLambda, tMax);
  const steps = kmSteps(times, status);
  const km_points = reportTimes.map(t => ({ t, S: round(evalKM(steps, t), 4) }));
  const nar_points = narTimes.map(t => ({ t, n: times.filter(x => x >= t - 1e-9).length }));
  const total_events = status.reduce((a, b) => a + b, 0);
  const median = trueMedian(times, status);
  return {
    arm: { arm_id, label, role, N: n, total_events, follow_up_max: tMax, km_points, nar_points },
    raw: { times, status },
    median,
  };
}

function trueMedian(times, status) {
  const steps = kmSteps(times, status);
  for (const s of steps) if (s.S < 0.5) return s.t; // first step below 0.5
  return null;
}

/**
 * @param {object} cfg {seed,nExp,nCtrl,medianCtrl,hr,shape,cLambda,tMax,reportTimes,narTimes}
 * @returns {trial, truth, raw}
 */
export function generate(cfg) {
  const c = {
    seed: 20260615, nExp: 600, nCtrl: 600,
    medianCtrl: 12, hr: 0.65, shape: 1, cLambda: 0.01, tMax: 36,
    reportTimes: Array.from({ length: 13 }, (_, i) => i * 3),
    narTimes: [0, 6, 12, 18, 24, 30, 36],
    ...(cfg || {}),
  };
  const k = c.shape;
  // control scale from median:  median = b*(ln2)^(1/k)  => b = median/(ln2)^(1/k)
  const bCtrl = c.medianCtrl / Math.pow(Math.log(2), 1 / k);
  // proportional hazards: experimental hazard = hr * ctrl hazard => scale b_exp = b_ctrl * hr^(-1/k)
  const bExp = bCtrl * Math.pow(c.hr, -1 / k);
  const medianExp = bExp * Math.pow(Math.log(2), 1 / k);

  const rng = mulberry32(c.seed);
  const exp = buildArm(rng, c, { n: c.nExp, b: bExp, role: 'experimental', arm_id: 'EXP', label: 'Experimental' });
  const ctrl = buildArm(rng, c, { n: c.nCtrl, b: bCtrl, role: 'comparator', arm_id: 'CTL', label: 'Comparator' });

  const trial = {
    nct_id: 'SYNTH-TRUTH-' + c.seed,
    source: 'truth-recovery DGP', time_unit: 'months',
    arms: [exp.arm, ctrl.arm],
    hr: { value: round(c.hr, 3), method: 'Cox', favors_arm_id: 'EXP' },
  };

  return {
    cfg: c,
    trial,
    truth: {
      medianExpTheory: round(medianExp, 2),
      medianCtrlTheory: round(c.medianCtrl, 2),
      logHRtheory: round(Math.log(c.hr), 4),
      hrTheory: round(c.hr, 3),
      eventsExp: exp.arm.total_events,
      eventsCtrl: ctrl.arm.total_events,
      medianExpSample: exp.median,
      medianCtrlSample: ctrl.median,
    },
    raw: { exp: exp.raw, ctrl: ctrl.raw },
  };
}
