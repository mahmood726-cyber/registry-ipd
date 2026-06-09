/*
 * Deterministic fixture generator. Writes test/fixtures/*.json from known ground truth
 * so every fixture is self-consistent and hand-checkable. Run: node test/gen_fixtures.js
 *
 * IMPORTANT: ground truth is generated with a seeded PRNG independent of the engine's
 * own RNG, so a passing test reflects genuine recovery, not a tautology.
 */
const fs = require('fs');
const path = require('path');
const OUT = path.join(__dirname, 'fixtures');
fs.mkdirSync(OUT, { recursive: true });

function mulberry32(seed) {
  let a = seed >>> 0;
  return function () { a |= 0; a = (a + 0x6D2B79F5) | 0; let t = Math.imul(a ^ (a >>> 15), 1 | a); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
}

// exact KM from an IPD array [{time,status}]
function km(ipd) {
  const rows = ipd.slice().sort((a, b) => a.time - b.time);
  const times = [...new Set(rows.map(r => r.time))].sort((a, b) => a - b);
  let n = rows.length, S = 1; const steps = [];
  for (const t of times) {
    const at = rows.filter(r => Math.abs(r.time - t) < 1e-12);
    const d = at.filter(r => r.status === 1).length, c = at.filter(r => r.status === 0).length;
    if (d > 0) S *= (1 - d / n);
    steps.push({ t, S, nRisk: n, d, c }); n -= (d + c);
  }
  return steps;
}
function nAtRisk(ipd, t) { return ipd.filter(r => r.time >= t - 1e-9).length; }
function median(steps) { for (const s of steps) if (s.S <= 0.5 + 1e-12) return s.t; return null; }

// ---------- fixture 1: exponential round-trip (Tier A, deterministic) ----------
function genExpArm(seed, N, lambda, cutoff) {
  const rng = mulberry32(seed); const ipd = [];
  for (let i = 0; i < N; i++) { const u = Math.max(1e-12, rng()); const T = -Math.log(u) / lambda; ipd.push({ time: Math.min(T, cutoff), status: T <= cutoff ? 1 : 0 }); }
  return ipd;
}
function tierAfromIPD(ipd, narTimes) {
  const steps = km(ipd);
  const km_points = steps.map(s => ({ t: round(s.t, 6), S: round(s.S, 6) }));
  const nar_points = narTimes.map(t => ({ t, n: nAtRisk(ipd, t) }));
  const total_events = ipd.filter(r => r.status === 1).length;
  return { km_points, nar_points, total_events, N: ipd.length, _true_median: median(steps) };
}
function round(x, d) { const f = Math.pow(10, d); return Math.round(x * f) / f; }

(function exp_known() {
  const ctl = genExpArm(101, 200, Math.log(2) / 12, 36); // median 12 mo
  const exp = genExpArm(202, 200, Math.log(2) / 18, 36); // median 18 mo
  const aCtl = tierAfromIPD(ctl, [0, 6, 12, 18, 24, 30]);
  const aExp = tierAfromIPD(exp, [0, 6, 12, 18, 24, 30]);
  const trial = {
    nct_id: 'SYNTH-EXP-A', source: 'synthetic', time_unit: 'months',
    arms: [
      { arm_id: 'ctl', label: 'Control', role: 'comparator', N: aCtl.N, total_events: aCtl.total_events, follow_up_max: 36, km_points: aCtl.km_points, nar_points: aCtl.nar_points },
      { arm_id: 'exp', label: 'Experimental', role: 'experimental', N: aExp.N, total_events: aExp.total_events, follow_up_max: 36, km_points: aExp.km_points, nar_points: aExp.nar_points }
    ],
    hr: { value: round((Math.log(2) / 18) / (Math.log(2) / 12), 4), ci_low: 0.55, ci_high: 0.98, method: 'Cox', favors_arm_id: 'exp' },
    _truth: { ctl_median: aCtl._true_median, exp_median: aExp._true_median, ctl_events: aCtl.total_events, exp_events: aExp.total_events }
  };
  write('fixture_exp_known.json', trial);
})();

// ---------- fixture 2: integer round-trip parity (Tier A) ----------
// Hand-specified integer death/censor schedule -> exact KM + NAR -> must recover d[] exactly.
(function guyot_roundtrip() {
  // single arm, explicit events at integer times
  const schedule = [ // {t, d, c}
    { t: 1, d: 2, c: 0 }, { t: 2, d: 3, c: 1 }, { t: 3, d: 1, c: 2 },
    { t: 4, d: 4, c: 0 }, { t: 5, d: 2, c: 3 }, { t: 6, d: 1, c: 1 }
  ];
  const ipd = [];
  for (const s of schedule) { for (let i = 0; i < s.d; i++) ipd.push({ time: s.t, status: 1 }); for (let i = 0; i < s.c; i++) ipd.push({ time: s.t, status: 0 }); }
  const N = ipd.length;
  const steps = km(ipd);
  const km_points = steps.map(s => ({ t: s.t, S: round(s.S, 8) }));
  const nar_points = [0, 3, 6].map(t => ({ t, n: nAtRisk(ipd, t) }));
  const total_events = schedule.reduce((a, s) => a + s.d, 0);
  const trial = {
    nct_id: 'SYNTH-GUYOT-RT', source: 'synthetic', time_unit: 'months',
    arms: [{ arm_id: 'a', label: 'Arm', role: 'experimental', N, total_events, follow_up_max: 6, km_points, nar_points }],
    hr: null,
    _truth: { d_schedule: schedule.map(s => s.d), total_events, N }
  };
  write('fixture_guyot_roundtrip.json', trial);
})();

// ---------- fixture 3: Tier B parametric (stochastic) ----------
(function weibull_tierB() {
  const ctlMed = 10, expMed = 16, hr = round((Math.log(2) / expMed) / (Math.log(2) / ctlMed), 4);
  const trial = {
    nct_id: 'SYNTH-TIERB', source: 'synthetic', time_unit: 'months',
    arms: [
      // event fractions must exceed 50% for the KM median to be observable (200/300, 165/300)
      { arm_id: 'ctl', label: 'Control', role: 'comparator', N: 300, total_events: 200, follow_up_max: 30, median: { value: ctlMed, ci_low: 8.5, ci_high: 11.8 }, km_points: [], nar_points: [] },
      { arm_id: 'exp', label: 'Experimental', role: 'experimental', N: 300, total_events: 165, follow_up_max: 30, median: { value: expMed, ci_low: 13.0, ci_high: 20.5 }, km_points: [], nar_points: [] }
    ],
    hr: { value: hr, ci_low: 0.50, ci_high: 0.78, method: 'Cox', favors_arm_id: 'exp' },
    _truth: { hr, logHR: Math.log(hr), exp_median: expMed, ctl_median: ctlMed }
  };
  write('fixture_weibull_tierB.json', trial);
})();

// ---------- fixture 4: Tier C fail-closed ----------
(function tierC() {
  write('fixture_tierC_failclosed.json', {
    nct_id: 'SYNTH-TIERC', source: 'synthetic', time_unit: 'months',
    arms: [
      { arm_id: 'ctl', label: 'Control', role: 'comparator', N: null, total_events: null, km_points: [], nar_points: [], median: null },
      { arm_id: 'exp', label: 'Experimental', role: 'experimental', N: null, total_events: null, km_points: [], nar_points: [], median: null }
    ],
    hr: { value: 0.82, ci_low: 0.66, ci_high: 1.02, method: 'Cox', favors_arm_id: 'exp' }
  });
})();

// ---------- fixture 5: edge cases ----------
(function edge() {
  // zero-event experimental arm (flat S=1, declining NAR from censoring) + one-sided CI
  // + ambiguous HR direction. Zero-event arm is still Tier-A-eligible: 3 flat anchors + NAR.
  const ctl = genExpArm(303, 100, Math.log(2) / 10, 24);
  const aCtl = tierAfromIPD(ctl, [0, 12, 24]);
  const expArm = {
    arm_id: 'exp', label: 'Experimental', role: 'experimental', N: 100, total_events: 0, follow_up_max: 24,
    km_points: [{ t: 0, S: 1 }, { t: 12, S: 1 }, { t: 24, S: 1 }],
    nar_points: [{ t: 0, n: 100 }, { t: 12, n: 60 }, { t: 24, n: 30 }] // pure censoring
  };
  write('fixture_edge.json', {
    nct_id: 'SYNTH-EDGE', source: 'synthetic', time_unit: 'months',
    arms: [
      { arm_id: 'ctl', label: 'Control', role: 'comparator', N: aCtl.N, total_events: aCtl.total_events, follow_up_max: 24, km_points: aCtl.km_points, nar_points: aCtl.nar_points },
      expArm
    ],
    hr: { value: 0.30, ci_low: null, ci_high: 0.70, method: 'Cox', one_sided: true, favors_arm_id: null } // ambiguous direction
  });
})();

function write(name, obj) { fs.writeFileSync(path.join(OUT, name), JSON.stringify(obj, null, 2)); console.log('wrote', name); }
console.log('fixtures generated in', OUT);
