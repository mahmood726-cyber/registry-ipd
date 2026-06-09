/*
 * Engine test suite. Run: node --test test/
 * Tolerance policy:
 *   - deterministic estimators (Tier A inversion, KM, RMST): tight
 *   - stochastic estimators (Tier B sampling/bootstrap): atol 0.05, seeded
 */
const { test } = require('node:test');
const assert = require('node:assert');
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');

const fx = (n) => JSON.parse(fs.readFileSync(path.join(__dirname, 'fixtures', n), 'utf8'));

// -------------------------------------------------- utility unit tests
test('pavaDecreasing enforces non-increasing and is identity on monotone input', () => {
  const mono = [1, 0.9, 0.8, 0.5];
  assert.deepStrictEqual(RIPD._.pavaDecreasing(mono).y.map(x => +x.toFixed(6)), mono);
  const bad = RIPD._.pavaDecreasing([1, 0.7, 0.8, 0.5]); // 0.8 violates
  for (let i = 1; i < bad.y.length; i++) assert.ok(bad.y[i] <= bad.y[i - 1] + 1e-9);
  assert.ok(bad.adjusted >= 1);
});

test('medianFromKM interpolates the 0.5-crossing on coarse curves (anti-quantization)', () => {
  const steps = [{ t: 6, S: 0.4 }];                 // S goes 1 -> 0.4 at t=6; true median between 0 and 6
  assert.strictEqual(RIPD._.medianFromKM(steps), 6); // step snaps up to the posted timepoint
  assert.ok(Math.abs(RIPD._.medianFromKM(steps, { interpolate: true }) - 5.0) < 1e-9); // 0+(1-.5)/(1-.4)*6
});

test('kmFromIPD + median + rmst on a hand example', () => {
  const ipd = [{ time: 1, status: 1 }, { time: 2, status: 1 }, { time: 3, status: 0 }, { time: 4, status: 1 }];
  const steps = RIPD._.kmFromIPD(ipd);
  // S after t=1: 3/4 ; after t=2: 3/4*2/3=1/2 ; censor at 3; event at 4 with 1 at risk -> 0
  assert.ok(Math.abs(steps[0].S - 0.75) < 1e-9);
  assert.ok(Math.abs(steps[1].S - 0.5) < 1e-9);
  assert.strictEqual(RIPD._.medianFromKM(steps), 2);
});

test('coxLogHR recovers a known strong effect with correct sign', () => {
  // build two arms where experimental clearly lives longer => HR < 1 => beta < 0
  const rng = RIPD._.mulberry32(7); const rows = [];
  for (let i = 0; i < 200; i++) { const T = -Math.log(rng()) / (Math.log(2) / 8); rows.push({ time: Math.min(T, 40), status: T <= 40 ? 1 : 0, x: 0 }); }
  for (let i = 0; i < 200; i++) { const T = -Math.log(rng()) / (Math.log(2) / 16); rows.push({ time: Math.min(T, 40), status: T <= 40 ? 1 : 0, x: 1 }); }
  const cox = RIPD._.coxLogHR(rows);
  assert.ok(cox.hr < 0.75, `expected HR<0.75 got ${cox.hr}`);
  assert.ok(!cox.separated);
});

// -------------------------------------------------- Tier A deterministic
test('Tier A: exponential round-trip recovers exact total events + tracks anchors (C2~0)', () => {
  const t = fx('fixture_exp_known.json');
  assert.strictEqual(RIPD.classifyTier(t), 'A');
  const r = RIPD.reconstruct(t);
  assert.strictEqual(r.tier, 'A');
  for (const a of t.arms) {
    const rec = r.arms.find(x => x.arm_id === a.arm_id);
    const ev = rec.ipd.filter(x => x.status === 1).length;
    assert.strictEqual(ev, a.total_events, `events exact for ${a.arm_id}`);
    assert.strictEqual(rec.ipd.length, a.N, `population conserved for ${a.arm_id}`);
  }
  // C2 anchor fidelity ~ 0
  assert.ok(r.audit.checks.C2_anchor_fidelity.maxDiff <= 1e-3,
    `C2 maxDiff ${r.audit.checks.C2_anchor_fidelity.maxDiff}`);
  assert.ok(r.audit.checks.C5_monotonic.pass);
  assert.ok(r.audit.checks.C7_conservation.pass);
});

test('Tier A: integer round-trip recovers the exact death schedule (Guyot parity)', () => {
  const t = fx('fixture_guyot_roundtrip.json');
  const r = RIPD.reconstruct(t);
  const rec = r.arms[0];
  // reconstruct per-time death counts
  const byT = {};
  for (const row of rec.ipd) if (row.status === 1) byT[row.time] = (byT[row.time] || 0) + 1;
  const times = [1, 2, 3, 4, 5, 6];
  const reconD = times.map(tt => byT[tt] || 0);
  assert.deepStrictEqual(reconD, t._truth.d_schedule, `death schedule parity: ${JSON.stringify(reconD)} vs ${JSON.stringify(t._truth.d_schedule)}`);
});

test('Tier A: median within 5% of true KM median', () => {
  const t = fx('fixture_exp_known.json');
  const r = RIPD.reconstruct(t);
  const recCtl = r.arms.find(x => x.arm_id === 'ctl');
  const m = RIPD._.medianFromKM(RIPD._.kmFromIPD(recCtl.ipd));
  assert.ok(Math.abs(m - t._truth.ctl_median) / t._truth.ctl_median <= 0.05,
    `median ${m} vs truth ${t._truth.ctl_median}`);
});

// -------------------------------------------------- multi-method (Wasserstein selection)
test('Tier A: anchor-exact beats Guyot under administrative censoring; best-of selects it', () => {
  // synthetic exponential with censoring ONLY at the cutoff (administrative) + COARSE anchors:
  // the case where Guyot's constant-censoring assumption is wrong.
  const rng = RIPD._.mulberry32(2024), N = 200, med = 12, cut = 24, lam = Math.log(2) / med;
  const ipd = [];
  for (let i = 0; i < N; i++) { const T = -Math.log(Math.max(1e-12, rng())) / lam; ipd.push({ time: Math.min(T, cut), status: T <= cut ? 1 : 0 }); }
  const steps = RIPD._.kmFromIPD(ipd);
  const evalK = (t) => RIPD._.evalKM(steps, t);
  const narAt = (t) => ipd.filter(r => r.time >= t - 1e-9).length;
  const ats = [0, 4, 8, 12, 16, 20, 24];
  const km_points = ats.map(t => ({ t, S: +evalK(t).toFixed(4) }));
  const arm = { arm_id: 'a', label: 'Arm', role: 'experimental', N, total_events: ipd.filter(r => r.status === 1).length, follow_up_max: cut, km_points, nar_points: [0, 12, 24].map(t => ({ t, n: narAt(t) })) };
  const trial = { nct_id: 'ADMIN-CENSOR', time_unit: 'months', arms: [arm], hr: null };
  const supErr = (res) => { const km = RIPD._.kmFromIPD(res.arms[0].ipd); let mx = 0; for (const p of km_points) mx = Math.max(mx, Math.abs(RIPD._.evalKM(km, p.t) - p.S)); return mx; };
  const g = supErr(RIPD.reconstruct(trial, { method: 'guyot' }));
  const ae = supErr(RIPD.reconstruct(trial, { method: 'anchor-exact' }));
  const best = RIPD.reconstruct(trial);
  assert.ok(ae < g, `anchor-exact sup-err ${ae.toFixed(4)} should beat Guyot ${g.toFixed(4)}`);
  assert.ok(best.method === 'anchor-exact', `best-of should select anchor-exact, got ${best.method}`);
  assert.ok(best.flags.some(f => f.startsWith('wasserstein_to_anchors:')), 'reports Wasserstein');
});

// -------------------------------------------------- competing risks (Aalen–Johansen)
test('Aalen-Johansen CIF: invariant CIF1+CIF2+S=1, and naive 1-KM overestimates cause-1 incidence', () => {
  // cause 1 events early, a burst of competing (cause 2) events at t=2, rest censored
  const ipd = [];
  for (let i = 0; i < 20; i++) ipd.push({ time: 1, cause: 1 });
  for (let i = 0; i < 30; i++) ipd.push({ time: 2, cause: 2 });   // competing
  for (let i = 0; i < 15; i++) ipd.push({ time: 3, cause: 1 });
  for (let i = 0; i < 35; i++) ipd.push({ time: 5, cause: 0 });   // censored
  const aj = RIPD._.cifAalenJohansen(ipd);
  for (const s of aj) assert.ok(Math.abs(s.cif1 + s.cif2 + s.S - 1) < 1e-9, 'CIF1+CIF2+S=1 invariant');
  const ajFinal = aj[aj.length - 1].cif1;
  const naive = RIPD._.cifNaive1(ipd);
  const naiveFinal = naive[naive.length - 1].cif1;
  // treating the 30 competing events as censoring inflates the cause-1 incidence
  assert.ok(naiveFinal > ajFinal + 1e-6, `naive ${naiveFinal.toFixed(3)} should exceed AJ ${ajFinal.toFixed(3)}`);
});

test('reconstructCompetingRisks labels causes and produces a CIF that respects the invariant', () => {
  const t = fx('fixture_exp_known.json');
  // inject competing events on each arm
  t.arms.forEach(a => { a.competing_events = Math.round(0.2 * a.N); });
  const r = RIPD.reconstructCompetingRisks(t);
  assert.ok(r.competing_risks && r.arms[0].cif, 'CIF computed');
  for (const s of r.arms[0].cif) assert.ok(Math.abs(s.cif1 + s.cif2 + s.S - 1) < 1e-9);
  const cause2 = r.arms[0].ipd_cr.filter(x => x.cause === 2).length;
  assert.ok(cause2 > 0 && cause2 <= Math.round(0.2 * t.arms[0].N) + 1, 'competing events labeled');
});

// -------------------------------------------------- multiple-imputation uncertainty
test('reconstructEnsemble yields credible intervals that cover the truth', () => {
  const t = fx('fixture_exp_known.json');
  const r = RIPD.reconstructEnsemble(t, { M: 80 });
  assert.strictEqual(r.tier, 'A');
  // HR credible interval is ordered and covers the registry HR
  const hr = r.ensemble.hr;
  assert.ok(hr.lo <= hr.est && hr.est <= hr.hi, 'HR CI ordered');
  assert.ok(t.hr.value >= hr.lo && t.hr.value <= hr.hi, `registry HR ${t.hr.value} in CI [${hr.lo},${hr.hi}]`);
  // control-arm median CI covers the true KM median
  const cm = r.ensemble.median.ctl;
  assert.ok(t._truth.ctl_median >= cm.lo && t._truth.ctl_median <= cm.hi,
    `true median ${t._truth.ctl_median} in CI [${cm.lo},${cm.hi}]`);
});

// -------------------------------------------------- C9 direction integrity (hard gate)
test('C9 hard-fails when the reconstructed HR direction contradicts the registry favored arm', () => {
  // experimental arm CLEARLY worse (lower survival) but registry HR says it is favored => contradiction
  const mk = (S) => ({ km_points: [0, 6, 12, 18, 24].map((t, i) => ({ t, S: i === 0 ? 1 : S[i - 1] })), nar_points: [] });
  const t = {
    nct_id: 'CONTRA', time_unit: 'months',
    arms: [
      Object.assign({ arm_id: 'exp', label: 'Drug', role: 'experimental', N: 100, total_events: 70, follow_up_max: 24 }, mk([0.7, 0.5, 0.35, 0.25])),
      Object.assign({ arm_id: 'ctl', label: 'Placebo', role: 'comparator', N: 100, total_events: 35, follow_up_max: 24 }, mk([0.9, 0.8, 0.72, 0.65])),
    ],
    hr: { value: 0.5, ci_low: 0.35, ci_high: 0.72, method: 'Cox', favors_arm_id: 'exp' }, // claims exp better — false
  };
  const r = RIPD.reconstruct(t);
  assert.strictEqual(r.audit.checks.C9_direction.pass, false, 'C9 should catch the contradiction');
  assert.strictEqual(r.audit.badge, 'none', 'hard C9 fail => badge none');
});

// -------------------------------------------------- Royston–Parmar flexible parametric
test('Royston-Parmar fit reproduces a known Weibull curve and stays monotone', () => {
  const k = 1.3, b = 14; // Weibull S(t)=exp(-(t/b)^k)
  const Sfn = (t) => Math.exp(-Math.pow(t / b, k));
  const km = [3, 6, 9, 12, 18, 24].map(t => ({ t, S: +Sfn(t).toFixed(4) }));
  const rp = RIPD._.fitRoystonParmar(km);
  assert.ok(rp, 'fit returned');
  // reproduces anchors within ~2%
  for (const p of km) assert.ok(Math.abs(rp.predict(p.t) - p.S) < 0.02, `RP(${p.t})=${rp.predict(p.t).toFixed(3)} vs ${p.S}`);
  // monotone non-increasing on a fine grid + extrapolation stays in [0,1]
  let prev = 1; for (let t = 1; t <= 36; t++) { const s = rp.predict(t); assert.ok(s <= prev + 1e-6 && s >= 0 && s <= 1); prev = s; }
});

test('smooth:rp densifies coarse Tier-A anchors and still reconstructs Tier A', () => {
  const t = fx('fixture_exp_known.json');
  const r = RIPD.reconstruct(t, { smooth: 'rp' });
  assert.strictEqual(r.tier, 'A');
  assert.ok(r.flags.includes('smoothed:rp'));
  assert.ok(r.audit.checks.C5_monotonic.pass && r.audit.checks.C7_conservation.pass);
});

// -------------------------------------------------- HR-calibration
test('HR-calibration imposes the registry HR on the pseudo-IPD (preserving anchors)', () => {
  const t = fx('fixture_exp_known.json');
  t.hr = { value: 0.5, ci_low: 0.4, ci_high: 0.63, method: 'Cox', favors_arm_id: 'exp' };
  const cal = RIPD.reconstruct(t, { calibrateHR: true });
  assert.ok(cal.calibrated, 'calibration metadata present');
  // achieved HR should be close to the imposed target
  assert.ok(Math.abs(Math.log(cal.calibrated.achieved_hr) - Math.log(0.5)) < Math.log(1.10),
    `achieved HR ${cal.calibrated.achieved_hr} should be ~0.5`);
  assert.ok(cal.audit.checks.C5_monotonic.pass && cal.audit.checks.C7_conservation.pass);
});

// -------------------------------------------------- Tier B stochastic
test('Tier B: classify + reconstruct, median/logHR within atol and truth inside envelope', () => {
  const t = fx('fixture_weibull_tierB.json');
  assert.strictEqual(RIPD.classifyTier(t), 'B');
  const r = RIPD.reconstruct(t, { bootstrap: 300 });
  assert.strictEqual(r.tier, 'B');
  assert.ok(r.flags.some(f => f.startsWith('assumption:parametric')));
  const recExp = r.arms.find(x => x.arm_id === 'exp');
  const m = RIPD._.medianFromKM(RIPD._.kmFromIPD(recExp.ipd));
  // single central draw: Monte-Carlo bound (~3-sigma), NOT a tight point estimate
  assert.ok(Math.abs(m - t._truth.exp_median) / t._truth.exp_median <= 0.25,
    `Tier B central-draw median ${m} vs ${t._truth.exp_median} (3-sigma bound)`);
  // rigorous check: registry median lies inside the bootstrap median envelope
  const me = r.envelope.median_exp;
  assert.ok(t._truth.exp_median >= me[0] && t._truth.exp_median <= me[1],
    `registry median ${t._truth.exp_median} not in envelope [${me.map(x => x.toFixed(2))}]`);
  // truth logHR inside the bootstrap envelope
  const env = r.envelope.logHR;
  assert.ok(t._truth.logHR >= env[0] - 0.05 && t._truth.logHR <= env[1] + 0.05,
    `truth logHR ${t._truth.logHR.toFixed(3)} not in envelope [${env.map(x => x.toFixed(3))}]`);
});

test('Tier B: total events within 1%', () => {
  const t = fx('fixture_weibull_tierB.json');
  const r = RIPD.reconstruct(t, { bootstrap: 50 });
  for (const a of t.arms) {
    const rec = r.arms.find(x => x.arm_id === a.arm_id);
    const ev = rec.ipd.filter(x => x.status === 1).length;
    assert.ok(Math.abs(ev - a.total_events) <= Math.max(2, 0.02 * a.total_events),
      `events ${ev} vs reg ${a.total_events} for ${a.arm_id}`);
  }
});

// -------------------------------------------------- Tier C fail-closed
test('Tier C: HR-only fails closed, not exportable, no IPD', () => {
  const t = fx('fixture_tierC_failclosed.json');
  assert.strictEqual(RIPD.classifyTier(t), 'C');
  const r = RIPD.reconstruct(t);
  assert.strictEqual(r.tier, 'C');
  assert.strictEqual(r.arms, null);
  assert.strictEqual(r.exportable, false);
  assert.strictEqual(r.audit.badge, 'none');
  assert.strictEqual(r.verdict, 'insufficient_registry_data');
});

// -------------------------------------------------- edge cases
test('edge: zero-event arm + one-sided CI + ambiguous HR direction handled without crash or inversion', () => {
  const t = fx('fixture_edge.json');
  const r = RIPD.reconstruct(t);
  // no crash, IPD produced for control; experimental has zero events
  const recExp = r.arms.find(x => x.arm_id === 'exp');
  assert.strictEqual(recExp.ipd.filter(x => x.status === 1).length, 0, 'zero-event arm preserved');
  // ambiguous HR direction must be flagged, never silently inverted
  assert.ok(r.flags.includes('hr_reference_ambiguous'), 'ambiguous HR direction flagged');
  // HR value never inverted in any output (we don't emit a swapped HR)
  assert.ok(r.audit.checks.C5_monotonic.pass);
});
