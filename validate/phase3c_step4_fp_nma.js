/*
 * PHASE 3c step 4 — the literal Jansen fractional-polynomial survival NMA, with reconstruction UQ.
 *
 * Closes the second "claimed, not yet proven" caveat of SYNTHESIS-VISION.md Sec 7: Jansen's (2011)
 * fractional-polynomial network meta-analysis of survival data models a TIME-VARYING log-HR(t) and, in
 * practice, ingests digitised pseudo-IPD as if exact -- ignoring reconstruction uncertainty. Sec 4i proved
 * the fix on a piecewise-exponential analogue (SurvivalNPHPooler); this step proves it on the LITERAL FP
 * parameterisation, reusing allmeta's own FP-NMA engine (HTA/src/engine/fpNMA.js, the FPNMAEngine class).
 *
 * The engine fits log(HR(t)) = d * basis(t, power) by inverse-variance weighted least squares over each
 * study's reported (timePoint, HR, se) tuples -- exactly the place to encode the reconstruction variance:
 * a reconstructed-curve study's late time-points are the LEAST identified (fewest at risk, most censoring;
 * the Sec 4i finding, time-resolved), so they carry a late-growing reconstruction bias AND a late-growing
 * reconstruction variance r(t)^2. Ignoring r(t) (naive: se = sampling only) equal-weights the biased late
 * points so the WLS fit is pulled off at late times; encoding it (honest: se = sqrt(sampling^2 + r(t)^2),
 * the Sec 5b weight w = 1/(var+delta^2)) down-weights them and recovers the curve.
 *
 * True effect: log(HR(t)) = d_true * log(t)  (power 0; a crossing/time-varying hazard). 6 IPD + 6
 * reconstructed-curve studies, A vs B, seeded Monte-Carlo. Reuses C:\Projects\allmeta. Run from repo root:
 *   node validate/phase3c_step4_fp_nma.js  ->  validate/phase3c_step4_results.json
 */
const fs = require('fs');
const path = require('path');
const { FPNMAEngine } = require('C:/Projects/allmeta/HTA/src/engine/fpNMA.js');

const OUT = path.join(__dirname, 'phase3c_step4_results.json');

const P_TRUE = 0;                 // FP power 0 -> basis(t)=log(t)
const D_TRUE = -0.5;              // log(HR(t)) = -0.5*log(t): HR>1 early, HR<1 late (time-varying)
const TS = [0.5, 1, 2, 4];        // reported time-points per study (t=4 is the least-identified late point)
const TMAX = 4;
const S_SAMP = 0.15;              // sampling SD of a reported logHR
const B_LATE = 0.45;             // reconstruction bias coefficient, concentrated late (Sec 4i, time-resolved)
const R_LATE = 0.5;              // reconstruction SD coefficient, grows late
const EVAL_LATE = 4;             // evaluate the fitted curve at this late time

const trueLogHR = (t) => D_TRUE * Math.log(Math.max(t, 1e-9));
const reconBias = (t) => B_LATE * Math.pow(t / TMAX, 2);     // grows quadratically, ~0 early, B_LATE at TMAX
const reconSD = (t) => R_LATE * (t / TMAX);                  // grows linearly with time

// seeded RNG (mulberry32 + Box-Muller) so the committed result is reproducible
function mulberry32(a) {
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
function makeNormal(seed) {
  const u = mulberry32(seed);
  return (mean = 0, sd = 1) => {
    const u1 = Math.max(u(), 1e-12), u2 = u();
    return mean + sd * Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
  };
}

// build the FPNMAEngine `data` array for one rep under a given mode
function buildData(rnd, mode) {
  const data = [];
  for (let j = 0; j < 12; j++) {
    const reconstructed = (j % 2 === 1);
    // reference arm A: HR=1 at every time, modest se
    data.push({ treatment: 'A', timePoints: [...TS], hazardRatios: TS.map(() => 1), ses: TS.map(() => S_SAMP) });
    // treatment arm B: reported logHR(t) with sampling noise (+ reconstruction bias/noise if reconstructed)
    const hrs = [], ses = [];
    for (const t of TS) {
      let lhr = trueLogHR(t) + rnd(0, S_SAMP);
      let se = S_SAMP;
      if (reconstructed && mode !== 'true') {
        lhr += reconBias(t) + rnd(0, reconSD(t));            // late-growing bias + noise
        se = (mode === 'honest') ? Math.sqrt(S_SAMP * S_SAMP + reconSD(t) * reconSD(t)) : S_SAMP;
      }
      hrs.push(Math.exp(lhr));
      ses.push(se);
    }
    data.push({ treatment: 'B' + (j), timePoints: [...TS], hazardRatios: hrs, ses, _b: 'B' });
    // NOTE: all treatment arms must map to the same comparison; relabel below
  }
  // FPNMAEngine groups by `treatment`; we want one B-vs-A comparison pooling all studies, so use 'B'
  for (const d of data) if (d._b === 'B') { d.treatment = 'B'; delete d._b; }
  return data;
}

function fitCoef(data) {
  const eng = new FPNMAEngine({ seed: 12345 });
  const m = eng.fit(data, { powers: [P_TRUE], order: 1, reference: 'A' });
  return m.treatmentEffects[0].coefficients[0];        // d in log(HR(t)) = d*log(t)
}

function run(reps = 400, seed = 20260612) {
  const rnd = makeNormal(seed);
  const acc = { true: [], naive: [], honest: [] };
  const rmseAcc = { true: [], naive: [], honest: [] };
  const gridTimes = [0.5, 1, 1.5, 2, 3, 4];
  for (let r = 0; r < reps; r++) {
    for (const mode of ['true', 'naive', 'honest']) {
      const coef = fitCoef(buildData(rnd, mode));
      acc[mode].push(coef * Math.log(EVAL_LATE) - trueLogHR(EVAL_LATE));   // late logHR bias
      let sse = 0;
      for (const t of gridTimes) { const e = coef * Math.log(t) - trueLogHR(t); sse += e * e; }
      rmseAcc[mode].push(Math.sqrt(sse / gridTimes.length));
    }
  }
  const mean = (a) => a.reduce((s, x) => s + x, 0) / (a.length || 1);
  const round = (x, n = 4) => Math.round(x * 10 ** n) / 10 ** n;
  const out = {
    config: {
      reps, network: 'A vs B, 6 IPD + 6 reconstructed-curve studies, literal FP power 0 (log t)',
      true_logHR_function: 'log(HR(t)) = -0.5 * log(t)', eval_late_time: EVAL_LATE,
      true_logHR_at_late: round(trueLogHR(EVAL_LATE)),
      late_recon_bias_at_tmax: B_LATE, late_recon_sd_at_tmax: R_LATE,
      engine: 'allmeta HTA/src/engine/fpNMA.js (FPNMAEngine)',
      note: 'Reconstructed studies carry a late-growing reconstruction bias + variance (the Sec 4i least-'
        + 'identified-late finding). naive ses ignore r(t) (equal weight -> biased late points pull the WLS '
        + 'FP fit); honest ses = sqrt(sampling^2 + r(t)^2) down-weight them (Sec 5b w = 1/(var+delta^2)).',
    },
  };
  for (const mode of ['true', 'naive', 'honest']) {
    out[mode] = { late_logHR_bias: round(mean(acc[mode])), curve_rmse: round(mean(rmseAcc[mode])) };
  }
  return out;
}

function main() {
  const out = run();
  fs.writeFileSync(OUT, JSON.stringify(out, null, 2));
  const c = out.config;
  console.log('=== Phase 3c step 4: literal Jansen FP survival NMA with reconstruction UQ (allmeta FPNMAEngine) ===\n');
  console.log(`  network ${c.network}`);
  console.log(`  truth ${c.true_logHR_function}; evaluate at late t=${c.eval_late_time} (true logHR ${c.true_logHR_at_late})`);
  console.log(`  late reconstruction bias ${c.late_recon_bias_at_tmax}, late reconstruction SD ${c.late_recon_sd_at_tmax}, ${c.reps} reps\n`);
  console.log(`  ${''.padEnd(28)}${'late logHR bias'.padEnd(18)}${'curve RMSE'}`);
  for (const [m, lab] of [['true', 'all-IPD (gold)'], ['naive', 'NAIVE (ignore r(t))'], ['honest', 'HONEST (encode r(t) in weight)']]) {
    console.log(`  ${lab.padEnd(28)}${String(out[m].late_logHR_bias).padEnd(18)}${out[m].curve_rmse}`);
  }
  console.log('\n  Reading: Jansen FP-NMA fits log(HR(t)) by inverse-variance weighted least squares. A reconstructed');
  console.log('  curve is least identified at LATE times, so its late log-HR points carry bias + extra variance.');
  console.log('  Ignoring that variance (naive) equal-weights the biased late points and pulls the fitted HR(t)');
  console.log('  curve off at late t; encoding it in the weight (honest) down-weights them and recovers the curve');
  console.log('  toward the all-IPD gold -- the Sec 4i lesson on the literal FP engine, not just its PWE analogue.');
  console.log(`\n  wrote ${path.relative(process.cwd(), OUT).replace(/\\/g, '/')}`);
  return out;
}

if (require.main === module) main();
module.exports = { run, trueLogHR, reconBias };
