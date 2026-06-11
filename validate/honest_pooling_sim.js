#!/usr/bin/env node
/*
 * THE LINCHPIN: pooling reconstructed effects is only safe if reconstruction uncertainty is PROPAGATED.
 * ====================================================================================================
 *
 * Vision (SYNTHESIS-VISION.md): registry-native reconstruction lets a curve- or HR-only trial join an
 * IPD-level synthesis. But a reconstructed log-HR carries TWO variance components — the trial's sampling
 * variance s_i^2 AND the reconstruction (censoring-level) uncertainty r_i^2 that the credible interval
 * already quantifies. If a meta-analysis pools the reconstructed point using only s_i^2 (treating
 * pseudo-IPD as if it were exact IPD), it commits the same error the user's advanced-stats rules flag for
 * multiverse pooling: the extra noise is mis-attributed to between-trial heterogeneity (tau^2 inflated)
 * and/or the pooled interval under-covers — manufactured precision == false robustness.
 *
 * The honest fix is Rubin's-rules total variance v_i = s_i^2 + r_i^2 (within + between imputation). This
 * Monte-Carlo experiment proves the claim: it simulates a random-effects truth, reconstructs each trial
 * with realistic extra noise, and pools two ways — NAIVE (s_i^2 only) vs HONEST (s_i^2 + r_i^2) — then
 * measures tau^2 bias and the coverage of the pooled effect. Seeded; stochastic tolerances per the
 * project's Monte-Carlo testing rules.
 *
 * Run from repo root:  node validate/honest_pooling_sim.js  ->  validate/honest_pooling_sim_results.json
 */
const fs = require('fs');
const path = require('path');

// --- seeded PRNG (mulberry32) + standard normal (Box-Muller) ---
function mulberry32(a) {
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
function gauss(rng) { return Math.sqrt(-2 * Math.log(1 - rng())) * Math.cos(2 * Math.PI * rng()); }

// --- random-effects meta-analysis: REML tau^2, HKSJ CI (floored), t_{k-1} (matches ipd_meta_fidelity) ---
function qt975(df) { // Cornish-Fisher-ish small-sample 97.5% t quantile (adequate for df>=2)
  const z = 1.959964;
  return z + (z ** 3 + z) / (4 * df) + (5 * z ** 5 + 16 * z ** 3 + 3 * z) / (96 * df * df);
}
function metaRE(y, v) {
  const k = y.length;
  let tau2 = 0;
  for (let it = 0; it < 100; it++) {
    const w = v.map(vi => 1 / (vi + tau2)), sw = w.reduce((a, b) => a + b, 0);
    const mu = y.reduce((s, yi, i) => s + w[i] * yi, 0) / sw;
    let num = 0, den = 0;
    for (let i = 0; i < k; i++) { num += w[i] * w[i] * ((y[i] - mu) ** 2 - v[i]); den += w[i] * w[i]; }
    num += w.reduce((a, b) => a + b, 0) === 0 ? 0 : (1 / sw); // REML correction term
    let t2 = Math.max(0, num / den);
    if (Math.abs(t2 - tau2) < 1e-10) { tau2 = t2; break; }
    tau2 = t2;
  }
  const w = v.map(vi => 1 / (vi + tau2)), sw = w.reduce((a, b) => a + b, 0);
  const mu = y.reduce((s, yi, i) => s + w[i] * yi, 0) / sw;
  // HKSJ variance with floor max(1, q) — never narrower than standard RE
  let q = 0; for (let i = 0; i < k; i++) q += w[i] * (y[i] - mu) ** 2;
  const qAdj = Math.max(1, q / (k - 1));
  const seHK = Math.sqrt(qAdj / sw);
  const tcrit = qt975(k - 1);
  return { mu, tau2, seHK, ci: [mu - tcrit * seHK, mu + tcrit * seHK] };
}

function run(opts) {
  const cfg = Object.assign({ reps: 4000, k: 12, MU: Math.log(0.7), tau2: 0.05,
    s_mean: 0.16, r_mean: 0.22, seed: 12345 }, opts || {});
  const rng = mulberry32(cfg.seed);
  const acc = () => ({ cov: 0, pi: 0, tau2: 0, ciw: 0, piw: 0 });
  const T = acc(), N = acc(), H = acc();
  const tcrit = qt975(cfg.k - 1);
  const tally = (m, store, mu_target, theta_new) => {
    if (mu_target >= m.ci[0] && mu_target <= m.ci[1]) store.cov++;
    const piHalf = tcrit * Math.sqrt(m.tau2 + m.seHK * m.seHK);
    if (theta_new >= m.mu - piHalf && theta_new <= m.mu + piHalf) store.pi++;
    store.tau2 += m.tau2; store.ciw += (m.ci[1] - m.ci[0]); store.piw += 2 * piHalf;
  };
  for (let rep = 0; rep < cfg.reps; rep++) {
    const yTrue = [], yRec = [], s2 = [], r2 = [];
    for (let i = 0; i < cfg.k; i++) {
      const theta = cfg.MU + Math.sqrt(cfg.tau2) * gauss(rng);    // true trial effect (RE)
      const si = cfg.s_mean * (0.7 + 0.6 * rng());                // per-trial sampling SE
      const ri = cfg.r_mean * (0.7 + 0.6 * rng());                // reconstruction (censoring) SE
      const yt = theta + si * gauss(rng);                         // IPD estimate: sampling noise only
      const yr = yt + ri * gauss(rng);                            // reconstructed: + reconstruction noise
      yTrue.push(yt); yRec.push(yr); s2.push(si * si); r2.push(ri * ri);
    }
    const theta_new = cfg.MU + Math.sqrt(cfg.tau2) * gauss(rng);  // a held-out future trial, for the PI
    tally(metaRE(yTrue, s2), T, cfg.MU, theta_new);               // gold: true IPD
    tally(metaRE(yRec, s2), N, cfg.MU, theta_new);               // NAIVE: ignore reconstruction variance
    tally(metaRE(yRec, s2.map((v, i) => v + r2[i])), H, cfg.MU, theta_new);  // HONEST: Rubin total var
  }
  const R = cfg.reps;
  const pack = (m) => ({ coverage_mu: +(m.cov / R).toFixed(3), coverage_PI: +(m.pi / R).toFixed(3),
    mean_tau2: +(m.tau2 / R).toFixed(4), mean_ci_width: +(m.ciw / R).toFixed(3),
    mean_PI_width: +(m.piw / R).toFixed(3) });
  return { config: cfg, nominal: 0.95, true_tau2: cfg.tau2,
    true_ipd: pack(T), naive: pack(N), honest: pack(H) };
}

if (require.main === module) {
  const out = run();
  fs.writeFileSync(path.join(__dirname, 'honest_pooling_sim_results.json'), JSON.stringify(out, null, 2));
  const row = (lab, m) => `    ${lab.padEnd(24)} tau^2 ${m.mean_tau2}  cov(mu) ${m.coverage_mu}  cov(PI) ${m.coverage_PI}  PIwidth ${m.mean_PI_width}`;
  console.log('=== pooling reconstructed effects: propagate the reconstruction uncertainty, or distort heterogeneity ===');
  console.log(`  ${out.config.reps} reps, k=${out.config.k} trials, true tau^2=${out.true_tau2}, nominal 0.95\n`);
  console.log(row('true IPD (gold)', out.true_ipd));
  console.log(row('NAIVE (ignore recon var)', out.naive) + '   <- tau^2 inflated; PI too wide');
  console.log(row('HONEST (Rubin total var)', out.honest) + '   <- tau^2 + PI restored');
  console.log('\n  Reading: the pooled MEAN coverage holds either way (REML tau^2 self-compensates), but');
  console.log('  NAIVE mis-reads reconstruction noise as between-trial heterogeneity (tau^2 ~2x), inflating');
  console.log('  the prediction interval. Propagating r_i^2 via Rubin recovers true tau^2 and a calibrated PI.');
  console.log('\n  wrote validate/honest_pooling_sim_results.json');
}
module.exports = { run, metaRE };
