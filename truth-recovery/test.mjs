// test.mjs — truth-recovery assertions for registry-ipd's pseudo-IPD reconstruction.
// Run: node truth-recovery/test.mjs   (exit 0 = all pass)
//
// Encodes the MEASURED recovery of the repo's OWN registry-native reconstruction
// (reconstruct + reconstructEnsemble + coxLogHR) against a known-truth survival DGP.
// Documents both the strong positives (near-exact HR/events/median when event
// counts are present) and the honest negatives (HR-CI under-coverage; curve-only
// attenuation when no event count is posted).

import { runReplicate, runStudy } from './harness.mjs';
import { generate } from './dgp-survival.mjs';
import RIPD from './engine.mjs';
const { reconstruct } = RIPD;
const { coxLogHR } = RIPD._;

let pass = 0, fail = 0;
function assert(name, cond, detail) {
  if (cond) { pass++; console.log('  PASS  ' + name + (detail ? '  [' + detail + ']' : '')); }
  else { fail++; console.log('  FAIL  ' + name + (detail ? '  [' + detail + ']' : '')); }
}

const cfg = { nExp: 600, nCtrl: 600, medianCtrl: 12, hr: 0.65, shape: 1, cLambda: 0.01, tMax: 36 };
const rep = runReplicate({ ...cfg, seed: 42 });
const study = runStudy(cfg, 40);

// --- 1. tier + method routing: rich registry data -> Tier A, Titman QP ---
assert('A1 rich registry trial routes to Tier A / QP method',
  rep.tier === 'A' && rep.method === 'qp', 'tier ' + rep.tier + ' method ' + rep.method);

// --- 2. EVENT-COUNT recovery (QP pins the posted count): within 1% ---
assert('A2 events(Exp) recovered within 1%', study.evExpRelErr < 0.01,
  'rel err ' + study.evExpRelErr);
assert('A3 events(Ctrl) recovered within 1%', study.evCtrlRelErr < 0.01,
  'rel err ' + study.evCtrlRelErr);

// --- 3. MEDIAN recovery: within 0.5 mo of true sample median ---
assert('A4 median recovered within 0.5 mo (both arms)',
  study.medExpAbsErr < 0.5 && study.medCtrlAbsErr < 0.5,
  'exp ' + study.medExpAbsErr + ' / ctrl ' + study.medCtrlAbsErr + ' mo');

// --- 4. log-HR recovery (POSITIVE): bias vs true IPD < 0.03 ---
// This is the headline: the censoring-informed reconstruction recovers the HR
// almost exactly, unlike a naive endpoint-placement digitizer.
assert('A5 reconstructed log-HR ~unbiased vs true (|bias| < 0.03)',
  Math.abs(study.biasVsTrueLogHR) < 0.03,
  'bias vs true = ' + study.biasVsTrueLogHR + ' (recHR ' + study.meanRecHR + ' vs trueHR ' + study.meanTrueHR + ', theory ' + study.theoryHR + ')');

// --- 5. HONEST NEGATIVE: ensemble HR-CI under-covers the truth ---
// Coverage measured ~0.78 over seeds (nominal 0.95). The point estimate is
// excellent but the credible interval is too narrow. Assert it is below 0.90 so
// this test FAILS LOUDLY if the CI width is ever changed without re-validation.
assert('A6 HR-CI under-covers true log-HR (coverage < 0.90, honest negative)',
  study.hrCICoverage !== null && study.hrCICoverage < 0.90,
  'coverage = ' + study.hrCICoverage + ' (nominal 0.95)');

// --- 6. HONEST NEGATIVE: curve-only (no event count) attenuates the HR ---
// Strip total_events -> engine falls back to anchor-exact; HR drifts toward null.
function recHRcurveOnly(seed) {
  const ds = generate({ ...cfg, seed });
  const r = reconstruct(ds.trial, { ignoreTotalEvents: true });
  const e = r.arms.find(a => a.role === 'experimental').ipd;
  const c = r.arms.find(a => a.role === 'comparator').ipd;
  const rows = e.map(x => ({ time: x.time, status: x.status, x: 1 }))
    .concat(c.map(x => ({ time: x.time, status: x.status, x: 0 })));
  return { method: r.method, hr: coxLogHR(rows).hr };
}
const co = recHRcurveOnly(42);
assert('A7 curve-only branch uses Guyot/anchor-exact (no QP) and attenuates toward null',
  co.method !== 'qp' && co.hr > rep.logHR.recHR,
  'curve-only method ' + co.method + ' HR ' + co.hr.toFixed(4) + ' vs QP HR ' + rep.logHR.recHR);

console.log('\n' + pass + ' passed, ' + fail + ' failed');
process.exit(fail === 0 ? 0 : 1);
