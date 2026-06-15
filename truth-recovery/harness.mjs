// harness.mjs — wire the repo's OWN registry-native reconstruction to the
// known-truth DGP and MEASURE recovery of median / events / log-HR + HR-CI coverage.
//
// Pipeline per replicate:
//   1. generate()          -> registry-shaped trial (km_points + nar + total_events) + truth
//   2. reconstruct()       [REPO CODE] -> pseudo-IPD per arm
//   3. coxLogHR()          [REPO CODE] on the pseudo-IPD -> reconstructed log-HR
//   4. reconstructEnsemble()[REPO CODE] -> HR/median/RMST credible intervals
//   5. compare against the known truth (theoretical + sample) using the SAME
//      Cox estimator on the TRUE individual data as the reference.
//
// Honest rule (advanced-stats): pseudo-IPD is NEVER exact. Report reconstruction
// error as measured; check it sits within tolerance and report CI coverage honestly.

import { generate } from './dgp-survival.mjs';
import RIPD from './engine.mjs';

const { reconstruct, reconstructEnsemble } = RIPD;
const { coxLogHR, kmFromIPD, medianFromKM } = RIPD._;

function armIpd(result, role) {
  const a = result.arms.find(x => x.role === role) ||
            result.arms.find(x => x.arm_id === (role === 'experimental' ? 'EXP' : 'CTL'));
  return a ? a.ipd : null;
}

function coxFromIpds(expIpd, ctlIpd) {
  const rows = expIpd.map(r => ({ time: r.time, status: r.status, x: 1 }))
    .concat(ctlIpd.map(r => ({ time: r.time, status: r.status, x: 0 })));
  return coxLogHR(rows);
}

function medianOfIpd(ipd) {
  const km = kmFromIPD(ipd);
  return medianFromKM(km, { interpolate: true });
}

export function runReplicate(cfg) {
  const ds = generate(cfg);
  const { trial, truth, raw } = ds;

  // --- REPO RECONSTRUCTION (registry-native) ---
  const rec = reconstruct(trial);
  const expIpd = armIpd(rec, 'experimental');
  const ctlIpd = armIpd(rec, 'comparator');

  const recEventsExp = expIpd.filter(r => r.status === 1).length;
  const recEventsCtrl = ctlIpd.filter(r => r.status === 1).length;
  const recMedExp = medianOfIpd(expIpd);
  const recMedCtrl = medianOfIpd(ctlIpd);
  const recCox = coxFromIpds(expIpd, ctlIpd);

  // --- TRUTH via SAME Cox estimator on true individual data ---
  const trueExp = raw.exp.times.map((t, i) => ({ time: t, status: raw.exp.status[i] }));
  const trueCtl = raw.ctrl.times.map((t, i) => ({ time: t, status: raw.ctrl.status[i] }));
  const trueCox = coxFromIpds(trueExp, trueCtl);

  // --- ENSEMBLE credible interval (repo's uncertainty quantification) ---
  const ens = reconstructEnsemble(trial, { pinEvents: true });
  const hrCI = ens.ensemble && ens.ensemble.hr ? ens.ensemble.hr : null;

  return {
    tier: rec.tier, method: rec.method,
    truth,
    median: { expRec: recMedExp, expTrue: truth.medianExpSample, ctrlRec: recMedCtrl, ctrlTrue: truth.medianCtrlSample },
    events: { expRec: recEventsExp, expTrue: truth.eventsExp, ctrlRec: recEventsCtrl, ctrlTrue: truth.eventsCtrl },
    logHR: {
      rec: recCox.beta, recHR: +recCox.hr.toFixed(4),
      true: trueCox.beta, trueHR: +trueCox.hr.toFixed(4),
      theory: truth.logHRtheory,
      ci: hrCI,
      covered: hrCI ? (Math.log(hrCI.lo) <= truth.logHRtheory && truth.logHRtheory <= Math.log(hrCI.hi)) : null,
    },
  };
}

export function runStudy(baseCfg, nRep = 60) {
  const reps = [];
  for (let r = 0; r < nRep; r++) reps.push(runReplicate({ ...baseCfg, seed: 700000 + r * 101 }));
  const mean = a => a.reduce((x, y) => x + y, 0) / a.length;

  const recLHR = reps.map(r => r.logHR.rec).filter(Number.isFinite);
  const trueLHR = reps.map(r => r.logHR.true).filter(Number.isFinite);
  const biasVsTrue = mean(reps.map(r => r.logHR.rec - r.logHR.true));
  const biasVsTheory = mean(recLHR) - reps[0].logHR.theory;
  const covArr = reps.map(r => r.logHR.covered).filter(v => v !== null);
  const coverage = covArr.length ? mean(covArr.map(v => (v ? 1 : 0))) : null;

  const medExpErr = mean(reps.map(r => Math.abs((r.median.expRec ?? NaN) - (r.median.expTrue ?? NaN))).filter(Number.isFinite));
  const medCtrlErr = mean(reps.map(r => Math.abs((r.median.ctrlRec ?? NaN) - (r.median.ctrlTrue ?? NaN))).filter(Number.isFinite));
  const evExpRel = mean(reps.map(r => Math.abs(r.events.expRec - r.events.expTrue) / r.events.expTrue));
  const evCtrlRel = mean(reps.map(r => Math.abs(r.events.ctrlRec - r.events.ctrlTrue) / r.events.ctrlTrue));

  return {
    nRep, tier: reps[0].tier, method: reps[0].method,
    meanRecHR: +Math.exp(mean(recLHR)).toFixed(4),
    meanTrueHR: +Math.exp(mean(trueLHR)).toFixed(4),
    theoryHR: +Math.exp(reps[0].logHR.theory).toFixed(4),
    biasVsTrueLogHR: +biasVsTrue.toFixed(4),
    biasVsTheoryLogHR: +biasVsTheory.toFixed(4),
    hrCICoverage: coverage == null ? null : +coverage.toFixed(3),
    medExpAbsErr: +medExpErr.toFixed(3),
    medCtrlAbsErr: +medCtrlErr.toFixed(3),
    evExpRelErr: +evExpRel.toFixed(4),
    evCtrlRelErr: +evCtrlRel.toFixed(4),
  };
}
