// Phase 4: the evidence-completeness atlas. Unit-tests the per-trial HR-identification classifier on
// crafted trials (deterministic, no cohort dependency), then asserts the committed atlas summary's
// invariants if realipd/evidence_atlas.json is present (cohort/ is gitignored, so the committed artifact
// is the thing to check, like coverage_census.json).
const { test } = require('node:test');
const assert = require('node:assert');
const fs = require('fs');
const path = require('path');
const { classify } = require('../validate/phase4_evidence_atlas.js');

const arm = (over) => Object.assign({ arm_id: 'a', N: 200, total_events: null, km_points: [], nar_points: [] }, over);
const km = (n) => Array.from({ length: n }, (_, i) => ({ t: i + 1, S: 0.9 - 0.05 * i }));

test('classifier places trials on the granularity manifold with the right identification class', () => {
  // numbers-at-risk present -> point-identified (Guyot)
  const nar = classify({ nct_id: 'N1', arms: [arm({ km_points: km(8), nar_points: [{ t: 1, n: 180 }] }),
    arm({ arm_id: 'b', km_points: km(8), nar_points: [{ t: 1, n: 175 }] })] });
  assert.strictEqual(nar.tier, 'point_nar');
  assert.ok(nar.infoScore > 0.95, 'NAR trial is essentially fully identified');

  // total events on every arm, no NAR -> Titman-QP point-identified
  const qp = classify({ nct_id: 'N2', arms: [arm({ total_events: 60, km_points: km(8) }),
    arm({ arm_id: 'b', total_events: 70, km_points: km(8) })] });
  assert.strictEqual(qp.tier, 'point_qp');
  assert.ok(qp.deltaRecon > 0 && qp.deltaRecon < 0.1);

  // a posted HR with CI -> logHR point-identified directly, no reconstruction penalty
  const hr = classify({ nct_id: 'N3', hr: { value: 0.7, ci_low: 0.55, ci_high: 0.9 },
    arms: [arm({}), arm({ arm_id: 'b' })] });
  assert.strictEqual(hr.tier, 'point_hr');
  assert.strictEqual(hr.deltaRecon, 0);
  assert.ok(Math.abs(hr.infoScore - 1) < 1e-9, 'posted HR is point-identified -> infoScore 1');

  // curve only, no lever -> PARTIALLY identified; sparse anchors widen the region (lower info score)
  const dense = classify({ nct_id: 'N4', arms: [arm({ km_points: km(8) }), arm({ arm_id: 'b', km_points: km(8) })] });
  const sparse = classify({ nct_id: 'N5', arms: [arm({ km_points: km(3) }), arm({ arm_id: 'b', km_points: km(3) })] });
  assert.strictEqual(dense.tier, 'partial_curve');
  assert.strictEqual(sparse.tier, 'partial_curve');
  assert.ok(sparse.deltaRecon > dense.deltaRecon, 'fewer anchors -> wider identification region');
  assert.ok(dense.infoScore < 1 && dense.infoScore > 0, 'curve-only is partially identified');

  // single arm / no curve & no HR -> not contributable to an HR synthesis
  const none = classify({ nct_id: 'N6', arms: [arm({})] });
  assert.strictEqual(none.tier, 'none');
  assert.strictEqual(none.contributable, false);
});

test('committed evidence atlas summary is internally consistent', (t) => {
  const p = path.join(__dirname, '..', 'realipd', 'evidence_atlas.json');
  if (!fs.existsSync(p)) { t.skip('realipd/evidence_atlas.json not generated'); return; }
  const s = JSON.parse(fs.readFileSync(p, 'utf8')).summary;
  const tc = s.tier_counts;
  // every trial lands in exactly one tier
  const point = tc.point_nar + tc.point_qp + tc.point_hr;
  assert.strictEqual(point + tc.partial_curve + tc.none, s.n_trials);
  // the scored point/partial split covers the contributable set (a tiered trial with no computable
  // sampling variance can be unscorable, so allow a small slack vs the raw tier counts)
  assert.strictEqual(s.point_identified_trials + s.partial_identified_trials, s.n_contributable_to_HR);
  assert.ok(point + tc.partial_curve - s.n_contributable_to_HR <= 3, 'few unscorable tiered trials');
  // trial-count and information-weighted point/partial each sum to ~100%
  assert.ok(Math.abs(s.point_identified_trial_pct + s.partial_identified_trial_pct - 100) < 0.2);
  assert.ok(Math.abs(s.point_identified_info_weighted_pct + s.partial_identified_info_weighted_pct - 100) < 0.2);
  // the binding finding: structured numbers-at-risk are essentially absent
  assert.ok(tc.point_nar <= 2, 'AACT structured NAR is ~0 (the binding limitation)');
  // a real corpus has a partially-identified majority by trial count -- the gap the pipeline fills
  assert.ok(s.partial_identified_trials > s.point_identified_trials);
  assert.ok(s.median_info_score > 0 && s.median_info_score <= 1);
});
