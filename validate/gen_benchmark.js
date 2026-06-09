/*
 * Generate a synthetic BENCHMARK COHORT to exercise the validation pipeline at scale and
 * preview the AACT-only vs digitization head-to-head BEFORE real AACT data lands.
 *
 * For each synthetic trial we know the ground-truth IPD, so we can build:
 *   - benchmark/aact/<id>.json       : exact registry KM anchors  (Tier A, no pixel error)
 *   - benchmark/digitized/<id>.json  : same anchors + pixel-like noise (simulated digitization)
 * Both reconstruct through the SAME engine; metrics compare each to the exact registry anchors.
 *
 * This is a METHODOLOGY demonstration on synthetic data — NOT a claim about real trials. The real
 * numbers come from running validate.js over harvested AACT trials once the snapshot is available.
 *
 * Run: node validate/gen_benchmark.js  then  node validate/validate.js benchmark/aact --digitized benchmark/digitized
 */
const fs = require('fs'), path = require('path');
const ROOT = path.join(__dirname, '..', 'benchmark');
const AACT = path.join(ROOT, 'aact'), DIG = path.join(ROOT, 'digitized');
fs.mkdirSync(AACT, { recursive: true }); fs.mkdirSync(DIG, { recursive: true });

function mulberry32(s){let a=s>>>0;return function(){a|=0;a=(a+0x6D2B79F5)|0;let t=Math.imul(a^(a>>>15),1|a);t=(t+Math.imul(t^(t>>>7),61|t))^t;return((t^(t>>>14))>>>0)/4294967296;};}
function km(ipd){const rows=ipd.slice().sort((a,b)=>a.time-b.time);const ts=[...new Set(rows.map(r=>r.time))].sort((a,b)=>a-b);let n=rows.length,S=1;const st=[];for(const t of ts){const at=rows.filter(r=>Math.abs(r.time-t)<1e-12);const d=at.filter(r=>r.status===1).length,c=at.filter(r=>r.status===0).length;if(d>0)S*=(1-d/n);st.push({t,S});n-=(d+c);}return st;}
function evalKM(st,t){let S=1;for(const s of st){if(s.t<=t+1e-9)S=s.S;else break;}return S;}
function nar(ipd,t){return ipd.filter(r=>r.time>=t-1e-9).length;}
function expArm(seed,N,med,cut){const rng=mulberry32(seed);const lam=Math.log(2)/med;const ipd=[];for(let i=0;i<N;i++){const u=Math.max(1e-12,rng());const T=-Math.log(u)/lam;ipd.push({time:Math.min(T,cut),status:T<=cut?1:0});}return ipd;}
function round(x,d){const f=Math.pow(10,d);return Math.round(x*f)/f;}

// vary trial parameters across a plausible oncology-like range
const SPECS = [];
for (let i = 0; i < 20; i++) {
  const ctlMed = 8 + (i % 5) * 3;            // 8..20
  const hr = 0.55 + (i % 4) * 0.1;           // 0.55..0.85
  const expMed = ctlMed / hr;
  const N = 150 + (i % 3) * 100;             // 150/250/350
  const cut = Math.max(expMed, ctlMed) * 1.6;
  SPECS.push({ id: `BENCH-${String(i + 1).padStart(2, '0')}`, ctlMed, expMed, hr: round(hr, 3), N, cut, seed: 1000 + i * 7 });
}

const anchorTs = (cut) => { const out = []; for (let t = 0; t <= cut; t += cut / 8) out.push(round(t, 2)); return out; };

for (const s of SPECS) {
  const ctl = expArm(s.seed, s.N, s.ctlMed, s.cut);
  const exp = expArm(s.seed + 1, s.N, s.expMed, s.cut);
  const ats = anchorTs(s.cut);
  const narTs = [0, round(s.cut / 2, 2), round(s.cut, 2)];
  function arm(ipd, role, label, id) {
    const st = km(ipd);
    return {
      arm_id: id, label, role, N: s.N,
      total_events: ipd.filter(r => r.status === 1).length, follow_up_max: round(s.cut, 2),
      km_points: ats.map(t => ({ t, S: round(evalKM(st, t), 4) })),
      nar_points: narTs.map(t => ({ t, n: nar(ipd, t) })),
    };
  }
  const exact = {
    nct_id: s.id, source_url: 'synthetic benchmark', time_unit: 'months',
    arms: [arm(exp, 'experimental', 'Drug', 'OG000'), arm(ctl, 'comparator', 'Placebo', 'OG001')],
    hr: { value: s.hr, ci_low: round(s.hr * 0.8, 3), ci_high: round(s.hr * 1.2, 3), method: 'Cox', favors_arm_id: 'OG000' },
  };
  fs.writeFileSync(path.join(AACT, s.id + '.json'), JSON.stringify(exact, null, 2));

  // digitized counterpart: identical EXCEPT km_points carry pixel-like noise (~+/-3% on S, +/-2% on t)
  const rng = mulberry32(s.seed + 555);
  const dig = JSON.parse(JSON.stringify(exact));
  for (const a of dig.arms) {
    a.km_points = a.km_points.map(p => ({
      t: round(Math.max(0, p.t + (rng() - 0.5) * 0.04 * s.cut), 3),
      S: round(Math.min(1, Math.max(0, p.S + (rng() - 0.5) * 0.06)), 4),
    }));
  }
  fs.writeFileSync(path.join(DIG, s.id + '.json'), JSON.stringify(dig, null, 2));
}
console.log(`wrote ${SPECS.length} trials to benchmark/aact and benchmark/digitized`);
