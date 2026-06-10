#!/usr/bin/env node
/*
 * PRODUCTION GALLERY — the method running on the REAL AACT trials it is built for (not the open-IPD
 * validation set). Scans the harvested cohort (cohort/*.json), reconstructs each Tier-A trial with the
 * shipped engine (Titman-QP default), and for trials that ALSO report a Cox HR scores the reconstructed
 * HR against that registry-reported HR. Emits a results JSON + a GALLERY.md table of diverse worked
 * examples (different conditions, badges) with full NCT provenance.
 *
 * This is "what does it produce", complementing "is it accurate" (VALIDATION.md). Run from repo root:
 *   node validate/gallery.js   ->   realipd/gallery_results.json  +  GALLERY.md
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');

const COHORT = path.join(__dirname, '..', 'cohort');
const files = fs.existsSync(COHORT) ? fs.readdirSync(COHORT).filter(f => f.endsWith('.json') && f.startsWith('NCT')) : [];
if (!files.length) { console.error('cohort/ not found or empty (gitignored; re-harvest with harvest/harvest_cohort.py)'); process.exit(2); }

const coxHR = (a, b) => Math.exp(RIPD._.coxLogHR(a.map(r => ({ time: r.time, status: r.status, x: 1 }))
  .concat(b.map(r => ({ time: r.time, status: r.status, x: 0 }))).slice()).beta);

const all = [];
for (const f of files) {
  let t; try { t = JSON.parse(fs.readFileSync(path.join(COHORT, f), 'utf8')); } catch { continue; }
  let r; try { r = RIPD.reconstruct(t); } catch { continue; }
  if (r.tier !== 'A' || !r.arms || r.arms.length !== 2) continue;
  const exp = r.arms.find(a => a.role === 'experimental') || r.arms[1], ctl = r.arms.find(a => a.role === 'comparator') || r.arms[0];
  const reconHR = coxHR(exp.ipd, ctl.ipd);
  const ev = exp.ipd.filter(x => x.status === 1).length + ctl.ipd.filter(x => x.status === 1).length;
  const anchors = Math.min(...(t.arms || []).map(a => (a.km_points || []).length));
  const regHR = t.hr && t.hr.value != null ? t.hr.value : null;
  const foldVsReg = regHR ? +Math.exp(Math.abs(Math.log(reconHR) - Math.log(regHR))).toFixed(2) : null;
  all.push({ nct: t.nct_id, condition: (t.condition || '').split(';')[0].trim().slice(0, 46), url: t.source_url,
    badge: r.audit && r.audit.badge, method: r.method, tier: r.tier, anchors,
    n_exp: exp.ipd.length, n_ctl: ctl.ipd.length, events: ev, recon_HR: +reconHR.toFixed(3), registry_HR: regHR, fold_vs_registry: foldVsReg });
}

const withHR = all.filter(a => a.registry_HR != null);
// pick a diverse gallery: trials with a registry HR, good fit, distinct conditions, mix of badges
const seen = new Set();
const gallery = withHR.filter(a => a.fold_vs_registry != null && a.fold_vs_registry < 1.6 && a.events >= 30)
  .sort((a, b) => a.fold_vs_registry - b.fold_vs_registry)
  .filter(a => { const key = a.condition.toLowerCase().slice(0, 12); if (seen.has(key)) return false; seen.add(key); return true; })
  .slice(0, 10);

const med = (xs) => { if (!xs.length) return null; const s = xs.slice().sort((a, b) => a - b); return +s[s.length >> 1].toFixed(3); };
const badgeCount = { gold: 0, silver: 0, bronze: 0, none: 0 };
for (const a of all) badgeCount[a.badge] = (badgeCount[a.badge] || 0) + 1;
const summary = {
  cohort_trials_reconstructed: all.length, with_registry_HR: withHR.length,
  median_fold_vs_registry_HR: med(withHR.map(a => a.fold_vs_registry).filter(Boolean)),
  badge_distribution: badgeCount, method_qp: all.filter(a => a.method === 'qp').length, method_other: all.filter(a => a.method !== 'qp').length,
  note: 'Real AACT trials reconstructed with the shipped engine (Titman-QP default). fold_vs_registry = '
    + 'reconstructed Cox HR vs the registry-reported HR (a coarse held-out truth). Provenance = NCT URL.',
};
fs.writeFileSync(path.join(__dirname, '..', 'realipd', 'gallery_results.json'), JSON.stringify({ summary, gallery, all_with_hr: withHR }, null, 2));

// GALLERY.md
const lines = [];
lines.push('# Production gallery — the method on real ClinicalTrials.gov trials', '');
lines.push('*The reconstruction running on the **real AACT trials it is built for** (not the open-IPD',
  'validation set). Each row: a real trial harvested from the 2026-06-01 AACT snapshot, reconstructed',
  'with the shipped engine (Titman-QP default), its reconstructed Cox HR scored against the',
  '**registry-reported HR** (a coarse held-out truth), with full NCT provenance. Reproduce:',
  '`node validate/gallery.js`. (Numbers from `realipd/gallery_results.json`.)*', '');
lines.push('## Cohort-wide', '');
lines.push('- **' + all.length + '** Tier-A trials reconstructed end-to-end; **' + withHR.length + '** also report a Cox HR.');
lines.push('- Median reconstructed-vs-registry HR fold-error: **' + summary.median_fold_vs_registry_HR + '**.');
lines.push('- Self-audit badges: ' + Object.entries(badgeCount).filter(([, v]) => v).map(([k, v]) => `${v} ${k}`).join(' · ') + '.', '');
lines.push('## Worked examples (diverse conditions, best-fit per condition)', '');
lines.push('| NCT | condition | N exp/ctl | events | anchors | badge | registry HR | reconstructed HR | fold |');
lines.push('|---|---|---|---|---|---|---|---|---|');
for (const g of gallery) lines.push(`| [${g.nct}](${g.url}) | ${g.condition} | ${g.n_exp}/${g.n_ctl} | ${g.events} | ${g.anchors} | ${g.badge} | ${g.registry_HR} | ${g.recon_HR} | ${g.fold_vs_registry} |`);
lines.push('', '*The HR is the hard estimand and the registry HR is itself coarse; read these as a triangulation',
  'check, not bit-exact recovery. RMST/median reconstruct far more tightly (see `VALIDATION.md`).*', '');
fs.writeFileSync(path.join(__dirname, '..', 'GALLERY.md'), lines.join('\n'));

console.log(JSON.stringify(summary, null, 2));
console.log('\ngallery (' + gallery.length + ' diverse examples):');
for (const g of gallery) console.log('  ' + g.nct.padEnd(12) + (g.badge || '').padEnd(7) + 'regHR ' + String(g.registry_HR).padStart(5) + ' recon ' + String(g.recon_HR).padStart(6) + ' fold ' + g.fold_vs_registry + '  ' + g.condition);
console.log('\nwrote realipd/gallery_results.json + GALLERY.md');
