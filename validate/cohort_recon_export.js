#!/usr/bin/env node
/*
 * Export reconstructed values for EVERY 2-arm Tier-A cohort trial -> realipd/cohort_recon.json.
 *
 * Step 1 of the full-cohort independent validation (extends the validation-grade PubMed check to all
 * reconstructed trials, including the ~220 where AACT posts no HR -- there the published HR is the ONLY
 * held-out truth). Emits per trial: reconstructed Cox HR, reconstructed per-arm medians, the registry HR
 * if any, condition, badge, anchors. harvest/cohort_pubmed.py then attaches PMIDs + curve endpoints,
 * fetches abstracts, and compares against the published HR + median.
 *
 * Run from repo root: node validate/cohort_recon_export.js
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');

const COHORT = path.join(__dirname, '..', 'cohort');
if (!fs.existsSync(COHORT)) { console.error('cohort/ missing (gitignored; re-harvest)'); process.exit(2); }
const files = fs.readdirSync(COHORT).filter(f => f.endsWith('.json') && f.startsWith('NCT'));

const coxHR = (a, b) => Math.exp(RIPD._.coxLogHR(a.map(r => ({ time: r.time, status: r.status, x: 1 }))
  .concat(b.map(r => ({ time: r.time, status: r.status, x: 0 }))).slice()).beta);
const armMedian = (ipd) => RIPD._.medianFromKM(RIPD._.kmFromIPD(ipd));

const out = [];
for (const f of files) {
  let t; try { t = JSON.parse(fs.readFileSync(path.join(COHORT, f), 'utf8')); } catch { continue; }
  if ((t.time_unit || 'months') !== 'months') continue;        // medians compared in months
  let r; try { r = RIPD.reconstruct(t); } catch { continue; }
  if (r.tier !== 'A' || !r.arms || r.arms.length !== 2) continue;
  const exp = r.arms.find(a => a.role === 'experimental') || r.arms[1];
  const ctl = r.arms.find(a => a.role === 'comparator') || r.arms[0];
  const meds = r.arms.map(a => armMedian(a.ipd)).filter(x => Number.isFinite(x) && x > 0);
  const ev = exp.ipd.filter(x => x.status === 1).length + ctl.ipd.filter(x => x.status === 1).length;
  out.push({
    nct: t.nct_id, condition: (t.condition || '').split(';')[0].trim().slice(0, 46),
    badge: r.audit && r.audit.badge, anchors: Math.min(...(t.arms || []).map(a => (a.km_points || []).length)),
    recon_HR: +coxHR(exp.ipd, ctl.ipd).toFixed(3),
    recon_medians: meds.length === 2 ? meds.map(x => +x.toFixed(2)) : null,
    registry_HR: (t.hr && t.hr.value != null) ? t.hr.value : null, events: ev,
  });
}
fs.writeFileSync(path.join(__dirname, '..', 'realipd', 'cohort_recon.json'), JSON.stringify(out, null, 1));
console.log(`exported ${out.length} reconstructed 2-arm Tier-A trials -> realipd/cohort_recon.json`);
console.log(`  with reconstructed medians: ${out.filter(r => r.recon_medians).length}; with a registry HR: ${out.filter(r => r.registry_HR != null).length}`);
