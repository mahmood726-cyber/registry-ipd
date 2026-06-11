#!/usr/bin/env node
/*
 * INDEPENDENT MEDIAN VALIDATION — reconstructed per-arm median survival vs the PUBLISHED median (PubMed).
 *
 * The HR is the reconstruction's hardest estimand; the MEDIAN is its tightest (~3% on the open gold
 * standard). This checks that tightest claim against an INDEPENDENT source: the per-arm medians stated in
 * the trial's own abstract (harvest/pubmed_medians.json). For each validation-grade trial with a clean
 * two-arm published pair, we reconstruct the per-arm medians and compare by SORTED magnitude (longer-to-
 * longer, shorter-to-shorter) so arm labelling never matters.
 *
 * Run from repo root: node validate/pubmed_median_validation.js  ->  realipd/pubmed_median_validation.json
 */
const fs = require('fs');
const path = require('path');
const RIPD = require('../src/engine.js');

const COHORT = path.join(__dirname, '..', 'cohort');
const RIPDdir = path.join(__dirname, '..', 'realipd');
const pub = JSON.parse(fs.readFileSync(path.join(RIPDdir, 'pubmed_medians.json'), 'utf8'));
if (!fs.existsSync(COHORT)) { console.error('cohort/ missing (gitignored; re-harvest)'); process.exit(2); }

const armMedian = (ipd) => RIPD._.medianFromKM(RIPD._.kmFromIPD(ipd));
const fold = (a, b) => Math.exp(Math.abs(Math.log(a) - Math.log(b)));
const med = (xs) => { if (!xs.length) return null; const s = xs.slice().sort((a, b) => a - b); return +s[s.length >> 1].toFixed(3); };

const rows = [];
const armFolds = [];
let nTrials = 0, nSkipUnit = 0, nSkipRecon = 0;
for (const nct in pub) {
  const p = pub[nct];
  if (!(p.medians && p.medians.length === 2 && p.n_numbers === 2)) continue;   // clean two-arm pairs only
  const fp = path.join(COHORT, `${nct}.json`);
  if (!fs.existsSync(fp)) continue;
  let t; try { t = JSON.parse(fs.readFileSync(fp, 'utf8')); } catch { continue; }
  if ((t.time_unit || 'months') !== 'months') { nSkipUnit++; continue; }        // medians are in months
  let r; try { r = RIPD.reconstruct(t); } catch { continue; }
  if (r.tier !== 'A' || !r.arms || r.arms.length !== 2) { nSkipRecon++; continue; }
  const recon = r.arms.map(a => armMedian(a.ipd)).filter(x => Number.isFinite(x) && x > 0);
  if (recon.length !== 2) { nSkipRecon++; continue; }                            // a median was undefined (heavy censoring)
  const ps = p.medians.slice().sort((a, b) => a - b);
  const rs = recon.slice().sort((a, b) => a - b);
  const f = [fold(rs[0], ps[0]), fold(rs[1], ps[1])];
  armFolds.push(...f);
  nTrials++;
  rows.push({ nct, pmid: p.pmid, published_medians: ps.map(x => +x.toFixed(1)),
    reconstructed_medians: rs.map(x => +x.toFixed(1)), arm_folds: f.map(x => +x.toFixed(3)) });
}

const within = (thr) => `${armFolds.filter(f => f < thr).length}/${armFolds.length}`;
const summary = {
  trials_with_clean_published_pair: Object.values(pub).filter(p => p.medians.length === 2 && p.n_numbers === 2).length,
  trials_validated: nTrials,
  arm_medians_compared: armFolds.length,
  median_arm_fold: med(armFolds),
  arm_medians_within_10pct: within(1.1),
  arm_medians_within_20pct: within(1.2),
  skipped_nonmonth_unit: nSkipUnit, skipped_reconstruction: nSkipRecon,
  note: 'Reconstructed per-arm median survival vs the trial’s PUBLISHED per-arm medians (PubMed '
    + 'abstract), matched by SORTED magnitude so arm labelling is irrelevant. Independent of the registry '
    + 'and of the HR; validates the reconstruction’s tightest quantity. fold = exp|log(recon/pub)|.',
};
fs.writeFileSync(path.join(RIPDdir, 'pubmed_median_validation.json'), JSON.stringify({ summary, rows }, null, 2));
console.log(JSON.stringify(summary, null, 2));
console.log(`\nIndependent median check: ${armFolds.length} arm-medians across ${nTrials} trials, median fold ${summary.median_arm_fold} (within 10%: ${summary.arm_medians_within_10pct})`);
console.log('wrote realipd/pubmed_median_validation.json');
