#!/usr/bin/env node
/*
 * Flexible cBioPortal cohort fetcher for NON-TCGA published cohorts (open API, no DUA).
 *
 * Generalises fetch_cbioportal.js to arbitrary survival + grouping fields, so real published cohorts
 * (METABRIC, MSKCC immunotherapy cohorts, etc.) can join the gold standard. Writes a CSV with
 * time / status (1=event) / group / sex columns.
 *
 * Run: node harvest/fetch_cbio_cohort.js <studyId> <survPrefix> <groupAttr>
 *   survPrefix: OS | PFS | DFS  (uses <prefix>_MONTHS + <prefix>_STATUS)
 *   groupAttr : a PATIENT clinical attribute id (e.g. ER_STATUS, GRADE, TUMOR_STAGE, TREATMENT)
 * Writes realipd/cbio2_<studyShort>_<group>.csv and prints the group distribution [N,events].
 */
const fs = require('fs');
const path = require('path');
const API = 'https://www.cbioportal.org/api';
const OUTDIR = path.join(__dirname, '..', 'realipd');

async function fetchClinical(study, type) {
  const url = `${API}/studies/${study}/clinical-data?clinicalDataType=${type}&projection=SUMMARY&pageSize=200000&pageNumber=0`;
  const res = await fetch(url, { headers: { accept: 'application/json' } });
  if (!res.ok) throw new Error(`HTTP ${res.status} for ${study} (${type})`);
  return res.json();
}
const csvCell = (x) => (x == null ? '' : String(x).replace(/[",\r\n]/g, ' '));

(async () => {
  const [study, surv = 'OS', group] = process.argv.slice(2);
  if (!study || !group) { console.error('usage: node harvest/fetch_cbio_cohort.js <studyId> <OS|PFS|DFS> <groupAttr>'); process.exit(2); }
  const tCol = `${surv}_MONTHS`, sCol = `${surv}_STATUS`;
  const want = new Set([tCol, sCol, group, 'SEX']);
  const [pat, samp] = await Promise.all([fetchClinical(study, 'PATIENT'), fetchClinical(study, 'SAMPLE').catch(() => [])]);
  const byPt = new Map();
  for (const r of pat) { if (!want.has(r.clinicalAttributeId)) continue; if (!byPt.has(r.patientId)) byPt.set(r.patientId, {}); byPt.get(r.patientId)[r.clinicalAttributeId] = r.value; }
  for (const r of samp) { if (r.clinicalAttributeId !== group) continue; const p = byPt.get(r.patientId); if (p && p[group] == null) p[group] = r.value; }
  const rows = [];
  for (const [pid, a] of byPt) {
    const t = parseFloat(a[tCol]); const st = a[sCol];
    if (!isFinite(t) || t <= 0 || !st) continue;
    rows.push({ patientId: pid, time: t, status: /^1/.test(st) ? 1 : 0, group: a[group] || '', sex: a.SEX || '' });
  }
  const head = ['patientId', 'time', 'status', 'group', 'sex'];
  const short = study.split('_')[0] + (study.includes('metabric') ? 'METABRIC' : '');
  const file = path.join(OUTDIR, `cbio2_${short}_${group.toLowerCase()}.csv`);
  fs.writeFileSync(file, [head.join(',')].concat(rows.map(r => head.map(h => csvCell(r[h])).join(','))).join('\n'));
  const dist = {}; for (const r of rows) { const g = r.group || '(blank)'; dist[g] = dist[g] || [0, 0]; dist[g][0]++; if (r.status === 1) dist[g][1]++; }
  console.log(`${study}  ${surv}  by ${group}  (${rows.length} patients)  -> ${path.basename(file)}`);
  console.log('  group -> [N, events]:', JSON.stringify(dist));
})();
