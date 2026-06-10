#!/usr/bin/env node
/*
 * Fetch real cancer-survival IPD from the public cBioPortal API (TCGA + published studies, open).
 *
 * For each study it pulls patient-level clinical data (overall survival + groupers), pivots to one
 * row per patient, derives clean binary arm columns, and writes a CSV the gold-standard harness can
 * read directly:  realipd/cbio_<study>.csv  with columns:
 *   patientId, time (OS months), status (1=deceased), sex, stage_group (early|late), subtype
 *
 * These are TRUE patient-level data used as a gold-standard target (the engine never sees them; we
 * build the registry-style coarse summary and reconstruct). Prognostic 2-arm contrasts (stage,
 * sex) — like the existing mgus2/aidssi-by-group datasets — not treatment RCTs.
 *
 * Open data, no credentials. Run: node harvest/fetch_cbioportal.js [study1 study2 ...]
 */
const fs = require('fs');
const path = require('path');

const API = 'https://www.cbioportal.org/api';
const OUTDIR = path.join(__dirname, '..', 'realipd');
const DEFAULT_STUDIES = [
  'luad_tcga_pan_can_atlas_2018', 'coadread_tcga_pan_can_atlas_2018',
  'stad_tcga_pan_can_atlas_2018', 'lihc_tcga_pan_can_atlas_2018',
  'kirc_tcga_pan_can_atlas_2018', 'hnsc_tcga_pan_can_atlas_2018',
  'skcm_tcga_pan_can_atlas_2018', 'blca_tcga_pan_can_atlas_2018',
];

const ATTRS = new Set(['OS_MONTHS', 'OS_STATUS', 'SEX', 'AJCC_PATHOLOGIC_TUMOR_STAGE', 'SUBTYPE']);

// map AJCC stage strings to early (I/II) vs late (III/IV). Order matters: IV, then III, then II/I.
// NOTE: do NOT strip non-letters first — that removes the word "STAGE" and the match silently fails
// (the original bug that made stage_group come back all-blank).
function stageGroup(v) {
  if (!v) return '';
  const s = v.toUpperCase();
  if (/STAGE\s*IV/.test(s)) return 'late';
  if (/STAGE\s*III/.test(s)) return 'late';
  if (/STAGE\s*II/.test(s) || /STAGE\s*I\b/.test(s) || /STAGE\s*I[^IVX]/.test(s)) return 'early';
  return '';
}
const csvCell = (x) => (x == null ? '' : String(x).replace(/[",\r\n]/g, ' '));

async function fetchClinical(study, type) {
  const url = `${API}/studies/${study}/clinical-data?clinicalDataType=${type}&projection=SUMMARY&pageSize=100000&pageNumber=0`;
  const res = await fetch(url, { headers: { accept: 'application/json' } });
  if (!res.ok) throw new Error(`HTTP ${res.status} for ${study} (${type})`);
  return res.json();
}

async function fetchStudy(study) {
  // OS + sex are PATIENT-level; AJCC stage is usually SAMPLE-level in pan-can studies.
  const [patRows, sampRows] = await Promise.all([
    fetchClinical(study, 'PATIENT'), fetchClinical(study, 'SAMPLE')]);
  const byPt = new Map();
  for (const r of patRows) {
    if (!ATTRS.has(r.clinicalAttributeId)) continue;
    if (!byPt.has(r.patientId)) byPt.set(r.patientId, {});
    byPt.get(r.patientId)[r.clinicalAttributeId] = r.value;
  }
  // merge sample-level stage onto the patient (first non-blank wins)
  for (const r of sampRows) {
    if (r.clinicalAttributeId !== 'AJCC_PATHOLOGIC_TUMOR_STAGE') continue;
    const p = byPt.get(r.patientId); if (p && !p.AJCC_PATHOLOGIC_TUMOR_STAGE) p.AJCC_PATHOLOGIC_TUMOR_STAGE = r.value;
  }
  const out = [];
  for (const [pid, a] of byPt) {
    const t = parseFloat(a.OS_MONTHS);
    const st = a.OS_STATUS; // "1:DECEASED" / "0:LIVING"
    if (!isFinite(t) || t <= 0 || !st) continue;
    const status = /^1/.test(st) ? 1 : 0;
    out.push({ patientId: pid, time: t, status, sex: a.SEX || '',
      stage_group: stageGroup(a.AJCC_PATHOLOGIC_TUMOR_STAGE), subtype: a.SUBTYPE || '' });
  }
  return out;
}

function writeCsv(study, rows) {
  const head = ['patientId', 'time', 'status', 'sex', 'stage_group', 'subtype'];
  const lines = [head.join(',')];
  for (const r of rows) lines.push(head.map(h => csvCell(r[h])).join(','));
  const file = path.join(OUTDIR, `cbio_${study.split('_')[0]}.csv`);
  fs.writeFileSync(file, lines.join('\n'));
  return file;
}

function splits(rows, col) {
  const m = {};
  for (const r of rows) { const k = r[col] || '(blank)'; m[k] = m[k] || [0, 0]; m[k][0]++; if (r.status === 1) m[k][1]++; }
  return m;
}

(async () => {
  const studies = process.argv.slice(2).length ? process.argv.slice(2) : DEFAULT_STUDIES;
  fs.mkdirSync(OUTDIR, { recursive: true });
  for (const study of studies) {
    try {
      const rows = await fetchStudy(study);
      const file = writeCsv(study, rows);
      const short = study.split('_')[0];
      console.log(`\n${short}  (${rows.length} patients)  -> ${path.basename(file)}`);
      console.log('  stage_group:', JSON.stringify(splits(rows, 'stage_group')));
      console.log('  sex        :', JSON.stringify(splits(rows, 'sex')));
    } catch (e) {
      console.error(`SKIP ${study}: ${e.message}`);
    }
  }
  console.log('\n(format: group -> [N, events])');
})();
