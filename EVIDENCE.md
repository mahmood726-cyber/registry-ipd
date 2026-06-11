# Evidence summary

*One-line headline for every validation in this repo, auto-generated from the committed result JSONs by
`node validate/evidence_summary.js`. Full detail and caveats in `VALIDATION.md`. Each row is reproducible
via the command shown. Engine: `src/engine.js` (Titman-QP default when a total-event count is posted).*

| Validation | Headline result | Reproduce |
|---|---|---|
| **Gold standard (true IPD)** | 51 datasets; ≥100/arm curve-only median fold 1.15 (17/29), Titman-QP **1.03** (28/29) | `node validate/goldstandard.js` |
| **Calibrated uncertainty** | 95% credible interval covers true HR 28/29 (97%) | `node validate/goldstandard_uncertainty.js` |
| **Uncertainty calibration** | nominal 50/80/90/95% → empirical 82.8/89.7/96.6/96.6% (conservative at narrow levels) | `node validate/uncertainty_calibration.js` |
| **Titman QP vs anchor-exact** | ≥100/arm: QP 1.044 (27/29) vs anchor-exact 1.148 (19/29) | `node validate/titman_qp.js` |
| **12-method benchmark** | best qp_roughness 1.034 … QP-L2 1.043 (near-optimal); classical Guyot/anchor-exact ~1.14 | `node validate/method_zoo.js` |
| **NAR fusion (curve + figure NAR)** | ≥100/arm: curve-only 1.153 → fusion **1.03** ≈ QP+events 1.035 | `node validate/nar_fusion.js` |
| **Advanced curve-only estimators** | heavily-censored TCGA: censor-to-tail 1.558, max-ent 1.551, Wasserstein-bary 1.51 (all fail ⇒ identifiability limit) | `node validate/advanced_estimators.js` |
| **Head-to-head vs digitization** | registry-exact vs digitized: equal-density near-tie, dense-digitized K=25 wins (anchor density > pixel noise) | `node validate/headtohead.js` |
| **Competing risks (Aalen–Johansen)** | reconstructed AJ CIF within ~1–6pp of truth; naive 1−KM overstates (aidssi +16pp) | `node validate/goldstandard_cr.js` |
| **Tier B (median+HR, no curve)** | exponential RMST ~6.6% median; fails 40–58% on non-exponential survival (Weibull fixes the tail) | `node validate/tierb_validation.js` |
| **IPD meta-analysis fidelity** | pooled HR true 2.497 vs reconstructed 2.643 (within 1.059); τ² 0.15→0.204 (conservative) | `node validate/ipd_meta_fidelity.js` |
| **AACT coverage census** | 76067 results-trials; 0 structured NAR rows; 288–514 post a reconstructable curve; only **77–112** post curve **+** HR (validation-grade, 21.8% of curve-posters) | `python harvest/census_full_aact.py` |
| **Production (real AACT trials)** | 250 trials reconstructed; median fold vs registry HR 1.13; only 27% post an event count | `node validate/gallery.js` |
| **Expanded held-out HR validation** | 57 validation-grade trials scored vs registry HR (was 30): curve-sourced 25 median fold 1.119, same-endpoint sibling 32 median fold 1.19 (endpoint-aware harvester HR recovery; OS-vs-PFS mismatches dropped) | `python harvest/backfill_validation_hr.py && node validate/gallery_expanded.js` |
| **Scale (registry-wide)** | 399 trials → 904 pairwise pseudo-IPD comparisons across 277 conditions; 886 exportable | `node validate/scale_run.js` |
| **Tier-B scale (median+HR, no curve)** | 1144 trials reconstruct (≈3× the Tier-A set); 53 self-audit silver/bronze, 1091 badge "none" ⇒ coverage tier, not fidelity | `python harvest/harvest_tierb.py && node validate/tierb_scale.js` |
| **When does the event count matter** | QP gap over curve-only mean ~0.186 fold but not predictable (Spearman vs censoring -0.017) ⇒ always prefer the event count | `node validate/censoring_stratified.js` |

**Tests:** `npm test` (engine + metric unit tests). **Smoke:** `python test/smoke_browser.py`.
Datasets under `realipd/` (open) and `cohort/` (AACT harvest) are gitignored; re-download from the script headers.
