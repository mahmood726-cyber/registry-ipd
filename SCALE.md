# Scale — registry-wide reconstruction across ClinicalTrials.gov

*The validated method run at scale: every pairwise arm comparison across the harvested AACT
Tier-A cohort (Titman-QP default, `validate/scale_run.js`). Metadata index at `realipd/scale_index.json`;
the pseudo-IPD itself is reproducible per trial (no patient rows committed). This quantifies what
ClinicalTrials.gov can yield as reconstructed survival evidence today.*

## At a glance

- **595** cohort trials scanned · **193** single-arm (no comparison) · **146** multi-arm.
- **399** trials yield ≥1 reconstructable comparison → **904** pairwise pseudo-IPD comparisons.
- Self-audit badges: 104 gold · 782 silver · 18 none · **886 exportable**.
- Reconstruction method: 124 qp · 12 anchor-exact · 768 guyot (QP fires only where a total-event count is posted).
- **277** distinct conditions — survival evidence well beyond oncology.

## Most-represented conditions

| condition | reconstructable trials |
|---|---|
| HIV Infections | 14 |
| Breast Cancer | 13 |
| Asthma | 10 |
| Rheumatoid Arthritis | 7 |
| Multiple Sclerosis | 6 |
| Melanoma | 6 |
| Multiple Myeloma | 5 |
| Atopic Dermatitis | 5 |
| Atrial Fibrillation | 4 |
| HIV Infection | 4 |
| Pulmonary Arterial Hypertension | 4 |
| Breast Neoplasms | 4 |
| Psoriasis | 4 |
| Hepatocellular Carcinoma | 4 |
| Ulcerative Colitis | 4 |

*The multi-arm trials are where scale compounds: a 4-arm trial is up to 6 pairwise comparisons.
The binding limit remains coverage (only a fraction of AACT posts a curve; see the census) and the
posted event count (only a minority enables the QP; see the production gallery).*
