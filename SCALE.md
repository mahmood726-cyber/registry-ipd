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

## Beyond the curve: the Tier-B population

The curve-based Tier-A cohort is **saturated** — 595 harvested trials already exceed the broad-census
count of 514 curve-posting trials, so re-harvesting Tier A does not widen the reconstructable set. The
remaining expansion lives in **Tier B**: trials that post a survival **median + a hazard ratio + arm N**
but no KM curve. Harvesting the full snapshot for this pattern
(`python harvest/harvest_tierb.py`) yields **1,144 Tier-B trials**, every one of which reconstructs via
the engine's parametric (exponential) path (`node validate/tierb_scale.js` → `realipd/tierb_scale.json`).

That roughly **triples** the reconstructable registry (≈400 Tier-A comparison-bearing trials → +1,144
Tier-B), but the fidelity is deliberately lower, and the self-audit says so out loud:

| Tier-B self-audit badge | trials | meaning |
|---|---:|---|
| silver / bronze | 53 | exponential PH happens to reconcile the posted median **and** HR |
| **none** | 1,091 | it cannot — the reconstructed exp median is forced to `cMed/HR`, disagreeing with the independently-posted exp median (the non-exponential common case) |

The **1,091 badge-`none`** split is the honest headline: Tier B is a **coverage/triage tier, not an
HR-recovery tier**. The HR is *imposed*, not recovered; the curve is *assumed* exponential; RMST carries
~7% error and worse on non-exponential survival (`VALIDATION.md` → Tier B). Use Tier-B pseudo-IPD for
scoping and feasibility, never as a fidelity claim — the curve-based Tier-A cohort remains the validated
deliverable.*
