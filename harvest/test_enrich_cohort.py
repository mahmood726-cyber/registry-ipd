"""Unit test for the cohort batch enrichment (harvest/enrich_cohort.py::enrich_rows)."""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import enrich_cohort as EC  # noqa: E402


def test_enrich_rows_promotes_confident_hr_and_attaches_median():
    rows = [
        {"nct": "NCT_A", "pmid": "1", "curve_endpoint": None},   # confident HR (AACT none) -> promoted
        {"nct": "NCT_B", "pmid": "2", "curve_endpoint": None},   # AACT already has HR -> cross-check only
        {"nct": "NCT_C", "pmid": "3", "curve_endpoint": None},   # nothing usable
        {"nct": "NCT_X", "pmid": "9", "curve_endpoint": None},   # no abstract -> skipped
    ]
    abstracts = {
        "1": ("Treatment reduced the risk of death (hazard ratio 0.62, 95% CI 0.50-0.77). "
              "Median overall survival was 14.0 months versus 9.0 months."),
        "2": "Treatment reduced death (hazard ratio 0.62, 95% CI 0.50-0.77).",
        "3": "The trial enrolled adults at 12 centres.",
    }
    trials = {
        "NCT_A": {"nct_id": "NCT_A", "arms": [{"N": 100}, {"N": 100}]},
        "NCT_B": {"nct_id": "NCT_B", "hr": {"value": 0.7, "method": "AACT"}, "arms": [{"N": 100}]},
        "NCT_C": {"nct_id": "NCT_C", "arms": [{"N": 50}]},
    }
    manifest, enriched = EC.enrich_rows(rows, abstracts, trials)
    assert manifest["considered"] == 3                 # NCT_X dropped (no abstract)
    assert manifest["hr_promoted"] == 1                # NCT_A
    assert manifest["hr_crosscheck"] == 1              # NCT_B (AACT already had an HR)
    assert manifest["median_attached"] == 1            # NCT_A
    assert manifest["enriched"] == 2                   # A (hr+median) and B (median? no) -> A and B
    # NCT_A got the HR promoted into the trial; NCT_B kept its AACT HR
    assert trials["NCT_A"]["hr"]["value"] == 0.62 and trials["NCT_A"]["hr"]["source"] == "pubmed_abstract"
    assert trials["NCT_B"]["hr"]["method"] == "AACT" and trials["NCT_B"]["hr_abstract"]["value"] == 0.62
    assert set(enriched) == {"NCT_A", "NCT_B"}
