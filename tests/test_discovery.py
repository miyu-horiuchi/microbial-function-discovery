import json
import unittest

from microbial_function_discovery.datasets import parse_labels_tsv
from microbial_function_discovery.discovery import annotate_discovery_candidates, export_safe_leads
from microbial_function_discovery.features import build_feature_matrix_from_annotation_tsv
from microbial_function_discovery.learning import train_baseline


ANNOTATIONS_TSV = """genome_id\tprotein_id\tdatabase\taccession\tname\tevalue
G1\tp1\tCAZy\tGH5\tcellulase\t1e-20
G2\tp1\tPfam\tPF00002\tnegative marker\t1e-20
G3\tp1\tCAZy\tGH5\tcandidate cellulase\t1e-20
"""

LABELS_TSV = """genome_id\tsplit\tfamily\tpanel\tlabel\tvalue
G1\ttrain\tFamilyA\tbiofuels_industrial\tcellulose_degradation\t1
G2\ttrain\tFamilyA\tbiofuels_industrial\tcellulose_degradation\t0
G1\ttrain\tFamilyA\tbiosafety\tpathogenicity_human\t0
G2\ttrain\tFamilyA\tbiosafety\tpathogenicity_human\t1
G3\ttest\tFamilyB\tbiofuels_industrial\tcellulose_degradation\t1
G3\ttest\tFamilyB\tbiosafety\tpathogenicity_human\t0
"""


class DiscoveryAnnotationTests(unittest.TestCase):
    def test_annotate_discovery_candidates_adds_taxonomy_evidence_and_safety(self):
        features = build_feature_matrix_from_annotation_tsv(ANNOTATIONS_TSV)
        labels = parse_labels_tsv(LABELS_TSV)
        model = train_baseline(features, labels)
        candidates = {
            "split": "test",
            "targets": [
                {
                    "target_key": "biofuels_industrial:cellulose_degradation",
                    "panel": "biofuels_industrial",
                    "label": "cellulose_degradation",
                    "source": "annotation",
                    "metrics": {"precision_at_10": 1.0},
                    "candidates": [
                        {
                            "genome_id": "G3",
                            "rank": 1,
                            "score": 0.9,
                            "truth": 1,
                        }
                    ],
                }
            ],
        }
        traits = [
            {
                "bacdive_id": "G3",
                "domain": "Bacteria",
                "family": "Candidateaceae",
                "genus": "Candidateus",
                "species": "Candidateus utilis",
                "biosafety_level": "BSL-1",
                "pathogenicity_human": False,
                "pathogenicity_animal": False,
            }
        ]
        accessions = [
            {
                "bacdive_id": "G3",
                "accession": "GCA_000000003",
                "assembly_level": "complete",
                "description": "Candidateus utilis assembly",
            }
        ]

        annotated = annotate_discovery_candidates(
            candidates,
            traits,
            accessions,
            model,
            features,
            labels,
            evidence_limit=3,
        )

        candidate = annotated["targets"][0]["candidates"][0]
        self.assertEqual(candidate["taxonomy"]["species"], "Candidateus utilis")
        self.assertEqual(candidate["genome_accession"]["accession"], "GCA_000000003")
        self.assertEqual(candidate["biosafety"]["known_biosafety_level"], "BSL-1")
        self.assertEqual(candidate["biosafety"]["risk_level"], "low")
        self.assertGreater(candidate["biosafety"]["predicted_pathogenicity_human"], 0.0)
        self.assertEqual(candidate["evidence"][0]["id"], "CAZy:GH5")
        json.dumps(annotated)

    def test_export_safe_leads_filters_for_validation_ready_candidates(self):
        annotated = {
            "targets": [
                {
                    "target_key": "biofuels_industrial:cellulose_degradation",
                    "panel": "biofuels_industrial",
                    "label": "cellulose_degradation",
                    "source": "annotation",
                    "metrics": {"precision_at_10": 0.9},
                    "candidates": [
                        {
                            "genome_id": "G3",
                            "rank": 1,
                            "score": 0.91,
                            "taxonomy": {"species": "Candidateus utilis", "genus": "Candidateus"},
                            "genome_accession": {"accession": "GCA_000000003"},
                            "biosafety": {"risk_level": "low", "flags": []},
                            "evidence": [{"id": "CAZy:GH5", "type": "annotation_feature"}],
                        },
                        {
                            "genome_id": "G4",
                            "rank": 2,
                            "score": 0.89,
                            "taxonomy": {"species": "Riskus pathogenus"},
                            "genome_accession": {"accession": "GCA_000000004"},
                            "biosafety": {"risk_level": "high", "flags": ["known_human_pathogen"]},
                            "evidence": [{"id": "CAZy:GH5", "type": "annotation_feature"}],
                        },
                    ],
                }
            ]
        }

        export = export_safe_leads(
            annotated,
            allowed_risks=["low", "unknown"],
            min_precision=0.8,
            precision_k=10,
            require_accession=True,
            require_evidence=True,
            max_leads=10,
        )

        self.assertEqual(export["n_leads"], 1)
        lead = export["leads"][0]
        self.assertEqual(lead["genome_id"], "G3")
        self.assertEqual(lead["species"], "Candidateus utilis")
        self.assertEqual(lead["accession"], "GCA_000000003")
        self.assertEqual(lead["evidence"], "CAZy:GH5")
        self.assertEqual(export["panels"]["biofuels_industrial"]["n_leads"], 1)
        json.dumps(export)


    def test_annotate_attaches_novelty_flag(self):
        features = build_feature_matrix_from_annotation_tsv(ANNOTATIONS_TSV)
        labels = parse_labels_tsv(LABELS_TSV)
        model = train_baseline(features, labels)
        candidates = {
            "split": "test",
            "targets": [
                {
                    "target_key": "biofuels_industrial:cellulose_degradation",
                    "panel": "biofuels_industrial",
                    "label": "cellulose_degradation",
                    "source": "annotation",
                    "metrics": {"precision_at_10": 1.0},
                    "candidates": [{"genome_id": "G3", "rank": 1, "score": 0.9}],
                }
            ],
        }
        annotated = annotate_discovery_candidates(candidates, [], [], model, features, labels)
        cand = annotated["targets"][0]["candidates"][0]
        self.assertIn("novelty", cand)
        self.assertIn(cand["novelty"]["level"], {"typical", "novel", "highly_novel", "unknown"})
        # G3's annotation set equals training genome G1's -> nearest distance 0 -> typical
        self.assertEqual(cand["novelty"]["level"], "typical")

    def test_annotate_novelty_unknown_for_missing_genome(self):
        features = build_feature_matrix_from_annotation_tsv(ANNOTATIONS_TSV)
        labels = parse_labels_tsv(LABELS_TSV)
        model = train_baseline(features, labels)
        candidates = {
            "split": "test",
            "targets": [
                {
                    "target_key": "biofuels_industrial:cellulose_degradation",
                    "panel": "biofuels_industrial",
                    "label": "cellulose_degradation",
                    "source": "annotation",
                    "metrics": {"precision_at_10": 1.0},
                    "candidates": [{"genome_id": "NOT_IN_MATRIX", "rank": 1, "score": 0.9}],
                }
            ],
        }
        annotated = annotate_discovery_candidates(candidates, [], [], model, features, labels)
        cand = annotated["targets"][0]["candidates"][0]
        self.assertEqual(cand["novelty"]["level"], "unknown")
        self.assertIsNone(cand["novelty"]["score"])


    def _annotated_for_reports(self):
        features = build_feature_matrix_from_annotation_tsv(ANNOTATIONS_TSV)
        labels = parse_labels_tsv(LABELS_TSV)
        model = train_baseline(features, labels)
        candidates = {
            "split": "test",
            "targets": [
                {
                    "target_key": "biofuels_industrial:cellulose_degradation",
                    "panel": "biofuels_industrial",
                    "label": "cellulose_degradation",
                    "source": "annotation",
                    "metrics": {"precision_at_10": 1.0},
                    "candidates": [{"genome_id": "G3", "rank": 1, "score": 0.9}],
                }
            ],
        }
        return annotate_discovery_candidates(candidates, [], [], model, features, labels)

    def test_candidate_report_has_novelty_column(self):
        from microbial_function_discovery.discovery import render_discovery_candidate_report
        report = render_discovery_candidate_report(self._annotated_for_reports())
        self.assertIn("Novelty", report)
        self.assertIn("representativeness", report.lower())

    def test_safe_leads_carry_novelty_level(self):
        from microbial_function_discovery.discovery import (
            export_safe_leads,
            render_safe_leads_report,
            SAFE_LEAD_FIELDS,
        )
        export = export_safe_leads(self._annotated_for_reports(), require_accession=False, require_evidence=False)
        self.assertIn("novelty_level", SAFE_LEAD_FIELDS)
        self.assertTrue(export["leads"], "expected at least one lead")
        self.assertIn("novelty_level", export["leads"][0])
        report = render_safe_leads_report(export)
        self.assertIn("Novelty", report)
        self.assertIn("representativeness", report.lower())


if __name__ == "__main__":
    unittest.main()
