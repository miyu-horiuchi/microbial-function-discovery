import json
import unittest

from microbial_function_discovery.datasets import parse_labels_tsv
from microbial_function_discovery.discovery import annotate_discovery_candidates
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


if __name__ == "__main__":
    unittest.main()
