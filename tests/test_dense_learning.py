import json
import unittest

from microbial_function_discovery.datasets import parse_labels_tsv
from microbial_function_discovery.dense_learning import (
    DenseFeatureMatrix,
    evaluate_dense_model,
    evaluate_dense_ranking,
    train_dense_baseline,
)


LABELS_TSV = """genome_id\tsplit\tfamily\tpanel\tlabel\tvalue
G1\ttrain\tFamilyA\tbiofuels_industrial\tthermophile\t1
G2\ttrain\tFamilyA\tbiofuels_industrial\tthermophile\t0
G3\ttest\tFamilyB\tbiofuels_industrial\tthermophile\t1
G4\ttest\tFamilyB\tbiofuels_industrial\tthermophile\t0
"""


class DenseLearningTests(unittest.TestCase):
    def test_train_dense_baseline_learns_float_signal(self):
        features = DenseFeatureMatrix(
            genome_ids=["G1", "G2", "G3", "G4"],
            feature_indexes=[0, 1],
            rows=[
                [2.0, 0.0],
                [-2.0, 0.0],
                [1.5, 0.1],
                [-1.5, -0.1],
            ],
        )
        labels = parse_labels_tsv(LABELS_TSV)

        model = train_dense_baseline(features, labels)

        self.assertGreater(
            model.predict_proba("biofuels_industrial:thermophile", "G3", features),
            model.predict_proba("biofuels_industrial:thermophile", "G4", features),
        )

    def test_dense_model_round_trips_json(self):
        features = DenseFeatureMatrix(
            genome_ids=["G1", "G2"],
            feature_indexes=[0],
            rows=[[1.0], [-1.0]],
        )
        labels = parse_labels_tsv(
            """genome_id\tsplit\tfamily\tpanel\tlabel\tvalue
G1\ttrain\tFamilyA\tbiofuels_industrial\tthermophile\t1
G2\ttrain\tFamilyA\tbiofuels_industrial\tthermophile\t0
"""
        )
        model = train_dense_baseline(features, labels)

        self.assertEqual(model, model.from_json(model.to_json()))

    def test_evaluate_dense_model_reports_accuracy(self):
        features = DenseFeatureMatrix(
            genome_ids=["G1", "G2", "G3", "G4"],
            feature_indexes=[0],
            rows=[[2.0], [-2.0], [1.5], [-1.5]],
        )
        labels = parse_labels_tsv(LABELS_TSV)
        model = train_dense_baseline(features, labels)

        report = evaluate_dense_model(model, features, labels, split="test")

        self.assertEqual(report["overall"]["n"], 2)
        self.assertEqual(report["overall"]["accuracy"], 1.0)
        json.dumps(report)

    def test_evaluate_dense_ranking_reports_precision_at_k(self):
        features = DenseFeatureMatrix(
            genome_ids=["G1", "G2", "G3", "G4"],
            feature_indexes=[0],
            rows=[[2.0], [-2.0], [1.5], [-1.5]],
        )
        labels = parse_labels_tsv(LABELS_TSV)
        model = train_dense_baseline(features, labels)

        report = evaluate_dense_ranking(
            model,
            features,
            labels,
            split="test",
            target_key="biofuels_industrial:thermophile",
            ks=[1, 2],
        )

        self.assertEqual(report["metrics"]["precision_at_1"], 1.0)
        self.assertEqual(report["ranked_candidates"][0]["genome_id"], "G3")


if __name__ == "__main__":
    unittest.main()
