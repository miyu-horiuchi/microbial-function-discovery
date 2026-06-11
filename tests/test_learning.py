import json
import unittest

from microbial_function_discovery.datasets import parse_labels_tsv
from microbial_function_discovery.features import build_feature_matrix_from_annotation_tsv
from microbial_function_discovery.learning import evaluate_model, predict_model, train_baseline
from microbial_function_discovery.validation import validate_prediction


ANNOTATIONS_TSV = """genome_id\tprotein_id\tdatabase\taccession\tname\tevalue
G1\tp1\tCAZy\tGH5\tglycoside hydrolase family 5 cellulase\t1e-40
G2\tp1\tPfam\tPF00001\thypothetical protein\t1e-5
G3\tp1\tCAZy\tGH5\tglycoside hydrolase family 5 cellulase\t1e-45
G4\tp1\tPfam\tPF00001\thypothetical protein\t1e-5
"""

LABELS_TSV = """genome_id\tsplit\tfamily\tpanel\tlabel\tvalue
G1\ttrain\tFamilyA\tbiofuels_industrial\tcellulose_degradation\t1
G2\ttrain\tFamilyA\tbiofuels_industrial\tcellulose_degradation\t0
G3\ttest\tFamilyB\tbiofuels_industrial\tcellulose_degradation\t1
G4\ttest\tFamilyB\tbiofuels_industrial\tcellulose_degradation\t0
"""


class LearningTests(unittest.TestCase):
    def test_train_baseline_learns_positive_feature(self):
        features = build_feature_matrix_from_annotation_tsv(ANNOTATIONS_TSV)
        labels = parse_labels_tsv(LABELS_TSV)

        model = train_baseline(features, labels)

        self.assertIn("biofuels_industrial:cellulose_degradation", model.targets)
        self.assertGreater(
            model.predict_proba("biofuels_industrial:cellulose_degradation", "G3", features),
            model.predict_proba("biofuels_industrial:cellulose_degradation", "G4", features),
        )

    def test_model_round_trips_json(self):
        features = build_feature_matrix_from_annotation_tsv(ANNOTATIONS_TSV)
        labels = parse_labels_tsv(LABELS_TSV)
        model = train_baseline(features, labels)

        self.assertEqual(model, model.from_json(model.to_json()))

    def test_evaluate_model_reports_accuracy(self):
        features = build_feature_matrix_from_annotation_tsv(ANNOTATIONS_TSV)
        labels = parse_labels_tsv(LABELS_TSV)
        model = train_baseline(features, labels)

        report = evaluate_model(model, features, labels, split="test")

        self.assertEqual(report["split"], "test")
        self.assertEqual(report["overall"]["n"], 2)
        self.assertEqual(report["overall"]["accuracy"], 1.0)
        json.dumps(report)

    def test_predict_model_returns_valid_prediction_json(self):
        features = build_feature_matrix_from_annotation_tsv(ANNOTATIONS_TSV)
        labels = parse_labels_tsv(LABELS_TSV)
        model = train_baseline(features, labels)

        prediction = predict_model(model, features, genome_id="G3")

        validate_prediction(prediction)
        self.assertEqual(prediction["genome_id"], "G3")
        self.assertEqual(prediction["functions"][0]["name"], "cellulose degradation")


if __name__ == "__main__":
    unittest.main()
