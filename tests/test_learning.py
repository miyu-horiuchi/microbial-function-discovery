import json
import unittest

from microbial_function_discovery.datasets import parse_labels_tsv
from microbial_function_discovery.features import build_feature_matrix_from_annotation_tsv
from microbial_function_discovery.learning import (
    BaselineModel,
    TargetModel,
    evaluate_model,
    evaluate_ranking,
    predict_model,
    rank_candidates,
    train_baseline,
)
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

    def test_train_baseline_ignores_labels_without_features(self):
        features = build_feature_matrix_from_annotation_tsv(ANNOTATIONS_TSV)
        labels = parse_labels_tsv(
            LABELS_TSV
            + "G_missing\ttrain\tFamilyA\tbiofuels_industrial\tcellulose_degradation\t1\n"
        )

        model = train_baseline(features, labels)

        self.assertIn("biofuels_industrial:cellulose_degradation", model.targets)

    def test_model_round_trips_json(self):
        features = build_feature_matrix_from_annotation_tsv(ANNOTATIONS_TSV)
        labels = parse_labels_tsv(LABELS_TSV)
        model = train_baseline(features, labels)

        self.assertEqual(model, model.from_json(model.to_json()))

    def test_predict_proba_handles_extreme_logits(self):
        features = build_feature_matrix_from_annotation_tsv(ANNOTATIONS_TSV)
        model = BaselineModel(
            feature_names=features.feature_names,
            targets={
                "biofuels_industrial:cellulose_degradation": TargetModel(
                    prior_log_odds=-1000.0,
                    feature_log_odds=[0.0 for _feature in features.feature_names],
                )
            },
        )

        self.assertEqual(
            model.predict_proba("biofuels_industrial:cellulose_degradation", "G1", features),
            0.0,
        )

    def test_evaluate_model_reports_accuracy(self):
        features = build_feature_matrix_from_annotation_tsv(ANNOTATIONS_TSV)
        labels = parse_labels_tsv(LABELS_TSV)
        model = train_baseline(features, labels)

        report = evaluate_model(model, features, labels, split="test")

        self.assertEqual(report["split"], "test")
        self.assertEqual(report["overall"]["n"], 2)
        self.assertEqual(report["overall"]["accuracy"], 1.0)
        json.dumps(report)

    def test_evaluate_model_ignores_labels_without_features(self):
        features = build_feature_matrix_from_annotation_tsv(ANNOTATIONS_TSV)
        labels = parse_labels_tsv(
            LABELS_TSV
            + "G_missing\ttest\tFamilyB\tbiofuels_industrial\tcellulose_degradation\t1\n"
        )
        model = train_baseline(features, labels)

        report = evaluate_model(model, features, labels, split="test")

        self.assertEqual(report["overall"]["n"], 2)

    def test_predict_model_returns_valid_prediction_json(self):
        features = build_feature_matrix_from_annotation_tsv(ANNOTATIONS_TSV)
        labels = parse_labels_tsv(LABELS_TSV)
        model = train_baseline(features, labels)

        prediction = predict_model(model, features, genome_id="G3")

        validate_prediction(prediction)
        self.assertEqual(prediction["genome_id"], "G3")
        self.assertEqual(prediction["functions"][0]["name"], "cellulose degradation")

    def test_rank_candidates_by_target(self):
        features = build_feature_matrix_from_annotation_tsv(ANNOTATIONS_TSV)
        labels = parse_labels_tsv(LABELS_TSV)
        model = train_baseline(features, labels)

        ranking = rank_candidates(model, features, target_key="biofuels_industrial:cellulose_degradation")

        self.assertEqual(ranking["mode"], "target")
        self.assertEqual(ranking["query"], "biofuels_industrial:cellulose_degradation")
        self.assertEqual(ranking["candidates"][0]["genome_id"], "G1")
        self.assertGreater(ranking["candidates"][0]["score"], ranking["candidates"][-1]["score"])

    def test_rank_candidates_by_panel(self):
        features = build_feature_matrix_from_annotation_tsv(ANNOTATIONS_TSV)
        labels = parse_labels_tsv(LABELS_TSV)
        model = train_baseline(features, labels)

        ranking = rank_candidates(model, features, panel="biofuels_industrial", limit=2)

        self.assertEqual(ranking["mode"], "panel")
        self.assertEqual(ranking["query"], "biofuels_industrial")
        self.assertEqual(len(ranking["candidates"]), 2)

    def test_evaluate_ranking_by_target_reports_precision_at_k(self):
        features = build_feature_matrix_from_annotation_tsv(ANNOTATIONS_TSV)
        labels = parse_labels_tsv(LABELS_TSV)
        model = train_baseline(features, labels)

        report = evaluate_ranking(
            model,
            features,
            labels,
            split="test",
            target_key="biofuels_industrial:cellulose_degradation",
            ks=[1, 2],
        )

        self.assertEqual(report["mode"], "target")
        self.assertEqual(report["query"], "biofuels_industrial:cellulose_degradation")
        self.assertEqual(report["split"], "test")
        self.assertEqual(report["n_labeled_candidates"], 2)
        self.assertEqual(report["n_positives"], 1)
        self.assertEqual(report["metrics"]["precision_at_1"], 1.0)
        self.assertEqual(report["metrics"]["recall_at_1"], 1.0)
        self.assertEqual(report["ranked_candidates"][0]["genome_id"], "G3")
        json.dumps(report)

    def test_evaluate_ranking_by_panel_treats_any_positive_panel_label_as_hit(self):
        features = build_feature_matrix_from_annotation_tsv(ANNOTATIONS_TSV)
        labels = parse_labels_tsv(LABELS_TSV)
        model = train_baseline(features, labels)

        report = evaluate_ranking(model, features, labels, split="test", panel="biofuels_industrial", ks=[1])

        self.assertEqual(report["mode"], "panel")
        self.assertEqual(report["query"], "biofuels_industrial")
        self.assertEqual(report["metrics"]["hits_at_1"], 1)


if __name__ == "__main__":
    unittest.main()
