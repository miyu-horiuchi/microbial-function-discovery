import json
import unittest

from microbial_function_discovery.datasets import parse_labels_tsv
from microbial_function_discovery.dense_learning import DenseFeatureMatrix, train_dense_baseline
from microbial_function_discovery.features import build_feature_matrix_from_annotation_tsv
from microbial_function_discovery.fusion import (
    evaluate_fusion_leaderboard,
    evaluate_fusion_ranking,
    rank_discovery_candidates,
)
from microbial_function_discovery.learning import train_baseline


ANNOTATIONS_TSV = """genome_id\tprotein_id\tdatabase\taccession\tname\tevalue
G1\tp1\tPfam\tPF00001\ttraining positive marker\t1e-20
G2\tp1\tPfam\tPF00002\ttraining negative marker\t1e-20
G3\tp1\tPfam\tPF00001\ttest false positive marker\t1e-20
G4\tp1\tPfam\tPF00002\ttest true positive dense marker\t1e-20
G5\tp1\tPfam\tPF00002\tdense missing candidate\t1e-20
"""

LABELS_TSV = """genome_id\tsplit\tfamily\tpanel\tlabel\tvalue
G1\ttrain\tFamilyA\tbiofuels_industrial\tthermophile\t1
G2\ttrain\tFamilyA\tbiofuels_industrial\tthermophile\t0
G3\ttest\tFamilyB\tbiofuels_industrial\tthermophile\t0
G4\ttest\tFamilyB\tbiofuels_industrial\tthermophile\t1
G5\ttest\tFamilyB\tbiofuels_industrial\tthermophile\t1
"""


class FusionTests(unittest.TestCase):
    def test_evaluate_fusion_ranking_selects_best_weight(self):
        annotation_features = build_feature_matrix_from_annotation_tsv(ANNOTATIONS_TSV)
        labels = parse_labels_tsv(LABELS_TSV)
        annotation_model = train_baseline(annotation_features, labels)
        dense_features = DenseFeatureMatrix(
            genome_ids=["G1", "G2", "G3", "G4"],
            feature_indexes=[0],
            rows=[[2.0], [-2.0], [-1.5], [1.5]],
        )
        dense_model = train_dense_baseline(dense_features, labels)

        report = evaluate_fusion_ranking(
            annotation_model,
            annotation_features,
            dense_model,
            dense_features,
            labels,
            split="test",
            target_key="biofuels_industrial:thermophile",
            ks=[1, 2],
            weights=[0.0, 1.0],
            select_k=1,
        )

        self.assertEqual(report["mode"], "target")
        self.assertEqual(report["best_weight"], 1.0)
        self.assertEqual(report["selection_metric"], "precision_at_1")
        self.assertEqual(report["best"]["metrics"]["precision_at_1"], 1.0)
        self.assertEqual(report["best"]["ranked_candidates"][0]["genome_id"], "G4")
        self.assertEqual(report["per_weight"][0]["weight"], 0.0)
        self.assertEqual(report["per_weight"][1]["weight"], 1.0)
        json.dumps(report)

    def test_evaluate_fusion_ranking_scores_only_common_labeled_genomes(self):
        annotation_features = build_feature_matrix_from_annotation_tsv(ANNOTATIONS_TSV)
        labels = parse_labels_tsv(LABELS_TSV)
        annotation_model = train_baseline(annotation_features, labels)
        dense_features = DenseFeatureMatrix(
            genome_ids=["G1", "G2", "G3", "G4"],
            feature_indexes=[0],
            rows=[[2.0], [-2.0], [-1.5], [1.5]],
        )
        dense_model = train_dense_baseline(dense_features, labels)

        report = evaluate_fusion_ranking(
            annotation_model,
            annotation_features,
            dense_model,
            dense_features,
            labels,
            split="test",
            target_key="biofuels_industrial:thermophile",
            ks=[1],
            weights=[1.0],
        )

        ranked_genomes = {row["genome_id"] for row in report["best"]["ranked_candidates"]}
        self.assertEqual(report["best"]["n_labeled_candidates"], 2)
        self.assertEqual(ranked_genomes, {"G3", "G4"})


LEADERBOARD_ANNOTATIONS_TSV = """genome_id\tprotein_id\tdatabase\taccession\tname\tevalue
G1\tp1\tPfam\tPF00001\ttraining positive marker\t1e-20
G2\tp1\tPfam\tPF00002\ttraining negative marker\t1e-20
G3\tp1\tPfam\tPF00001\tvalidation false positive marker\t1e-20
G4\tp1\tPfam\tPF00002\tvalidation true positive dense marker\t1e-20
G5\tp1\tPfam\tPF00001\ttest false positive marker\t1e-20
G6\tp1\tPfam\tPF00002\ttest true positive dense marker\t1e-20
"""

LEADERBOARD_LABELS_TSV = """genome_id\tsplit\tfamily\tpanel\tlabel\tvalue
G1\ttrain\tFamilyA\tbiofuels_industrial\tthermophile\t1
G2\ttrain\tFamilyA\tbiofuels_industrial\tthermophile\t0
G3\tval\tFamilyB\tbiofuels_industrial\tthermophile\t0
G4\tval\tFamilyB\tbiofuels_industrial\tthermophile\t1
G5\ttest\tFamilyC\tbiofuels_industrial\tthermophile\t0
G6\ttest\tFamilyC\tbiofuels_industrial\tthermophile\t1
"""


class FusionLeaderboardTests(unittest.TestCase):
    def test_evaluate_fusion_leaderboard_selects_on_validation_and_reports_test(self):
        annotation_features = build_feature_matrix_from_annotation_tsv(LEADERBOARD_ANNOTATIONS_TSV)
        labels = parse_labels_tsv(LEADERBOARD_LABELS_TSV)
        annotation_model = train_baseline(annotation_features, labels)
        dense_features = DenseFeatureMatrix(
            genome_ids=["G1", "G2", "G3", "G4", "G5", "G6"],
            feature_indexes=[0],
            rows=[[2.0], [-2.0], [-1.5], [1.5], [-1.4], [1.4]],
        )
        dense_model = train_dense_baseline(dense_features, labels)

        leaderboard = evaluate_fusion_leaderboard(
            annotation_model,
            annotation_features,
            {"dense": (dense_model, dense_features)},
            labels,
            validation_split="val",
            test_split="test",
            ks=[1, 2],
            weights=[0.0, 1.0],
            select_k=1,
        )

        row = leaderboard["targets"][0]
        self.assertEqual(row["target_key"], "biofuels_industrial:thermophile")
        self.assertEqual(row["best_source"], "dense")
        self.assertEqual(row["best_weight"], 1.0)
        self.assertEqual(row["selected_on_split"], "val")
        self.assertEqual(row["validation_metrics"]["precision_at_1"], 1.0)
        self.assertEqual(row["test_metrics"]["precision_at_1"], 1.0)
        self.assertTrue(row["beats_annotation"])
        json.dumps(leaderboard)

    def test_rank_discovery_candidates_uses_leaderboard_best_source(self):
        annotation_features = build_feature_matrix_from_annotation_tsv(LEADERBOARD_ANNOTATIONS_TSV)
        labels = parse_labels_tsv(LEADERBOARD_LABELS_TSV)
        annotation_model = train_baseline(annotation_features, labels)
        dense_features = DenseFeatureMatrix(
            genome_ids=["G1", "G2", "G3", "G4", "G5", "G6"],
            feature_indexes=[0],
            rows=[[2.0], [-2.0], [-1.5], [1.5], [-1.4], [1.4]],
        )
        dense_model = train_dense_baseline(dense_features, labels)
        leaderboard = evaluate_fusion_leaderboard(
            annotation_model,
            annotation_features,
            {"dense": (dense_model, dense_features)},
            labels,
            validation_split="val",
            test_split="test",
            ks=[1],
            weights=[0.0, 1.0],
            select_k=1,
        )

        candidates = rank_discovery_candidates(
            leaderboard,
            annotation_model,
            annotation_features,
            {"dense": (dense_model, dense_features)},
            labels,
            split="test",
            limit_per_target=1,
            min_precision=0.5,
            precision_k=1,
        )

        self.assertEqual(candidates["targets"][0]["target_key"], "biofuels_industrial:thermophile")
        self.assertEqual(candidates["targets"][0]["candidates"][0]["genome_id"], "G6")
        self.assertEqual(candidates["panels"]["biofuels_industrial"]["n_targets"], 1)
        json.dumps(candidates)


if __name__ == "__main__":
    unittest.main()
