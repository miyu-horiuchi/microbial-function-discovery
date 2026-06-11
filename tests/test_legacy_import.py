import unittest

from microbial_function_discovery.legacy_import import (
    dense_feature_matrix_from_legacy_arrays,
    feature_matrix_from_legacy_arrays,
    labels_tsv_from_legacy_tables,
    records_from_legacy_tables,
)


TRAIT_ROWS = [
    {
        "bacdive_id": 1,
        "family": "FamilyA",
        "temperature_class": "thermophile",
        "pathogenicity_human": True,
        "carbon_utilization": {"cellulose": True, "glucose": False, "acetate": None},
    },
    {
        "bacdive_id": 2,
        "family": "FamilyB",
        "temperature_class": "mesophile",
        "pathogenicity_human": False,
        "carbon_utilization": {"cellulose": False, "glucose": True, "acetate": None},
    },
]

SPLIT_ROWS = [
    {"bacdive_id": 1, "family": "FamilyA", "family_split": "train"},
    {"bacdive_id": 2, "family": "FamilyB", "family_split": "test"},
]


class LegacyImportTests(unittest.TestCase):
    def test_records_from_legacy_tables_exports_binary_and_one_vs_rest_labels(self):
        records = records_from_legacy_tables(TRAIT_ROWS, SPLIT_ROWS, split_column="family_split")

        by_key = {(record.genome_id, record.panel, record.label): record for record in records}
        self.assertEqual(by_key[("1", "biofuels_industrial", "temperature_class__thermophile")].value, 1)
        self.assertEqual(by_key[("2", "biofuels_industrial", "temperature_class__thermophile")].value, 0)
        self.assertEqual(by_key[("1", "biosafety", "pathogenicity_human")].value, 1)
        self.assertEqual(by_key[("2", "biosafety", "pathogenicity_human")].value, 0)
        self.assertEqual(by_key[("1", "biofuels_industrial", "carbon_utilization__cellulose")].value, 1)
        self.assertEqual(by_key[("2", "biofuels_industrial", "carbon_utilization__cellulose")].value, 0)

    def test_labels_tsv_from_legacy_tables_round_trips_to_required_format(self):
        text = labels_tsv_from_legacy_tables(TRAIT_ROWS, SPLIT_ROWS, split_column="family_split")

        self.assertTrue(text.startswith("genome_id\tsplit\tfamily\tpanel\tlabel\tvalue\n"))
        self.assertIn("1\ttrain\tFamilyA\tbiosafety\tpathogenicity_human\t1\n", text)

    def test_feature_matrix_from_legacy_arrays_selects_most_prevalent_features(self):
        matrix = feature_matrix_from_legacy_arrays(
            bacdive_ids=[1, 2, 3],
            feature_names=["rare", "common", "middle"],
            rows=[
                [0, 1, 1],
                [1, 1, 0],
                [0, 1, 0],
            ],
            max_features=2,
            feature_prefix="eggNOG",
        )

        self.assertEqual(matrix.genome_ids, ["1", "2", "3"])
        self.assertEqual(matrix.feature_names, ["eggNOG:common", "eggNOG:middle"])
        self.assertEqual(matrix.rows, [[1, 1], [1, 0], [1, 0]])

    def test_dense_feature_matrix_from_legacy_arrays_selects_high_variance_dimensions(self):
        matrix = dense_feature_matrix_from_legacy_arrays(
            bacdive_ids=[1, 2, 3],
            rows=[
                [0.1, 0.0, 2.0],
                [0.2, 10.0, 2.0],
                [0.3, -10.0, 2.0],
            ],
            max_features=2,
            feature_prefix="ESM2",
        )

        self.assertEqual(matrix.genome_ids, ["1", "2", "3"])
        self.assertEqual(matrix.feature_names, ["ESM2:dim_0000_gt_median", "ESM2:dim_0001_gt_median"])
        self.assertEqual(matrix.rows, [[0, 0], [0, 1], [1, 0]])


if __name__ == "__main__":
    unittest.main()
