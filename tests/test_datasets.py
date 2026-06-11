import unittest

from microbial_function_discovery.datasets import FamilyLeakageError, parse_labels_tsv, validate_family_holdout


LABELS_TSV = """genome_id\tsplit\tfamily\tpanel\tlabel\tvalue
G1\ttrain\tFamilyA\tbiofuels_industrial\tcellulose_degradation\t1
G2\ttest\tFamilyB\tbiofuels_industrial\tcellulose_degradation\t0
"""


class DatasetTests(unittest.TestCase):
    def test_parse_labels_tsv(self):
        records = parse_labels_tsv(LABELS_TSV)

        self.assertEqual(records[0].genome_id, "G1")
        self.assertEqual(records[0].target_key, "biofuels_industrial:cellulose_degradation")
        self.assertEqual(records[0].value, 1)

    def test_validate_family_holdout_allows_disjoint_families(self):
        validate_family_holdout(parse_labels_tsv(LABELS_TSV))

    def test_validate_family_holdout_rejects_family_leakage(self):
        records = parse_labels_tsv(
            """genome_id\tsplit\tfamily\tpanel\tlabel\tvalue
G1\ttrain\tFamilyA\tbiofuels_industrial\tcellulose_degradation\t1
G2\ttest\tFamilyA\tbiofuels_industrial\tcellulose_degradation\t0
"""
        )

        with self.assertRaisesRegex(FamilyLeakageError, "FamilyA"):
            validate_family_holdout(records)


if __name__ == "__main__":
    unittest.main()
