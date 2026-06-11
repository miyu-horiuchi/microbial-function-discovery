import unittest

from microbial_function_discovery.features import build_feature_matrix_from_annotation_tsv


ANNOTATIONS_TSV = """genome_id\tprotein_id\tdatabase\taccession\tname\tevalue
G1\tp1\tCAZy\tGH5\tglycoside hydrolase family 5 cellulase\t1e-40
G1\tp2\tKEGG\tK02588\tnitrogenase iron protein nifH\t1e-50
G2\tp1\tPfam\tPF00001\thypothetical protein\t1e-5
"""


class FeatureTests(unittest.TestCase):
    def test_build_feature_matrix_from_annotation_tsv(self):
        matrix = build_feature_matrix_from_annotation_tsv(ANNOTATIONS_TSV)

        self.assertEqual(matrix.genome_ids, ["G1", "G2"])
        self.assertEqual(matrix.feature_names, ["CAZy:GH5", "KEGG:K02588", "Pfam:PF00001"])
        self.assertEqual(matrix.rows, [[1, 1, 0], [0, 0, 1]])

    def test_feature_matrix_round_trips_json(self):
        matrix = build_feature_matrix_from_annotation_tsv(ANNOTATIONS_TSV)

        self.assertEqual(matrix, matrix.from_json(matrix.to_json()))


if __name__ == "__main__":
    unittest.main()
