import unittest

from microbial_function_discovery.annotations import AnnotationHit, parse_annotation_hits_tsv
from microbial_function_discovery.baseline import predict_from_annotation_hits
from microbial_function_discovery.validation import validate_prediction


ANNOTATION_TSV = """protein_id\tdatabase\taccession\tname\tevalue
p1\tCAZy\tGH5\tglycoside hydrolase family 5 cellulase\t1e-40
p2\tKEGG\tK02588\tnitrogenase iron protein nifH\t1e-50
p3\tCARD\tblaTEM\tbeta-lactamase\t1e-20
p4\tVFDB\tVF0001\tadhesin virulence factor\t1e-12
"""


class AnnotationHitTests(unittest.TestCase):
    def test_parse_annotation_hits_tsv(self):
        hits = parse_annotation_hits_tsv(ANNOTATION_TSV)

        self.assertEqual(
            hits[0],
            AnnotationHit(
                protein_id="p1",
                database="CAZy",
                accession="GH5",
                name="glycoside hydrolase family 5 cellulase",
                evalue=1e-40,
            ),
        )

    def test_parse_annotation_hits_rejects_missing_columns(self):
        with self.assertRaisesRegex(ValueError, "missing required columns"):
            parse_annotation_hits_tsv("protein_id\tdatabase\np1\tCAZy\n")

    def test_predict_from_annotation_hits_returns_valid_prediction(self):
        prediction = predict_from_annotation_hits(parse_annotation_hits_tsv(ANNOTATION_TSV), genome_id="bin_001")

        validate_prediction(prediction)
        self.assertEqual(prediction["genome_id"], "bin_001")

    def test_predict_from_annotation_hits_scores_functions(self):
        prediction = predict_from_annotation_hits(parse_annotation_hits_tsv(ANNOTATION_TSV), genome_id="bin_001")
        function_names = {function["name"] for function in prediction["functions"]}

        self.assertIn("cellulose degradation", function_names)
        self.assertIn("nitrogen fixation", function_names)
        self.assertIn("antimicrobial resistance", function_names)
        self.assertGreaterEqual(prediction["biosafety"]["pathogenicity_risk"], 0.5)


if __name__ == "__main__":
    unittest.main()
