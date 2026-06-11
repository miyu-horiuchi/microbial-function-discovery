import unittest

from microbial_function_discovery.baseline import predict_from_fasta
from microbial_function_discovery.validation import validate_prediction


SAMPLE_FASTA = """>bin_001_protein_1 GH5 glycoside hydrolase cellulase
MKKLAA
>bin_001_protein_2 nitrogenase iron protein nifH
MSSGTT
>bin_001_protein_3 beta-lactamase antimicrobial resistance protein
MEEEKK
>bin_001_protein_4 adhesin virulence factor
MCCCAA
"""


class BaselinePredictorTests(unittest.TestCase):
    def test_predict_from_fasta_returns_valid_prediction_contract(self):
        prediction = predict_from_fasta(SAMPLE_FASTA, genome_id="bin_001")

        validate_prediction(prediction)
        self.assertEqual(prediction["genome_id"], "bin_001")

    def test_predict_from_fasta_scores_functions_from_evidence(self):
        prediction = predict_from_fasta(SAMPLE_FASTA, genome_id="bin_001")
        function_names = {function["name"] for function in prediction["functions"]}

        self.assertIn("cellulose degradation", function_names)
        self.assertIn("nitrogen fixation", function_names)
        self.assertIn("antimicrobial resistance", function_names)

    def test_biosafety_risk_increases_from_amr_and_virulence_evidence(self):
        prediction = predict_from_fasta(SAMPLE_FASTA, genome_id="bin_001")

        self.assertGreaterEqual(prediction["biosafety"]["amr_risk"], 0.5)
        self.assertGreaterEqual(prediction["biosafety"]["pathogenicity_risk"], 0.5)
        self.assertTrue(prediction["biosafety"]["warnings"])


if __name__ == "__main__":
    unittest.main()
