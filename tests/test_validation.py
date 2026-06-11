import json
import unittest
from pathlib import Path

from microbial_function_discovery.validation import PredictionValidationError, validate_prediction


class PredictionValidationTests(unittest.TestCase):
    def load_example(self):
        return json.loads(Path("examples/prediction.example.json").read_text())

    def test_example_prediction_is_valid(self):
        validate_prediction(self.load_example())

    def test_invalid_application_area_is_rejected(self):
        prediction = self.load_example()
        prediction["application_scores"][0]["area"] = "unknown"

        with self.assertRaisesRegex(PredictionValidationError, r"application_scores\[0\].area"):
            validate_prediction(prediction)

    def test_out_of_range_function_score_is_rejected(self):
        prediction = self.load_example()
        prediction["functions"][0]["score"] = 1.2

        with self.assertRaisesRegex(PredictionValidationError, r"functions\[0\].score"):
            validate_prediction(prediction)

    def test_missing_evidence_id_is_rejected(self):
        prediction = self.load_example()
        del prediction["functions"][0]["evidence"][0]["id"]

        with self.assertRaisesRegex(PredictionValidationError, r"functions\[0\].evidence\[0\].id"):
            validate_prediction(prediction)


if __name__ == "__main__":
    unittest.main()
