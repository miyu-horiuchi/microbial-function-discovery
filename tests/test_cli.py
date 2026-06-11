import subprocess
import sys
import unittest
import json


class CliTests(unittest.TestCase):
    def test_panels_command_lists_application_areas(self):
        result = subprocess.run(
            [sys.executable, "-m", "microbial_function_discovery.cli", "panels"],
            check=True,
            capture_output=True,
            text=True,
        )

        self.assertIn("environmental_terraforming", result.stdout)
        self.assertIn("Biofuels / Industrial Enzymes", result.stdout)

    def test_validate_command_accepts_example(self):
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "microbial_function_discovery.cli",
                "validate",
                "examples/prediction.example.json",
            ],
            check=True,
            capture_output=True,
            text=True,
        )

        self.assertIn("valid", result.stdout)

    def test_predict_command_emits_prediction_json(self):
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "microbial_function_discovery.cli",
                "predict",
                "examples/useful_functions.faa",
                "--genome-id",
                "candidate_001",
            ],
            check=True,
            capture_output=True,
            text=True,
        )

        prediction = json.loads(result.stdout)
        self.assertEqual(prediction["genome_id"], "candidate_001")
        self.assertIn("application_scores", prediction)


if __name__ == "__main__":
    unittest.main()
