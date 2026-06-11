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

    def test_predict_annotations_command_emits_prediction_json(self):
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "microbial_function_discovery.cli",
                "predict-annotations",
                "examples/annotation_hits.tsv",
                "--genome-id",
                "candidate_001",
            ],
            check=True,
            capture_output=True,
            text=True,
        )

        prediction = json.loads(result.stdout)
        self.assertEqual(prediction["genome_id"], "candidate_001")
        function_names = {function["name"] for function in prediction["functions"]}
        self.assertIn("cellulose degradation", function_names)

    def test_import_eggnog_command_outputs_annotation_tsv(self):
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "microbial_function_discovery.cli",
                "import-eggnog",
                "examples/eggnog_mapper.emapper.annotations",
            ],
            check=True,
            capture_output=True,
            text=True,
        )

        self.assertIn("protein_id\tdatabase\taccession\tname\tevalue", result.stdout)
        self.assertIn("p1\tCAZy\tGH5", result.stdout)

    def test_import_domtblout_command_outputs_annotation_tsv(self):
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "microbial_function_discovery.cli",
                "import-domtblout",
                "examples/pfam.domtblout",
                "--database",
                "Pfam",
            ],
            check=True,
            capture_output=True,
            text=True,
        )

        self.assertIn("p1\tPfam\tPF00150.20", result.stdout)


if __name__ == "__main__":
    unittest.main()
