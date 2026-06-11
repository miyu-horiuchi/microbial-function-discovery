import subprocess
import sys
import unittest
import json
import stat
import tempfile
import textwrap
from pathlib import Path


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

    def test_run_eggnog_command_outputs_annotation_tsv(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fasta = tmp_path / "proteins.faa"
            fasta.write_text(">p1\nMKK\n")
            fake = _write_fake_executable(
                tmp_path / "fake_emapper.py",
                """
                import pathlib
                import sys

                args = sys.argv
                out_dir = pathlib.Path(args[args.index("--output_dir") + 1])
                prefix = args[args.index("-o") + 1]
                out_dir.mkdir(parents=True, exist_ok=True)
                (out_dir / f"{prefix}.emapper.annotations").write_text(
                    "#query\\tevalue\\tDescription\\tKEGG_ko\\tCAZy\\tPFAMs\\n"
                    "p1\\t1e-40\\tglycoside hydrolase family 5 cellulase\\tko:K01179\\tGH5\\tPF00150\\n"
                )
                """,
            )

            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "microbial_function_discovery.cli",
                    "run-eggnog",
                    str(fasta),
                    "--output-dir",
                    str(tmp_path),
                    "--executable",
                    str(fake),
                ],
                check=True,
                capture_output=True,
                text=True,
            )

        self.assertIn("p1\tCAZy\tGH5", result.stdout)

    def test_run_domtblout_command_outputs_annotation_tsv(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fasta = tmp_path / "proteins.faa"
            fasta.write_text(">p1\nMKK\n")
            hmm = tmp_path / "Pfam-A.hmm"
            hmm.write_text("HMMER3/f\n")
            out = tmp_path / "pfam.domtblout"
            fake = _write_fake_executable(
                tmp_path / "fake_hmmscan.py",
                """
                import pathlib
                import sys

                args = sys.argv
                out = pathlib.Path(args[args.index("--domtblout") + 1])
                out.write_text(
                    "# domtblout\\n"
                    "GH5.hmm PF00150.20 300 p1 - 320 1e-40 180.0 0.0 1 1 1e-42 1e-40 180.0 0.0 5 290 10 300 8 305 0.98 glycoside hydrolase family 5 cellulase\\n"
                )
                """,
            )

            result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "microbial_function_discovery.cli",
                    "run-domtblout",
                    str(fasta),
                    "--hmm",
                    str(hmm),
                    "--out",
                    str(out),
                    "--database",
                    "Pfam",
                    "--executable",
                    str(fake),
                ],
                check=True,
                capture_output=True,
                text=True,
            )

        self.assertIn("p1\tPfam\tPF00150.20", result.stdout)

    def test_validate_splits_command_accepts_family_holdout_labels(self):
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "microbial_function_discovery.cli",
                "validate-splits",
                "examples/benchmark_labels.tsv",
            ],
            check=True,
            capture_output=True,
            text=True,
        )

        self.assertIn("family holdout valid", result.stdout)

    def test_build_features_train_and_evaluate_commands(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            features_path = tmp_path / "features.json"
            model_path = tmp_path / "model.json"
            report_path = tmp_path / "report.json"

            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "microbial_function_discovery.cli",
                    "build-features",
                    "examples/multi_genome_annotation_hits.tsv",
                    "--out",
                    str(features_path),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "microbial_function_discovery.cli",
                    "train-baseline",
                    str(features_path),
                    "examples/benchmark_labels.tsv",
                    "--out",
                    str(model_path),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "microbial_function_discovery.cli",
                    "evaluate",
                    str(model_path),
                    str(features_path),
                    "examples/benchmark_labels.tsv",
                    "--split",
                    "test",
                    "--out",
                    str(report_path),
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            report = json.loads(report_path.read_text())
            self.assertEqual(report["overall"]["accuracy"], 1.0)

            prediction_result = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "microbial_function_discovery.cli",
                    "predict-baseline",
                    str(model_path),
                    str(features_path),
                    "--genome-id",
                    "G3",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            prediction = json.loads(prediction_result.stdout)
            self.assertEqual(prediction["genome_id"], "G3")
            self.assertIn("functions", prediction)


def _write_fake_executable(path: Path, body: str) -> Path:
    script = "#!/usr/bin/env python3\n" + textwrap.dedent(body).strip() + "\n"
    path.write_text(script)
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


if __name__ == "__main__":
    unittest.main()
