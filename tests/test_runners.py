import os
import stat
import tempfile
import textwrap
import unittest
from pathlib import Path

from microbial_function_discovery.runners import run_eggnog_mapper, run_hmmer_domtblout


class RunnerTests(unittest.TestCase):
    def test_run_eggnog_mapper_skips_existing_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            existing = tmp_path / "eggnog.emapper.annotations"
            existing.write_text("#query\tevalue\tDescription\n")

            output = run_eggnog_mapper(
                input_fasta=tmp_path / "proteins.faa",
                output_dir=tmp_path,
                executable="missing-emapper",
            )

            self.assertEqual(output, existing)

    def test_run_eggnog_mapper_invokes_executable(self):
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
                (out_dir / f"{prefix}.emapper.annotations").write_text("#query\\tevalue\\tDescription\\n")
                """,
            )

            output = run_eggnog_mapper(input_fasta=fasta, output_dir=tmp_path, executable=str(fake))

            self.assertTrue(output.exists())

    def test_run_hmmer_domtblout_invokes_executable(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fasta = tmp_path / "proteins.faa"
            fasta.write_text(">p1\nMKK\n")
            hmm = tmp_path / "Pfam-A.hmm"
            hmm.write_text("HMMER3/f\n")
            output = tmp_path / "pfam.domtblout"
            fake = _write_fake_executable(
                tmp_path / "fake_hmmscan.py",
                """
                import pathlib
                import sys

                args = sys.argv
                out = pathlib.Path(args[args.index("--domtblout") + 1])
                out.write_text("# domtblout\\n")
                """,
            )

            result = run_hmmer_domtblout(
                input_fasta=fasta,
                hmm_database=hmm,
                output_path=output,
                executable=str(fake),
                database="Pfam",
            )

            self.assertEqual(result, output)
            self.assertTrue(output.exists())


def _write_fake_executable(path: Path, body: str) -> Path:
    script = "#!/usr/bin/env python3\n" + textwrap.dedent(body).strip() + "\n"
    path.write_text(script)
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


if __name__ == "__main__":
    unittest.main()
