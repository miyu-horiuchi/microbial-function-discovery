"""Wrappers for external CPU annotation tools."""

from __future__ import annotations

import subprocess
from pathlib import Path


def run_eggnog_mapper(
    *,
    input_fasta: Path,
    output_dir: Path,
    executable: str = "emapper.py",
    data_dir: Path | None = None,
    cpus: int = 1,
    output_prefix: str = "eggnog",
    force: bool = False,
) -> Path:
    """Run eggNOG-mapper on a protein FASTA file.

    Returns the expected ``*.emapper.annotations`` path. If the output already
    exists and ``force`` is false, the command is skipped so prior CPU work is
    reused.
    """

    output_dir = Path(output_dir)
    output_path = output_dir / f"{output_prefix}.emapper.annotations"
    if _usable_existing(output_path) and not force:
        return output_path

    _require_existing(input_fasta, "input FASTA")
    output_dir.mkdir(parents=True, exist_ok=True)

    command = [
        executable,
        "-i",
        str(input_fasta),
        "--itype",
        "proteins",
        "-m",
        "diamond",
        "--cpu",
        str(cpus),
        "--output_dir",
        str(output_dir),
        "-o",
        output_prefix,
        "--override",
    ]
    if data_dir is not None:
        command.extend(["--data_dir", str(data_dir)])

    subprocess.run(command, check=True)
    if not _usable_existing(output_path):
        raise RuntimeError(f"eggNOG-mapper did not produce expected output: {output_path}")
    return output_path


def run_hmmer_domtblout(
    *,
    input_fasta: Path,
    hmm_database: Path,
    output_path: Path,
    database: str,
    executable: str = "hmmscan",
    cpus: int = 1,
    evalue: float = 1e-5,
    force: bool = False,
) -> Path:
    """Run HMMER hmmscan and write domain table output.

    ``database`` is kept in the signature for caller clarity and CLI symmetry;
    normalization uses it when parsing the produced file.
    """

    del database
    output_path = Path(output_path)
    if _usable_existing(output_path) and not force:
        return output_path

    _require_existing(input_fasta, "input FASTA")
    _require_existing(hmm_database, "HMM database")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    command = [
        executable,
        "--cpu",
        str(cpus),
        "-E",
        str(evalue),
        "--domtblout",
        str(output_path),
        str(hmm_database),
        str(input_fasta),
    ]
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL)
    if not _usable_existing(output_path):
        raise RuntimeError(f"HMMER did not produce expected output: {output_path}")
    return output_path


def _usable_existing(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


def _require_existing(path: Path, label: str) -> None:
    if not Path(path).exists():
        raise FileNotFoundError(f"{label} not found: {path}")
