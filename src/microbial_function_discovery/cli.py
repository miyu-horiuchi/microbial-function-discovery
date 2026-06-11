"""Command line interface for microbial function discovery scaffolding."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from microbial_function_discovery.annotations import parse_annotation_hits_tsv
from microbial_function_discovery.baseline import predict_from_annotation_hits, predict_from_fasta
from microbial_function_discovery.importers import format_annotation_hits_tsv, parse_eggnog_mapper, parse_hmmer_domtblout
from microbial_function_discovery.panels import list_panels
from microbial_function_discovery.validation import PredictionValidationError, validate_prediction


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mfd",
        description="Microbial function discovery utilities.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("panels", help="List application panels.")

    validate_parser = subparsers.add_parser("validate", help="Validate a prediction JSON file.")
    validate_parser.add_argument("path", type=Path, help="Path to prediction JSON.")

    predict_parser = subparsers.add_parser("predict", help="Run the baseline FASTA-to-prediction model.")
    predict_parser.add_argument("path", type=Path, help="Path to genome/protein FASTA.")
    predict_parser.add_argument("--genome-id", default=None, help="Genome id to place in the prediction JSON.")

    predict_annotations_parser = subparsers.add_parser(
        "predict-annotations",
        help="Run the baseline annotation-hit-to-prediction model.",
    )
    predict_annotations_parser.add_argument("path", type=Path, help="Path to annotation-hit TSV.")
    predict_annotations_parser.add_argument(
        "--genome-id",
        default=None,
        help="Genome id to place in the prediction JSON.",
    )

    import_eggnog_parser = subparsers.add_parser(
        "import-eggnog",
        help="Convert eggNOG-mapper annotations to normalized annotation-hit TSV.",
    )
    import_eggnog_parser.add_argument("path", type=Path, help="Path to .emapper.annotations file.")

    import_domtblout_parser = subparsers.add_parser(
        "import-domtblout",
        help="Convert HMMER --domtblout output to normalized annotation-hit TSV.",
    )
    import_domtblout_parser.add_argument("path", type=Path, help="Path to domtblout file.")
    import_domtblout_parser.add_argument("--database", required=True, help="Database name, e.g. Pfam or dbCAN.")

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "panels":
        return _panels()
    if args.command == "validate":
        return _validate(args.path)
    if args.command == "predict":
        return _predict(args.path, args.genome_id)
    if args.command == "predict-annotations":
        return _predict_annotations(args.path, args.genome_id)
    if args.command == "import-eggnog":
        return _import_eggnog(args.path)
    if args.command == "import-domtblout":
        return _import_domtblout(args.path, args.database)

    parser.error(f"unknown command: {args.command}")
    return 2


def _panels() -> int:
    for panel in list_panels():
        print(f"{panel.key}\t{panel.name}")
    return 0


def _validate(path: Path) -> int:
    try:
        prediction = json.loads(path.read_text())
        validate_prediction(prediction)
    except FileNotFoundError:
        print(f"{path}: file not found", file=sys.stderr)
        return 1
    except json.JSONDecodeError as exc:
        print(f"{path}: invalid JSON: {exc}", file=sys.stderr)
        return 1
    except PredictionValidationError as exc:
        print(f"{path}: invalid prediction: {exc}", file=sys.stderr)
        return 1

    print(f"{path}: valid")
    return 0


def _predict(path: Path, genome_id: str | None) -> int:
    try:
        fasta_text = path.read_text()
        prediction = predict_from_fasta(fasta_text, genome_id=genome_id or path.stem)
        validate_prediction(prediction)
    except FileNotFoundError:
        print(f"{path}: file not found", file=sys.stderr)
        return 1
    except ValueError as exc:
        print(f"{path}: invalid FASTA: {exc}", file=sys.stderr)
        return 1
    except PredictionValidationError as exc:
        print(f"{path}: invalid prediction generated: {exc}", file=sys.stderr)
        return 1

    print(json.dumps(prediction, indent=2, sort_keys=True))
    return 0


def _predict_annotations(path: Path, genome_id: str | None) -> int:
    try:
        annotation_text = path.read_text()
        hits = parse_annotation_hits_tsv(annotation_text)
        prediction = predict_from_annotation_hits(hits, genome_id=genome_id or path.stem)
        validate_prediction(prediction)
    except FileNotFoundError:
        print(f"{path}: file not found", file=sys.stderr)
        return 1
    except ValueError as exc:
        print(f"{path}: invalid annotation TSV: {exc}", file=sys.stderr)
        return 1
    except PredictionValidationError as exc:
        print(f"{path}: invalid prediction generated: {exc}", file=sys.stderr)
        return 1

    print(json.dumps(prediction, indent=2, sort_keys=True))
    return 0


def _import_eggnog(path: Path) -> int:
    try:
        hits = parse_eggnog_mapper(path.read_text())
    except FileNotFoundError:
        print(f"{path}: file not found", file=sys.stderr)
        return 1
    except ValueError as exc:
        print(f"{path}: invalid eggNOG annotations: {exc}", file=sys.stderr)
        return 1

    print(format_annotation_hits_tsv(hits), end="")
    return 0


def _import_domtblout(path: Path, database: str) -> int:
    try:
        hits = parse_hmmer_domtblout(path.read_text(), database=database)
    except FileNotFoundError:
        print(f"{path}: file not found", file=sys.stderr)
        return 1
    except ValueError as exc:
        print(f"{path}: invalid domtblout: {exc}", file=sys.stderr)
        return 1

    print(format_annotation_hits_tsv(hits), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
