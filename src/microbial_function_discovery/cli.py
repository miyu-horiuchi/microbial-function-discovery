"""Command line interface for microbial function discovery scaffolding."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from microbial_function_discovery.baseline import predict_from_fasta
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


if __name__ == "__main__":
    raise SystemExit(main())
