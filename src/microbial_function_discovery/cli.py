"""Command line interface for microbial function discovery scaffolding."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from microbial_function_discovery.annotations import parse_annotation_hits_tsv
from microbial_function_discovery.baseline import predict_from_annotation_hits, predict_from_fasta
from microbial_function_discovery.datasets import FamilyLeakageError, parse_labels_tsv, validate_family_holdout
from microbial_function_discovery.dense_learning import (
    DenseBaselineModel,
    DenseFeatureMatrix,
    evaluate_dense_model,
    evaluate_dense_ranking,
    load_dense_npz_feature_matrix,
    train_dense_baseline,
)
from microbial_function_discovery.discovery import (
    annotate_discovery_candidates,
    render_discovery_candidate_report,
    to_json,
)
from microbial_function_discovery.features import FeatureMatrix, build_feature_matrix_from_annotation_tsv
from microbial_function_discovery.fusion import (
    evaluate_fusion_leaderboard,
    evaluate_fusion_ranking,
    rank_discovery_candidates,
)
from microbial_function_discovery.importers import format_annotation_hits_tsv, parse_eggnog_mapper, parse_hmmer_domtblout
from microbial_function_discovery.legacy_import import (
    labels_tsv_from_legacy_tables,
    load_legacy_dense_npz_feature_matrix,
    load_legacy_npz_feature_matrix,
    load_table_records,
)
from microbial_function_discovery.learning import (
    BaselineModel,
    evaluate_model,
    evaluate_ranking,
    predict_model,
    rank_candidates,
    train_baseline,
)
from microbial_function_discovery.panels import list_panels
from microbial_function_discovery.runners import run_eggnog_mapper, run_hmmer_domtblout
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

    run_eggnog_parser = subparsers.add_parser(
        "run-eggnog",
        help="Run eggNOG-mapper and print normalized annotation-hit TSV.",
    )
    run_eggnog_parser.add_argument("fasta", type=Path, help="Path to protein FASTA.")
    run_eggnog_parser.add_argument("--output-dir", type=Path, required=True, help="Directory for eggNOG output.")
    run_eggnog_parser.add_argument("--executable", default="emapper.py", help="Path to emapper.py.")
    run_eggnog_parser.add_argument("--data-dir", type=Path, default=None, help="Path to eggNOG data directory.")
    run_eggnog_parser.add_argument("--cpus", type=int, default=1, help="CPU threads for eggNOG-mapper.")
    run_eggnog_parser.add_argument("--force", action="store_true", help="Rerun even if output already exists.")

    run_domtblout_parser = subparsers.add_parser(
        "run-domtblout",
        help="Run HMMER hmmscan and print normalized annotation-hit TSV.",
    )
    run_domtblout_parser.add_argument("fasta", type=Path, help="Path to protein FASTA.")
    run_domtblout_parser.add_argument("--hmm", type=Path, required=True, help="Path to HMM database.")
    run_domtblout_parser.add_argument("--out", type=Path, required=True, help="Path for domtblout output.")
    run_domtblout_parser.add_argument("--database", required=True, help="Database name, e.g. Pfam or dbCAN.")
    run_domtblout_parser.add_argument("--executable", default="hmmscan", help="Path to hmmscan.")
    run_domtblout_parser.add_argument("--cpus", type=int, default=1, help="CPU threads for HMMER.")
    run_domtblout_parser.add_argument("--evalue", type=float, default=1e-5, help="HMMER E-value cutoff.")
    run_domtblout_parser.add_argument("--force", action="store_true", help="Rerun even if output already exists.")

    validate_splits_parser = subparsers.add_parser(
        "validate-splits",
        help="Validate benchmark label splits for family leakage.",
    )
    validate_splits_parser.add_argument("labels", type=Path, help="Path to benchmark labels TSV.")

    build_features_parser = subparsers.add_parser(
        "build-features",
        help="Build binary annotation feature matrix JSON.",
    )
    build_features_parser.add_argument("annotations", type=Path, help="Path to multi-genome annotation TSV.")
    build_features_parser.add_argument("--out", type=Path, required=True, help="Output feature matrix JSON.")

    import_legacy_labels_parser = subparsers.add_parser(
        "import-legacy-labels",
        help="Convert microbe-foundation traits/splits tables to benchmark label TSV.",
    )
    import_legacy_labels_parser.add_argument("traits", type=Path, help="Legacy traits table, CSV/TSV/parquet.")
    import_legacy_labels_parser.add_argument("splits", type=Path, help="Legacy splits table, CSV/TSV/parquet.")
    import_legacy_labels_parser.add_argument("--out", type=Path, required=True, help="Output benchmark labels TSV.")
    import_legacy_labels_parser.add_argument("--split-column", default="family_split", help="Split column to use.")
    import_legacy_labels_parser.add_argument(
        "--max-multilabel-classes",
        type=int,
        default=50,
        help="Maximum observed classes per legacy multilabel trait.",
    )

    import_legacy_features_parser = subparsers.add_parser(
        "import-legacy-eggnog-features",
        help="Convert cached microbe-foundation eggNOG NPZ features to feature matrix JSON.",
    )
    import_legacy_features_parser.add_argument("npz", type=Path, help="Legacy eggNOG feature NPZ.")
    import_legacy_features_parser.add_argument("vocab", type=Path, help="Legacy eggNOG vocabulary JSON.")
    import_legacy_features_parser.add_argument("--out", type=Path, required=True, help="Output feature matrix JSON.")
    import_legacy_features_parser.add_argument(
        "--labels",
        type=Path,
        default=None,
        help="Optional benchmark labels TSV used to restrict genomes.",
    )
    import_legacy_features_parser.add_argument("--max-features", type=int, default=1000, help="Top features to keep.")
    import_legacy_features_parser.add_argument("--min-prevalence", type=int, default=1, help="Minimum feature prevalence.")
    import_legacy_features_parser.add_argument("--feature-prefix", default="eggNOG", help="Feature namespace prefix.")

    import_legacy_dense_parser = subparsers.add_parser(
        "import-legacy-dense-features",
        help="Convert cached dense embedding NPZ features to median-binarized feature matrix JSON.",
    )
    import_legacy_dense_parser.add_argument("npz", type=Path, help="Legacy dense embedding NPZ.")
    import_legacy_dense_parser.add_argument("--out", type=Path, required=True, help="Output feature matrix JSON.")
    import_legacy_dense_parser.add_argument(
        "--labels",
        type=Path,
        default=None,
        help="Optional benchmark labels TSV used to restrict genomes.",
    )
    import_legacy_dense_parser.add_argument("--max-features", type=int, default=1000, help="High-variance features to keep.")
    import_legacy_dense_parser.add_argument("--feature-prefix", default="embedding", help="Feature namespace prefix.")

    train_parser = subparsers.add_parser(
        "train-baseline",
        help="Train a no-GPU baseline model from features and labels.",
    )
    train_parser.add_argument("features", type=Path, help="Path to feature matrix JSON.")
    train_parser.add_argument("labels", type=Path, help="Path to benchmark labels TSV.")
    train_parser.add_argument("--out", type=Path, required=True, help="Output model JSON.")
    train_parser.add_argument("--train-split", default="train", help="Split to train on.")

    train_dense_parser = subparsers.add_parser(
        "train-dense-baseline",
        help="Train a direct dense-embedding baseline model from a legacy NPZ.",
    )
    train_dense_parser.add_argument("npz", type=Path, help="Dense feature NPZ with bacdive_ids/features.")
    train_dense_parser.add_argument("labels", type=Path, help="Path to benchmark labels TSV.")
    train_dense_parser.add_argument("--out", type=Path, required=True, help="Output dense model JSON.")
    train_dense_parser.add_argument("--train-split", default="train", help="Split to train on.")
    train_dense_parser.add_argument("--max-features", type=int, default=None, help="Optional high-variance dimensions to keep.")

    evaluate_parser = subparsers.add_parser(
        "evaluate",
        help="Evaluate a trained baseline model.",
    )
    evaluate_parser.add_argument("model", type=Path, help="Path to model JSON.")
    evaluate_parser.add_argument("features", type=Path, help="Path to feature matrix JSON.")
    evaluate_parser.add_argument("labels", type=Path, help="Path to benchmark labels TSV.")
    evaluate_parser.add_argument("--split", default="test", help="Split to evaluate.")
    evaluate_parser.add_argument("--out", type=Path, default=None, help="Optional output report JSON.")

    evaluate_dense_parser = subparsers.add_parser(
        "evaluate-dense",
        help="Evaluate a direct dense-embedding baseline model.",
    )
    evaluate_dense_parser.add_argument("model", type=Path, help="Path to dense model JSON.")
    evaluate_dense_parser.add_argument("npz", type=Path, help="Dense feature NPZ with bacdive_ids/features.")
    evaluate_dense_parser.add_argument("labels", type=Path, help="Path to benchmark labels TSV.")
    evaluate_dense_parser.add_argument("--split", default="test", help="Split to evaluate.")
    evaluate_dense_parser.add_argument("--out", type=Path, default=None, help="Optional output report JSON.")

    evaluate_ranking_parser = subparsers.add_parser(
        "evaluate-ranking",
        help="Evaluate top-k candidate ranking against held-out labels.",
    )
    evaluate_ranking_parser.add_argument("model", type=Path, help="Path to model JSON.")
    evaluate_ranking_parser.add_argument("features", type=Path, help="Path to feature matrix JSON.")
    evaluate_ranking_parser.add_argument("labels", type=Path, help="Path to benchmark labels TSV.")
    evaluate_ranking_parser.add_argument("--split", default="test", help="Split to evaluate.")
    ranking_group = evaluate_ranking_parser.add_mutually_exclusive_group(required=True)
    ranking_group.add_argument("--target", help="Target key, e.g. biofuels_industrial:cellulose_degradation.")
    ranking_group.add_argument("--panel", help="Application panel, e.g. biofuels_industrial.")
    evaluate_ranking_parser.add_argument(
        "--k",
        type=int,
        action="append",
        default=None,
        help="Top-k cutoff to evaluate. Repeat for multiple cutoffs.",
    )
    evaluate_ranking_parser.add_argument("--out", type=Path, default=None, help="Optional output report JSON.")

    evaluate_dense_ranking_parser = subparsers.add_parser(
        "evaluate-dense-ranking",
        help="Evaluate top-k candidate ranking for a dense-embedding model.",
    )
    evaluate_dense_ranking_parser.add_argument("model", type=Path, help="Path to dense model JSON.")
    evaluate_dense_ranking_parser.add_argument("npz", type=Path, help="Dense feature NPZ with bacdive_ids/features.")
    evaluate_dense_ranking_parser.add_argument("labels", type=Path, help="Path to benchmark labels TSV.")
    evaluate_dense_ranking_parser.add_argument("--split", default="test", help="Split to evaluate.")
    dense_ranking_group = evaluate_dense_ranking_parser.add_mutually_exclusive_group(required=True)
    dense_ranking_group.add_argument("--target", help="Target key, e.g. biosafety:pathogenicity_human.")
    dense_ranking_group.add_argument("--panel", help="Application panel, e.g. biosafety.")
    evaluate_dense_ranking_parser.add_argument(
        "--k",
        type=int,
        action="append",
        default=None,
        help="Top-k cutoff to evaluate. Repeat for multiple cutoffs.",
    )
    evaluate_dense_ranking_parser.add_argument("--out", type=Path, default=None, help="Optional output report JSON.")

    evaluate_fusion_ranking_parser = subparsers.add_parser(
        "evaluate-fusion-ranking",
        help="Evaluate top-k ranking after fusing annotation and dense model scores.",
    )
    evaluate_fusion_ranking_parser.add_argument("annotation_model", type=Path, help="Path to annotation model JSON.")
    evaluate_fusion_ranking_parser.add_argument("annotation_features", type=Path, help="Path to annotation feature matrix JSON.")
    evaluate_fusion_ranking_parser.add_argument("dense_model", type=Path, help="Path to dense model JSON.")
    evaluate_fusion_ranking_parser.add_argument("labels", type=Path, help="Path to benchmark labels TSV.")
    dense_source_group = evaluate_fusion_ranking_parser.add_mutually_exclusive_group(required=True)
    dense_source_group.add_argument("--dense-npz", type=Path, help="Dense feature NPZ with bacdive_ids/features.")
    dense_source_group.add_argument("--dense-feature-json", help="Dense feature matrix JSON, used by lightweight tests.")
    evaluate_fusion_ranking_parser.add_argument("--split", default="test", help="Split to evaluate.")
    fusion_group = evaluate_fusion_ranking_parser.add_mutually_exclusive_group(required=True)
    fusion_group.add_argument("--target", help="Target key, e.g. biosafety:pathogenicity_human.")
    fusion_group.add_argument("--panel", help="Application panel, e.g. biosafety.")
    evaluate_fusion_ranking_parser.add_argument(
        "--k",
        type=int,
        action="append",
        default=None,
        help="Top-k cutoff to evaluate. Repeat for multiple cutoffs.",
    )
    evaluate_fusion_ranking_parser.add_argument(
        "--weight",
        type=float,
        action="append",
        default=None,
        help="Dense score weight to test. Repeat for a grid.",
    )
    evaluate_fusion_ranking_parser.add_argument(
        "--select-k",
        type=int,
        default=None,
        help="Precision@k metric used to select the best weight.",
    )
    evaluate_fusion_ranking_parser.add_argument("--out", type=Path, default=None, help="Optional output report JSON.")

    leaderboard_parser = subparsers.add_parser(
        "evaluate-fusion-leaderboard",
        help="Select best annotation/dense/fusion source per target on validation and report test metrics.",
    )
    leaderboard_parser.add_argument("annotation_model", type=Path, help="Path to annotation model JSON.")
    leaderboard_parser.add_argument("annotation_features", type=Path, help="Path to annotation feature matrix JSON.")
    leaderboard_parser.add_argument("labels", type=Path, help="Path to benchmark labels TSV.")
    leaderboard_parser.add_argument(
        "--dense-source",
        nargs=3,
        action="append",
        metavar=("NAME", "MODEL", "NPZ"),
        default=None,
        help="Dense source name, dense model JSON, and dense NPZ. Repeat for multiple sources.",
    )
    leaderboard_parser.add_argument(
        "--dense-feature-json-source",
        nargs=3,
        action="append",
        metavar=("NAME", "MODEL", "JSON"),
        default=None,
        help="Dense source name, dense model JSON, and dense feature JSON.",
    )
    leaderboard_parser.add_argument("--validation-split", default="val", help="Split used to select source/weight.")
    leaderboard_parser.add_argument("--test-split", default="test", help="Split used for final held-out metrics.")
    leaderboard_parser.add_argument(
        "--k",
        type=int,
        action="append",
        default=None,
        help="Top-k cutoff to evaluate. Repeat for multiple cutoffs.",
    )
    leaderboard_parser.add_argument(
        "--weight",
        type=float,
        action="append",
        default=None,
        help="Dense score weight to test for fusion. Repeat for a grid.",
    )
    leaderboard_parser.add_argument("--select-k", type=int, default=None, help="Precision@k metric used for selection.")
    leaderboard_parser.add_argument("--out", type=Path, default=None, help="Optional output leaderboard JSON.")

    discovery_parser = subparsers.add_parser(
        "rank-discovery-candidates",
        help="Generate candidate rankings for targets selected in a fusion leaderboard.",
    )
    discovery_parser.add_argument("annotation_model", type=Path, help="Path to annotation model JSON.")
    discovery_parser.add_argument("annotation_features", type=Path, help="Path to annotation feature matrix JSON.")
    discovery_parser.add_argument("labels", type=Path, help="Path to benchmark labels TSV.")
    discovery_parser.add_argument("leaderboard", type=Path, help="Path to target leaderboard JSON.")
    discovery_parser.add_argument(
        "--dense-source",
        nargs=3,
        action="append",
        metavar=("NAME", "MODEL", "NPZ"),
        default=None,
        help="Dense source name, dense model JSON, and dense NPZ. Repeat for multiple sources.",
    )
    discovery_parser.add_argument("--split", default="test", help="Split to rank candidates from.")
    discovery_parser.add_argument("--limit-per-target", type=int, default=10, help="Maximum candidates per target.")
    discovery_parser.add_argument("--min-precision", type=float, default=0.0, help="Minimum leaderboard precision to include.")
    discovery_parser.add_argument("--precision-k", type=int, default=10, help="Precision@k field used for filtering.")
    discovery_parser.add_argument("--max-targets", type=int, default=None, help="Maximum targets to include.")
    discovery_parser.add_argument("--out", type=Path, default=None, help="Optional output candidates JSON.")

    annotate_candidates_parser = subparsers.add_parser(
        "annotate-discovery-candidates",
        help="Attach taxonomy, accessions, evidence, and biosafety flags to discovery candidate rankings.",
    )
    annotate_candidates_parser.add_argument("candidates", type=Path, help="Path to discovery candidates JSON.")
    annotate_candidates_parser.add_argument("traits", type=Path, help="Legacy traits table, CSV/TSV/parquet.")
    annotate_candidates_parser.add_argument("model", type=Path, help="Path to annotation model JSON.")
    annotate_candidates_parser.add_argument("features", type=Path, help="Path to annotation feature matrix JSON.")
    annotate_candidates_parser.add_argument("labels", type=Path, help="Path to benchmark labels TSV.")
    annotate_candidates_parser.add_argument("--accessions", type=Path, default=None, help="Optional genome accessions TSV/CSV/parquet.")
    annotate_candidates_parser.add_argument("--evidence-limit", type=int, default=5, help="Maximum evidence rows per candidate.")
    annotate_candidates_parser.add_argument("--out", type=Path, required=True, help="Output annotated candidate JSON.")
    annotate_candidates_parser.add_argument("--report-out", type=Path, default=None, help="Optional Markdown report output.")

    predict_baseline_parser = subparsers.add_parser(
        "predict-baseline",
        help="Emit product-style prediction JSON from a trained baseline model.",
    )
    predict_baseline_parser.add_argument("model", type=Path, help="Path to model JSON.")
    predict_baseline_parser.add_argument("features", type=Path, help="Path to feature matrix JSON.")
    predict_baseline_parser.add_argument("--genome-id", required=True, help="Genome id to predict.")

    rank_parser = subparsers.add_parser(
        "rank-candidates",
        help="Rank genomes by a trained baseline target or application panel.",
    )
    rank_parser.add_argument("model", type=Path, help="Path to model JSON.")
    rank_parser.add_argument("features", type=Path, help="Path to feature matrix JSON.")
    group = rank_parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--target", help="Target key, e.g. biofuels_industrial:cellulose_degradation.")
    group.add_argument("--panel", help="Application panel, e.g. biofuels_industrial.")
    rank_parser.add_argument("--limit", type=int, default=None, help="Maximum candidates to return.")

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
    if args.command == "run-eggnog":
        return _run_eggnog(args)
    if args.command == "run-domtblout":
        return _run_domtblout(args)
    if args.command == "validate-splits":
        return _validate_splits(args.labels)
    if args.command == "build-features":
        return _build_features(args.annotations, args.out)
    if args.command == "import-legacy-labels":
        return _import_legacy_labels(
            args.traits,
            args.splits,
            args.out,
            args.split_column,
            args.max_multilabel_classes,
        )
    if args.command == "import-legacy-eggnog-features":
        return _import_legacy_eggnog_features(
            args.npz,
            args.vocab,
            args.out,
            args.labels,
            args.max_features,
            args.min_prevalence,
            args.feature_prefix,
        )
    if args.command == "import-legacy-dense-features":
        return _import_legacy_dense_features(
            args.npz,
            args.out,
            args.labels,
            args.max_features,
            args.feature_prefix,
        )
    if args.command == "train-baseline":
        return _train_baseline(args.features, args.labels, args.out, args.train_split)
    if args.command == "train-dense-baseline":
        return _train_dense_baseline(args.npz, args.labels, args.out, args.train_split, args.max_features)
    if args.command == "evaluate":
        return _evaluate(args.model, args.features, args.labels, args.split, args.out)
    if args.command == "evaluate-dense":
        return _evaluate_dense(args.model, args.npz, args.labels, args.split, args.out)
    if args.command == "evaluate-ranking":
        return _evaluate_ranking(args.model, args.features, args.labels, args.split, args.target, args.panel, args.k, args.out)
    if args.command == "evaluate-dense-ranking":
        return _evaluate_dense_ranking(args.model, args.npz, args.labels, args.split, args.target, args.panel, args.k, args.out)
    if args.command == "evaluate-fusion-ranking":
        return _evaluate_fusion_ranking(
            args.annotation_model,
            args.annotation_features,
            args.dense_model,
            args.labels,
            args.dense_npz,
            args.dense_feature_json,
            args.split,
            args.target,
            args.panel,
            args.k,
            args.weight,
            args.select_k,
            args.out,
        )
    if args.command == "evaluate-fusion-leaderboard":
        return _evaluate_fusion_leaderboard(
            args.annotation_model,
            args.annotation_features,
            args.labels,
            args.dense_source,
            args.dense_feature_json_source,
            args.validation_split,
            args.test_split,
            args.k,
            args.weight,
            args.select_k,
            args.out,
        )
    if args.command == "rank-discovery-candidates":
        return _rank_discovery_candidates(
            args.annotation_model,
            args.annotation_features,
            args.labels,
            args.leaderboard,
            args.dense_source,
            args.split,
            args.limit_per_target,
            args.min_precision,
            args.precision_k,
            args.max_targets,
            args.out,
        )
    if args.command == "annotate-discovery-candidates":
        return _annotate_discovery_candidates(
            args.candidates,
            args.traits,
            args.model,
            args.features,
            args.labels,
            args.accessions,
            args.evidence_limit,
            args.out,
            args.report_out,
        )
    if args.command == "predict-baseline":
        return _predict_baseline(args.model, args.features, args.genome_id)
    if args.command == "rank-candidates":
        return _rank_candidates(args.model, args.features, args.target, args.panel, args.limit)

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


def _run_eggnog(args: argparse.Namespace) -> int:
    try:
        output_path = run_eggnog_mapper(
            input_fasta=args.fasta,
            output_dir=args.output_dir,
            executable=args.executable,
            data_dir=args.data_dir,
            cpus=args.cpus,
            force=args.force,
        )
        hits = parse_eggnog_mapper(output_path.read_text())
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    except ValueError as exc:
        print(f"invalid eggNOG annotations: {exc}", file=sys.stderr)
        return 1
    except RuntimeError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    print(format_annotation_hits_tsv(hits), end="")
    return 0


def _run_domtblout(args: argparse.Namespace) -> int:
    try:
        output_path = run_hmmer_domtblout(
            input_fasta=args.fasta,
            hmm_database=args.hmm,
            output_path=args.out,
            database=args.database,
            executable=args.executable,
            cpus=args.cpus,
            evalue=args.evalue,
            force=args.force,
        )
        hits = parse_hmmer_domtblout(output_path.read_text(), database=args.database)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    except ValueError as exc:
        print(f"invalid domtblout: {exc}", file=sys.stderr)
        return 1
    except RuntimeError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    print(format_annotation_hits_tsv(hits), end="")
    return 0


def _validate_splits(labels_path: Path) -> int:
    try:
        labels = parse_labels_tsv(labels_path.read_text())
        validate_family_holdout(labels)
    except FileNotFoundError:
        print(f"{labels_path}: file not found", file=sys.stderr)
        return 1
    except (ValueError, FamilyLeakageError) as exc:
        print(f"{labels_path}: invalid split file: {exc}", file=sys.stderr)
        return 1

    print(f"{labels_path}: family holdout valid")
    return 0


def _build_features(annotations_path: Path, out_path: Path) -> int:
    try:
        features = build_feature_matrix_from_annotation_tsv(annotations_path.read_text())
    except FileNotFoundError:
        print(f"{annotations_path}: file not found", file=sys.stderr)
        return 1
    except ValueError as exc:
        print(f"{annotations_path}: invalid annotation TSV: {exc}", file=sys.stderr)
        return 1

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(features.to_json())
    print(f"wrote {len(features.genome_ids)} genomes x {len(features.feature_names)} features to {out_path}")
    return 0


def _import_legacy_labels(
    traits_path: Path,
    splits_path: Path,
    out_path: Path,
    split_column: str,
    max_multilabel_classes: int,
) -> int:
    try:
        trait_rows = load_table_records(traits_path)
        split_rows = load_table_records(splits_path)
        labels_text = labels_tsv_from_legacy_tables(
            trait_rows,
            split_rows,
            split_column=split_column,
            max_multilabel_classes=max_multilabel_classes,
        )
        records = parse_labels_tsv(labels_text)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    except (RuntimeError, ValueError) as exc:
        print(f"cannot import legacy labels: {exc}", file=sys.stderr)
        return 1

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(labels_text)
    targets = {record.target_key for record in records}
    genomes = {record.genome_id for record in records}
    print(f"wrote {len(records)} labels for {len(genomes)} genomes and {len(targets)} targets to {out_path}")
    return 0


def _import_legacy_eggnog_features(
    npz_path: Path,
    vocab_path: Path,
    out_path: Path,
    labels_path: Path | None,
    max_features: int,
    min_prevalence: int,
    feature_prefix: str,
) -> int:
    try:
        keep_genome_ids = None
        if labels_path is not None:
            keep_genome_ids = {record.genome_id for record in parse_labels_tsv(labels_path.read_text())}
        features = load_legacy_npz_feature_matrix(
            npz_path,
            vocab_path,
            max_features=max_features,
            min_prevalence=min_prevalence,
            feature_prefix=feature_prefix,
            keep_genome_ids=keep_genome_ids,
        )
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    except (RuntimeError, ValueError, KeyError) as exc:
        print(f"cannot import legacy eggNOG features: {exc}", file=sys.stderr)
        return 1

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(features.to_json())
    print(f"wrote {len(features.genome_ids)} genomes x {len(features.feature_names)} features to {out_path}")
    return 0


def _import_legacy_dense_features(
    npz_path: Path,
    out_path: Path,
    labels_path: Path | None,
    max_features: int,
    feature_prefix: str,
) -> int:
    try:
        keep_genome_ids = None
        if labels_path is not None:
            keep_genome_ids = {record.genome_id for record in parse_labels_tsv(labels_path.read_text())}
        features = load_legacy_dense_npz_feature_matrix(
            npz_path,
            max_features=max_features,
            feature_prefix=feature_prefix,
            keep_genome_ids=keep_genome_ids,
        )
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    except (RuntimeError, ValueError, KeyError) as exc:
        print(f"cannot import legacy dense features: {exc}", file=sys.stderr)
        return 1

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(features.to_json())
    print(f"wrote {len(features.genome_ids)} genomes x {len(features.feature_names)} features to {out_path}")
    return 0


def _train_baseline(features_path: Path, labels_path: Path, out_path: Path, train_split: str) -> int:
    try:
        features = FeatureMatrix.from_json(features_path.read_text())
        labels = parse_labels_tsv(labels_path.read_text())
        validate_family_holdout(labels)
        model = train_baseline(features, labels, train_split=train_split)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    except (ValueError, FamilyLeakageError, KeyError) as exc:
        print(f"cannot train baseline: {exc}", file=sys.stderr)
        return 1

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(model.to_json())
    print(f"wrote {len(model.targets)} target models to {out_path}")
    return 0


def _train_dense_baseline(
    npz_path: Path,
    labels_path: Path,
    out_path: Path,
    train_split: str,
    max_features: int | None,
) -> int:
    try:
        labels = parse_labels_tsv(labels_path.read_text())
        keep_genome_ids = {record.genome_id for record in labels}
        features = load_dense_npz_feature_matrix(npz_path, keep_genome_ids=keep_genome_ids, max_features=max_features)
        model = train_dense_baseline(features, labels, train_split=train_split)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    except (RuntimeError, ValueError, KeyError) as exc:
        print(f"cannot train dense baseline: {exc}", file=sys.stderr)
        return 1

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(model.to_json())
    print(
        f"wrote {len(model.targets)} dense target models using "
        f"{len(model.feature_indexes)} dimensions to {out_path}"
    )
    return 0


def _evaluate(model_path: Path, features_path: Path, labels_path: Path, split: str, out_path: Path | None) -> int:
    try:
        model = BaselineModel.from_json(model_path.read_text())
        features = FeatureMatrix.from_json(features_path.read_text())
        labels = parse_labels_tsv(labels_path.read_text())
        report = evaluate_model(model, features, labels, split=split)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    except (ValueError, KeyError) as exc:
        print(f"cannot evaluate baseline: {exc}", file=sys.stderr)
        return 1

    text = json.dumps(report, indent=2, sort_keys=True)
    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(text)
        print(f"wrote evaluation report to {out_path}")
    else:
        print(text)
    return 0


def _evaluate_dense(
    model_path: Path,
    npz_path: Path,
    labels_path: Path,
    split: str,
    out_path: Path | None,
) -> int:
    try:
        model = DenseBaselineModel.from_json(model_path.read_text())
        labels = parse_labels_tsv(labels_path.read_text())
        keep_genome_ids = {record.genome_id for record in labels}
        features = load_dense_npz_feature_matrix(
            npz_path,
            keep_genome_ids=keep_genome_ids,
            feature_indexes=model.feature_indexes,
        )
        report = evaluate_dense_model(model, features, labels, split=split)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    except (RuntimeError, ValueError, KeyError) as exc:
        print(f"cannot evaluate dense baseline: {exc}", file=sys.stderr)
        return 1

    text = json.dumps(report, indent=2, sort_keys=True)
    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(text)
        print(f"wrote dense evaluation report to {out_path}")
    else:
        print(text)
    return 0


def _evaluate_ranking(
    model_path: Path,
    features_path: Path,
    labels_path: Path,
    split: str,
    target_key: str | None,
    panel: str | None,
    ks: list[int] | None,
    out_path: Path | None,
) -> int:
    try:
        model = BaselineModel.from_json(model_path.read_text())
        features = FeatureMatrix.from_json(features_path.read_text())
        labels = parse_labels_tsv(labels_path.read_text())
        report = evaluate_ranking(
            model,
            features,
            labels,
            split=split,
            target_key=target_key,
            panel=panel,
            ks=ks,
        )
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    except (ValueError, KeyError) as exc:
        print(f"cannot evaluate ranking: {exc}", file=sys.stderr)
        return 1

    text = json.dumps(report, indent=2, sort_keys=True)
    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(text)
        print(f"wrote ranking evaluation report to {out_path}")
    else:
        print(text)
    return 0


def _evaluate_dense_ranking(
    model_path: Path,
    npz_path: Path,
    labels_path: Path,
    split: str,
    target_key: str | None,
    panel: str | None,
    ks: list[int] | None,
    out_path: Path | None,
) -> int:
    try:
        model = DenseBaselineModel.from_json(model_path.read_text())
        labels = parse_labels_tsv(labels_path.read_text())
        keep_genome_ids = {record.genome_id for record in labels}
        features = load_dense_npz_feature_matrix(
            npz_path,
            keep_genome_ids=keep_genome_ids,
            feature_indexes=model.feature_indexes,
        )
        report = evaluate_dense_ranking(
            model,
            features,
            labels,
            split=split,
            target_key=target_key,
            panel=panel,
            ks=ks,
        )
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    except (RuntimeError, ValueError, KeyError) as exc:
        print(f"cannot evaluate dense ranking: {exc}", file=sys.stderr)
        return 1

    text = json.dumps(report, indent=2, sort_keys=True)
    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(text)
        print(f"wrote dense ranking evaluation report to {out_path}")
    else:
        print(text)
    return 0


def _evaluate_fusion_ranking(
    annotation_model_path: Path,
    annotation_features_path: Path,
    dense_model_path: Path,
    labels_path: Path,
    dense_npz_path: Path | None,
    dense_feature_json: str | None,
    split: str,
    target_key: str | None,
    panel: str | None,
    ks: list[int] | None,
    weights: list[float] | None,
    select_k: int | None,
    out_path: Path | None,
) -> int:
    try:
        annotation_model = BaselineModel.from_json(annotation_model_path.read_text())
        annotation_features = FeatureMatrix.from_json(annotation_features_path.read_text())
        dense_model = DenseBaselineModel.from_json(dense_model_path.read_text())
        labels = parse_labels_tsv(labels_path.read_text())
        if dense_npz_path is not None:
            keep_genome_ids = {record.genome_id for record in labels}
            dense_features = load_dense_npz_feature_matrix(
                dense_npz_path,
                keep_genome_ids=keep_genome_ids,
                feature_indexes=dense_model.feature_indexes,
            )
        else:
            dense_features = _dense_features_from_json(dense_feature_json or "")
        report = evaluate_fusion_ranking(
            annotation_model,
            annotation_features,
            dense_model,
            dense_features,
            labels,
            split=split,
            target_key=target_key,
            panel=panel,
            ks=ks,
            weights=weights,
            select_k=select_k,
        )
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    except (RuntimeError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print(f"cannot evaluate fusion ranking: {exc}", file=sys.stderr)
        return 1

    text = json.dumps(report, indent=2, sort_keys=True)
    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(text)
        print(f"wrote fusion ranking evaluation report to {out_path}")
    else:
        print(text)
    return 0


def _evaluate_fusion_leaderboard(
    annotation_model_path: Path,
    annotation_features_path: Path,
    labels_path: Path,
    dense_source_specs: list[list[str]] | None,
    dense_json_source_specs: list[list[str]] | None,
    validation_split: str,
    test_split: str,
    ks: list[int] | None,
    weights: list[float] | None,
    select_k: int | None,
    out_path: Path | None,
) -> int:
    try:
        annotation_model = BaselineModel.from_json(annotation_model_path.read_text())
        annotation_features = FeatureMatrix.from_json(annotation_features_path.read_text())
        labels = parse_labels_tsv(labels_path.read_text())
        dense_sources = _load_dense_sources(labels, dense_source_specs, dense_json_source_specs)
        report = evaluate_fusion_leaderboard(
            annotation_model,
            annotation_features,
            dense_sources,
            labels,
            validation_split=validation_split,
            test_split=test_split,
            ks=ks,
            weights=weights,
            select_k=select_k,
        )
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    except (RuntimeError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print(f"cannot evaluate fusion leaderboard: {exc}", file=sys.stderr)
        return 1

    _write_or_print_json(report, out_path, f"wrote fusion leaderboard to {out_path}" if out_path else "")
    return 0


def _rank_discovery_candidates(
    annotation_model_path: Path,
    annotation_features_path: Path,
    labels_path: Path,
    leaderboard_path: Path,
    dense_source_specs: list[list[str]] | None,
    split: str,
    limit_per_target: int,
    min_precision: float,
    precision_k: int,
    max_targets: int | None,
    out_path: Path | None,
) -> int:
    try:
        annotation_model = BaselineModel.from_json(annotation_model_path.read_text())
        annotation_features = FeatureMatrix.from_json(annotation_features_path.read_text())
        labels = parse_labels_tsv(labels_path.read_text())
        leaderboard = json.loads(leaderboard_path.read_text())
        dense_sources = _load_dense_sources(labels, dense_source_specs, None)
        report = rank_discovery_candidates(
            leaderboard,
            annotation_model,
            annotation_features,
            dense_sources,
            labels,
            split=split,
            limit_per_target=limit_per_target,
            min_precision=min_precision,
            precision_k=precision_k,
            max_targets=max_targets,
        )
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    except (RuntimeError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print(f"cannot rank discovery candidates: {exc}", file=sys.stderr)
        return 1

    _write_or_print_json(report, out_path, f"wrote discovery candidate rankings to {out_path}" if out_path else "")
    return 0


def _annotate_discovery_candidates(
    candidates_path: Path,
    traits_path: Path,
    model_path: Path,
    features_path: Path,
    labels_path: Path,
    accessions_path: Path | None,
    evidence_limit: int,
    out_path: Path,
    report_out_path: Path | None,
) -> int:
    try:
        candidates = json.loads(candidates_path.read_text())
        trait_rows = load_table_records(traits_path)
        accession_rows = load_table_records(accessions_path) if accessions_path is not None else []
        model = BaselineModel.from_json(model_path.read_text())
        features = FeatureMatrix.from_json(features_path.read_text())
        labels = parse_labels_tsv(labels_path.read_text())
        annotated = annotate_discovery_candidates(
            candidates,
            trait_rows,
            accession_rows,
            model,
            features,
            labels,
            evidence_limit=evidence_limit,
        )
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    except (RuntimeError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print(f"cannot annotate discovery candidates: {exc}", file=sys.stderr)
        return 1

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(to_json(annotated))
    print(f"wrote annotated discovery candidates to {out_path}")
    if report_out_path is not None:
        report_out_path.parent.mkdir(parents=True, exist_ok=True)
        report_out_path.write_text(render_discovery_candidate_report(annotated))
        print(f"wrote discovery candidate report to {report_out_path}")
    return 0


def _dense_features_from_json(text: str) -> DenseFeatureMatrix:
    data = json.loads(text)
    return DenseFeatureMatrix(
        genome_ids=[str(value) for value in data["genome_ids"]],
        feature_indexes=[int(value) for value in data["feature_indexes"]],
        rows=[[float(value) for value in row] for row in data["rows"]],
    )


def _load_dense_sources(
    labels: list,
    dense_source_specs: list[list[str]] | None,
    dense_json_source_specs: list[list[str]] | None,
) -> dict[str, tuple[DenseBaselineModel, DenseFeatureMatrix]]:
    dense_sources: dict[str, tuple[DenseBaselineModel, DenseFeatureMatrix]] = {}
    keep_genome_ids = {record.genome_id for record in labels}
    for name, model_path, npz_path in dense_source_specs or []:
        dense_model = DenseBaselineModel.from_json(Path(model_path).read_text())
        dense_features = load_dense_npz_feature_matrix(
            Path(npz_path),
            keep_genome_ids=keep_genome_ids,
            feature_indexes=dense_model.feature_indexes,
        )
        dense_sources[name] = (dense_model, dense_features)
    for name, model_path, feature_json in dense_json_source_specs or []:
        dense_model = DenseBaselineModel.from_json(Path(model_path).read_text())
        dense_sources[name] = (dense_model, _dense_features_from_json(feature_json))
    if not dense_sources:
        raise ValueError("at least one dense source is required")
    return dense_sources


def _write_or_print_json(report: dict, out_path: Path | None, message: str) -> None:
    text = json.dumps(report, indent=2, sort_keys=True)
    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(text)
        print(message)
    else:
        print(text)


def _predict_baseline(model_path: Path, features_path: Path, genome_id: str) -> int:
    try:
        model = BaselineModel.from_json(model_path.read_text())
        features = FeatureMatrix.from_json(features_path.read_text())
        prediction = predict_model(model, features, genome_id=genome_id)
        validate_prediction(prediction)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    except (ValueError, KeyError, PredictionValidationError) as exc:
        print(f"cannot predict with baseline: {exc}", file=sys.stderr)
        return 1

    print(json.dumps(prediction, indent=2, sort_keys=True))
    return 0


def _rank_candidates(
    model_path: Path,
    features_path: Path,
    target_key: str | None,
    panel: str | None,
    limit: int | None,
) -> int:
    try:
        model = BaselineModel.from_json(model_path.read_text())
        features = FeatureMatrix.from_json(features_path.read_text())
        ranking = rank_candidates(model, features, target_key=target_key, panel=panel, limit=limit)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    except (ValueError, KeyError) as exc:
        print(f"cannot rank candidates: {exc}", file=sys.stderr)
        return 1

    print(json.dumps(ranking, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
