"""Late-fusion ranking for annotation and dense embedding baselines."""

from __future__ import annotations

from typing import Any

from microbial_function_discovery.datasets import LabelRecord
from microbial_function_discovery.dense_learning import DenseBaselineModel, DenseFeatureMatrix
from microbial_function_discovery.features import FeatureMatrix
from microbial_function_discovery.learning import BaselineModel
from microbial_function_discovery.panels import APPLICATION_AREAS


def evaluate_fusion_ranking(
    annotation_model: BaselineModel,
    annotation_features: FeatureMatrix,
    dense_model: DenseBaselineModel,
    dense_features: DenseFeatureMatrix,
    labels: list[LabelRecord],
    *,
    split: str,
    target_key: str | None = None,
    panel: str | None = None,
    ks: list[int] | None = None,
    weights: list[float] | None = None,
    select_k: int | None = None,
) -> dict[str, Any]:
    """Grid-search dense-score weights for late-fused top-k ranking."""

    if (target_key is None) == (panel is None):
        raise ValueError("provide exactly one of target_key or panel")
    if target_key is not None:
        _require_target(annotation_model.targets, target_key, "annotation model")
        _require_target(dense_model.targets, target_key, "dense model")
    if panel is not None:
        if panel not in APPLICATION_AREAS:
            raise ValueError(f"unknown panel: {panel}")
        _require_panel_targets(annotation_model.targets, dense_model.targets, panel)

    cutoffs = _normalize_cutoffs(ks or [1, 5, 10])
    dense_weights = _normalize_weights(weights or [0.0, 0.25, 0.5, 0.75, 1.0])
    selection_k = select_k if select_k is not None else (10 if 10 in cutoffs else min(cutoffs))
    if selection_k not in cutoffs:
        cutoffs = _normalize_cutoffs([*cutoffs, selection_k])
    selection_metric = f"precision_at_{selection_k}"
    common_genomes = set(annotation_features.genome_ids) & set(dense_features.genome_ids)
    truth_by_genome = _ranking_truth(
        labels,
        split=split,
        target_key=target_key,
        panel=panel,
        feature_genomes=common_genomes,
    )
    if not truth_by_genome:
        raise ValueError("no labeled candidates found in both feature matrices")

    per_weight = [
        _ranking_report_for_weight(
            dense_weight,
            annotation_model,
            annotation_features,
            dense_model,
            dense_features,
            truth_by_genome,
            target_key=target_key,
            panel=panel,
            cutoffs=cutoffs,
        )
        for dense_weight in dense_weights
    ]
    best = max(
        per_weight,
        key=lambda report: (
            report["metrics"][selection_metric],
            report["metrics"].get(f"recall_at_{selection_k}", 0.0),
            -report["weight"],
        ),
    )

    return {
        "mode": "target" if target_key is not None else "panel",
        "query": target_key if target_key is not None else panel,
        "split": split,
        "weights": dense_weights,
        "selection_metric": selection_metric,
        "best_weight": best["weight"],
        "best": _without_weight(best),
        "per_weight": [
            {
                "weight": report["weight"],
                "n_labeled_candidates": report["n_labeled_candidates"],
                "n_positives": report["n_positives"],
                "metrics": report["metrics"],
            }
            for report in per_weight
        ],
    }


def _ranking_report_for_weight(
    dense_weight: float,
    annotation_model: BaselineModel,
    annotation_features: FeatureMatrix,
    dense_model: DenseBaselineModel,
    dense_features: DenseFeatureMatrix,
    truth_by_genome: dict[str, int],
    *,
    target_key: str | None,
    panel: str | None,
    cutoffs: list[int],
) -> dict[str, Any]:
    annotation_weight = 1.0 - dense_weight
    ranked = []
    for genome_id, truth in truth_by_genome.items():
        annotation_score = _annotation_score(
            annotation_model,
            annotation_features,
            genome_id,
            target_key=target_key,
            panel=panel,
        )
        dense_score = _dense_score(
            dense_model,
            dense_features,
            genome_id,
            target_key=target_key,
            panel=panel,
        )
        fused_score = annotation_weight * annotation_score + dense_weight * dense_score
        ranked.append(
            {
                "genome_id": genome_id,
                "truth": truth,
                "score": round(fused_score, 6),
                "annotation_score": round(annotation_score, 6),
                "dense_score": round(dense_score, 6),
            }
        )

    ranked.sort(key=lambda row: (-row["score"], row["genome_id"]))
    for rank, row in enumerate(ranked, start=1):
        row["rank"] = rank

    n_positives = sum(row["truth"] for row in ranked)
    metrics: dict[str, float | int] = {}
    for k in cutoffs:
        top_k = ranked[:k]
        hits = sum(row["truth"] for row in top_k)
        denominator = min(k, len(ranked))
        metrics[f"hits_at_{k}"] = hits
        metrics[f"precision_at_{k}"] = hits / denominator if denominator else 0.0
        metrics[f"recall_at_{k}"] = hits / n_positives if n_positives else 0.0

    return {
        "weight": dense_weight,
        "n_labeled_candidates": len(ranked),
        "n_positives": n_positives,
        "metrics": metrics,
        "ranked_candidates": ranked,
    }


def _ranking_truth(
    labels: list[LabelRecord],
    *,
    split: str,
    target_key: str | None,
    panel: str | None,
    feature_genomes: set[str],
) -> dict[str, int]:
    truth_by_genome: dict[str, int] = {}
    for record in labels:
        if record.split != split or record.genome_id not in feature_genomes:
            continue
        if target_key is not None and record.target_key != target_key:
            continue
        if panel is not None and record.panel != panel:
            continue
        truth_by_genome[record.genome_id] = max(truth_by_genome.get(record.genome_id, 0), record.value)
    return truth_by_genome


def _annotation_score(
    model: BaselineModel,
    features: FeatureMatrix,
    genome_id: str,
    *,
    target_key: str | None,
    panel: str | None,
) -> float:
    if target_key is not None:
        return model.predict_proba(target_key, genome_id, features)
    target_keys = [key for key in model.targets if key.startswith(f"{panel}:")]
    return max(model.predict_proba(key, genome_id, features) for key in target_keys)


def _dense_score(
    model: DenseBaselineModel,
    features: DenseFeatureMatrix,
    genome_id: str,
    *,
    target_key: str | None,
    panel: str | None,
) -> float:
    if target_key is not None:
        return model.predict_proba(target_key, genome_id, features)
    target_keys = [key for key in model.targets if key.startswith(f"{panel}:")]
    return max(model.predict_proba(key, genome_id, features) for key in target_keys)


def _require_target(targets: dict[str, Any], target_key: str, owner: str) -> None:
    if target_key not in targets:
        raise KeyError(f"target not found in {owner}: {target_key}")


def _require_panel_targets(annotation_targets: dict[str, Any], dense_targets: dict[str, Any], panel: str) -> None:
    annotation_panel_targets = {key for key in annotation_targets if key.startswith(f"{panel}:")}
    dense_panel_targets = {key for key in dense_targets if key.startswith(f"{panel}:")}
    if not annotation_panel_targets:
        raise KeyError(f"panel has no targets in annotation model: {panel}")
    if not dense_panel_targets:
        raise KeyError(f"panel has no targets in dense model: {panel}")


def _normalize_cutoffs(ks: list[int]) -> list[int]:
    cutoffs = sorted(set(ks))
    if not cutoffs or any(k <= 0 for k in cutoffs):
        raise ValueError("ranking cutoffs must be positive integers")
    return cutoffs


def _normalize_weights(weights: list[float]) -> list[float]:
    normalized = sorted(set(float(weight) for weight in weights))
    if not normalized or any(weight < 0.0 or weight > 1.0 for weight in normalized):
        raise ValueError("fusion weights must be between 0 and 1")
    return normalized


def _without_weight(report: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in report.items() if key != "weight"}
