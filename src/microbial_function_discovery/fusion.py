"""Late-fusion ranking for annotation and dense embedding baselines."""

from __future__ import annotations

from typing import Any

from microbial_function_discovery.datasets import LabelRecord
from microbial_function_discovery.dense_learning import (
    DenseBaselineModel,
    DenseFeatureMatrix,
    evaluate_dense_ranking,
)
from microbial_function_discovery.features import FeatureMatrix
from microbial_function_discovery.learning import BaselineModel, evaluate_ranking
from microbial_function_discovery.panels import APPLICATION_AREAS


DenseSourceMap = dict[str, tuple[DenseBaselineModel, DenseFeatureMatrix]]


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


def evaluate_fusion_leaderboard(
    annotation_model: BaselineModel,
    annotation_features: FeatureMatrix,
    dense_sources: DenseSourceMap,
    labels: list[LabelRecord],
    *,
    validation_split: str = "val",
    test_split: str = "test",
    ks: list[int] | None = None,
    weights: list[float] | None = None,
    select_k: int | None = None,
) -> dict[str, Any]:
    """Select the best source per target on validation and report held-out test metrics."""

    cutoffs = _normalize_cutoffs(ks or [10, 50])
    selection_k = select_k if select_k is not None else (10 if 10 in cutoffs else min(cutoffs))
    if selection_k not in cutoffs:
        cutoffs = _normalize_cutoffs([*cutoffs, selection_k])
    selection_metric = f"precision_at_{selection_k}"
    target_rows = []
    for target_key in _target_keys(labels):
        selected = _select_target_source(
            target_key,
            annotation_model,
            annotation_features,
            dense_sources,
            labels,
            validation_split=validation_split,
            test_split=test_split,
            cutoffs=cutoffs,
            weights=weights,
            select_k=selection_k,
            selection_metric=selection_metric,
        )
        if selected is not None:
            target_rows.append(selected)

    target_rows.sort(
        key=lambda row: (
            -float(row["test_metrics"].get(selection_metric, 0.0)),
            -float(row["test_metrics"].get(f"recall_at_{selection_k}", 0.0)),
            row["target_key"],
        )
    )
    return {
        "validation_split": validation_split,
        "test_split": test_split,
        "selection_metric": selection_metric,
        "cutoffs": cutoffs,
        "weights": _normalize_weights(weights or [0.0, 0.25, 0.5, 0.75, 1.0]),
        "dense_sources": sorted(dense_sources),
        "n_targets": len(target_rows),
        "targets": target_rows,
        "panels": _panel_summary(target_rows, selection_metric),
    }


def rank_discovery_candidates(
    leaderboard: dict[str, Any],
    annotation_model: BaselineModel,
    annotation_features: FeatureMatrix,
    dense_sources: DenseSourceMap,
    labels: list[LabelRecord],
    *,
    split: str = "test",
    limit_per_target: int = 10,
    min_precision: float = 0.0,
    precision_k: int = 10,
    max_targets: int | None = None,
) -> dict[str, Any]:
    """Emit candidate rankings for targets that pass a leaderboard precision threshold."""

    precision_metric = f"precision_at_{precision_k}"
    eligible = [
        row
        for row in leaderboard.get("targets", [])
        if float(row.get("test_metrics", {}).get(precision_metric, 0.0)) >= min_precision
    ]
    eligible.sort(
        key=lambda row: (
            -float(row.get("test_metrics", {}).get(precision_metric, 0.0)),
            row["target_key"],
        )
    )
    if max_targets is not None:
        eligible = eligible[:max_targets]

    target_reports = [
        _rank_from_leaderboard_row(
            row,
            annotation_model,
            annotation_features,
            dense_sources,
            labels,
            split=split,
            limit=limit_per_target,
        )
        for row in eligible
    ]
    return {
        "split": split,
        "precision_metric": precision_metric,
        "min_precision": min_precision,
        "limit_per_target": limit_per_target,
        "n_targets": len(target_reports),
        "targets": target_reports,
        "panels": _candidate_panel_summary(target_reports),
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


def _select_target_source(
    target_key: str,
    annotation_model: BaselineModel,
    annotation_features: FeatureMatrix,
    dense_sources: DenseSourceMap,
    labels: list[LabelRecord],
    *,
    validation_split: str,
    test_split: str,
    cutoffs: list[int],
    weights: list[float] | None,
    select_k: int,
    selection_metric: str,
) -> dict[str, Any] | None:
    candidates: list[dict[str, Any]] = []
    if target_key in annotation_model.targets:
        annotation_candidate = _source_candidate(
            "annotation",
            None,
            _try_annotation_ranking(annotation_model, annotation_features, labels, split=validation_split, target_key=target_key, ks=cutoffs),
            _try_annotation_ranking(annotation_model, annotation_features, labels, split=test_split, target_key=target_key, ks=cutoffs),
            complexity=0,
            weight=0.0,
        )
        if annotation_candidate is not None:
            candidates.append(annotation_candidate)

    for source_name, (dense_model, dense_features) in dense_sources.items():
        if target_key in dense_model.targets:
            dense_candidate = _source_candidate(
                source_name,
                source_name,
                _try_dense_ranking(dense_model, dense_features, labels, split=validation_split, target_key=target_key, ks=cutoffs),
                _try_dense_ranking(dense_model, dense_features, labels, split=test_split, target_key=target_key, ks=cutoffs),
                complexity=1,
                weight=1.0,
            )
            if dense_candidate is not None:
                candidates.append(dense_candidate)
        if target_key in annotation_model.targets and target_key in dense_model.targets:
            validation_fusion = _try_fusion_ranking(
                annotation_model,
                annotation_features,
                dense_model,
                dense_features,
                labels,
                split=validation_split,
                target_key=target_key,
                ks=cutoffs,
                weights=weights,
                select_k=select_k,
            )
            if validation_fusion is None:
                continue
            best_weight = validation_fusion["best_weight"]
            test_fusion = _try_fusion_ranking(
                annotation_model,
                annotation_features,
                dense_model,
                dense_features,
                labels,
                split=test_split,
                target_key=target_key,
                ks=cutoffs,
                weights=[best_weight],
                select_k=select_k,
            )
            fusion_candidate = _source_candidate(
                f"fusion:{source_name}",
                source_name,
                validation_fusion["best"],
                test_fusion["best"] if test_fusion is not None else None,
                complexity=2,
                weight=best_weight,
            )
            if fusion_candidate is not None:
                candidates.append(fusion_candidate)

    if not candidates:
        return None
    best = max(
        candidates,
        key=lambda row: (
            float(row["validation_metrics"].get(selection_metric, 0.0)),
            float(row["validation_metrics"].get(f"recall_at_{select_k}", 0.0)),
            -int(row["complexity"]),
            row["best_source"],
        ),
    )
    annotation_test_metrics = next(
        (row["test_metrics"] for row in candidates if row["best_source"] == "annotation"),
        {},
    )
    panel, label = _split_target_key(target_key)
    return {
        "target_key": target_key,
        "panel": panel,
        "label": label,
        "best_source": best["best_source"],
        "dense_source": best["dense_source"],
        "best_weight": best["best_weight"],
        "selected_on_split": validation_split,
        "selection_metric": selection_metric,
        "validation_metrics": best["validation_metrics"],
        "test_metrics": best["test_metrics"],
        "validation_n_labeled_candidates": best["validation_n_labeled_candidates"],
        "test_n_labeled_candidates": best["test_n_labeled_candidates"],
        "test_n_positives": best["test_n_positives"],
        "annotation_test_metrics": annotation_test_metrics,
        "beats_annotation": _beats_annotation(best["test_metrics"], annotation_test_metrics, selection_metric),
        "source_metrics": [
            {
                "source": row["best_source"],
                "dense_source": row["dense_source"],
                "weight": row["best_weight"],
                "validation_metrics": row["validation_metrics"],
                "test_metrics": row["test_metrics"],
            }
            for row in sorted(candidates, key=lambda item: item["best_source"])
        ],
    }


def _source_candidate(
    source: str,
    dense_source: str | None,
    validation_report: dict[str, Any] | None,
    test_report: dict[str, Any] | None,
    *,
    complexity: int,
    weight: float | None = None,
) -> dict[str, Any] | None:
    if validation_report is None or test_report is None:
        return None
    return {
        "best_source": source,
        "dense_source": dense_source,
        "best_weight": weight,
        "complexity": complexity,
        "validation_metrics": validation_report["metrics"],
        "test_metrics": test_report["metrics"],
        "validation_n_labeled_candidates": validation_report["n_labeled_candidates"],
        "test_n_labeled_candidates": test_report["n_labeled_candidates"],
        "test_n_positives": test_report["n_positives"],
    }


def _try_annotation_ranking(
    model: BaselineModel,
    features: FeatureMatrix,
    labels: list[LabelRecord],
    *,
    split: str,
    target_key: str,
    ks: list[int],
) -> dict[str, Any] | None:
    try:
        report = evaluate_ranking(model, features, labels, split=split, target_key=target_key, ks=ks)
    except (KeyError, ValueError):
        return None
    return report if report["n_labeled_candidates"] else None


def _try_dense_ranking(
    model: DenseBaselineModel,
    features: DenseFeatureMatrix,
    labels: list[LabelRecord],
    *,
    split: str,
    target_key: str,
    ks: list[int],
) -> dict[str, Any] | None:
    try:
        report = evaluate_dense_ranking(model, features, labels, split=split, target_key=target_key, ks=ks)
    except (KeyError, ValueError):
        return None
    return report if report["n_labeled_candidates"] else None


def _try_fusion_ranking(
    annotation_model: BaselineModel,
    annotation_features: FeatureMatrix,
    dense_model: DenseBaselineModel,
    dense_features: DenseFeatureMatrix,
    labels: list[LabelRecord],
    *,
    split: str,
    target_key: str,
    ks: list[int],
    weights: list[float] | None,
    select_k: int,
) -> dict[str, Any] | None:
    try:
        return evaluate_fusion_ranking(
            annotation_model,
            annotation_features,
            dense_model,
            dense_features,
            labels,
            split=split,
            target_key=target_key,
            ks=ks,
            weights=weights,
            select_k=select_k,
        )
    except (KeyError, ValueError):
        return None


def _rank_from_leaderboard_row(
    row: dict[str, Any],
    annotation_model: BaselineModel,
    annotation_features: FeatureMatrix,
    dense_sources: DenseSourceMap,
    labels: list[LabelRecord],
    *,
    split: str,
    limit: int,
) -> dict[str, Any]:
    target_key = row["target_key"]
    source = row["best_source"]
    if source == "annotation":
        ranking = evaluate_ranking(annotation_model, annotation_features, labels, split=split, target_key=target_key, ks=[limit])
    elif source.startswith("fusion:"):
        source_name = row["dense_source"]
        dense_model, dense_features = dense_sources[source_name]
        ranking = evaluate_fusion_ranking(
            annotation_model,
            annotation_features,
            dense_model,
            dense_features,
            labels,
            split=split,
            target_key=target_key,
            ks=[limit],
            weights=[row["best_weight"]],
        )["best"]
    else:
        dense_model, dense_features = dense_sources[source]
        ranking = evaluate_dense_ranking(dense_model, dense_features, labels, split=split, target_key=target_key, ks=[limit])
    return {
        "target_key": target_key,
        "panel": row["panel"],
        "label": row["label"],
        "source": source,
        "dense_source": row.get("dense_source"),
        "weight": row.get("best_weight"),
        "metrics": row["test_metrics"],
        "candidates": ranking["ranked_candidates"][:limit],
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


def _target_keys(labels: list[LabelRecord]) -> list[str]:
    return sorted({record.target_key for record in labels})


def _split_target_key(target_key: str) -> tuple[str, str]:
    if ":" not in target_key:
        return "unknown", target_key
    panel, label = target_key.split(":", 1)
    return panel, label


def _beats_annotation(
    metrics: dict[str, Any],
    annotation_metrics: dict[str, Any],
    selection_metric: str,
) -> bool:
    if not annotation_metrics:
        return False
    return float(metrics.get(selection_metric, 0.0)) > float(annotation_metrics.get(selection_metric, 0.0))


def _panel_summary(rows: list[dict[str, Any]], selection_metric: str) -> dict[str, dict[str, Any]]:
    panels: dict[str, dict[str, Any]] = {}
    for row in rows:
        panel = row["panel"]
        summary = panels.setdefault(
            panel,
            {
                "n_targets": 0,
                "n_beating_annotation": 0,
                "best_targets": [],
            },
        )
        summary["n_targets"] += 1
        summary["n_beating_annotation"] += int(row["beats_annotation"])
        summary["best_targets"].append(
            {
                "target_key": row["target_key"],
                "best_source": row["best_source"],
                "precision": row["test_metrics"].get(selection_metric, 0.0),
            }
        )
    for summary in panels.values():
        summary["best_targets"] = sorted(
            summary["best_targets"],
            key=lambda item: (-float(item["precision"]), item["target_key"]),
        )[:10]
    return panels


def _candidate_panel_summary(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    panels: dict[str, dict[str, Any]] = {}
    for row in rows:
        summary = panels.setdefault(row["panel"], {"n_targets": 0, "targets": []})
        summary["n_targets"] += 1
        summary["targets"].append(row["target_key"])
    return panels


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
