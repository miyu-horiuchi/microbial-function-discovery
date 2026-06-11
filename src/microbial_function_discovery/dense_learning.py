"""Dense embedding baselines for cached microbial genome feature arrays."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from microbial_function_discovery.datasets import LabelRecord
from microbial_function_discovery.panels import APPLICATION_AREAS


@dataclass(frozen=True)
class DenseFeatureMatrix:
    genome_ids: list[str]
    feature_indexes: list[int]
    rows: list[list[float]]
    _row_index: dict[str, int] = field(default_factory=dict, init=False, compare=False, repr=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "_row_index", {genome_id: index for index, genome_id in enumerate(self.genome_ids)})

    def row_for(self, genome_id: str) -> list[float]:
        try:
            return self.rows[self._row_index[genome_id]]
        except KeyError as exc:
            raise KeyError(f"genome id not found in dense feature matrix: {genome_id}") from exc


@dataclass(frozen=True)
class DenseTargetModel:
    bias: float
    weights: list[float]
    threshold: float = 0.5


@dataclass(frozen=True)
class DenseBaselineModel:
    feature_indexes: list[int]
    feature_mean: list[float]
    feature_scale: list[float]
    targets: dict[str, DenseTargetModel]

    def predict_proba(self, target_key: str, genome_id: str, features: DenseFeatureMatrix) -> float:
        target = self.targets[target_key]
        row = _standardize(features.row_for(genome_id), self.feature_mean, self.feature_scale)
        logit = target.bias + sum(value * weight for value, weight in zip(row, target.weights))
        return _sigmoid(logit)

    def to_json(self) -> str:
        return json.dumps(
            {
                "feature_indexes": self.feature_indexes,
                "feature_mean": self.feature_mean,
                "feature_scale": self.feature_scale,
                "targets": {
                    key: {
                        "bias": target.bias,
                        "weights": target.weights,
                        "threshold": target.threshold,
                    }
                    for key, target in self.targets.items()
                },
            },
            indent=2,
            sort_keys=True,
        )

    @classmethod
    def from_json(cls, text: str) -> "DenseBaselineModel":
        data = json.loads(text)
        return cls(
            feature_indexes=[int(value) for value in data["feature_indexes"]],
            feature_mean=[float(value) for value in data["feature_mean"]],
            feature_scale=[float(value) for value in data["feature_scale"]],
            targets={
                key: DenseTargetModel(
                    bias=float(value["bias"]),
                    weights=[float(weight) for weight in value["weights"]],
                    threshold=float(value["threshold"]),
                )
                for key, value in data["targets"].items()
            },
        )


def load_dense_npz_feature_matrix(
    npz_path: Path,
    *,
    keep_genome_ids: set[str] | None = None,
    max_features: int | None = None,
    feature_indexes: list[int] | None = None,
) -> DenseFeatureMatrix:
    """Load dense legacy NPZ features without binarizing them."""

    try:
        import numpy as np  # type: ignore[import-not-found]
    except ImportError as exc:
        raise RuntimeError("dense NPZ loading requires optional dependency: numpy") from exc

    data = np.load(npz_path, allow_pickle=True)
    genome_ids = [_id(value) for value in data["bacdive_ids"].tolist()]
    matrix = data["features"]
    kept_row_indexes = [
        index
        for index, genome_id in enumerate(genome_ids)
        if genome_id and (keep_genome_ids is None or genome_id in keep_genome_ids)
    ]
    if not kept_row_indexes:
        raise ValueError("dense feature matrix contains no selected genomes")
    selected_rows = matrix[kept_row_indexes]
    if feature_indexes is None:
        n_features = selected_rows.shape[1]
        feature_indexes = list(range(n_features))
        if max_features is not None and max_features < n_features:
            variances = selected_rows.var(axis=0)
            ranked = sorted(enumerate(variances.tolist()), key=lambda item: (-item[1], item[0]))
            feature_indexes = sorted(index for index, _variance in ranked[:max_features])
    else:
        feature_indexes = [int(index) for index in feature_indexes]
    return DenseFeatureMatrix(
        genome_ids=[genome_ids[index] for index in kept_row_indexes],
        feature_indexes=feature_indexes,
        rows=selected_rows[:, feature_indexes].astype(float).tolist(),
    )


def train_dense_baseline(
    features: DenseFeatureMatrix,
    labels: list[LabelRecord],
    *,
    train_split: str = "train",
) -> DenseBaselineModel:
    """Train one linear dense-feature classifier per target."""

    feature_genomes = set(features.genome_ids)
    train_labels = [record for record in labels if record.split == train_split and record.genome_id in feature_genomes]
    if not train_labels:
        raise ValueError(f"no labels found for split: {train_split}")

    mean, scale = _feature_mean_scale(features.rows)
    standardized_rows = {
        genome_id: _standardize(row, mean, scale)
        for genome_id, row in zip(features.genome_ids, features.rows)
    }
    by_target: dict[str, list[LabelRecord]] = {}
    for record in train_labels:
        by_target.setdefault(record.target_key, []).append(record)

    targets = {
        target_key: _train_dense_target(standardized_rows, target_records, len(mean))
        for target_key, target_records in sorted(by_target.items())
    }
    return DenseBaselineModel(
        feature_indexes=features.feature_indexes,
        feature_mean=mean,
        feature_scale=scale,
        targets=targets,
    )


def evaluate_dense_model(
    model: DenseBaselineModel,
    features: DenseFeatureMatrix,
    labels: list[LabelRecord],
    *,
    split: str,
) -> dict[str, Any]:
    feature_genomes = set(features.genome_ids)
    records = [
        record
        for record in labels
        if record.split == split and record.target_key in model.targets and record.genome_id in feature_genomes
    ]
    by_target: dict[str, list[dict[str, Any]]] = {}
    for record in records:
        proba = model.predict_proba(record.target_key, record.genome_id, features)
        pred = int(proba >= model.targets[record.target_key].threshold)
        by_target.setdefault(record.target_key, []).append(
            {
                "genome_id": record.genome_id,
                "truth": record.value,
                "prediction": pred,
                "score": round(proba, 6),
            }
        )

    correct_total = 0
    n_total = 0
    target_reports = {}
    for target_key, rows in sorted(by_target.items()):
        correct = sum(int(row["truth"] == row["prediction"]) for row in rows)
        n = len(rows)
        correct_total += correct
        n_total += n
        target_reports[target_key] = {
            "n": n,
            "accuracy": correct / n if n else 0.0,
            "predictions": rows,
        }

    return {
        "split": split,
        "overall": {
            "n": n_total,
            "accuracy": correct_total / n_total if n_total else 0.0,
        },
        "targets": target_reports,
    }


def evaluate_dense_ranking(
    model: DenseBaselineModel,
    features: DenseFeatureMatrix,
    labels: list[LabelRecord],
    *,
    split: str,
    target_key: str | None = None,
    panel: str | None = None,
    ks: list[int] | None = None,
) -> dict[str, Any]:
    if (target_key is None) == (panel is None):
        raise ValueError("provide exactly one of target_key or panel")
    if target_key is not None and target_key not in model.targets:
        raise KeyError(f"target not found in model: {target_key}")
    if panel is not None and panel not in APPLICATION_AREAS:
        raise ValueError(f"unknown panel: {panel}")

    truth_by_genome = _ranking_truth(
        labels,
        split=split,
        target_key=target_key,
        panel=panel,
        feature_genomes=set(features.genome_ids),
    )
    ranked = []
    for genome_id, truth in truth_by_genome.items():
        score = _ranking_score(model, features, genome_id, target_key=target_key, panel=panel)
        ranked.append({"genome_id": genome_id, "truth": truth, "score": round(score, 6)})
    ranked.sort(key=lambda row: (-row["score"], row["genome_id"]))
    for rank, row in enumerate(ranked, start=1):
        row["rank"] = rank

    n_positives = sum(row["truth"] for row in ranked)
    metrics: dict[str, float | int] = {}
    for k in _normalize_cutoffs(ks or [1, 5, 10]):
        top_k = ranked[:k]
        hits = sum(row["truth"] for row in top_k)
        denominator = min(k, len(ranked))
        metrics[f"hits_at_{k}"] = hits
        metrics[f"precision_at_{k}"] = hits / denominator if denominator else 0.0
        metrics[f"recall_at_{k}"] = hits / n_positives if n_positives else 0.0

    return {
        "mode": "target" if target_key is not None else "panel",
        "query": target_key if target_key is not None else panel,
        "split": split,
        "n_labeled_candidates": len(ranked),
        "n_positives": n_positives,
        "metrics": metrics,
        "ranked_candidates": ranked,
    }


def _train_dense_target(
    standardized_rows: dict[str, list[float]],
    labels: list[LabelRecord],
    n_features: int,
) -> DenseTargetModel:
    positives = [record for record in labels if record.value == 1]
    negatives = [record for record in labels if record.value == 0]
    alpha = 1.0
    prior = (len(positives) + alpha) / (len(labels) + 2 * alpha)
    pos_mean = _mean_row([standardized_rows[record.genome_id] for record in positives], n_features)
    neg_mean = _mean_row([standardized_rows[record.genome_id] for record in negatives], n_features)
    weights = [pos - neg for pos, neg in zip(pos_mean, neg_mean)]
    bias = _logit(prior) - 0.5 * (
        sum(value * value for value in pos_mean) - sum(value * value for value in neg_mean)
    )
    return DenseTargetModel(bias=bias, weights=weights)


def _feature_mean_scale(rows: list[list[float]]) -> tuple[list[float], list[float]]:
    n_rows = len(rows)
    n_features = len(rows[0]) if rows else 0
    mean = [0.0] * n_features
    for row in rows:
        for index, value in enumerate(row):
            mean[index] += value
    mean = [value / n_rows for value in mean]
    variance = [0.0] * n_features
    for row in rows:
        for index, value in enumerate(row):
            delta = value - mean[index]
            variance[index] += delta * delta
    scale = [math.sqrt(value / n_rows) or 1.0 for value in variance]
    return mean, scale


def _standardize(row: list[float], mean: list[float], scale: list[float]) -> list[float]:
    return [(value - center) / spread for value, center, spread in zip(row, mean, scale)]


def _mean_row(rows: list[list[float]], n_features: int) -> list[float]:
    if not rows:
        return [0.0] * n_features
    mean = [0.0] * n_features
    for row in rows:
        for index, value in enumerate(row):
            mean[index] += value
    return [value / len(rows) for value in mean]


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


def _ranking_score(
    model: DenseBaselineModel,
    features: DenseFeatureMatrix,
    genome_id: str,
    *,
    target_key: str | None,
    panel: str | None,
) -> float:
    if target_key is not None:
        return model.predict_proba(target_key, genome_id, features)
    panel_targets = [key for key in model.targets if key.startswith(f"{panel}:")]
    if not panel_targets:
        return 0.05
    return max(model.predict_proba(key, genome_id, features) for key in panel_targets)


def _normalize_cutoffs(ks: list[int]) -> list[int]:
    cutoffs = sorted(set(ks))
    if not cutoffs or any(k <= 0 for k in cutoffs):
        raise ValueError("ranking cutoffs must be positive integers")
    return cutoffs


def _id(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value).strip()


def _logit(p: float) -> float:
    return math.log(p / (1 - p))


def _sigmoid(logit: float) -> float:
    if logit >= 0:
        z = math.exp(-logit)
        return 1 / (1 + z)
    z = math.exp(logit)
    return z / (1 + z)
