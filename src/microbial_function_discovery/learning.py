"""Simple learned baselines for annotation-derived feature matrices."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from typing import Any

from microbial_function_discovery.datasets import LabelRecord
from microbial_function_discovery.features import FeatureMatrix


@dataclass(frozen=True)
class TargetModel:
    prior_log_odds: float
    feature_log_odds: list[float]
    threshold: float = 0.5


@dataclass(frozen=True)
class BaselineModel:
    feature_names: list[str]
    targets: dict[str, TargetModel]

    def predict_proba(self, target_key: str, genome_id: str, features: FeatureMatrix) -> float:
        target = self.targets[target_key]
        row = _align_row(features, genome_id, self.feature_names)
        logit = target.prior_log_odds + sum(value * weight for value, weight in zip(row, target.feature_log_odds))
        return 1 / (1 + math.exp(-logit))

    def predict_label(self, target_key: str, genome_id: str, features: FeatureMatrix) -> int:
        return int(self.predict_proba(target_key, genome_id, features) >= self.targets[target_key].threshold)

    def to_json(self) -> str:
        return json.dumps(
            {
                "feature_names": self.feature_names,
                "targets": {
                    key: {
                        "prior_log_odds": target.prior_log_odds,
                        "feature_log_odds": target.feature_log_odds,
                        "threshold": target.threshold,
                    }
                    for key, target in self.targets.items()
                },
            },
            indent=2,
            sort_keys=True,
        )

    @classmethod
    def from_json(cls, text: str) -> "BaselineModel":
        data = json.loads(text)
        return cls(
            feature_names=list(data["feature_names"]),
            targets={
                key: TargetModel(
                    prior_log_odds=float(value["prior_log_odds"]),
                    feature_log_odds=[float(weight) for weight in value["feature_log_odds"]],
                    threshold=float(value["threshold"]),
                )
                for key, value in data["targets"].items()
            },
        )

    def to_dict(self) -> dict[str, Any]:
        return json.loads(self.to_json())


def train_baseline(features: FeatureMatrix, labels: list[LabelRecord], *, train_split: str = "train") -> BaselineModel:
    """Train one simple Bernoulli-style linear classifier per target."""

    train_labels = [record for record in labels if record.split == train_split]
    if not train_labels:
        raise ValueError(f"no labels found for split: {train_split}")

    by_target: dict[str, list[LabelRecord]] = {}
    for record in train_labels:
        by_target.setdefault(record.target_key, []).append(record)

    targets = {
        target_key: _train_target_model(features, target_records)
        for target_key, target_records in sorted(by_target.items())
    }
    return BaselineModel(feature_names=features.feature_names, targets=targets)


def evaluate_model(
    model: BaselineModel,
    features: FeatureMatrix,
    labels: list[LabelRecord],
    *,
    split: str,
) -> dict[str, Any]:
    """Evaluate a baseline model on one split."""

    records = [record for record in labels if record.split == split and record.target_key in model.targets]
    by_target: dict[str, list[dict[str, Any]]] = {}
    for record in records:
        proba = model.predict_proba(record.target_key, record.genome_id, features)
        pred = int(proba >= 0.5)
        by_target.setdefault(record.target_key, []).append(
            {
                "genome_id": record.genome_id,
                "truth": record.value,
                "prediction": pred,
                "score": round(proba, 6),
            }
        )

    target_reports = {}
    correct_total = 0
    n_total = 0
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


def _train_target_model(features: FeatureMatrix, labels: list[LabelRecord]) -> TargetModel:
    positives = [record for record in labels if record.value == 1]
    negatives = [record for record in labels if record.value == 0]
    alpha = 1.0
    prior = (len(positives) + alpha) / (len(labels) + 2 * alpha)
    weights = []
    for feature_index in range(len(features.feature_names)):
        pos_with = sum(features.row_for(record.genome_id)[feature_index] for record in positives)
        neg_with = sum(features.row_for(record.genome_id)[feature_index] for record in negatives)
        p_feature_pos = (pos_with + alpha) / (len(positives) + 2 * alpha)
        p_feature_neg = (neg_with + alpha) / (len(negatives) + 2 * alpha)
        weights.append(_logit(p_feature_pos) - _logit(p_feature_neg))
    return TargetModel(prior_log_odds=_logit(prior), feature_log_odds=weights)


def _align_row(features: FeatureMatrix, genome_id: str, model_feature_names: list[str]) -> list[int]:
    source = dict(zip(features.feature_names, features.row_for(genome_id)))
    return [source.get(feature_name, 0) for feature_name in model_feature_names]


def _logit(p: float) -> float:
    return math.log(p / (1 - p))
