"""Simple learned baselines for annotation-derived feature matrices."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from typing import Any

from microbial_function_discovery.datasets import LabelRecord
from microbial_function_discovery.features import FeatureMatrix
from microbial_function_discovery.panels import APPLICATION_AREAS


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
        return _sigmoid(logit)

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

    feature_genomes = set(features.genome_ids)
    train_labels = [record for record in labels if record.split == train_split and record.genome_id in feature_genomes]
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

    feature_genomes = set(features.genome_ids)
    records = [
        record
        for record in labels
        if record.split == split and record.target_key in model.targets and record.genome_id in feature_genomes
    ]
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


def evaluate_ranking(
    model: BaselineModel,
    features: FeatureMatrix,
    labels: list[LabelRecord],
    *,
    split: str,
    target_key: str | None = None,
    panel: str | None = None,
    ks: list[int] | None = None,
) -> dict[str, Any]:
    """Evaluate whether top-ranked candidates recover held-out positives."""

    if (target_key is None) == (panel is None):
        raise ValueError("provide exactly one of target_key or panel")
    if target_key is not None and target_key not in model.targets:
        raise KeyError(f"target not found in model: {target_key}")
    if panel is not None and panel not in APPLICATION_AREAS:
        raise ValueError(f"unknown panel: {panel}")

    cutoffs = _normalize_cutoffs(ks or [1, 5, 10])
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
        ranked.append(
            {
                "genome_id": genome_id,
                "truth": truth,
                "score": round(score, 6),
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
        "mode": "target" if target_key is not None else "panel",
        "query": target_key if target_key is not None else panel,
        "split": split,
        "n_labeled_candidates": len(ranked),
        "n_positives": n_positives,
        "metrics": metrics,
        "ranked_candidates": ranked,
    }


def predict_model(model: BaselineModel, features: FeatureMatrix, *, genome_id: str) -> dict[str, Any]:
    """Emit product-style prediction JSON from a trained baseline model."""

    scored_targets = []
    for target_key in sorted(model.targets):
        panel, label = _split_target_key(target_key)
        score = model.predict_proba(target_key, genome_id, features)
        scored_targets.append(
            {
                "target_key": target_key,
                "panel": panel,
                "label": label,
                "score": score,
                "evidence": _active_evidence(model, features, genome_id, target_key),
            }
        )

    functions = [
        {
            "name": item["label"].replace("_", " "),
            "score": round(item["score"], 6),
            "confidence": round(_confidence_from_score(item["score"]), 6),
            "evidence": item["evidence"],
        }
        for item in sorted(scored_targets, key=lambda row: row["score"], reverse=True)
    ]
    if not functions:
        functions = [
            {
                "name": "unresolved function",
                "score": 0.05,
                "confidence": 0.2,
                "evidence": [
                    {
                        "type": "database_hit",
                        "id": "no_trained_targets",
                        "annotation": "No trained target models are available.",
                        "database": "baseline",
                        "weight": 0.05,
                    }
                ],
            }
        ]

    return {
        "genome_id": genome_id,
        "application_scores": _application_scores_from_targets(scored_targets),
        "functions": functions,
        "biosafety": _biosafety_from_targets(scored_targets),
        "novelty": {
            "nearest_training_family": "unknown",
            "generalization_risk": "unknown",
            "notes": [
                "Prediction uses a no-GPU annotation-feature baseline; calibrated taxonomic novelty is not available yet."
            ],
        },
    }


def rank_candidates(
    model: BaselineModel,
    features: FeatureMatrix,
    *,
    target_key: str | None = None,
    panel: str | None = None,
    limit: int | None = None,
) -> dict[str, Any]:
    """Rank candidate genomes by one target or by best score in a panel."""

    if (target_key is None) == (panel is None):
        raise ValueError("provide exactly one of target_key or panel")
    if target_key is not None and target_key not in model.targets:
        raise KeyError(f"target not found in model: {target_key}")
    if panel is not None and panel not in APPLICATION_AREAS:
        raise ValueError(f"unknown panel: {panel}")

    candidates = []
    for genome_id in features.genome_ids:
        if target_key is not None:
            score = model.predict_proba(target_key, genome_id, features)
            panel_key, label = _split_target_key(target_key)
            evidence = _active_evidence(model, features, genome_id, target_key)
            target_scores = [
                {
                    "target_key": target_key,
                    "label": label.replace("_", " "),
                    "score": round(score, 6),
                }
            ]
        else:
            panel_targets = [key for key in model.targets if key.startswith(f"{panel}:")]
            if not panel_targets:
                score = 0.05
                panel_key = panel or "unknown"
                evidence = []
                target_scores = []
            else:
                target_rows = [
                    {
                        "target_key": key,
                        "label": _split_target_key(key)[1].replace("_", " "),
                        "score": model.predict_proba(key, genome_id, features),
                    }
                    for key in panel_targets
                ]
                best = max(target_rows, key=lambda row: row["score"])
                score = best["score"]
                panel_key = panel or "unknown"
                evidence = _active_evidence(model, features, genome_id, best["target_key"])
                target_scores = [
                    {**row, "score": round(row["score"], 6)}
                    for row in sorted(target_rows, key=lambda row: row["score"], reverse=True)
                ]

        candidates.append(
            {
                "genome_id": genome_id,
                "score": round(score, 6),
                "panel": panel_key,
                "target_scores": target_scores,
                "evidence": evidence,
            }
        )

    ranked = sorted(candidates, key=lambda row: (-row["score"], row["genome_id"]))
    if limit is not None:
        ranked = ranked[:limit]
    return {
        "mode": "target" if target_key is not None else "panel",
        "query": target_key if target_key is not None else panel,
        "n_candidates": len(candidates),
        "candidates": ranked,
    }


def _train_target_model(features: FeatureMatrix, labels: list[LabelRecord]) -> TargetModel:
    positives = [record for record in labels if record.value == 1]
    negatives = [record for record in labels if record.value == 0]
    alpha = 1.0
    prior = (len(positives) + alpha) / (len(labels) + 2 * alpha)
    pos_counts = _feature_counts(features, positives)
    neg_counts = _feature_counts(features, negatives)
    weights = []
    for feature_index in range(len(features.feature_names)):
        pos_with = pos_counts[feature_index]
        neg_with = neg_counts[feature_index]
        p_feature_pos = (pos_with + alpha) / (len(positives) + 2 * alpha)
        p_feature_neg = (neg_with + alpha) / (len(negatives) + 2 * alpha)
        weights.append(_logit(p_feature_pos) - _logit(p_feature_neg))
    return TargetModel(prior_log_odds=_logit(prior), feature_log_odds=weights)


def _feature_counts(features: FeatureMatrix, labels: list[LabelRecord]) -> list[int]:
    counts = [0] * len(features.feature_names)
    for record in labels:
        for feature_index, value in enumerate(features.row_for(record.genome_id)):
            counts[feature_index] += value
    return counts


def _normalize_cutoffs(ks: list[int]) -> list[int]:
    cutoffs = sorted(set(ks))
    if not cutoffs or any(k <= 0 for k in cutoffs):
        raise ValueError("ranking cutoffs must be positive integers")
    return cutoffs


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
        if record.split != split:
            continue
        if record.genome_id not in feature_genomes:
            continue
        if target_key is not None and record.target_key != target_key:
            continue
        if panel is not None and record.panel != panel:
            continue
        truth_by_genome[record.genome_id] = max(truth_by_genome.get(record.genome_id, 0), record.value)
    return truth_by_genome


def _ranking_score(
    model: BaselineModel,
    features: FeatureMatrix,
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


def _active_evidence(
    model: BaselineModel,
    features: FeatureMatrix,
    genome_id: str,
    target_key: str,
) -> list[dict[str, Any]]:
    row = _align_row(features, genome_id, model.feature_names)
    weights = model.targets[target_key].feature_log_odds
    evidence = []
    for feature_name, value, weight in sorted(zip(model.feature_names, row, weights), key=lambda item: item[2], reverse=True):
        if not value or weight <= 0:
            continue
        database, accession = _split_feature_name(feature_name)
        evidence.append(
            {
                "type": "database_hit",
                "id": feature_name,
                "annotation": f"active feature {feature_name} supports {target_key}",
                "database": database,
                "weight": round(min(1.0, weight / 2), 6),
            }
        )
    if evidence:
        return evidence[:5]
    return [
        {
            "type": "database_hit",
            "id": "no_positive_active_features",
            "annotation": f"No positive active annotation features support {target_key}.",
            "database": "baseline",
            "weight": 0.05,
        }
    ]


def _application_scores_from_targets(scored_targets: list[dict[str, Any]]) -> list[dict[str, Any]]:
    scores = []
    for panel in APPLICATION_AREAS:
        panel_items = [item for item in scored_targets if item["panel"] == panel]
        if panel_items:
            best = max(panel_items, key=lambda item: item["score"])
            score = best["score"]
            confidence = _confidence_from_score(score)
            top_functions = [
                item["label"].replace("_", " ")
                for item in sorted(panel_items, key=lambda item: item["score"], reverse=True)[:3]
            ]
        else:
            score = 0.05
            confidence = 0.2
            top_functions = []
        scores.append(
            {
                "area": panel,
                "score": round(score, 6),
                "confidence": round(confidence, 6),
                "top_functions": top_functions,
            }
        )
    return scores


def _biosafety_from_targets(scored_targets: list[dict[str, Any]]) -> dict[str, Any]:
    biosafety_items = [item for item in scored_targets if item["panel"] == "biosafety"]
    pathogenicity_scores = [
        item["score"]
        for item in biosafety_items
        if "pathogen" in item["label"] or "virulence" in item["label"] or "toxin" in item["label"]
    ]
    amr_scores = [
        item["score"]
        for item in biosafety_items
        if "amr" in item["label"] or "resistance" in item["label"] or "antimicrobial" in item["label"]
    ]
    pathogenicity_risk = max(pathogenicity_scores, default=0.05)
    amr_risk = max(amr_scores, default=0.05)
    warnings = []
    if pathogenicity_risk >= 0.5:
        warnings.append("Learned baseline predicts elevated pathogenicity or virulence risk.")
    if amr_risk >= 0.5:
        warnings.append("Learned baseline predicts elevated antimicrobial resistance risk.")
    return {
        "pathogenicity_risk": round(pathogenicity_risk, 6),
        "amr_risk": round(amr_risk, 6),
        "warnings": warnings,
    }


def _confidence_from_score(score: float) -> float:
    return min(0.95, 0.2 + abs(score - 0.5) * 1.5)


def _split_target_key(target_key: str) -> tuple[str, str]:
    panel, label = target_key.split(":", maxsplit=1)
    return panel, label


def _split_feature_name(feature_name: str) -> tuple[str, str]:
    if ":" not in feature_name:
        return "annotation", feature_name
    return feature_name.split(":", maxsplit=1)


def _align_row(features: FeatureMatrix, genome_id: str, model_feature_names: list[str]) -> list[int]:
    source = dict(zip(features.feature_names, features.row_for(genome_id)))
    return [source.get(feature_name, 0) for feature_name in model_feature_names]


def _logit(p: float) -> float:
    return math.log(p / (1 - p))


def _sigmoid(logit: float) -> float:
    if logit >= 0:
        z = math.exp(-logit)
        return 1 / (1 + z)
    z = math.exp(logit)
    return z / (1 + z)
