"""Import helpers for reusing artifacts from the earlier microbe-foundation repo."""

from __future__ import annotations

import csv
import io
import math
import re
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

from microbial_function_discovery.datasets import LabelRecord
from microbial_function_discovery.features import FeatureMatrix


MULTICLASS_TARGETS = {
    "gram_stain": ("biosafety", ["positive", "negative", "variable"]),
    "cell_shape": ("food_fermentation_agriculture", ["coccus", "rod", "vibrio", "spiral", "filament", "pleomorphic", "other"]),
    "oxygen_tolerance": (
        "environmental_terraforming",
        [
            "obligate_aerobe",
            "facultative_aerobe",
            "microaerophile",
            "aerotolerant",
            "facultative_anaerobe",
            "obligate_anaerobe",
        ],
    ),
    "halophily": ("environmental_terraforming", ["non_halophile", "halotolerant", "halophile", "extreme_halophile"]),
    "temperature_class": (
        "biofuels_industrial",
        ["psychrophile", "psychrotroph", "mesophile", "thermophile", "hyperthermophile"],
    ),
    "ph_class": ("biofuels_industrial", ["acidophile", "neutrophile", "alkaliphile"]),
    "biosafety_level": ("biosafety", ["BSL-1", "BSL-2", "BSL-3", "BSL-4"]),
}

BINARY_TARGETS = {
    "motility": ("environmental_terraforming", "motility"),
    "sporulation": ("environmental_terraforming", "sporulation"),
    "pigmentation": ("food_fermentation_agriculture", "pigmentation"),
    "catalase": ("environmental_terraforming", "catalase_activity"),
    "cytochrome_oxidase": ("environmental_terraforming", "cytochrome_oxidase_activity"),
    "pathogenicity_human": ("biosafety", "pathogenicity_human"),
    "pathogenicity_animal": ("biosafety", "pathogenicity_animal"),
}

MULTILABEL_TARGETS = {
    "carbon_utilization": "biofuels_industrial",
    "metabolite_production": "food_fermentation_agriculture",
    "amr_phenotype": "biosafety",
}


def load_table_records(path: Path) -> list[dict[str, Any]]:
    """Load CSV/TSV directly, or parquet through optional pandas/pyarrow."""

    suffix = path.suffix.lower()
    if suffix in {".tsv", ".txt"}:
        return _load_delimited(path, delimiter="\t")
    if suffix == ".csv":
        return _load_delimited(path, delimiter=",")
    if suffix == ".parquet":
        try:
            import pandas as pd  # type: ignore[import-not-found]
        except ImportError as exc:
            raise RuntimeError("parquet import requires optional dependencies: pandas and pyarrow") from exc
        frame = pd.read_parquet(path)
        return frame.to_dict(orient="records")
    raise ValueError(f"unsupported table format: {path.suffix}")


def labels_tsv_from_legacy_tables(
    trait_rows: Iterable[dict[str, Any]],
    split_rows: Iterable[dict[str, Any]],
    *,
    split_column: str = "family_split",
    max_multilabel_classes: int = 50,
) -> str:
    """Convert legacy BacDive trait/split rows to benchmark label TSV text."""

    records = records_from_legacy_tables(
        trait_rows,
        split_rows,
        split_column=split_column,
        max_multilabel_classes=max_multilabel_classes,
    )
    out = io.StringIO()
    writer = csv.DictWriter(
        out,
        fieldnames=["genome_id", "split", "family", "panel", "label", "value"],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    for record in records:
        writer.writerow(
            {
                "genome_id": record.genome_id,
                "split": record.split,
                "family": record.family,
                "panel": record.panel,
                "label": record.label,
                "value": record.value,
            }
        )
    return out.getvalue()


def records_from_legacy_tables(
    trait_rows: Iterable[dict[str, Any]],
    split_rows: Iterable[dict[str, Any]],
    *,
    split_column: str = "family_split",
    max_multilabel_classes: int = 50,
) -> list[LabelRecord]:
    """Convert parsed BacDive trait rows into binary target records."""

    traits = list(trait_rows)
    split_by_id = {_id(row.get("bacdive_id")): row for row in split_rows if _id(row.get("bacdive_id"))}
    multilabel_classes = _select_multilabel_classes(traits, max_multilabel_classes=max_multilabel_classes)
    records: list[LabelRecord] = []

    for row in traits:
        genome_id = _id(row.get("bacdive_id"))
        if not genome_id or genome_id not in split_by_id:
            continue
        split_row = split_by_id[genome_id]
        split = _string(split_row.get(split_column))
        if not split:
            continue
        family = _string(split_row.get("family")) or _string(row.get("family")) or "unknown_family"

        for trait, (panel, label) in BINARY_TARGETS.items():
            value = _coerce_binary(row.get(trait))
            if value is not None:
                records.append(LabelRecord(genome_id, split, family, panel, label, value))

        for trait, (panel, classes) in MULTICLASS_TARGETS.items():
            observed = _string(row.get(trait))
            if not observed:
                continue
            for class_name in classes:
                label = f"{_slug(trait)}__{_slug(class_name)}"
                records.append(LabelRecord(genome_id, split, family, panel, label, int(observed == class_name)))

        for trait, panel in MULTILABEL_TARGETS.items():
            values = _as_observed_mapping(row.get(trait))
            if not values:
                continue
            for class_name in multilabel_classes.get(trait, []):
                if class_name not in values:
                    continue
                label = f"{_slug(trait)}__{_slug(class_name)}"
                records.append(LabelRecord(genome_id, split, family, panel, label, values[class_name]))

    return sorted(records, key=lambda record: (record.genome_id, record.panel, record.label))


def feature_matrix_from_legacy_arrays(
    *,
    bacdive_ids: Iterable[Any],
    feature_names: Iterable[Any],
    rows: Iterable[Iterable[Any]],
    max_features: int | None = 1000,
    min_prevalence: int = 1,
    feature_prefix: str = "eggNOG",
    keep_genome_ids: set[str] | None = None,
) -> FeatureMatrix:
    """Convert a legacy dense feature matrix into the scaffold FeatureMatrix."""

    genome_ids = [_id(value) for value in bacdive_ids]
    names = [str(name) for name in feature_names]
    dense_rows = [[1 if _numeric(cell) > 0 else 0 for cell in row] for row in rows]

    kept_row_indexes = [
        index
        for index, genome_id in enumerate(genome_ids)
        if genome_id and (keep_genome_ids is None or genome_id in keep_genome_ids)
    ]
    if not kept_row_indexes:
        raise ValueError("legacy feature matrix contains no selected genomes")

    prevalences = []
    for feature_index, name in enumerate(names):
        prevalence = sum(dense_rows[row_index][feature_index] for row_index in kept_row_indexes)
        if prevalence >= min_prevalence:
            prevalences.append((feature_index, name, prevalence))
    prevalences.sort(key=lambda item: (-item[2], item[1]))
    if max_features is not None:
        prevalences = prevalences[:max_features]
    selected = sorted(prevalences, key=lambda item: item[1])
    selected_indexes = [item[0] for item in selected]

    return FeatureMatrix(
        genome_ids=[genome_ids[index] for index in kept_row_indexes],
        feature_names=[f"{feature_prefix}:{names[index]}" for index in selected_indexes],
        rows=[[dense_rows[row_index][feature_index] for feature_index in selected_indexes] for row_index in kept_row_indexes],
    )


def load_legacy_npz_feature_matrix(
    npz_path: Path,
    vocab_path: Path,
    *,
    max_features: int | None = 1000,
    min_prevalence: int = 1,
    feature_prefix: str = "eggNOG",
    keep_genome_ids: set[str] | None = None,
) -> FeatureMatrix:
    """Load old `.npz` feature arrays and convert them to FeatureMatrix."""

    try:
        import json
        import numpy as np  # type: ignore[import-not-found]
    except ImportError as exc:
        raise RuntimeError("legacy NPZ import requires optional dependency: numpy") from exc

    data = np.load(npz_path, allow_pickle=True)
    vocab_data = json.loads(vocab_path.read_text())
    feature_names = vocab_data["vocab"] if isinstance(vocab_data, dict) and "vocab" in vocab_data else vocab_data
    bacdive_ids = [_id(value) for value in data["bacdive_ids"].tolist()]
    names = [str(name) for name in feature_names]
    matrix = data["features"]
    kept_row_indexes = [
        index
        for index, genome_id in enumerate(bacdive_ids)
        if genome_id and (keep_genome_ids is None or genome_id in keep_genome_ids)
    ]
    if not kept_row_indexes:
        raise ValueError("legacy feature matrix contains no selected genomes")
    selected_rows = matrix[kept_row_indexes]
    prevalences = (selected_rows > 0).sum(axis=0)
    ranked = [
        (index, names[index], int(prevalence))
        for index, prevalence in enumerate(prevalences.tolist())
        if int(prevalence) >= min_prevalence
    ]
    ranked.sort(key=lambda item: (-item[2], item[1]))
    if max_features is not None:
        ranked = ranked[:max_features]
    selected = sorted(ranked, key=lambda item: item[1])
    selected_indexes = [index for index, _name, _prevalence in selected]
    selected_matrix = (selected_rows[:, selected_indexes] > 0).astype(int)
    return FeatureMatrix(
        genome_ids=[bacdive_ids[index] for index in kept_row_indexes],
        feature_names=[f"{feature_prefix}:{names[index]}" for index in selected_indexes],
        rows=selected_matrix.tolist(),
    )


def _load_delimited(path: Path, *, delimiter: str) -> list[dict[str, Any]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def _select_multilabel_classes(
    rows: list[dict[str, Any]],
    *,
    max_multilabel_classes: int,
) -> dict[str, list[str]]:
    selected = {}
    for trait in MULTILABEL_TARGETS:
        counts: Counter[str] = Counter()
        for row in rows:
            for class_name, value in _as_observed_mapping(row.get(trait)).items():
                if value in {0, 1}:
                    counts[class_name] += 1
        selected[trait] = [name for name, _count in counts.most_common(max_multilabel_classes)]
    return selected


def _as_observed_mapping(value: Any) -> dict[str, int]:
    if not isinstance(value, dict):
        return {}
    out = {}
    for key, raw in value.items():
        label_value = _coerce_binary(raw)
        if label_value is not None:
            out[str(key)] = label_value
    return out


def _coerce_binary(value: Any) -> int | None:
    if _is_missing(value):
        return None
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if value in {0, 1}:
            return int(value)
        return None
    text = str(value).strip().lower()
    if text in {"true", "yes", "positive", "+", "1", "r", "resistant"}:
        return 1
    if text in {"false", "no", "negative", "-", "0", "s", "susceptible"}:
        return 0
    return None


def _id(value: Any) -> str:
    if _is_missing(value):
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value).strip()


def _string(value: Any) -> str:
    if _is_missing(value):
        return ""
    return str(value).strip()


def _numeric(value: Any) -> float:
    if _is_missing(value):
        return 0.0
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def _is_missing(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    return False


def _slug(value: str) -> str:
    slug = re.sub(r"[^a-zA-Z0-9]+", "_", str(value).strip().lower()).strip("_")
    return slug or "unknown"
