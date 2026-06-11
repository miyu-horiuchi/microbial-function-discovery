"""Feature matrix construction from normalized annotation hits."""

from __future__ import annotations

import csv
import io
import json
from dataclasses import dataclass
from typing import Any


FEATURE_COLUMNS = {"genome_id", "protein_id", "database", "accession", "name", "evalue"}


@dataclass(frozen=True)
class FeatureMatrix:
    genome_ids: list[str]
    feature_names: list[str]
    rows: list[list[int]]

    def to_json(self) -> str:
        return json.dumps(
            {
                "genome_ids": self.genome_ids,
                "feature_names": self.feature_names,
                "rows": self.rows,
            },
            indent=2,
            sort_keys=True,
        )

    @classmethod
    def from_json(cls, text: str) -> "FeatureMatrix":
        data = json.loads(text)
        return cls(
            genome_ids=list(data["genome_ids"]),
            feature_names=list(data["feature_names"]),
            rows=[list(map(int, row)) for row in data["rows"]],
        )

    def row_for(self, genome_id: str) -> list[int]:
        try:
            index = self.genome_ids.index(genome_id)
        except ValueError as exc:
            raise KeyError(f"genome id not found in feature matrix: {genome_id}") from exc
        return self.rows[index]

    def to_dict(self) -> dict[str, Any]:
        return {
            "genome_ids": self.genome_ids,
            "feature_names": self.feature_names,
            "rows": self.rows,
        }


def build_feature_matrix_from_annotation_tsv(text: str) -> FeatureMatrix:
    """Build binary database/accession features from a multi-genome annotation TSV."""

    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    if reader.fieldnames is None:
        raise ValueError("annotation TSV is empty")
    missing = FEATURE_COLUMNS - set(reader.fieldnames)
    if missing:
        raise ValueError(f"missing required columns: {', '.join(sorted(missing))}")

    genome_features: dict[str, set[str]] = {}
    for row_index, row in enumerate(reader, start=2):
        genome_id = _require_cell(row, "genome_id", row_index)
        database = _require_cell(row, "database", row_index)
        accession = _require_cell(row, "accession", row_index)
        genome_features.setdefault(genome_id, set()).add(f"{database}:{accession}")

    if not genome_features:
        raise ValueError("annotation TSV contains no hits")

    genome_ids = sorted(genome_features)
    feature_names = sorted({feature for features in genome_features.values() for feature in features})
    rows = [[1 if feature in genome_features[genome_id] else 0 for feature in feature_names] for genome_id in genome_ids]
    return FeatureMatrix(genome_ids=genome_ids, feature_names=feature_names, rows=rows)


def _require_cell(row: dict[str, str | None], column: str, row_index: int) -> str:
    value = row.get(column)
    if value is None or not value.strip():
        raise ValueError(f"line {row_index}: {column} is required")
    return value.strip()
