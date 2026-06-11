"""Benchmark label parsing and split validation."""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass


LABEL_COLUMNS = {"genome_id", "split", "family", "panel", "label", "value"}


class FamilyLeakageError(ValueError):
    """Raised when the same family appears in multiple splits."""


@dataclass(frozen=True)
class LabelRecord:
    genome_id: str
    split: str
    family: str
    panel: str
    label: str
    value: int

    @property
    def target_key(self) -> str:
        return f"{self.panel}:{self.label}"


def parse_labels_tsv(text: str) -> list[LabelRecord]:
    """Parse benchmark labels from TSV text."""

    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    if reader.fieldnames is None:
        raise ValueError("labels TSV is empty")
    missing = LABEL_COLUMNS - set(reader.fieldnames)
    if missing:
        raise ValueError(f"missing required columns: {', '.join(sorted(missing))}")

    records: list[LabelRecord] = []
    for row_index, row in enumerate(reader, start=2):
        value_text = _require_cell(row, "value", row_index)
        if value_text not in {"0", "1"}:
            raise ValueError(f"line {row_index}: value must be 0 or 1")
        records.append(
            LabelRecord(
                genome_id=_require_cell(row, "genome_id", row_index),
                split=_require_cell(row, "split", row_index),
                family=_require_cell(row, "family", row_index),
                panel=_require_cell(row, "panel", row_index),
                label=_require_cell(row, "label", row_index),
                value=int(value_text),
            )
        )
    if not records:
        raise ValueError("labels TSV contains no records")
    return records


def validate_family_holdout(records: list[LabelRecord]) -> None:
    """Ensure each family belongs to exactly one split."""

    family_to_splits: dict[str, set[str]] = {}
    for record in records:
        family_to_splits.setdefault(record.family, set()).add(record.split)

    leaked = {family: splits for family, splits in family_to_splits.items() if len(splits) > 1}
    if leaked:
        details = ", ".join(f"{family}={sorted(splits)}" for family, splits in sorted(leaked.items()))
        raise FamilyLeakageError(f"families appear in multiple splits: {details}")


def _require_cell(row: dict[str, str | None], column: str, row_index: int) -> str:
    value = row.get(column)
    if value is None or not value.strip():
        raise ValueError(f"line {row_index}: {column} is required")
    return value.strip()
