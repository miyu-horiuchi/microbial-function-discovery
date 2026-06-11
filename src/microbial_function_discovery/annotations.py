"""Structured annotation-hit parsing."""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass


REQUIRED_COLUMNS = {"protein_id", "database", "accession", "name", "evalue"}


@dataclass(frozen=True)
class AnnotationHit:
    """One functional annotation hit for a protein."""

    protein_id: str
    database: str
    accession: str
    name: str
    evalue: float | None = None


def parse_annotation_hits_tsv(text: str) -> list[AnnotationHit]:
    """Parse a tab-separated annotation-hit table.

    Required columns: protein_id, database, accession, name, evalue.
    """

    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    if reader.fieldnames is None:
        raise ValueError("annotation TSV is empty")
    missing = REQUIRED_COLUMNS - set(reader.fieldnames)
    if missing:
        raise ValueError(f"missing required columns: {', '.join(sorted(missing))}")

    hits: list[AnnotationHit] = []
    for row_index, row in enumerate(reader, start=2):
        protein_id = _require_cell(row, "protein_id", row_index)
        database = _require_cell(row, "database", row_index)
        accession = _require_cell(row, "accession", row_index)
        name = _require_cell(row, "name", row_index)
        evalue_text = _require_cell(row, "evalue", row_index)
        try:
            evalue = float(evalue_text)
        except ValueError as exc:
            raise ValueError(f"line {row_index}: evalue must be numeric") from exc
        hits.append(
            AnnotationHit(
                protein_id=protein_id,
                database=database,
                accession=accession,
                name=name,
                evalue=evalue,
            )
        )

    if not hits:
        raise ValueError("annotation TSV contains no hits")
    return hits


def _require_cell(row: dict[str, str | None], column: str, row_index: int) -> str:
    value = row.get(column)
    if value is None or not value.strip():
        raise ValueError(f"line {row_index}: {column} is required")
    return value.strip()
