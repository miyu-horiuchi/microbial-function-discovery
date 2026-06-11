"""Utilities for browsing exported microbial discovery leads."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path
from typing import Any


NUMERIC_FIELDS = {"score", "target_precision"}
INTEGER_FIELDS = {"rank"}
QUERY_FIELDS = (
    "species",
    "genus",
    "family",
    "genome_id",
    "accession",
    "target_key",
    "label",
    "evidence",
)


def load_leads(path: str | Path) -> list[dict[str, Any]]:
    """Load an exported lead TSV while preserving IDs as strings."""

    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [_coerce_row(row) for row in reader]


def option_values(leads: list[dict[str, Any]], field: str) -> list[str]:
    """Return sorted non-empty option values for a lead field."""

    return sorted({str(lead.get(field, "")) for lead in leads if lead.get(field, "") != ""})


def summarize_leads(leads: list[dict[str, Any]]) -> dict[str, Any]:
    """Summarize a lead set for the demo header."""

    return {
        "n_leads": len(leads),
        "n_panels": len(option_values(leads, "panel")),
        "n_species": len(option_values(leads, "species")),
        "risk_counts": dict(sorted(Counter(str(lead.get("risk_level", "")) for lead in leads).items())),
        "panel_counts": dict(sorted(Counter(str(lead.get("panel", "")) for lead in leads).items())),
    }


def filter_leads(
    leads: list[dict[str, Any]],
    *,
    panel: str | None = "All",
    risk: str | None = "All",
    source: str | None = "All",
    query: str | None = "",
    min_precision: float = 0.0,
) -> list[dict[str, Any]]:
    """Filter and rank leads for review."""

    query_normalized = (query or "").strip().casefold()
    selected = []
    for lead in leads:
        if not _matches_option(lead, "panel", panel):
            continue
        if not _matches_option(lead, "risk_level", risk):
            continue
        if not _matches_option(lead, "source", source):
            continue
        if float(lead.get("target_precision", 0.0) or 0.0) < min_precision:
            continue
        if query_normalized and not _matches_query(lead, query_normalized):
            continue
        selected.append(lead)

    return sorted(
        selected,
        key=lambda lead: (
            -float(lead.get("target_precision", 0.0) or 0.0),
            -float(lead.get("score", 0.0) or 0.0),
            int(lead.get("rank", 999999) or 999999),
            str(lead.get("species", "")),
            str(lead.get("genome_id", "")),
        ),
    )


def table_rows(leads: list[dict[str, Any]]) -> list[list[Any]]:
    """Convert lead dictionaries into stable table rows for Gradio."""

    columns = table_columns()
    return [[lead.get(column, "") for column in columns] for lead in leads]


def table_columns() -> list[str]:
    """Columns shown in the demo table."""

    return [
        "panel",
        "target_key",
        "species",
        "accession",
        "target_precision",
        "score",
        "risk_level",
        "source",
        "evidence",
    ]


def lead_detail(lead: dict[str, Any] | None) -> str:
    """Render one lead as a compact Markdown evidence card."""

    if not lead:
        return "No lead selected."

    flags = str(lead.get("biosafety_flags", "") or "none")
    return "\n".join(
        [
            f"### {lead.get('species', 'Unknown species')}",
            "",
            f"- Genome: `{lead.get('genome_id', '')}`",
            f"- Accession: `{lead.get('accession', '')}`",
            f"- Target: `{lead.get('target_key', '')}`",
            f"- Panel: `{lead.get('panel', '')}`",
            f"- Source: `{lead.get('source', '')}`",
            f"- Precision: `{lead.get('target_precision', '')}`",
            f"- Score: `{lead.get('score', '')}`",
            f"- Risk: `{lead.get('risk_level', '')}`",
            f"- Biosafety flags: `{flags}`",
            f"- Evidence: `{lead.get('evidence', '')}`",
        ]
    )


def summary_markdown(leads: list[dict[str, Any]], dataset_name: str) -> str:
    """Render a compact summary for the current filtered result set."""

    summary = summarize_leads(leads)
    risk = ", ".join(f"{key}: {value}" for key, value in summary["risk_counts"].items()) or "none"
    return (
        f"**{dataset_name}**  \n"
        f"{summary['n_leads']} leads across {summary['n_panels']} panels and "
        f"{summary['n_species']} species.  \n"
        f"Risk mix: {risk}."
    )


def _coerce_row(row: dict[str, str]) -> dict[str, Any]:
    coerced: dict[str, Any] = {}
    for key, value in row.items():
        if key in NUMERIC_FIELDS:
            coerced[key] = float(value) if value not in ("", None) else 0.0
        elif key in INTEGER_FIELDS:
            coerced[key] = int(value) if value not in ("", None) else 0
        else:
            coerced[key] = value
    return coerced


def _matches_option(lead: dict[str, Any], field: str, value: str | None) -> bool:
    if value in (None, "", "All"):
        return True
    return str(lead.get(field, "")) == value


def _matches_query(lead: dict[str, Any], query: str) -> bool:
    return any(query in str(lead.get(field, "")).casefold() for field in QUERY_FIELDS)

