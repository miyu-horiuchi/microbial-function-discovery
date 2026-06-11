"""FASTA parsing helpers."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class FastaRecord:
    """One FASTA record."""

    identifier: str
    description: str
    sequence: str


def parse_fasta(text: str) -> list[FastaRecord]:
    """Parse FASTA text into records.

    The parser accepts protein or nucleotide-like sequences. It removes common
    separators and stop markers so tiny examples and copied FASTA snippets work
    without extra preprocessing.
    """

    records: list[FastaRecord] = []
    current_id = "sequence_1"
    current_description = ""
    sequence_parts: list[str] = []
    saw_header = False

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if sequence_parts:
                records.append(
                    FastaRecord(
                        identifier=current_id,
                        description=current_description,
                        sequence=_clean_sequence("".join(sequence_parts)),
                    )
                )
            current_id, current_description = _parse_header(line[1:].strip(), len(records) + 1)
            sequence_parts = []
            saw_header = True
            continue
        sequence_parts.append(line)

    if sequence_parts:
        records.append(
            FastaRecord(
                identifier=current_id if saw_header else "sequence_1",
                description=current_description if saw_header else "",
                sequence=_clean_sequence("".join(sequence_parts)),
            )
        )

    records = [record for record in records if record.sequence]
    if not records:
        raise ValueError("no FASTA records found")
    return records


def _parse_header(header: str, fallback_index: int) -> tuple[str, str]:
    if not header:
        return f"sequence_{fallback_index}", ""
    parts = header.split(maxsplit=1)
    identifier = parts[0]
    description = parts[1] if len(parts) == 2 else ""
    return identifier, description


def _clean_sequence(sequence: str) -> str:
    return "".join(char for char in sequence.upper() if char.isalpha())
