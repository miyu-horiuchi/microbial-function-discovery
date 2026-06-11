"""Import real annotation tool outputs into normalized annotation hits."""

from __future__ import annotations

import csv
import io

from microbial_function_discovery.annotations import AnnotationHit


def parse_eggnog_mapper(text: str) -> list[AnnotationHit]:
    """Parse eggNOG-mapper annotation TSV text.

    Supports files with a commented header line such as
    ``#query ... Description ... KEGG_ko ... CAZy ... PFAMs``.
    """

    lines = [line for line in text.splitlines() if line.strip()]
    header_line = next((line for line in lines if line.startswith("#") and "\t" in line), None)
    if header_line is None:
        raise ValueError("eggNOG mapper output must include a tab-separated header")

    header = [column.lstrip("#") for column in header_line.split("\t")]
    data_lines = [line for line in lines if not line.startswith("#")]
    reader = csv.DictReader(io.StringIO("\n".join(data_lines)), delimiter="\t", fieldnames=header)

    hits: list[AnnotationHit] = []
    for row in reader:
        protein_id = _first_present(row, "query", "query_name")
        description = _clean_missing(row.get("Description")) or _clean_missing(row.get("Preferred_name")) or "eggNOG annotation"
        evalue = _parse_evalue(row.get("evalue"))

        for ko in _split_multi_value(row.get("KEGG_ko")):
            hits.append(
                AnnotationHit(
                    protein_id=protein_id,
                    database="KEGG",
                    accession=ko.replace("ko:", ""),
                    name=description,
                    evalue=evalue,
                )
            )
        for cazy in _split_multi_value(row.get("CAZy")):
            hits.append(
                AnnotationHit(
                    protein_id=protein_id,
                    database="CAZy",
                    accession=cazy,
                    name=description,
                    evalue=evalue,
                )
            )
        for pfam in _split_multi_value(row.get("PFAMs")):
            hits.append(
                AnnotationHit(
                    protein_id=protein_id,
                    database="Pfam",
                    accession=pfam,
                    name=description,
                    evalue=evalue,
                )
            )

    if not hits:
        raise ValueError("eggNOG mapper output contains no importable KEGG, CAZy, or Pfam hits")
    return hits


def parse_hmmer_domtblout(text: str, *, database: str) -> list[AnnotationHit]:
    """Parse HMMER ``--domtblout`` output.

    The normalized hit uses the query name as ``protein_id`` and the target
    accession/name as the imported domain/database evidence.
    """

    hits: list[AnnotationHit] = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(maxsplit=22)
        if len(parts) < 22:
            raise ValueError("domtblout line has fewer than 22 required fields")
        target_name = parts[0]
        target_accession = parts[1]
        query_name = parts[3]
        evalue = _parse_evalue(parts[6])
        description = parts[22] if len(parts) > 22 else target_name
        hits.append(
            AnnotationHit(
                protein_id=query_name,
                database=database,
                accession=target_accession if target_accession != "-" else target_name,
                name=description,
                evalue=evalue,
            )
        )

    if not hits:
        raise ValueError("domtblout contains no hits")
    return hits


def format_annotation_hits_tsv(hits: list[AnnotationHit]) -> str:
    """Format normalized annotation hits as the project TSV contract."""

    output = io.StringIO()
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(["protein_id", "database", "accession", "name", "evalue"])
    for hit in hits:
        writer.writerow(
            [
                hit.protein_id,
                hit.database,
                hit.accession,
                hit.name,
                "" if hit.evalue is None else f"{hit.evalue:.6g}",
            ]
        )
    return output.getvalue()


def _first_present(row: dict[str, str | None], *keys: str) -> str:
    for key in keys:
        value = _clean_missing(row.get(key))
        if value:
            return value
    raise ValueError(f"missing required field: {'/'.join(keys)}")


def _clean_missing(value: str | None) -> str | None:
    if value is None:
        return None
    value = value.strip()
    if not value or value == "-":
        return None
    return value


def _split_multi_value(value: str | None) -> list[str]:
    cleaned = _clean_missing(value)
    if cleaned is None:
        return []
    parts = []
    for chunk in cleaned.replace(",", ";").split(";"):
        chunk = chunk.strip()
        if chunk and chunk != "-":
            parts.append(chunk)
    return parts


def _parse_evalue(value: str | None) -> float | None:
    cleaned = _clean_missing(value)
    if cleaned is None:
        return None
    return float(cleaned)
