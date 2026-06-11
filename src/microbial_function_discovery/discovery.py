"""Discovery candidate annotation and reporting helpers."""

from __future__ import annotations

import json
import csv
import io
from typing import Any

from microbial_function_discovery.datasets import LabelRecord
from microbial_function_discovery.features import FeatureMatrix
from microbial_function_discovery.learning import BaselineModel
from microbial_function_discovery.novelty import AnnotationNoveltyReference


TAXONOMY_FIELDS = ["domain", "phylum", "class", "order", "family", "genus", "species", "type_strain"]
BIOSAFETY_TARGETS = {
    "pathogenicity_human": "biosafety:pathogenicity_human",
    "pathogenicity_animal": "biosafety:pathogenicity_animal",
}
SAFE_LEAD_FIELDS = [
    "panel",
    "target_key",
    "label",
    "source",
    "target_precision",
    "rank",
    "genome_id",
    "species",
    "genus",
    "family",
    "accession",
    "score",
    "risk_level",
    "biosafety_flags",
    "evidence",
    "novelty_level",
]


def annotate_discovery_candidates(
    candidate_report: dict[str, Any],
    trait_rows: list[dict[str, Any]],
    accession_rows: list[dict[str, Any]],
    annotation_model: BaselineModel,
    annotation_features: FeatureMatrix,
    labels: list[LabelRecord],
    *,
    evidence_limit: int = 5,
) -> dict[str, Any]:
    """Attach taxonomy, accessions, annotation evidence, and biosafety flags to rankings."""

    traits_by_id = {_id(row.get("bacdive_id")): row for row in trait_rows if _id(row.get("bacdive_id"))}
    accessions_by_id = {
        _id(row.get("bacdive_id")): row
        for row in accession_rows
        if _id(row.get("bacdive_id"))
    }
    novelty_ref = _build_novelty_reference(labels, annotation_features)
    present_ids = set(annotation_features.genome_ids)
    annotated_targets = []
    for target in candidate_report.get("targets", []):
        target_key = str(target.get("target_key", ""))
        annotated_candidates = []
        for candidate in target.get("candidates", []):
            genome_id = _id(candidate.get("genome_id"))
            trait_row = traits_by_id.get(genome_id, {})
            accession_row = accessions_by_id.get(genome_id, {})
            annotated_candidates.append(
                {
                    **candidate,
                    "genome_id": genome_id,
                    "taxonomy": _taxonomy(trait_row),
                    "genome_accession": _accession(accession_row),
                    "biosafety": _biosafety(
                        genome_id,
                        trait_row,
                        annotation_model,
                        annotation_features,
                    ),
                    "novelty": _novelty(genome_id, novelty_ref, present_ids, annotation_features),
                    "evidence": _candidate_evidence(
                        genome_id,
                        target_key,
                        candidate,
                        annotation_model,
                        annotation_features,
                        limit=evidence_limit,
                    ),
                }
            )
        annotated_targets.append({**target, "candidates": annotated_candidates})

    return {
        **candidate_report,
        "annotation": {
            "trait_records": len(trait_rows),
            "accession_records": len(accession_rows),
            "evidence_limit": evidence_limit,
            "biosafety_targets": sorted(BIOSAFETY_TARGETS.values()),
        },
        "targets": annotated_targets,
        "panels": _panel_summary(annotated_targets),
    }


def render_discovery_candidate_report(annotated_report: dict[str, Any], *, max_targets: int = 25) -> str:
    """Render a compact Markdown report from annotated discovery candidates."""

    risk_counts = _risk_counts(annotated_report)
    lines = [
        "# Discovery Candidate Report",
        "",
        "This report annotates ranked held-out candidates with taxonomy, genome accessions, evidence, and biosafety flags.",
        "",
        "## Summary",
        "",
        f"- Candidate target sets: {annotated_report.get('n_targets', len(annotated_report.get('targets', [])))}",
        f"- Split: {annotated_report.get('split', 'unknown')}",
        f"- Minimum precision filter: {annotated_report.get('min_precision', 'n/a')}",
        f"- Biosafety risk counts: {_risk_counts_text(risk_counts)}",
        "",
        "## Top Candidates",
        "",
        "| Target | Rank | Genome | Species | Source | Score | Safety | Novelty | Evidence |",
        "|---|---:|---|---|---|---:|---|---|---|",
    ]
    for target in annotated_report.get("targets", [])[:max_targets]:
        for candidate in target.get("candidates", [])[:1]:
            taxonomy = candidate.get("taxonomy", {})
            evidence = candidate.get("evidence", [])
            evidence_text = evidence[0]["id"] if evidence else "none"
            lines.append(
                "| "
                + " | ".join(
                    [
                        f"`{target.get('target_key', '')}`",
                        str(candidate.get("rank", "")),
                        str(candidate.get("genome_id", "")),
                        _markdown_cell(taxonomy.get("species") or taxonomy.get("genus") or "unknown"),
                        str(target.get("source", "")),
                        _format_score(candidate.get("score")),
                        str(candidate.get("biosafety", {}).get("risk_level", "unknown")),
                        str(candidate.get("novelty", {}).get("level", "unknown")),
                        _markdown_cell(evidence_text),
                    ]
                )
                + " |"
            )

    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- These are ranked benchmark candidates, not wet-lab validated recommendations.",
            "- Safety flags combine known BacDive fields with annotation-baseline biosafety scores when available.",
            "- Evidence lists active annotation features with positive target weights; dense-only leads may have limited feature evidence.",
            "- Novelty = how unlike the training set a candidate is (annotation-feature Jaccard distance); it is a representativeness signal, NOT a confidence or correctness estimate.",
            "",
        ]
    )
    return "\n".join(lines)


def export_safe_leads(
    annotated_report: dict[str, Any],
    *,
    allowed_risks: list[str] | None = None,
    min_precision: float = 0.5,
    precision_k: int = 10,
    require_accession: bool = True,
    require_evidence: bool = True,
    max_leads: int | None = None,
) -> dict[str, Any]:
    """Flatten annotated candidate rankings into validation-ready lead rows."""

    risks = set(allowed_risks or ["low", "unknown"])
    precision_metric = f"precision_at_{precision_k}"
    leads = []
    for target in annotated_report.get("targets", []):
        target_precision = float(target.get("metrics", {}).get(precision_metric, 0.0))
        if target_precision < min_precision:
            continue
        for candidate in target.get("candidates", []):
            lead = _lead_row(target, candidate, target_precision)
            if lead["risk_level"] not in risks:
                continue
            if require_accession and not lead["accession"]:
                continue
            if require_evidence and not lead["evidence"]:
                continue
            leads.append(lead)

    leads.sort(
        key=lambda row: (
            -float(row["target_precision"]),
            int(row["rank"] or 0),
            row["panel"],
            row["target_key"],
            row["genome_id"],
        )
    )
    if max_leads is not None:
        leads = leads[:max_leads]
    return {
        "precision_metric": precision_metric,
        "min_precision": min_precision,
        "allowed_risks": sorted(risks),
        "require_accession": require_accession,
        "require_evidence": require_evidence,
        "n_leads": len(leads),
        "leads": leads,
        "panels": _lead_panel_summary(leads),
    }


def safe_leads_delimited(export: dict[str, Any], *, delimiter: str = "\t") -> str:
    """Render safe leads as TSV or CSV text."""

    out = io.StringIO()
    writer = csv.DictWriter(out, fieldnames=SAFE_LEAD_FIELDS, delimiter=delimiter, lineterminator="\n")
    writer.writeheader()
    for lead in export.get("leads", []):
        writer.writerow({field: lead.get(field, "") for field in SAFE_LEAD_FIELDS})
    return out.getvalue()


def render_safe_leads_report(export: dict[str, Any], *, max_leads: int = 10) -> str:
    """Render a product-style shortlist report for validation-ready leads."""

    lines = [
        "# Top Validation-Ready Microbial Leads",
        "",
        "This shortlist filters annotated model candidates by target precision, biosafety risk, accession availability, and evidence availability.",
        "",
        "## Summary",
        "",
        f"- Leads exported: {export.get('n_leads', 0)}",
        f"- Precision metric: {export.get('precision_metric', 'unknown')}",
        f"- Minimum precision: {export.get('min_precision', 'n/a')}",
        f"- Allowed risks: {', '.join(export.get('allowed_risks', []))}",
        "",
        "## Leads",
        "",
        "| Panel | Target | Genome | Species | Accession | Score | Risk | Novelty | Evidence |",
        "|---|---|---|---|---|---:|---|---|---|",
    ]
    for lead in export.get("leads", [])[:max_leads]:
        lines.append(
            "| "
            + " | ".join(
                [
                    _markdown_cell(lead.get("panel")),
                    f"`{lead.get('target_key', '')}`",
                    _markdown_cell(lead.get("genome_id")),
                    _markdown_cell(lead.get("species") or lead.get("genus") or "unknown"),
                    _markdown_cell(lead.get("accession")),
                    _format_score(lead.get("score")),
                    _markdown_cell(lead.get("risk_level")),
                    _markdown_cell(lead.get("novelty_level")),
                    _markdown_cell(lead.get("evidence")),
                ]
            )
            + " |"
        )

    lines.extend(["", "## Leads by Panel", ""])
    leads_by_panel: dict[str, list[dict[str, Any]]] = {}
    for lead in export.get("leads", []):
        leads_by_panel.setdefault(str(lead.get("panel", "unknown")), []).append(lead)
    for panel, leads in sorted(leads_by_panel.items()):
        lines.extend([f"### {panel}", ""])
        for lead in leads[:max_leads]:
            lines.append(
                f"- `{lead.get('target_key', '')}`: {lead.get('species') or lead.get('genus') or 'unknown'} "
                f"({lead.get('genome_id')}, {lead.get('accession')}) "
                f"score={_format_score(lead.get('score'))}, risk={lead.get('risk_level')}"
            )
        lines.append("")
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- These leads are benchmark-derived shortlists for review, not recommendations for release or deployment.",
            "- Wet-lab triage should verify biosafety, taxonomy, culturing feasibility, and functional assay design.",
            "",
        ]
    )
    return "\n".join(lines)


def parse_lead_table(text: str, *, delimiter: str = "\t") -> list[dict[str, Any]]:
    """Parse an exported lead table into typed lead rows."""

    reader = csv.DictReader(io.StringIO(text), delimiter=delimiter)
    return [_coerce_lead_table_row(row) for row in reader]


def select_validation_packet_leads(
    leads: list[dict[str, Any]],
    *,
    limit: int = 20,
    min_precision: float = 0.0,
    allowed_risks: list[str] | None = None,
    include_panels: list[str] | None = None,
    exclude_panels: list[str] | None = None,
) -> list[dict[str, Any]]:
    """Select lead rows for wet-lab validation packet rendering."""

    risks = set(allowed_risks or [])
    included = set(include_panels or [])
    excluded = set(exclude_panels or [])
    selected = []
    for lead in leads:
        panel = str(lead.get("panel", ""))
        risk = str(lead.get("risk_level", ""))
        if included and panel not in included:
            continue
        if excluded and panel in excluded:
            continue
        if risks and risk not in risks:
            continue
        if float(lead.get("target_precision", 0.0) or 0.0) < min_precision:
            continue
        selected.append(lead)

    selected.sort(
        key=lambda row: (
            -float(row.get("target_precision", 0.0) or 0.0),
            int(row.get("rank", 999999) or 999999),
            -float(row.get("score", 0.0) or 0.0),
            str(row.get("panel", "")),
            str(row.get("target_key", "")),
            str(row.get("genome_id", "")),
        )
    )
    return selected[: max(limit, 0)]


def render_validation_packets(
    leads: list[dict[str, Any]],
    *,
    title: str = "Wet-Lab Validation Packets",
) -> str:
    """Render selected lead rows as wet-lab validation packets."""

    lines = [
        f"# {title}",
        "",
        "These packets convert model-ranked microbial leads into review units for partner triage.",
        "They are not release or deployment recommendations.",
        "",
        "## Summary",
        "",
        f"- Packets: {len(leads)}",
        f"- Risk mix: {_lead_risk_counts_text(leads)}",
        f"- Panels: {_lead_panels_text(leads)}",
        "",
    ]
    for index, lead in enumerate(leads, start=1):
        species = lead.get("species") or lead.get("genus") or "Unknown microbe"
        label = str(lead.get("label", ""))
        evidence = _evidence_items(lead)
        flags = str(lead.get("biosafety_flags", "") or "none")
        lines.extend(
            [
                f"## Packet {index}: {species}",
                "",
                f"- Genome ID: `{lead.get('genome_id', '')}`",
                f"- Accession: `{lead.get('accession', '')}`",
                f"- Candidate function: {_readable_label(label)}",
                f"- Target key: `{lead.get('target_key', '')}`",
                f"- Application panel: `{lead.get('panel', '')}`",
                f"- Why it matters: {_why_it_matters(lead)}",
                f"- Model support: source `{lead.get('source', '')}`, target precision `{lead.get('target_precision', '')}`, rank `{lead.get('rank', '')}`, score `{lead.get('score', '')}`",
                f"- Evidence to review: {', '.join(f'`{item}`' for item in evidence) if evidence else 'none'}",
                f"- Biosafety review: risk `{lead.get('risk_level', 'unknown')}`, flags `{flags}`",
                f"- First validation step: {_first_validation_step(lead)}",
                "",
            ]
        )
    lines.extend(
        [
            "## Operating Notes",
            "",
            "- Confirm taxonomy, accession metadata, and strain availability before experimental planning.",
            "- Re-check biosafety and AMR/pathogenicity signals before culturing or partner handoff.",
            "- Treat model evidence as prioritization support; require orthogonal assay confirmation.",
            "",
        ]
    )
    return "\n".join(lines)


def _candidate_evidence(
    genome_id: str,
    target_key: str,
    candidate: dict[str, Any],
    model: BaselineModel,
    features: FeatureMatrix,
    *,
    limit: int,
) -> list[dict[str, Any]]:
    evidence = _active_annotation_evidence(model, features, genome_id, target_key, limit=limit)
    evidence.append(
        {
            "type": "model_score",
            "id": "candidate_score",
            "description": "Final discovery ranking score.",
            "score": _round(candidate.get("score")),
        }
    )
    if "annotation_score" in candidate:
        evidence.append(
            {
                "type": "model_score",
                "id": "annotation_score",
                "description": "Annotation model component score.",
                "score": _round(candidate.get("annotation_score")),
            }
        )
    if "dense_score" in candidate:
        evidence.append(
            {
                "type": "model_score",
                "id": "dense_score",
                "description": "Dense embedding model component score.",
                "score": _round(candidate.get("dense_score")),
            }
        )
    return evidence[: max(limit, 1)]


def _build_novelty_reference(labels: list[LabelRecord], features: FeatureMatrix):
    """Reference = split=="train" labeled genomes (fall back to all labels)."""
    train_ids = [rec.genome_id for rec in labels if getattr(rec, "split", "") == "train"]
    if not train_ids:
        train_ids = [rec.genome_id for rec in labels]
    try:
        return AnnotationNoveltyReference().fit(features, train_ids)
    except ValueError:
        return None


def _novelty(genome_id: str, novelty_ref, present_ids: set, features: FeatureMatrix) -> dict[str, Any]:
    if novelty_ref is None or genome_id not in present_ids:
        return {"score": None, "level": "unknown", "ref_percentile": None}
    score = novelty_ref.score(features.row_for(genome_id))
    return {
        "score": round(score, 6),
        "level": novelty_ref.level(score),
        "ref_percentile": round(novelty_ref.ref_percentile(score), 6),
    }


def _lead_row(target: dict[str, Any], candidate: dict[str, Any], target_precision: float) -> dict[str, Any]:
    taxonomy = candidate.get("taxonomy", {})
    accession = candidate.get("genome_accession", {})
    biosafety = candidate.get("biosafety", {})
    return {
        "panel": str(target.get("panel", "")),
        "target_key": str(target.get("target_key", "")),
        "label": str(target.get("label", "")),
        "source": str(target.get("source", "")),
        "target_precision": round(target_precision, 6),
        "rank": candidate.get("rank", ""),
        "genome_id": str(candidate.get("genome_id", "")),
        "species": _clean_value(taxonomy.get("species")) or "",
        "genus": _clean_value(taxonomy.get("genus")) or "",
        "family": _clean_value(taxonomy.get("family")) or "",
        "accession": _clean_value(accession.get("accession")) or "",
        "score": _round(candidate.get("score")),
        "risk_level": str(biosafety.get("risk_level", "unknown")),
        "biosafety_flags": ",".join(str(flag) for flag in biosafety.get("flags", [])),
        "evidence": ";".join(str(item.get("id", "")) for item in candidate.get("evidence", []) if item.get("id")),
        "novelty_level": str(candidate.get("novelty", {}).get("level", "unknown")),
    }


def _lead_panel_summary(leads: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    panels: dict[str, dict[str, Any]] = {}
    for lead in leads:
        summary = panels.setdefault(lead["panel"], {"n_leads": 0, "targets": []})
        summary["n_leads"] += 1
        if lead["target_key"] not in summary["targets"]:
            summary["targets"].append(lead["target_key"])
    return panels


def _coerce_lead_table_row(row: dict[str, str]) -> dict[str, Any]:
    coerced: dict[str, Any] = {}
    for key, value in row.items():
        if key in {"target_precision", "score"}:
            coerced[key] = float(value) if value not in ("", None) else 0.0
        elif key == "rank":
            coerced[key] = int(value) if value not in ("", None) else 0
        else:
            coerced[key] = value
    return coerced


def _lead_risk_counts_text(leads: list[dict[str, Any]]) -> str:
    counts: dict[str, int] = {}
    for lead in leads:
        risk = str(lead.get("risk_level", "unknown") or "unknown")
        counts[risk] = counts.get(risk, 0) + 1
    return _risk_counts_text(counts)


def _lead_panels_text(leads: list[dict[str, Any]]) -> str:
    panels = sorted({str(lead.get("panel", "unknown") or "unknown") for lead in leads})
    return ", ".join(panels) if panels else "none"


def _evidence_items(lead: dict[str, Any]) -> list[str]:
    return [item for item in str(lead.get("evidence", "") or "").split(";") if item]


def _readable_label(label: str) -> str:
    return label.replace("__", ": ").replace("_", " ")


def _why_it_matters(lead: dict[str, Any]) -> str:
    panel = str(lead.get("panel", ""))
    label = _readable_label(str(lead.get("label", "")))
    if panel == "biofuels_industrial":
        return f"{label} can prioritize strains or enzymes for feedstock conversion and industrial biocatalysis."
    if panel == "environmental_terraforming":
        return f"{label} can prioritize organisms for environmental stress response, remediation, or closed-system bioprocessing assays."
    if panel == "food_fermentation_agriculture":
        return f"{label} can prioritize microbes for fermentation, metabolite production, food, or agriculture screens."
    if panel == "therapeutics_antimicrobials":
        return f"{label} can prioritize microbial products or interactions for therapeutic and antimicrobial discovery."
    if panel == "biosafety":
        return f"{label} is useful for safety characterization and go/no-go triage before application testing."
    return f"{label} is a targetable function for follow-up microbial screening."


def _first_validation_step(lead: dict[str, Any]) -> str:
    label = str(lead.get("label", ""))
    panel = str(lead.get("panel", ""))
    if "carbon_utilization" in label:
        substrate = label.split("__")[-1].replace("_", " ")
        return f"Run a growth or activity assay with {substrate} as the target carbon substrate."
    if "metabolite_production" in label:
        metabolite = label.split("__")[-1].replace("_", " ")
        return f"Run targeted {metabolite} quantification under matched culture or enrichment conditions."
    if "catalase" in label:
        return "Run catalase and peroxide-stress assays, then compare activity against close relatives."
    if panel == "biosafety":
        return "Confirm the predicted safety phenotype with strain-level metadata and standard biosafety screens."
    return "Design a targeted phenotype assay and confirm the supporting annotation evidence with an orthogonal method."


def _risk_counts(annotated_report: dict[str, Any]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for target in annotated_report.get("targets", []):
        for candidate in target.get("candidates", []):
            risk = candidate.get("biosafety", {}).get("risk_level", "unknown")
            counts[risk] = counts.get(risk, 0) + 1
    return counts


def _risk_counts_text(counts: dict[str, int]) -> str:
    if not counts:
        return "none"
    return ", ".join(f"{key}={counts[key]}" for key in sorted(counts))


def _active_annotation_evidence(
    model: BaselineModel,
    features: FeatureMatrix,
    genome_id: str,
    target_key: str,
    *,
    limit: int,
) -> list[dict[str, Any]]:
    if target_key not in model.targets:
        return []
    try:
        row = _align_row(features, genome_id, model.feature_names)
    except KeyError:
        return []
    weights = model.targets[target_key].feature_log_odds
    evidence = []
    for feature_name, value, weight in sorted(zip(model.feature_names, row, weights), key=lambda item: item[2], reverse=True):
        if not value or weight <= 0:
            continue
        evidence.append(
            {
                "type": "annotation_feature",
                "id": feature_name,
                "description": f"Active feature supports {target_key}.",
                "weight": round(weight, 6),
            }
        )
    return evidence[:limit]


def _biosafety(
    genome_id: str,
    trait_row: dict[str, Any],
    model: BaselineModel,
    features: FeatureMatrix,
) -> dict[str, Any]:
    known_bsl = _string(trait_row.get("biosafety_level"))
    known_human = _boolish(trait_row.get("pathogenicity_human"))
    known_animal = _boolish(trait_row.get("pathogenicity_animal"))
    human_score = _predict_optional(model, features, genome_id, BIOSAFETY_TARGETS["pathogenicity_human"])
    animal_score = _predict_optional(model, features, genome_id, BIOSAFETY_TARGETS["pathogenicity_animal"])
    amr_score = _max_prefix_score(model, features, genome_id, "biosafety:amr_phenotype__")
    risk_level = _risk_level(known_bsl, known_human, known_animal, human_score, animal_score, amr_score)
    flags = []
    if known_human:
        flags.append("known_human_pathogen")
    if known_animal:
        flags.append("known_animal_pathogen")
    if known_bsl and known_bsl not in {"BSL-1", "1"}:
        flags.append(f"known_{known_bsl.lower()}")
    if human_score is not None and human_score >= 0.5:
        flags.append("predicted_human_pathogenicity")
    if amr_score is not None and amr_score >= 0.5:
        flags.append("predicted_amr_signal")
    return {
        "risk_level": risk_level,
        "flags": flags,
        "known_biosafety_level": known_bsl,
        "known_pathogenicity_human": known_human,
        "known_pathogenicity_animal": known_animal,
        "predicted_pathogenicity_human": _round(human_score),
        "predicted_pathogenicity_animal": _round(animal_score),
        "predicted_amr_max": _round(amr_score),
    }


def _risk_level(
    known_bsl: str | None,
    known_human: bool | None,
    known_animal: bool | None,
    human_score: float | None,
    animal_score: float | None,
    amr_score: float | None,
) -> str:
    if known_human or known_animal or known_bsl in {"BSL-3", "BSL-4", "3", "4"}:
        return "high"
    if known_bsl == "BSL-2" or any(score is not None and score >= 0.5 for score in [human_score, animal_score, amr_score]):
        return "moderate"
    if known_bsl == "BSL-1" or known_bsl == "1":
        return "low"
    return "unknown"


def _predict_optional(model: BaselineModel, features: FeatureMatrix, genome_id: str, target_key: str) -> float | None:
    if target_key not in model.targets:
        return None
    try:
        return model.predict_proba(target_key, genome_id, features)
    except KeyError:
        return None


def _max_prefix_score(model: BaselineModel, features: FeatureMatrix, genome_id: str, prefix: str) -> float | None:
    scores = [
        score
        for target_key in sorted(model.targets)
        if target_key.startswith(prefix)
        for score in [_predict_optional(model, features, genome_id, target_key)]
        if score is not None
    ]
    return max(scores) if scores else None


def _taxonomy(row: dict[str, Any]) -> dict[str, Any]:
    return {field: _clean_value(row.get(field)) for field in TAXONOMY_FIELDS}


def _accession(row: dict[str, Any]) -> dict[str, Any]:
    return {
        "accession": _clean_value(row.get("accession")),
        "assembly_level": _clean_value(row.get("assembly_level")),
        "description": _clean_value(row.get("description")),
    }


def _panel_summary(targets: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    panels: dict[str, dict[str, Any]] = {}
    for target in targets:
        panel = str(target.get("panel", "unknown"))
        summary = panels.setdefault(panel, {"n_targets": 0, "n_candidates": 0, "risk_counts": {}})
        summary["n_targets"] += 1
        summary["n_candidates"] += len(target.get("candidates", []))
        for candidate in target.get("candidates", []):
            risk = candidate.get("biosafety", {}).get("risk_level", "unknown")
            summary["risk_counts"][risk] = summary["risk_counts"].get(risk, 0) + 1
    return panels


def _align_row(features: FeatureMatrix, genome_id: str, feature_names: list[str]) -> list[int]:
    source_row = features.row_for(genome_id)
    index_by_name = {name: index for index, name in enumerate(features.feature_names)}
    return [source_row[index_by_name[name]] if name in index_by_name else 0 for name in feature_names]


def _id(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value).strip()


def _clean_value(value: Any) -> Any:
    if value is None:
        return None
    if isinstance(value, float) and value != value:
        return None
    if isinstance(value, (dict, list)):
        return value
    text = str(value).strip()
    return text if text else None


def _string(value: Any) -> str | None:
    clean = _clean_value(value)
    return str(clean) if clean is not None else None


def _boolish(value: Any) -> bool | None:
    if value is None:
        return None
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"true", "1", "yes", "y"}:
        return True
    if text in {"false", "0", "no", "n"}:
        return False
    return None


def _round(value: Any) -> float | None:
    if value is None:
        return None
    try:
        return round(float(value), 6)
    except (TypeError, ValueError):
        return None


def _format_score(value: Any) -> str:
    rounded = _round(value)
    return "" if rounded is None else f"{rounded:.3f}"


def _markdown_cell(value: Any) -> str:
    text = "" if value is None else str(value)
    return text.replace("|", "\\|").replace("\n", " ")


def to_json(data: dict[str, Any]) -> str:
    return json.dumps(data, indent=2, sort_keys=True)
