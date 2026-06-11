"""Discovery candidate annotation and reporting helpers."""

from __future__ import annotations

import json
from typing import Any

from microbial_function_discovery.datasets import LabelRecord
from microbial_function_discovery.features import FeatureMatrix
from microbial_function_discovery.learning import BaselineModel


TAXONOMY_FIELDS = ["domain", "phylum", "class", "order", "family", "genus", "species", "type_strain"]
BIOSAFETY_TARGETS = {
    "pathogenicity_human": "biosafety:pathogenicity_human",
    "pathogenicity_animal": "biosafety:pathogenicity_animal",
}


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
        "| Target | Rank | Genome | Species | Source | Score | Safety | Evidence |",
        "|---|---:|---|---|---|---:|---|---|",
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
