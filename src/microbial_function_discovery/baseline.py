"""First-pass evidence baseline for microbial useful-function discovery."""

from __future__ import annotations

from dataclasses import dataclass

from microbial_function_discovery.annotations import AnnotationHit
from microbial_function_discovery.fasta import FastaRecord, parse_fasta
from microbial_function_discovery.panels import APPLICATION_AREAS


@dataclass(frozen=True)
class FunctionRule:
    name: str
    panel: str
    keywords: tuple[str, ...]
    evidence_database: str
    evidence_type: str = "protein"
    biosafety_kind: str | None = None


FUNCTION_RULES: tuple[FunctionRule, ...] = (
    FunctionRule(
        name="cellulose degradation",
        panel="biofuels_industrial",
        keywords=("cellulase", "glycoside hydrolase", "gh5", "gh6", "gh9", "endoglucanase"),
        evidence_database="CAZy",
    ),
    FunctionRule(
        name="xylan degradation",
        panel="biofuels_industrial",
        keywords=("xylanase", "gh10", "gh11", "beta-xylosidase"),
        evidence_database="CAZy",
    ),
    FunctionRule(
        name="lignin degradation",
        panel="biofuels_industrial",
        keywords=("laccase", "peroxidase", "lignin"),
        evidence_database="KEGG",
    ),
    FunctionRule(
        name="nitrogen fixation",
        panel="environmental_terraforming",
        keywords=("nitrogenase", "nifh", "nifd", "nifk"),
        evidence_database="KEGG",
    ),
    FunctionRule(
        name="carbon fixation",
        panel="environmental_terraforming",
        keywords=("rubisco", "rbcl", "rbcs", "carbon monoxide dehydrogenase"),
        evidence_database="KEGG",
    ),
    FunctionRule(
        name="sulfur metabolism",
        panel="environmental_terraforming",
        keywords=("sulfite reductase", "sulfate adenylyltransferase", "sox", "dsra", "dsrb"),
        evidence_database="KEGG",
    ),
    FunctionRule(
        name="biosynthetic gene cluster potential",
        panel="therapeutics_antimicrobials",
        keywords=("polyketide synthase", "nonribosomal peptide", "nrps", "pks", "terpene synthase"),
        evidence_database="antiSMASH",
    ),
    FunctionRule(
        name="antimicrobial peptide potential",
        panel="therapeutics_antimicrobials",
        keywords=("bacteriocin", "lantibiotic", "microcin"),
        evidence_database="BAGEL",
    ),
    FunctionRule(
        name="fermentation traits",
        panel="food_fermentation_agriculture",
        keywords=("lactate dehydrogenase", "alcohol dehydrogenase", "acetate kinase", "pyruvate decarboxylase"),
        evidence_database="KEGG",
    ),
    FunctionRule(
        name="plant growth promotion",
        panel="food_fermentation_agriculture",
        keywords=("indole-3-acetic acid", "phosphate solubilization", "siderophore", "acc deaminase"),
        evidence_database="MetaCyc",
    ),
    FunctionRule(
        name="antimicrobial resistance",
        panel="biosafety",
        keywords=("beta-lactamase", "antimicrobial resistance", "multidrug efflux", "aminoglycoside resistance"),
        evidence_database="CARD",
        biosafety_kind="amr",
    ),
    FunctionRule(
        name="virulence-associated factors",
        panel="biosafety",
        keywords=("virulence", "toxin", "adhesin", "invasin", "type iii secretion", "hemolysin"),
        evidence_database="VFDB",
        biosafety_kind="pathogenicity",
    ),
)


def predict_from_fasta(text: str, *, genome_id: str = "input_genome") -> dict:
    """Predict useful functions from protein FASTA descriptions.

    This is a transparent baseline, not the final learned foundation model. It
    turns recognizable annotations into the same JSON contract the future model
    will emit.
    """

    records = parse_fasta(text)
    function_hits = _score_rules(records)

    return _prediction_from_function_hits(
        function_hits,
        genome_id=genome_id,
        novelty_note="Baseline prediction uses annotation keywords only; learned cross-clade calibration is not available yet.",
    )


def predict_from_annotation_hits(hits: list[AnnotationHit], *, genome_id: str = "input_genome") -> dict:
    """Predict useful functions from structured annotation hits."""

    function_hits = _score_annotation_rules(hits)
    return _prediction_from_function_hits(
        function_hits,
        genome_id=genome_id,
        novelty_note="Baseline prediction uses structured annotation hits; learned cross-clade calibration is not available yet.",
    )


def _prediction_from_function_hits(function_hits: list[dict], *, genome_id: str, novelty_note: str) -> dict:

    functions = [_function_to_json(hit) for hit in function_hits]
    application_scores = _application_scores(function_hits)
    biosafety = _biosafety(function_hits)

    if not functions:
        functions = [
            {
                "name": "unresolved function",
                "score": 0.05,
                "confidence": 0.2,
                "evidence": [
                    {
                        "type": "database_hit",
                        "id": "no_rule_hits",
                        "annotation": "No baseline keyword evidence matched this input.",
                        "database": "baseline",
                        "weight": 0.05,
                    }
                ],
            }
        ]

    return {
        "genome_id": genome_id,
        "application_scores": application_scores,
        "functions": functions,
        "biosafety": biosafety,
        "novelty": {
            "nearest_training_family": "unknown",
            "generalization_risk": "unknown",
            "notes": [
                novelty_note
            ],
        },
    }


def _score_rules(records: list[FastaRecord]) -> list[dict]:
    function_hits = [_score_rule(rule, records) for rule in FUNCTION_RULES]
    return [hit for hit in function_hits if hit is not None]


def _score_rule(rule: FunctionRule, records: list[FastaRecord]) -> dict | None:
    evidence = []
    for record in records:
        haystack = f"{record.identifier} {record.description} {record.sequence}".lower()
        matched = [keyword for keyword in rule.keywords if keyword in haystack]
        if not matched:
            continue
        evidence.append(
            {
                "type": rule.evidence_type,
                "id": record.identifier,
                "annotation": f"{rule.name}: matched {', '.join(matched)}",
                "database": rule.evidence_database,
                "weight": round(min(1.0, 0.35 + 0.15 * len(matched)), 3),
            }
        )

    if not evidence:
        return None

    score = min(0.95, 0.45 + 0.18 * len(evidence))
    confidence = min(0.9, 0.35 + 0.12 * len(evidence))
    return {
        "rule": rule,
        "score": score,
        "confidence": confidence,
        "evidence": evidence,
    }


def _score_annotation_rules(hits: list[AnnotationHit]) -> list[dict]:
    function_hits = [_score_annotation_rule(rule, hits) for rule in FUNCTION_RULES]
    return [hit for hit in function_hits if hit is not None]


def _score_annotation_rule(rule: FunctionRule, hits: list[AnnotationHit]) -> dict | None:
    evidence = []
    for hit in hits:
        haystack = f"{hit.database} {hit.accession} {hit.name}".lower()
        matched = [keyword for keyword in rule.keywords if keyword in haystack]
        if not matched:
            continue
        evidence.append(
            {
                "type": "database_hit",
                "id": hit.protein_id,
                "annotation": f"{rule.name}: matched {', '.join(matched)} ({hit.database}:{hit.accession})",
                "database": hit.database,
                "weight": _annotation_weight(hit, matched),
            }
        )

    if not evidence:
        return None

    score = min(0.97, 0.5 + 0.2 * len(evidence))
    confidence = min(0.92, 0.4 + 0.12 * len(evidence))
    return {
        "rule": rule,
        "score": score,
        "confidence": confidence,
        "evidence": evidence,
    }


def _annotation_weight(hit: AnnotationHit, matched: list[str]) -> float:
    evalue_bonus = 0.2 if hit.evalue is not None and hit.evalue <= 1e-20 else 0.0
    return round(min(1.0, 0.4 + 0.12 * len(matched) + evalue_bonus), 3)


def _function_to_json(hit: dict) -> dict:
    rule: FunctionRule = hit["rule"]
    return {
        "name": rule.name,
        "score": round(hit["score"], 3),
        "confidence": round(hit["confidence"], 3),
        "evidence": hit["evidence"],
    }


def _application_scores(hits: list[dict]) -> list[dict]:
    by_panel = {panel: [] for panel in APPLICATION_AREAS}
    for hit in hits:
        rule: FunctionRule = hit["rule"]
        by_panel[rule.panel].append(hit)

    scores = []
    for panel in APPLICATION_AREAS:
        panel_hits = by_panel[panel]
        if panel_hits:
            score = min(0.95, sum(hit["score"] for hit in panel_hits) / len(panel_hits) + 0.05 * (len(panel_hits) - 1))
            confidence = min(0.9, sum(hit["confidence"] for hit in panel_hits) / len(panel_hits))
            top_functions = [hit["rule"].name for hit in sorted(panel_hits, key=lambda item: item["score"], reverse=True)[:3]]
        else:
            score = 0.05
            confidence = 0.2
            top_functions = []
        scores.append(
            {
                "area": panel,
                "score": round(score, 3),
                "confidence": round(confidence, 3),
                "top_functions": top_functions,
            }
        )
    return scores


def _biosafety(hits: list[dict]) -> dict:
    amr_hits = [hit for hit in hits if hit["rule"].biosafety_kind == "amr"]
    pathogenicity_hits = [hit for hit in hits if hit["rule"].biosafety_kind == "pathogenicity"]
    amr_risk = min(0.95, 0.1 + 0.45 * len(amr_hits))
    pathogenicity_risk = min(0.95, 0.1 + 0.45 * len(pathogenicity_hits))

    warnings = []
    if amr_hits:
        warnings.append("AMR-associated evidence detected by the baseline keyword model.")
    if pathogenicity_hits:
        warnings.append("Virulence-associated evidence detected by the baseline keyword model.")

    return {
        "pathogenicity_risk": round(pathogenicity_risk, 3),
        "amr_risk": round(amr_risk, 3),
        "warnings": warnings,
    }
