"""Validation for microbial function discovery prediction outputs."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from microbial_function_discovery.panels import APPLICATION_AREAS


class PredictionValidationError(ValueError):
    """Raised when a prediction object violates the output contract."""


TOP_LEVEL_KEYS = {"genome_id", "application_scores", "functions", "biosafety", "novelty"}
EVIDENCE_TYPES = {"protein", "gene", "domain", "pathway", "database_hit"}
NOVELTY_RISKS = {"low", "medium", "high", "unknown"}


def validate_prediction(prediction: Mapping[str, Any]) -> None:
    """Validate a prediction object.

    The project ships a JSON Schema for external tools. This focused runtime validator
    keeps the CLI dependency-free and produces path-specific error messages.
    """

    if not isinstance(prediction, Mapping):
        raise PredictionValidationError("prediction must be an object")

    for key in sorted(TOP_LEVEL_KEYS):
        _require_key(prediction, key, key)

    _require_string(prediction["genome_id"], "genome_id")
    _validate_application_scores(prediction["application_scores"])
    _validate_functions(prediction["functions"])
    _validate_biosafety(prediction["biosafety"])
    _validate_novelty(prediction["novelty"])


def _validate_application_scores(value: Any) -> None:
    items = _require_sequence(value, "application_scores")
    for index, item in enumerate(items):
        path = f"application_scores[{index}]"
        obj = _require_mapping(item, path)
        _require_key(obj, "area", f"{path}.area")
        _require_key(obj, "score", f"{path}.score")
        _require_key(obj, "confidence", f"{path}.confidence")
        if obj["area"] not in APPLICATION_AREAS:
            raise PredictionValidationError(f"{path}.area must be one of {', '.join(APPLICATION_AREAS)}")
        _require_probability(obj["score"], f"{path}.score")
        _require_probability(obj["confidence"], f"{path}.confidence")
        if "top_functions" in obj:
            for fn_index, fn_value in enumerate(_require_sequence(obj["top_functions"], f"{path}.top_functions")):
                _require_string(fn_value, f"{path}.top_functions[{fn_index}]")


def _validate_functions(value: Any) -> None:
    items = _require_sequence(value, "functions")
    for index, item in enumerate(items):
        path = f"functions[{index}]"
        obj = _require_mapping(item, path)
        _require_key(obj, "name", f"{path}.name")
        _require_key(obj, "score", f"{path}.score")
        _require_key(obj, "evidence", f"{path}.evidence")
        _require_string(obj["name"], f"{path}.name")
        _require_probability(obj["score"], f"{path}.score")
        if "confidence" in obj:
            _require_probability(obj["confidence"], f"{path}.confidence")
        _validate_evidence(obj["evidence"], path)


def _validate_evidence(value: Any, function_path: str) -> None:
    items = _require_sequence(value, f"{function_path}.evidence")
    for index, item in enumerate(items):
        path = f"{function_path}.evidence[{index}]"
        obj = _require_mapping(item, path)
        _require_key(obj, "type", f"{path}.type")
        _require_key(obj, "id", f"{path}.id")
        if obj["type"] not in EVIDENCE_TYPES:
            raise PredictionValidationError(f"{path}.type must be one of {', '.join(sorted(EVIDENCE_TYPES))}")
        _require_string(obj["id"], f"{path}.id")
        for optional_key in ("annotation", "database", "status"):
            if optional_key in obj:
                _require_string(obj[optional_key], f"{path}.{optional_key}")
        if "weight" in obj:
            _require_probability(obj["weight"], f"{path}.weight")


def _validate_biosafety(value: Any) -> None:
    obj = _require_mapping(value, "biosafety")
    _require_key(obj, "pathogenicity_risk", "biosafety.pathogenicity_risk")
    _require_key(obj, "amr_risk", "biosafety.amr_risk")
    _require_key(obj, "warnings", "biosafety.warnings")
    _require_probability(obj["pathogenicity_risk"], "biosafety.pathogenicity_risk")
    _require_probability(obj["amr_risk"], "biosafety.amr_risk")
    for index, warning in enumerate(_require_sequence(obj["warnings"], "biosafety.warnings")):
        _require_string(warning, f"biosafety.warnings[{index}]")


def _validate_novelty(value: Any) -> None:
    obj = _require_mapping(value, "novelty")
    _require_key(obj, "generalization_risk", "novelty.generalization_risk")
    if obj["generalization_risk"] not in NOVELTY_RISKS:
        raise PredictionValidationError(
            f"novelty.generalization_risk must be one of {', '.join(sorted(NOVELTY_RISKS))}"
        )
    for optional_key in ("nearest_training_family", "nearest_training_genus"):
        if optional_key in obj:
            _require_string(obj[optional_key], f"novelty.{optional_key}")
    if "notes" in obj:
        for index, note in enumerate(_require_sequence(obj["notes"], "novelty.notes")):
            _require_string(note, f"novelty.notes[{index}]")


def _require_key(obj: Mapping[str, Any], key: str, path: str) -> None:
    if key not in obj:
        raise PredictionValidationError(f"{path} is required")


def _require_mapping(value: Any, path: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise PredictionValidationError(f"{path} must be an object")
    return value


def _require_sequence(value: Any, path: str) -> Sequence[Any]:
    if isinstance(value, str) or not isinstance(value, Sequence):
        raise PredictionValidationError(f"{path} must be an array")
    return value


def _require_string(value: Any, path: str) -> str:
    if not isinstance(value, str) or not value:
        raise PredictionValidationError(f"{path} must be a non-empty string")
    return value


def _require_probability(value: Any, path: str) -> float:
    if not isinstance(value, (int, float)) or isinstance(value, bool) or not 0 <= value <= 1:
        raise PredictionValidationError(f"{path} must be a number from 0 to 1")
    return float(value)
