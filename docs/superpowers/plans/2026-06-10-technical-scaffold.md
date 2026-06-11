# Technical Scaffold Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn the repository from strategy documents into a small working Python project that can define application panels and validate prediction outputs.

**Architecture:** Add a standard-library Python package under `src/` with focused modules for panel definitions, prediction validation, and CLI commands. Keep the model/training code out of this first slice; the scaffold establishes stable contracts that later model code can emit.

**Tech Stack:** Python 3.11+, standard library, `unittest`, JSON files.

---

## File Structure

- Create `pyproject.toml`: package metadata and console script.
- Create `src/microbial_function_discovery/__init__.py`: package exports.
- Create `src/microbial_function_discovery/panels.py`: application panel registry and allowed application areas.
- Create `src/microbial_function_discovery/validation.py`: validation for prediction JSON objects.
- Create `src/microbial_function_discovery/cli.py`: `panels` and `validate` commands.
- Create `examples/prediction.example.json`: valid prediction output example.
- Create `tests/test_panels.py`: panel registry behavior.
- Create `tests/test_validation.py`: valid and invalid prediction output behavior.
- Create `tests/test_cli.py`: CLI command smoke tests.

## Task 1: Panel Registry

**Files:**
- Create: `src/microbial_function_discovery/panels.py`
- Create: `tests/test_panels.py`

- [ ] **Step 1: Write failing tests**

```python
from microbial_function_discovery.panels import APPLICATION_AREAS, get_panel, list_panels


def test_all_schema_application_areas_have_panels():
    assert set(APPLICATION_AREAS) == {
        "environmental_terraforming",
        "therapeutics_antimicrobials",
        "biofuels_industrial",
        "food_fermentation_agriculture",
        "biosafety",
    }
    assert {panel.key for panel in list_panels()} == set(APPLICATION_AREAS)


def test_get_panel_returns_targets_and_evidence_sources():
    panel = get_panel("biofuels_industrial")
    assert "cellulose degradation" in panel.example_targets
    assert "CAZy" in panel.evidence_sources
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m unittest tests.test_panels -v`
Expected: FAIL because `microbial_function_discovery.panels` does not exist.

- [ ] **Step 3: Implement panel registry**

Create a frozen `Panel` dataclass, `APPLICATION_AREAS`, `PANELS`, `list_panels()`, and `get_panel()`.

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH=src python3 -m unittest tests.test_panels -v`
Expected: PASS.

## Task 2: Prediction Validation

**Files:**
- Create: `src/microbial_function_discovery/validation.py`
- Create: `examples/prediction.example.json`
- Create: `tests/test_validation.py`

- [ ] **Step 1: Write failing tests**

```python
import json
import unittest
from pathlib import Path

from microbial_function_discovery.validation import PredictionValidationError, validate_prediction


class PredictionValidationTests(unittest.TestCase):
    def test_example_prediction_is_valid(self):
        prediction = json.loads(Path("examples/prediction.example.json").read_text())
        validate_prediction(prediction)

    def test_invalid_application_area_is_rejected(self):
        prediction = json.loads(Path("examples/prediction.example.json").read_text())
        prediction["application_scores"][0]["area"] = "unknown"
        with self.assertRaisesRegex(PredictionValidationError, "application_scores\\[0\\].area"):
            validate_prediction(prediction)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m unittest tests.test_validation -v`
Expected: FAIL because `validation.py` and the example JSON do not exist.

- [ ] **Step 3: Implement validator and example**

Implement a focused validator for required top-level keys, application areas, numeric score ranges, biosafety ranges, function evidence, and novelty risk values.

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH=src python3 -m unittest tests.test_validation -v`
Expected: PASS.

## Task 3: CLI

**Files:**
- Create: `src/microbial_function_discovery/cli.py`
- Create: `tests/test_cli.py`
- Modify: `pyproject.toml`

- [ ] **Step 1: Write failing tests**

```python
import subprocess
import sys
import unittest


class CliTests(unittest.TestCase):
    def test_panels_command_lists_application_areas(self):
        result = subprocess.run(
            [sys.executable, "-m", "microbial_function_discovery.cli", "panels"],
            check=True,
            capture_output=True,
            text=True,
        )
        assert "environmental_terraforming" in result.stdout

    def test_validate_command_accepts_example(self):
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "microbial_function_discovery.cli",
                "validate",
                "examples/prediction.example.json",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        assert "valid" in result.stdout
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=src python3 -m unittest tests.test_cli -v`
Expected: FAIL because CLI does not exist.

- [ ] **Step 3: Implement CLI and package metadata**

Add `argparse` commands:

- `panels`: print panel keys and display names.
- `validate <path>`: load JSON and run `validate_prediction`.

- [ ] **Step 4: Run all tests**

Run: `PYTHONPATH=src python3 -m unittest -v`
Expected: PASS.

- [ ] **Step 5: Commit and push**

```bash
git add .
git commit -m "Add initial Python scaffold"
git push
```
