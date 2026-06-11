# Annotation-Space Novelty Flag Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Attach a per-candidate representativeness/novelty flag (annotation-feature Jaccard distance from the training genomes) to discovery candidates and show it in both report renderers, framed explicitly as novelty — not confidence.

**Architecture:** A pure-stdlib `AnnotationNoveltyReference` (in a new `novelty.py`) scores a candidate's binary annotation-feature set by mean Jaccard distance to its k nearest training-genome sets, with a reference-calibrated threshold for a typical/novel/highly_novel level. `annotate_discovery_candidates` builds the reference from `split=="train"` labels and attaches a `novelty` dict per candidate; the two report renderers add a Novelty column plus a caveat footnote.

**Tech Stack:** Python standard library only (the repo's core `dependencies = []`; numpy is an optional `legacy` extra). Tests are `unittest`-style, run with `python -m pytest`.

**Spec:** `docs/superpowers/specs/2026-06-11-annotation-novelty-flag-design.md`

---

### Task 1: novelty.py — Jaccard novelty reference

**Files:**
- Create: `src/microbial_function_discovery/novelty.py`
- Test: `tests/test_novelty.py`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_novelty.py`:

```python
import unittest

from microbial_function_discovery.features import FeatureMatrix
from microbial_function_discovery.novelty import AnnotationNoveltyReference, jaccard_distance


class JaccardTests(unittest.TestCase):
    def test_identical_sets_distance_zero(self):
        self.assertEqual(jaccard_distance(frozenset({0, 1}), frozenset({0, 1})), 0.0)

    def test_disjoint_sets_distance_one(self):
        self.assertEqual(jaccard_distance(frozenset({0}), frozenset({1})), 1.0)

    def test_both_empty_distance_zero(self):
        self.assertEqual(jaccard_distance(frozenset(), frozenset()), 0.0)

    def test_partial_overlap(self):
        # |{0,1} ∩ {1,2}| / |{0,1,2}| = 1/3 -> distance 2/3
        self.assertAlmostEqual(jaccard_distance(frozenset({0, 1}), frozenset({1, 2})), 2.0 / 3.0)


def _matrix():
    # 3 features; three reference genomes + room to score new vectors
    return FeatureMatrix(
        genome_ids=["R1", "R2", "R3"],
        feature_names=["f0", "f1", "f2"],
        rows=[[1, 1, 0], [1, 1, 0], [1, 0, 0]],
    )


class NoveltyReferenceTests(unittest.TestCase):
    def test_fit_requires_two_reference_genomes(self):
        m = FeatureMatrix(genome_ids=["R1"], feature_names=["f0"], rows=[[1]])
        with self.assertRaises(ValueError):
            AnnotationNoveltyReference().fit(m, ["R1"])

    def test_vector_matching_reference_scores_low(self):
        ref = AnnotationNoveltyReference(k=2).fit(_matrix(), ["R1", "R2", "R3"])
        # identical to R1/R2 -> nearest distance 0
        self.assertEqual(ref.score([1, 1, 0]), 0.0)

    def test_disjoint_vector_scores_high_and_flags_novel(self):
        ref = AnnotationNoveltyReference(k=2).fit(_matrix(), ["R1", "R2", "R3"])
        score = ref.score([0, 0, 1])  # shares nothing with any reference
        self.assertGreater(score, ref.threshold_)
        self.assertIn(ref.level(score), {"novel", "highly_novel"})

    def test_level_typical_below_threshold(self):
        ref = AnnotationNoveltyReference(k=2).fit(_matrix(), ["R1", "R2", "R3"])
        self.assertEqual(ref.level(ref.score([1, 1, 0])), "typical")

    def test_ref_percentile_in_unit_interval(self):
        ref = AnnotationNoveltyReference(k=2).fit(_matrix(), ["R1", "R2", "R3"])
        p = ref.ref_percentile(ref.score([0, 0, 1]))
        self.assertGreaterEqual(p, 0.0)
        self.assertLessEqual(p, 1.0)

    def test_ignores_ids_absent_from_matrix(self):
        ref = AnnotationNoveltyReference(k=2).fit(_matrix(), ["R1", "R2", "R3", "GHOST"])
        self.assertEqual(len(ref._ref_sets), 3)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_novelty.py -q`
Expected: FAIL with `ModuleNotFoundError: microbial_function_discovery.novelty`

- [ ] **Step 3: Implement novelty.py**

Create `src/microbial_function_discovery/novelty.py`:

```python
"""Annotation-space novelty (representativeness) scoring.

Novelty = how unlike the training genomes a candidate is, measured in the binary
annotation feature space the model's predictions use (Jaccard distance over the set
of features each genome possesses).

This is NOT a confidence/correctness signal. A microbe-foundation finding established
that out-of-distribution distance does not predict per-genome model error in general,
so this flag reports representativeness only.

Pure standard library: the repo's core dependencies are empty.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

from microbial_function_discovery.features import FeatureMatrix


def jaccard_distance(a: frozenset, b: frozenset) -> float:
    """1 - |a ∩ b| / |a ∪ b|. Two empty sets are identical -> distance 0.0."""
    if not a and not b:
        return 0.0
    union = len(a | b)
    return 1.0 - len(a & b) / union


def _row_to_set(row: list[int]) -> frozenset:
    return frozenset(index for index, value in enumerate(row) if value)


def _quantile(sorted_values: list[float], q: float) -> float:
    if not sorted_values:
        return 0.0
    pos = q * (len(sorted_values) - 1)
    lo = int(pos)
    hi = min(lo + 1, len(sorted_values) - 1)
    frac = pos - lo
    return sorted_values[lo] * (1.0 - frac) + sorted_values[hi] * frac


@dataclass
class AnnotationNoveltyReference:
    """Mean k-NN Jaccard distance to the training-genome feature sets."""

    k: int = 10
    threshold_quantile: float = 0.95
    extreme_quantile: float = 0.99
    _ref_sets: list = field(default_factory=list, repr=False)
    k_: int = 0
    threshold_: float = 0.0
    extreme_threshold_: float = 0.0
    _ref_scores: list = field(default_factory=list, repr=False)

    def fit(self, feature_matrix: FeatureMatrix, reference_genome_ids: Iterable[str]) -> "AnnotationNoveltyReference":
        present = set(feature_matrix.genome_ids)
        ref_ids = [g for g in reference_genome_ids if g in present]
        sets = [_row_to_set(feature_matrix.row_for(g)) for g in ref_ids]
        if len(sets) < 2:
            raise ValueError(f"need at least 2 reference genomes with features; got {len(sets)}")
        self._ref_sets = sets
        # Adapt k down to the available neighbour count (cf. EuclideanBackend min(k, n)).
        self.k_ = max(1, min(self.k, len(sets) - 1))

        scores = []
        for i, s in enumerate(sets):
            dists = sorted(jaccard_distance(s, t) for j, t in enumerate(sets) if j != i)
            scores.append(sum(dists[: self.k_]) / self.k_)
        scores.sort()
        self._ref_scores = scores
        self.threshold_ = _quantile(scores, self.threshold_quantile)
        self.extreme_threshold_ = _quantile(scores, self.extreme_quantile)
        return self

    def score(self, feature_vector: list[int]) -> float:
        s = _row_to_set(feature_vector)
        dists = sorted(jaccard_distance(s, t) for t in self._ref_sets)
        kk = min(self.k_, len(dists))
        return sum(dists[:kk]) / kk

    def level(self, score: float) -> str:
        if score >= self.extreme_threshold_ and self.extreme_threshold_ > self.threshold_:
            return "highly_novel"
        if score > self.threshold_:
            return "novel"
        return "typical"

    def ref_percentile(self, score: float) -> float:
        below = sum(1 for x in self._ref_scores if x < score)
        return below / len(self._ref_scores)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/test_novelty.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/microbial_function_discovery/novelty.py tests/test_novelty.py
git commit -m "feat(novelty): annotation-space Jaccard novelty reference"
```

---

### Task 2: Attach novelty to annotated candidates

**Files:**
- Modify: `src/microbial_function_discovery/discovery.py`
- Test: `tests/test_discovery.py`

- [ ] **Step 1: Add the failing test**

Append this test method inside the existing `DiscoveryAnnotationTests` class in `tests/test_discovery.py` (it reuses the module-level `ANNOTATIONS_TSV` / `LABELS_TSV` fixtures already defined there):

```python
    def test_annotate_attaches_novelty_flag(self):
        features = build_feature_matrix_from_annotation_tsv(ANNOTATIONS_TSV)
        labels = parse_labels_tsv(LABELS_TSV)
        model = train_baseline(features, labels)
        candidates = {
            "split": "test",
            "targets": [
                {
                    "target_key": "biofuels_industrial:cellulose_degradation",
                    "panel": "biofuels_industrial",
                    "label": "cellulose_degradation",
                    "source": "annotation",
                    "metrics": {"precision_at_10": 1.0},
                    "candidates": [{"genome_id": "G3", "rank": 1, "score": 0.9}],
                }
            ],
        }
        annotated = annotate_discovery_candidates(candidates, [], [], model, features, labels)
        cand = annotated["targets"][0]["candidates"][0]
        self.assertIn("novelty", cand)
        self.assertIn(cand["novelty"]["level"], {"typical", "novel", "highly_novel", "unknown"})
        # G3's annotation set equals training genome G1's -> nearest distance 0 -> typical
        self.assertEqual(cand["novelty"]["level"], "typical")

    def test_annotate_novelty_unknown_for_missing_genome(self):
        features = build_feature_matrix_from_annotation_tsv(ANNOTATIONS_TSV)
        labels = parse_labels_tsv(LABELS_TSV)
        model = train_baseline(features, labels)
        candidates = {
            "split": "test",
            "targets": [
                {
                    "target_key": "biofuels_industrial:cellulose_degradation",
                    "panel": "biofuels_industrial",
                    "label": "cellulose_degradation",
                    "source": "annotation",
                    "metrics": {"precision_at_10": 1.0},
                    "candidates": [{"genome_id": "NOT_IN_MATRIX", "rank": 1, "score": 0.9}],
                }
            ],
        }
        annotated = annotate_discovery_candidates(candidates, [], [], model, features, labels)
        cand = annotated["targets"][0]["candidates"][0]
        self.assertEqual(cand["novelty"]["level"], "unknown")
        self.assertIsNone(cand["novelty"]["score"])
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_discovery.py -q -k novelty`
Expected: FAIL (`KeyError: 'novelty'`)

- [ ] **Step 3: Add the import to discovery.py**

In `src/microbial_function_discovery/discovery.py`, add to the import block near the top (after the existing `from microbial_function_discovery.learning import BaselineModel` line):

```python
from microbial_function_discovery.novelty import AnnotationNoveltyReference
```

- [ ] **Step 4: Add the `_novelty` helper to discovery.py**

Add this function near the other private helpers (e.g. just above `def _lead_row`):

```python
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
```

- [ ] **Step 5: Wire it into `annotate_discovery_candidates`**

In `annotate_discovery_candidates`, immediately after the `accessions_by_id = {...}` dict comprehension and before `annotated_targets = []`, add:

```python
    novelty_ref = _build_novelty_reference(labels, annotation_features)
    present_ids = set(annotation_features.genome_ids)
```

Then inside the candidate loop, add the `"novelty"` key to the appended candidate dict (alongside `"biosafety"` and `"evidence"`):

```python
                    "novelty": _novelty(genome_id, novelty_ref, present_ids, annotation_features),
```

- [ ] **Step 6: Run tests to verify they pass**

Run: `python -m pytest tests/test_discovery.py -q -k novelty`
Expected: PASS (2 passed)

- [ ] **Step 7: Commit**

```bash
git add src/microbial_function_discovery/discovery.py tests/test_discovery.py
git commit -m "feat(novelty): attach novelty flag to annotated candidates"
```

---

### Task 3: Surface novelty in leads + both reports

**Files:**
- Modify: `src/microbial_function_discovery/discovery.py`
- Test: `tests/test_discovery.py`

- [ ] **Step 1: Add the failing test**

Append inside `DiscoveryAnnotationTests` in `tests/test_discovery.py`:

```python
    def _annotated_for_reports(self):
        features = build_feature_matrix_from_annotation_tsv(ANNOTATIONS_TSV)
        labels = parse_labels_tsv(LABELS_TSV)
        model = train_baseline(features, labels)
        candidates = {
            "split": "test",
            "targets": [
                {
                    "target_key": "biofuels_industrial:cellulose_degradation",
                    "panel": "biofuels_industrial",
                    "label": "cellulose_degradation",
                    "source": "annotation",
                    "metrics": {"precision_at_10": 1.0},
                    "candidates": [{"genome_id": "G3", "rank": 1, "score": 0.9}],
                }
            ],
        }
        return annotate_discovery_candidates(candidates, [], [], model, features, labels)

    def test_candidate_report_has_novelty_column(self):
        from microbial_function_discovery.discovery import render_discovery_candidate_report
        report = render_discovery_candidate_report(self._annotated_for_reports())
        self.assertIn("Novelty", report)
        self.assertIn("representativeness", report.lower())

    def test_safe_leads_carry_novelty_level(self):
        from microbial_function_discovery.discovery import (
            export_safe_leads,
            render_safe_leads_report,
            SAFE_LEAD_FIELDS,
        )
        export = export_safe_leads(self._annotated_for_reports(), require_accession=False, require_evidence=False)
        self.assertIn("novelty_level", SAFE_LEAD_FIELDS)
        self.assertTrue(export["leads"], "expected at least one lead")
        self.assertIn("novelty_level", export["leads"][0])
        report = render_safe_leads_report(export)
        self.assertIn("Novelty", report)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_discovery.py -q -k "novelty_column or novelty_level"`
Expected: FAIL (`Novelty` not in report / `novelty_level` not in `SAFE_LEAD_FIELDS`)

- [ ] **Step 3: Add `novelty_level` to `SAFE_LEAD_FIELDS`**

In `discovery.py`, append `"novelty_level"` to the `SAFE_LEAD_FIELDS` list (after `"evidence"`):

```python
    "evidence",
    "novelty_level",
```

- [ ] **Step 4: Carry novelty into `_lead_row`**

In `_lead_row`, add this key to the returned dict (after the `"evidence"` entry):

```python
        "novelty_level": str(candidate.get("novelty", {}).get("level", "unknown")),
```

- [ ] **Step 5: Add the Novelty column to `render_discovery_candidate_report`**

Change the table header and separator lines (currently ending `... | Safety | Evidence |`) to include Novelty before Evidence:

```python
        "| Target | Rank | Genome | Species | Source | Score | Safety | Novelty | Evidence |",
        "|---|---:|---|---|---|---:|---|---|---|",
```

In the per-candidate row list (the `" | ".join([...])` block), insert the novelty cell right before the evidence cell (`_markdown_cell(evidence_text)`):

```python
                        str(candidate.get("biosafety", {}).get("risk_level", "unknown")),
                        str(candidate.get("novelty", {}).get("level", "unknown")),
                        _markdown_cell(evidence_text),
```

Add a caveat bullet to the `## Notes` list (after the Evidence note):

```python
            "- Novelty = how unlike the training set a candidate is (annotation-feature Jaccard distance); it is a representativeness signal, NOT a confidence or correctness estimate.",
```

- [ ] **Step 6: Add the Novelty column to `render_safe_leads_report`**

Change the Leads table header and separator (currently `... | Risk | Evidence |`) to:

```python
        "| Panel | Target | Genome | Species | Accession | Score | Risk | Novelty | Evidence |",
        "|---|---|---|---|---|---:|---|---|---|",
```

In the per-lead `" | ".join([...])` block, insert the novelty cell right before the evidence cell (`_markdown_cell(lead.get("evidence"))`):

```python
                    _markdown_cell(lead.get("risk_level")),
                    _markdown_cell(lead.get("novelty_level")),
                    _markdown_cell(lead.get("evidence")),
```

- [ ] **Step 7: Run tests to verify they pass**

Run: `python -m pytest tests/test_discovery.py -q -k "novelty"`
Expected: PASS

- [ ] **Step 8: Commit**

```bash
git add src/microbial_function_discovery/discovery.py tests/test_discovery.py
git commit -m "feat(novelty): show novelty in leads export and both reports"
```

---

### Task 4: Full-suite verification

**Files:**
- Test: (whole repo)

- [ ] **Step 1: Run the full suite**

Run: `python -m pytest -q`
Expected: PASS — all pre-existing tests plus the new novelty tests green. (If `pytest` is unavailable, use `python -m unittest discover -s tests -q`.)

- [ ] **Step 2: Smoke-check the rendered output by eye (optional)**

Run:
```bash
python - <<'PY'
from microbial_function_discovery.datasets import parse_labels_tsv
from microbial_function_discovery.features import build_feature_matrix_from_annotation_tsv
from microbial_function_discovery.learning import train_baseline
from microbial_function_discovery.discovery import annotate_discovery_candidates, render_discovery_candidate_report
ann = "genome_id\tprotein_id\tdatabase\taccession\tname\tevalue\nG1\tp1\tCAZy\tGH5\tx\t1e-20\nG2\tp1\tPfam\tPF00002\ty\t1e-20\nG3\tp1\tCAZy\tGH5\tz\t1e-20\n"
lab = "genome_id\tsplit\tfamily\tpanel\tlabel\tvalue\nG1\ttrain\tA\tbiofuels_industrial\tcellulose_degradation\t1\nG2\ttrain\tA\tbiofuels_industrial\tcellulose_degradation\t0\nG3\ttest\tB\tbiofuels_industrial\tcellulose_degradation\t1\n"
f = build_feature_matrix_from_annotation_tsv(ann); L = parse_labels_tsv(lab); m = train_baseline(f, L)
rep = {"split":"test","targets":[{"target_key":"biofuels_industrial:cellulose_degradation","panel":"biofuels_industrial","label":"cellulose_degradation","source":"annotation","metrics":{"precision_at_10":1.0},"candidates":[{"genome_id":"G3","rank":1,"score":0.9}]}]}
print(render_discovery_candidate_report(annotate_discovery_candidates(rep, [], [], m, f, L)))
PY
```
Expected: a Markdown table containing a `Novelty` column and the representativeness caveat note.

- [ ] **Step 3: Commit (only if any incidental fixes were needed)**

```bash
git add -A && git commit -m "test(novelty): verify full suite green" || echo "nothing to commit"
```

---

## Self-review notes

- **Spec coverage:** Jaccard metric + pure-stdlib (Task 1), `AnnotationNoveltyReference.fit/score/level/ref_percentile` (Task 1), reference = train labels with fallback (Task 2 `_build_novelty_reference`), per-candidate `novelty` dict with unknown-for-missing (Task 2), Novelty column + caveat footnote in both reports (Task 3), leads export field (Task 3), tests for all (Tasks 1–3), full-suite check (Task 4). The spec's "raise on too-small reference" is intentionally relaxed to "raise only below 2 + adapt k" so it works on real and fixture-sized references — documented in `fit`.
- **Type consistency:** `AnnotationNoveltyReference(k, threshold_quantile, extreme_quantile)` with `fit/score(list[int])/level(float)/ref_percentile(float)`; `_novelty` returns `{score,level,ref_percentile}`; `_lead_row` adds `novelty_level`; `SAFE_LEAD_FIELDS` includes `novelty_level`. Names match across tasks.
- **No placeholders:** every code step is complete.
