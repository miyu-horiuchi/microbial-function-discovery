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
