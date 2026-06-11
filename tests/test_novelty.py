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
