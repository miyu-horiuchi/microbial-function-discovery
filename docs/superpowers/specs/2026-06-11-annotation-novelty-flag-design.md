# Annotation-Space Novelty Flag — Design

**Date:** 2026-06-11
**Status:** Approved
**Home:** `src/microbial_function_discovery/novelty.py` (microbial-function-discovery repo)

## Motivation

Discovery candidates are ranked held-out genomes proposed for wet-lab validation. A
reviewer needs to know **how unlike the training set** each candidate is — a genome in
genuinely novel sequence/function space is a different kind of bet than one resembling
well-characterized organisms.

This flag reports **representativeness / novelty**, computed in the same annotation
feature space the model's predictions use.

### Explicit non-goal: this is NOT a confidence score

A companion finding in the microbe-foundation repo established that per-genome OOD
*distance* does **not** predict per-genome model *error* in general (it is null for most
traits and anti-correlated for imbalanced traits like pathogenicity). Therefore this flag
must be framed and documented as **"how novel / how unlike training"**, never as
"how trustworthy / how likely correct". Reports must state this caveat.

## Scope

Add one per-candidate `novelty` field to annotated discovery candidates, computed from
the existing binary annotation `FeatureMatrix`, plus a Novelty column in the two report
renderers. Self-contained to this repo; no dependency on microbe-foundation or ESM2.

## Why Jaccard (not Euclidean)

The annotation features are **binary presence/absence** vectors (`FeatureMatrix.rows`
are `list[int]` of 0/1). The right distance is **Jaccard** = 1 − |A∩B| / |A∪B|, treating
each genome as the set of feature indices it possesses. Jaccard ignores joint-absence,
which dominates sparse binary profiles and would wash out a Euclidean/Hamming distance.
This mirrors the validated ESM2 monitor's k-NN-distance structure, swapping
Euclidean → Jaccard for the binary space.

Alternatives considered and rejected: Euclidean/Hamming (joint-absence inflation),
cosine (workable but less interpretable for sets).

## Constraint: pure standard library

The repo's core `dependencies = []` (numpy is only an optional `legacy` extra). The core
modules are pure-Python. `novelty.py` MUST be pure stdlib — Jaccard over Python `set`s
needs no numpy.

## Architecture — `novelty.py`

One focused unit:

### `AnnotationNoveltyReference`

State after fit: the reference feature-sets, `k`, and a calibrated `threshold_`.

- `fit(feature_matrix: FeatureMatrix, reference_genome_ids: Iterable[str], *, k: int = 10, threshold_quantile: float = 0.95) -> AnnotationNoveltyReference`
  - Convert each reference genome's row to a `frozenset` of indices where the value is
    truthy. Store the list of reference sets.
  - Calibrate `threshold_` = the `threshold_quantile` quantile of each reference genome's
    own novelty score (mean Jaccard distance to its k nearest *other* reference genomes).
  - Raise `ValueError` if the reference is empty or smaller than k+1.
- `score(feature_vector: list[int]) -> float`
  - Convert to a feature-set; return the mean Jaccard distance to the k nearest reference
    sets. Jaccard distance of two empty sets is defined as 0.0 (identical emptiness).
- `level(score: float) -> str`
  - `"typical"` if score ≤ `threshold_`; `"highly_novel"` if score ≥ the 99th-percentile
    reference score (`extreme_threshold_`, also stored at fit); else `"novel"`.
- `ref_percentile(score: float) -> float`
  - Fraction of reference self-scores below `score` (0..1), for an interpretable readout.

Helper (module-level, testable in isolation):
- `jaccard_distance(a: frozenset, b: frozenset) -> float`

## Reference set

Reference = `genome_id`s of `LabelRecord`s with `split == "train"` (the genomes the
baseline model learned from). If no records carry `split == "train"`, fall back to all
labeled genome_ids and note it. Candidates (held-out) are scored against this reference.

## Data flow / integration into `discovery.py`

`annotate_discovery_candidates` already attaches `taxonomy`, `genome_accession`,
`biosafety`, `evidence` per candidate from `annotation_features` + `labels`. Add:

1. Build one `AnnotationNoveltyReference` from `(annotation_features, train genome_ids)`
   before the candidate loop.
2. Per candidate, attach:
   ```json
   "novelty": {"score": 0.62, "level": "novel", "ref_percentile": 0.88}
   ```
   If the candidate's `genome_id` is absent from the feature matrix, attach
   `{"score": null, "level": "unknown", "ref_percentile": null}` — never fabricate.

## Reports

- `render_discovery_candidate_report` and `render_safe_leads_report` gain a **Novelty**
  column showing `level` (and score). A footnote states: *"Novelty = how unlike the
  training set a candidate is (annotation-feature Jaccard distance); it is a
  representativeness signal, NOT a confidence or correctness estimate."*

## Error handling

- Candidate genome missing from feature matrix → `level: "unknown"`, null score.
- Empty / too-small reference → `ValueError` at fit.
- Genome with no annotation features (empty set) → scored normally (Jaccard vs reference);
  typically scores high novelty, which is correct.

## Testing

- `jaccard_distance`: known set pairs; empty/empty → 0.0; disjoint → 1.0.
- `score` monotonic: a vector equal to a reference row → ~0; a disjoint vector → ~1.
- Calibration: reference genomes mostly score ≤ `threshold_`; a synthetic disjoint
  profile is flagged `novel`/`highly_novel`.
- `fit` raises on empty / sub-k reference.
- `annotate_discovery_candidates`: attaches a well-formed `novelty` dict; missing genome
  → `unknown`; field shape stable.
- Report renderers include the Novelty column and the caveat footnote.

## Success criterion

Every annotated candidate carries a `novelty` field; reports show the Novelty column with
the representativeness caveat; all behavior is covered by tests; no new runtime
dependency is introduced (pure stdlib).
