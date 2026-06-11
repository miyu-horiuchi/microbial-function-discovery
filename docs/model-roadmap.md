# Model Roadmap

## North Star

Build a microbial genome foundation model that predicts useful biological
functions from sequenced but uncultured microbes, ranks candidate organisms for
applications, and explains predictions with gene/pathway evidence.

## Phase 1: Benchmark and Schema

Deliverables:

- application panel definitions
- prediction output schema
- family-held-out benchmark protocol
- baseline list
- evidence validation rules

No new training is required in this phase. The goal is to make the model target
precise.

## Phase 2: Multi-Panel Dataset

Deliverables:

- genome-indexed label matrix
- application panel metadata
- per-label provenance
- family-held-out splits
- baseline feature matrices
- evidence tables for function databases

Candidate label/evidence sources:

- BacDive and MediaDive for traits and cultivation links
- eggNOG, KEGG, COG, Pfam, TIGRFAM for functional annotations
- CAZy for carbohydrate-active enzymes
- VFDB for virulence factors
- CARD / AMRFinderPlus for AMR
- antiSMASH / MIBiG for biosynthetic gene clusters
- GTDB / MGnify / IMG/M for genome metadata and uncultured diversity

## Phase 3: Shared Encoder

Recommended first architecture:

- protein set encoder over per-protein protein-language-model embeddings
- domain/pathway token stream for interpretable biological evidence
- set attention or Perceiver-style latent bottleneck for gene combinations
- application-specific prediction heads
- calibrated novelty and uncertainty output

The first model should compare against simple annotation and mean-pool baselines.
The goal is not to beat every method on every panel immediately; the goal is to
show credible discovery value on at least one high-value panel while preserving
biological evidence.

## Phase 4: Discovery Engine

Deliverables:

- upload one genome or batch candidate set
- choose target application
- rank candidate organisms
- show predicted functions and confidence
- show evidence genes/pathways
- show biosafety and novelty caveats
- export JSON/CSV reports

## Phase 5: Wet-Lab Feedback Loop

Deliverables:

- candidate selection for experiments
- result ingestion
- active-learning acquisition function
- model retraining with new validation outcomes

The long-term flywheel is: predict useful organisms, test them, learn from the
results, and improve discovery in regions of microbial space where labels are
sparse.

