# Microbial Function Discovery

**Search microbial dark matter for useful biology.**

Microbial Function Discovery is a foundation-model project for predicting useful
functions from microbial genomes, including metagenome-assembled genomes (MAGs),
single-amplified genomes (SAGs), isolates, and protein FASTA files. The goal is
to identify uncultured microbes worth validating for environmental, therapeutic,
industrial, food, agriculture, and bioenergy applications.

## Vision

Most microbial life has never been cultured, but its genomes are already being
sequenced. This project turns those genomes into ranked, evidence-backed
discovery candidates:

- What can this microbe do?
- Which applications might it be useful for?
- Which genes, proteins, domains, or pathways support the prediction?
- Is it safe enough to pursue?
- How trustworthy is the prediction for a taxonomically novel organism?

## Product Contract

Input:

- assembled genome FASTA
- MAG or SAG FASTA
- predicted protein FASTA
- batch candidate sets from metagenomic studies

Output:

- predicted functions
- application-specific usefulness scores
- biosafety and pathogenicity risk
- confidence and taxonomic novelty warnings
- evidence genes, proteins, domains, and pathways
- ranked candidates for wet-lab validation

## Application Panels

| Panel | Example targets |
|---|---|
| Environmental / terraforming | carbon fixation, nitrogen cycling, sulfur metabolism, metal reduction, plastic degradation, desiccation/radiation/salinity tolerance |
| Therapeutics / antimicrobials | biosynthetic gene clusters, antimicrobial peptides, pathogen suppression, microbiome-relevant functions, toxin and virulence risk |
| Biofuels / industrial enzymes | cellulose/xylan/lignin degradation, lipid accumulation, fermentation pathways, thermostable enzymes, acid/salt-stable enzymes |
| Food / fermentation / agriculture | fermentation traits, flavor/aroma metabolism, probiotic-relevant functions, plant growth promotion, nitrogen fixation, spoilage and safety flags |
| Biosafety | pathogenicity, virulence factors, AMR phenotype, toxin genes, novelty caveats |

## Model Direction

The project uses one shared microbial genome encoder with application-specific
prediction heads:

```text
genome / MAG / SAG
  -> protein calling
  -> protein language model embeddings
  -> domain / pathway / gene-family evidence tokens
  -> shared genome encoder
  -> application heads
  -> discovery ranking layer
```

The model should predict useful functions while exposing the evidence behind
each prediction. The first technical challenge is robust generalization to
taxonomically novel organisms, not just high accuracy on close relatives.

## Why This Is Different

Existing genome annotation tools often report what known genes are present.
This project aims to rank organisms by useful application potential while
combining learned representations, pathway evidence, biosafety filtering,
uncertainty, and cross-clade evaluation.

The benchmark is designed around the hard case: organisms unlike those seen in
training.

## Repository Status

This repository is the new standalone surface for the broader microbial
function-discovery model. It starts with:

- product and model definition
- benchmark protocol
- output schema
- roadmap

The earlier `microbe-foundation` research repository contains the first wedge:
BacDive trait prediction, per-protein attention pooling, pathogenicity
attribution, and the predictability-gradient paper.

## Initial Documents

- [Benchmark protocol](docs/benchmark.md)
- [Discovery candidate report](docs/discovery-candidate-report.md)
- [Legacy baseline comparison](docs/legacy-baseline-comparison.md)
- [Model roadmap](docs/model-roadmap.md)
- [Product leads report](docs/product-leads-report.md)
- [Prediction output schema](schemas/prediction.schema.json)
- [Safe leads report](docs/safe-leads-report.md)
- [Target leaderboard](docs/target-leaderboard.md)

## Developer Quickstart

Run the scaffold without installing dependencies:

```bash
PYTHONPATH=src python3 -m microbial_function_discovery.cli panels
PYTHONPATH=src python3 -m microbial_function_discovery.cli predict examples/useful_functions.faa --genome-id candidate_001
PYTHONPATH=src python3 -m microbial_function_discovery.cli import-eggnog examples/eggnog_mapper.emapper.annotations
PYTHONPATH=src python3 -m microbial_function_discovery.cli import-domtblout examples/pfam.domtblout --database Pfam
PYTHONPATH=src python3 -m microbial_function_discovery.cli run-eggnog examples/useful_functions.faa --output-dir outputs/eggnog
PYTHONPATH=src python3 -m microbial_function_discovery.cli run-domtblout examples/useful_functions.faa --hmm /path/to/Pfam-A.hmm --out outputs/pfam.domtblout --database Pfam
PYTHONPATH=src python3 -m microbial_function_discovery.cli predict-annotations examples/annotation_hits.tsv --genome-id candidate_001
PYTHONPATH=src python3 -m microbial_function_discovery.cli validate-splits examples/benchmark_labels.tsv
PYTHONPATH=src python3 -m microbial_function_discovery.cli build-features examples/multi_genome_annotation_hits.tsv --out outputs/features.json
PYTHONPATH=src python3 -m microbial_function_discovery.cli train-baseline outputs/features.json examples/benchmark_labels.tsv --out outputs/model.json
PYTHONPATH=src python3 -m microbial_function_discovery.cli evaluate outputs/model.json outputs/features.json examples/benchmark_labels.tsv --split test
PYTHONPATH=src python3 -m microbial_function_discovery.cli evaluate-ranking outputs/model.json outputs/features.json examples/benchmark_labels.tsv --split test --target biofuels_industrial:cellulose_degradation --k 1 --k 10
PYTHONPATH=src python3 -m microbial_function_discovery.cli predict-baseline outputs/model.json outputs/features.json --genome-id G3
PYTHONPATH=src python3 -m microbial_function_discovery.cli rank-candidates outputs/model.json outputs/features.json --target biofuels_industrial:cellulose_degradation --limit 10
PYTHONPATH=src python3 -m microbial_function_discovery.cli validate examples/prediction.example.json
PYTHONPATH=src python3 -m unittest discover -s tests -v
```

Install the local CLI:

```bash
python3 -m pip install -e .
mfd panels
mfd predict examples/useful_functions.faa --genome-id candidate_001
mfd import-eggnog examples/eggnog_mapper.emapper.annotations
mfd import-domtblout examples/pfam.domtblout --database Pfam
mfd run-eggnog examples/useful_functions.faa --output-dir outputs/eggnog
mfd run-domtblout examples/useful_functions.faa --hmm /path/to/Pfam-A.hmm --out outputs/pfam.domtblout --database Pfam
mfd predict-annotations examples/annotation_hits.tsv --genome-id candidate_001
mfd validate-splits examples/benchmark_labels.tsv
mfd build-features examples/multi_genome_annotation_hits.tsv --out outputs/features.json
mfd train-baseline outputs/features.json examples/benchmark_labels.tsv --out outputs/model.json
mfd evaluate outputs/model.json outputs/features.json examples/benchmark_labels.tsv --split test
mfd evaluate-ranking outputs/model.json outputs/features.json examples/benchmark_labels.tsv --split test --target biofuels_industrial:cellulose_degradation --k 1 --k 10
mfd predict-baseline outputs/model.json outputs/features.json --genome-id G3
mfd rank-candidates outputs/model.json outputs/features.json --target biofuels_industrial:cellulose_degradation --limit 10
mfd validate examples/prediction.example.json
```

Reuse the earlier `microbe-foundation` local artifacts with the optional legacy
dependencies:

```bash
python3 -m pip install -e ".[legacy]"
mfd import-legacy-labels \
  /Users/miyuhoriuchi/microbe-foundation/data/traits.parquet \
  /Users/miyuhoriuchi/microbe-foundation/data/splits.parquet \
  --out data/legacy/labels.tsv \
  --split-column family_split \
  --max-multilabel-classes 25
mfd validate-splits data/legacy/labels.tsv
mfd import-legacy-eggnog-features \
  /Users/miyuhoriuchi/microbe-foundation/data/eggnog_features_5851.npz \
  /Users/miyuhoriuchi/microbe-foundation/data/eggnog_vocab.json \
  --labels data/legacy/labels.tsv \
  --out data/legacy/features.1000.json \
  --max-features 1000 \
  --min-prevalence 5
mfd train-baseline data/legacy/features.1000.json data/legacy/labels.tsv --out data/legacy/model.1000.json
mfd evaluate data/legacy/model.1000.json data/legacy/features.1000.json data/legacy/labels.tsv --split test --out data/legacy/eval.test.1000.json
mfd evaluate-ranking data/legacy/model.1000.json data/legacy/features.1000.json data/legacy/labels.tsv --split test --target biofuels_industrial:temperature_class__thermophile --k 10 --k 50
mfd train-dense-baseline \
  /Users/miyuhoriuchi/microbe-foundation/data/esm2_features.npz \
  data/legacy/labels.tsv \
  --out data/legacy/dense_model.esm2.json \
  --max-features 640
mfd evaluate-dense \
  data/legacy/dense_model.esm2.json \
  /Users/miyuhoriuchi/microbe-foundation/data/esm2_features.npz \
  data/legacy/labels.tsv \
  --split test \
  --out data/legacy/dense_eval.test.esm2.json
mfd evaluate-dense-ranking \
  data/legacy/dense_model.esm2.json \
  /Users/miyuhoriuchi/microbe-foundation/data/esm2_features.npz \
  data/legacy/labels.tsv \
  --split test \
  --target biofuels_industrial:temperature_class__thermophile \
  --k 10 --k 50
```

The current `predict` command is a transparent annotation-keyword baseline. It
is useful for exercising the product contract now and will be replaced by
learned genome encoders as benchmark datasets come online.

For real annotation pipelines, use `predict-annotations` with a TSV containing:

```text
protein_id	database	accession	name	evalue
p1	CAZy	GH5	glycoside hydrolase family 5 cellulase	1e-40
p2	KEGG	K02588	nitrogenase iron protein nifH	1e-50
p3	CARD	blaTEM	beta-lactamase	1e-20
```

You can also convert common annotation tool outputs into this normalized TSV:

```bash
mfd import-eggnog examples/eggnog_mapper.emapper.annotations > annotation_hits.tsv
mfd import-domtblout examples/pfam.domtblout --database Pfam >> annotation_hits.tsv
mfd predict-annotations annotation_hits.tsv --genome-id candidate_001
```

To run the external CPU tools directly:

```bash
mfd run-eggnog proteins.faa --output-dir outputs/eggnog --data-dir /path/to/eggnog_data > annotation_hits.tsv
mfd run-domtblout proteins.faa --hmm /path/to/Pfam-A.hmm --out outputs/pfam.domtblout --database Pfam >> annotation_hits.tsv
mfd predict-annotations annotation_hits.tsv --genome-id candidate_001
```

The runner commands skip existing non-empty output files unless `--force` is
provided. This is intentionally CPU-first and cache-friendly, reusing the
annotation strategy from the earlier `microbe-foundation` work instead of
starting with expensive GPU embedding runs.

## No-GPU Benchmark Baseline

The first learned baseline uses binary annotation-accession features and a
simple per-label classifier:

```bash
mfd validate-splits examples/benchmark_labels.tsv
mfd build-features examples/multi_genome_annotation_hits.tsv --out outputs/features.json
mfd train-baseline outputs/features.json examples/benchmark_labels.tsv --out outputs/model.json
mfd evaluate outputs/model.json outputs/features.json examples/benchmark_labels.tsv --split test
mfd evaluate-ranking outputs/model.json outputs/features.json examples/benchmark_labels.tsv --split test --target biofuels_industrial:cellulose_degradation --k 1 --k 10
mfd evaluate-ranking outputs/model.json outputs/features.json examples/benchmark_labels.tsv --split test --panel biofuels_industrial --k 1 --k 10
mfd predict-baseline outputs/model.json outputs/features.json --genome-id G3
mfd rank-candidates outputs/model.json outputs/features.json --target biofuels_industrial:cellulose_degradation --limit 10
mfd rank-candidates outputs/model.json outputs/features.json --panel biofuels_industrial --limit 10
```

Use `evaluate` for per-label classification accuracy and `evaluate-ranking`
for the discovery question: are true positives enriched near the top of the
candidate list? The ranking report includes `hits_at_k`, `precision_at_k`,
`recall_at_k`, and the labeled held-out candidates in ranked order.

Label TSV format:

```text
genome_id	split	family	panel	label	value
G1	train	FamilyA	biofuels_industrial	cellulose_degradation	1
G3	test	FamilyB	biofuels_industrial	cellulose_degradation	1
```

Feature input is the normalized multi-genome annotation TSV:

```text
genome_id	protein_id	database	accession	name	evalue
G1	p1	CAZy	GH5	glycoside hydrolase family 5 cellulase	1e-40
```

## Current Real-Data Baseline

Using local artifacts from the earlier `microbe-foundation` repo, the importer
currently produces:

- `data/legacy/labels.tsv`: 1,258,841 labels, 60,611 genomes, 114 targets
- `data/legacy/features.1000.json`: 5,838 genomes x 1,000 eggNOG features
- `data/legacy/model.1000.json`: 107 trained target models

Initial family-held-out test result:

- overall accuracy: 0.743 across 39,509 evaluated labels
- thermophile ranking: precision@10 = 0.50, precision@50 = 0.16
- human pathogenicity ranking: precision@10 = 0.00, precision@50 = 0.16

Cached ESM2, BacFormer, ESM2+BacFormer, and hybrid-v3 feature files have also
been evaluated with both median-binarized and direct dense heads. Current
comparison: eggNOG remains stronger for broad held-out classification, but
direct BacFormer and ESM2+BacFormer heads beat eggNOG on thermophile ranking.
See [Legacy baseline comparison](docs/legacy-baseline-comparison.md).

Target-specific late fusion is now implemented with `evaluate-fusion-ranking`.
On the priority targets, eggNOG + BacFormer and eggNOG + ESM2+BacFormer fusion
reach thermophile precision@10 = 0.70 and precision@50 = 0.82, while eggNOG +
Hybrid v3 fusion preserves the strongest human pathogenicity precision@10 =
0.40. The useful-function ranking layer should therefore select feature sources
per target instead of assuming one global model wins every application.

The all-target validation-selected leaderboard evaluates 103 targets and finds
37 where dense or fusion beats annotation. It also emits top held-out discovery
candidate rankings for 25 high-precision targets. See
[Target leaderboard](docs/target-leaderboard.md) and
[Discovery candidate report](docs/discovery-candidate-report.md).
The lead exporter also produces a strict low/unknown-risk shortlist in
[Safe leads report](docs/safe-leads-report.md) and a broader application review
shortlist in [Product leads report](docs/product-leads-report.md).

This is a CPU baseline and a data-integration check, not the final foundation
model. The next model step is calibration: validate the selected fusion weights
on additional useful-function targets and only advance to a shared genome
encoder for targets where embeddings show lift over annotation-only features.
