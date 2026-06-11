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
- [Model roadmap](docs/model-roadmap.md)
- [Prediction output schema](schemas/prediction.schema.json)

## Developer Quickstart

Run the scaffold without installing dependencies:

```bash
PYTHONPATH=src python3 -m microbial_function_discovery.cli panels
PYTHONPATH=src python3 -m microbial_function_discovery.cli predict examples/useful_functions.faa --genome-id candidate_001
PYTHONPATH=src python3 -m microbial_function_discovery.cli predict-annotations examples/annotation_hits.tsv --genome-id candidate_001
PYTHONPATH=src python3 -m microbial_function_discovery.cli validate examples/prediction.example.json
PYTHONPATH=src python3 -m unittest discover -s tests -v
```

Install the local CLI:

```bash
python3 -m pip install -e .
mfd panels
mfd predict examples/useful_functions.faa --genome-id candidate_001
mfd predict-annotations examples/annotation_hits.tsv --genome-id candidate_001
mfd validate examples/prediction.example.json
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
