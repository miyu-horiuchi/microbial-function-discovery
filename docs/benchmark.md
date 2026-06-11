# Function Discovery Benchmark

The benchmark measures whether a model can discover useful microbial functions
in organisms outside the taxonomic comfort zone of training data.

## Primary Question

Given a genome, MAG, SAG, or protein FASTA from a microbe that may never have
been cultured, can the model predict useful functions, application potential,
biosafety risk, and supporting genes/pathways?

## Evaluation Protocol

The primary split is **family-held-out**. Whole taxonomic families must not span
train, validation, and test sets.

Challenge subsets should include stricter class- or phylum-held-out evaluation
when label coverage supports it.

Required leakage controls:

- no target label fields as features
- no raw family/genus/species strings as predictive features
- close homolog removal for gene-level function claims
- exact train/validation/test genome identifiers released with each benchmark
- taxonomy-majority baseline reported for every label panel

## Metrics

| Metric | Purpose |
|---|---|
| AUROC / macro-F1 / accuracy | Basic per-function prediction quality |
| precision@k | Whether top-ranked candidates are useful for discovery workflows |
| calibration error | Whether confidence is trustworthy |
| novelty-stratified performance | Whether the model fails gracefully on distant clades |
| evidence enrichment | Whether highlighted genes/pathways match known function databases |
| ablation sensitivity | Whether predictions depend on highlighted evidence |

The scaffold reports top-k discovery metrics with `mfd evaluate-ranking`.
Target-level ranking evaluates one function label, while panel-level ranking
treats any positive label inside the application panel as a useful hit.

## Baselines

Every report should include:

- taxonomy-majority baseline
- nearest-neighbor by genome or protein similarity
- eggNOG / KEGG / Pfam feature baseline
- ESM-2 mean-pool baseline
- per-protein attention baseline
- public benchmark comparisons where labels overlap

## Application Panels

The benchmark is organized into panels so the shared model can be evaluated
across useful-function discovery tasks without forcing every label source into
one metric.

### Environmental / Terraforming

Targets include carbon fixation, nitrogen fixation, nitrification,
denitrification, sulfur metabolism, methane metabolism, stress tolerance,
metal reduction, plastic degradation, and survival in extreme conditions.

### Therapeutics / Antimicrobials

Targets include biosynthetic gene cluster potential, antimicrobial peptide
potential, pathogen suppression, microbiome-relevant functions, toxin genes,
virulence factors, and AMR risk.

### Biofuels / Industrial Enzymes

Targets include cellulose, xylan, lignin, starch, lipid, and fermentation
pathways, plus enzyme stability proxies such as thermotolerance, acid tolerance,
and halotolerance.

### Food / Fermentation / Agriculture

Targets include fermentation functions, flavor/aroma metabolism, spoilage risk,
plant growth promotion, nitrogen fixation, phosphate solubilization, probiotic
traits, and biosafety filters.

### Biosafety

Targets include pathogenicity, virulence factors, antimicrobial resistance,
toxins, biosafety level proxies, and confidence caveats under novelty.

## Evidence Validation

Predictions should include genes, proteins, domains, or pathways that support
the result. Evidence is evaluated separately from prediction accuracy.

Examples:

- VFDB enrichment for virulence predictions
- CARD or AMRFinderPlus enrichment for AMR predictions
- CAZy enrichment for carbohydrate-active enzyme predictions
- KEGG/MetaCyc module completeness for metabolic predictions
- antiSMASH/MIBiG evidence for biosynthetic potential

## Output Contract

Each model should emit one JSON object per genome following
`schemas/prediction.schema.json`.
