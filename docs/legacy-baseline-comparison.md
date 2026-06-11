# Legacy Baseline Comparison

This report compares CPU-friendly feature sources imported from the earlier
`microbe-foundation` repo on the same BacDive-derived benchmark format.

All runs use:

- labels: `data/legacy/labels.tsv`
- split: `family_split`, evaluated on `test`
- model: the scaffold no-GPU one-vs-rest baseline
- dense embeddings: high-variance dimensions median-binarized before training

## Results

| Feature source | Test labels evaluated | Targets evaluated | Overall accuracy | Thermophile P@10 | Thermophile P@50 | Human pathogenicity P@10 | Human pathogenicity P@50 |
|---|---:|---:|---:|---:|---:|---:|---:|
| eggNOG top 1,000 | 39,509 | 101 | 0.7433 | 0.50 | 0.16 | 0.00 | 0.16 |
| ESM2 640-dim | 100,150 | 107 | 0.5543 | 0.50 | 0.18 | 0.10 | 0.12 |
| BacFormer 960-dim | 38,333 | 104 | 0.5612 | 0.30 | 0.14 | 0.20 | 0.30 |
| ESM2 + BacFormer | 37,309 | 104 | 0.5612 | 0.30 | 0.14 | 0.20 | 0.16 |
| Hybrid v3 | 28,679 | 101 | 0.5464 | 0.20 | 0.16 | 0.00 | 0.10 |

## Interpretation

The annotation baseline is still the strongest default for broad held-out trait
classification. ESM2 has wider genome coverage but lower overall accuracy under
the current median-binarized baseline. BacFormer is weaker overall but improves
human pathogenicity precision@50, suggesting embeddings may help specific
safety/risk tasks once modeled with a better dense-feature learner.

The immediate model lesson is not "start GPU training." It is:

- keep eggNOG as the required CPU baseline
- add a dense-feature learner before judging embeddings
- evaluate each application target separately, because broad panel metrics can
  hide target-level failures
- only move to expensive GPU foundation-model training after embeddings beat
  the eggNOG baseline on priority targets

## Reproduction Commands

```bash
mfd import-legacy-dense-features /Users/miyuhoriuchi/microbe-foundation/data/esm2_features.npz \
  --labels data/legacy/labels.tsv \
  --out data/legacy/features.esm2.json \
  --max-features 640 \
  --feature-prefix ESM2
mfd train-baseline data/legacy/features.esm2.json data/legacy/labels.tsv --out data/legacy/model.esm2.json
mfd evaluate data/legacy/model.esm2.json data/legacy/features.esm2.json data/legacy/labels.tsv \
  --split test \
  --out data/legacy/eval.test.esm2.json
mfd evaluate-ranking data/legacy/model.esm2.json data/legacy/features.esm2.json data/legacy/labels.tsv \
  --split test \
  --target biofuels_industrial:temperature_class__thermophile \
  --k 10 --k 50
```
