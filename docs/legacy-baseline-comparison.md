# Legacy Baseline Comparison

This report compares CPU-friendly feature sources imported from the earlier
`microbe-foundation` repo on the same BacDive-derived benchmark format.

All runs use:

- labels: `data/legacy/labels.tsv`
- split: `family_split`, evaluated on `test`
- model: the scaffold no-GPU one-vs-rest baseline
- median-bin embeddings: high-variance dimensions median-binarized before training
- direct dense embeddings: raw float dimensions with standardized linear target heads

## Results

| Feature source | Test labels evaluated | Targets evaluated | Overall accuracy | Thermophile P@10 | Thermophile P@50 | Human pathogenicity P@10 | Human pathogenicity P@50 |
|---|---:|---:|---:|---:|---:|---:|---:|
| eggNOG top 1,000 | 39,509 | 101 | 0.7433 | 0.50 | 0.16 | 0.00 | 0.16 |
| ESM2 median-bin | 100,150 | 107 | 0.5543 | 0.50 | 0.18 | 0.10 | 0.12 |
| BacFormer median-bin | 38,333 | 104 | 0.5612 | 0.30 | 0.14 | 0.20 | 0.30 |
| ESM2 direct dense | 100,150 | 107 | 0.6415 | 0.50 | 0.32 | 0.10 | 0.12 |
| BacFormer direct dense | 38,333 | 104 | 0.7282 | 0.70 | 0.50 | 0.00 | 0.14 |
| ESM2 + BacFormer direct dense | 37,309 | 104 | 0.7177 | 0.70 | 0.52 | 0.20 | 0.28 |
| Hybrid v3 direct dense | 28,679 | 101 | 0.7138 | 0.20 | 0.26 | 0.40 | 0.16 |

## Interpretation

The annotation baseline is still the strongest default for broad held-out trait
classification. The direct dense learner substantially improves embedding
results over median binning: BacFormer and ESM2+BacFormer now beat eggNOG on
thermophile ranking, while Hybrid v3 has the best human pathogenicity
precision@10. Embeddings are therefore useful for target-level discovery
ranking, but they do not yet beat eggNOG for broad trait classification.

The immediate model lesson is:

- keep eggNOG as the required CPU baseline
- use direct dense embedding heads for priority target ranking
- evaluate each application target separately, because broad panel metrics can
  hide target-level failures
- move toward a shared genome encoder only for targets where embeddings already
  show lift over annotation-only features

## Reproduction Commands

```bash
mfd train-dense-baseline /Users/miyuhoriuchi/microbe-foundation/data/esm2_features.npz \
  data/legacy/labels.tsv \
  --out data/legacy/dense_model.esm2.json \
  --max-features 640
mfd evaluate-dense data/legacy/dense_model.esm2.json \
  /Users/miyuhoriuchi/microbe-foundation/data/esm2_features.npz \
  data/legacy/labels.tsv \
  --split test \
  --out data/legacy/dense_eval.test.esm2.json
mfd evaluate-dense-ranking data/legacy/dense_model.esm2.json \
  /Users/miyuhoriuchi/microbe-foundation/data/esm2_features.npz \
  data/legacy/labels.tsv \
  --split test \
  --target biofuels_industrial:temperature_class__thermophile \
  --k 10 --k 50
```
