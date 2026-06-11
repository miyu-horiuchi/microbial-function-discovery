# Target Leaderboard

This report summarizes validation-selected source choices for the legacy
BacDive-derived useful-function benchmark.

All target winners were selected on `val` using `precision_at_10`, then reported
on held-out `test`. The leaderboard compares:

- eggNOG annotation baseline
- BacFormer direct dense head
- ESM2+BacFormer direct dense head
- Hybrid v3 direct dense head
- target-specific late fusion of eggNOG plus each dense source

Raw outputs:

- `data/legacy/target_leaderboard.json`
- `data/legacy/discovery_candidates.top_targets.json`
- `data/legacy/discovery_candidates.annotated.json`
- `data/legacy/safe_leads.tsv`
- `data/legacy/product_leads_review.tsv`
- [Discovery candidate report](discovery-candidate-report.md)
- [Safe leads report](safe-leads-report.md)
- [Product leads report](product-leads-report.md)

## Summary

| Result | Count |
|---|---:|
| Targets evaluated | 103 |
| Targets where dense/fusion beats annotation | 37 |
| Top-target discovery candidate sets emitted | 25 |

## Winning Sources

| Source | Targets won |
|---|---:|
| eggNOG + ESM2+BacFormer fusion | 36 |
| eggNOG annotation | 30 |
| eggNOG + Hybrid v3 fusion | 15 |
| ESM2+BacFormer direct dense | 8 |
| Hybrid v3 direct dense | 7 |
| eggNOG + BacFormer fusion | 4 |
| BacFormer direct dense | 3 |

## Panel Coverage

| Panel | Targets evaluated | Targets beating annotation |
|---|---:|---:|
| Biofuels / industrial | 33 | 15 |
| Biosafety | 34 | 12 |
| Environmental / terraforming | 14 | 5 |
| Food / fermentation / agriculture | 22 | 5 |

## Top Held-Out Targets

| Target | Best source | Weight | Test P@10 | Test P@50 | Test positives | Beats annotation |
|---|---|---:|---:|---:|---:|---|
| `biosafety:amr_phenotype__oxacillin` | annotation | 0.0 | 1.00 | 1.00 | 3 | no |
| `food_fermentation_agriculture:metabolite_production__butyrate` | annotation | 0.0 | 1.00 | 1.00 | 4 | no |
| `food_fermentation_agriculture:metabolite_production__carbon_dioxide` | annotation | 0.0 | 1.00 | 1.00 | 3 | no |
| `food_fermentation_agriculture:metabolite_production__dihydrogen` | annotation | 0.0 | 1.00 | 1.00 | 5 | no |
| `food_fermentation_agriculture:metabolite_production__formate` | annotation | 0.0 | 1.00 | 1.00 | 3 | no |
| `food_fermentation_agriculture:metabolite_production__propionate` | annotation | 0.0 | 1.00 | 1.00 | 4 | no |
| `food_fermentation_agriculture:metabolite_production__acetate` | annotation | 0.0 | 1.00 | 1.00 | 12 | no |
| `biofuels_industrial:carbon_utilization__trehalose` | fusion:hybrid_esm2_bacformer | 0.25 | 1.00 | 0.74 | 58 | yes |
| `biosafety:gram_stain__positive` | fusion:hybrid_esm2_bacformer | 0.25 | 1.00 | 0.98 | 103 | yes |
| `biosafety:gram_stain__negative` | fusion:hybrid_esm2_bacformer | 1.0 | 1.00 | 0.94 | 163 | yes |
| `environmental_terraforming:catalase_activity` | fusion:hybrid_esm2_bacformer | 0.0 | 1.00 | 0.84 | 179 | no |
| `biosafety:biosafety_level__bsl_1` | hybrid_esm2_bacformer | 1.0 | 1.00 | 1.00 | 827 | no |
| `environmental_terraforming:sporulation` | fusion:hybrid_v3 | 0.25 | 0.90 | 0.80 | 46 | yes |
| `biofuels_industrial:carbon_utilization__d_mannose` | fusion:bacformer | 0.25 | 0.90 | 0.58 | 61 | yes |
| `biofuels_industrial:carbon_utilization__cellobiose` | fusion:hybrid_esm2_bacformer | 0.25 | 0.90 | 0.74 | 69 | yes |
| `food_fermentation_agriculture:pigmentation` | fusion:bacformer | 0.5 | 0.90 | 0.84 | 79 | yes |
| `environmental_terraforming:cytochrome_oxidase_activity` | fusion:hybrid_esm2_bacformer | 0.0 | 0.90 | 0.86 | 117 | yes |
| `environmental_terraforming:oxygen_tolerance__obligate_aerobe` | fusion:hybrid_esm2_bacformer | 0.25 | 0.90 | 0.92 | 304 | no |
| `biofuels_industrial:temperature_class__mesophile` | fusion:hybrid_esm2_bacformer | 0.25 | 0.90 | 0.98 | 390 | yes |
| `biosafety:amr_phenotype__streptomycin` | fusion:hybrid_esm2_bacformer | 0.0 | 0.83 | 0.83 | 5 | yes |

## Candidate Ranking Artifact

`data/legacy/discovery_candidates.top_targets.json` contains top-10 held-out
candidates for 25 target/source choices with test `precision_at_10 >= 0.5`.
Example first-ranked candidates:

| Target | Source | Test P@10 | Top genome | Top score |
|---|---|---:|---:|---:|
| `biofuels_industrial:carbon_utilization__trehalose` | fusion:hybrid_esm2_bacformer | 1.00 | 1268 | 1.00 |
| `biosafety:amr_phenotype__oxacillin` | annotation | 1.00 | 131334 | 1.00 |
| `biosafety:biosafety_level__bsl_1` | hybrid_esm2_bacformer | 1.00 | 10173 | 1.00 |
| `biosafety:gram_stain__negative` | fusion:hybrid_esm2_bacformer | 1.00 | 10423 | 1.00 |
| `biosafety:gram_stain__positive` | fusion:hybrid_esm2_bacformer | 1.00 | 10173 | 1.00 |
| `environmental_terraforming:catalase_activity` | fusion:hybrid_esm2_bacformer | 1.00 | 10173 | 1.00 |
| `food_fermentation_agriculture:metabolite_production__acetate` | annotation | 1.00 | 132618 | 1.00 |
| `biofuels_industrial:carbon_utilization__cellobiose` | fusion:hybrid_esm2_bacformer | 0.90 | 10173 | 1.00 |
| `biofuels_industrial:carbon_utilization__d_mannose` | fusion:bacformer | 0.90 | 130271 | 1.00 |
| `biofuels_industrial:temperature_class__mesophile` | fusion:hybrid_esm2_bacformer | 0.90 | 10173 | 1.00 |

These genome IDs are ranked held-out benchmark candidates, not wet-lab validated
recommendations. The next product step is to attach taxonomy, source metadata,
biosafety filters, and evidence features before treating any row as an
experimental lead.

That product step is now started in
[Discovery candidate report](discovery-candidate-report.md), which annotates the
top candidates with BacDive taxonomy, genome accessions, evidence features, and
biosafety flags.

The strict low/unknown-risk shortlist is in [Safe leads report](safe-leads-report.md).
The broader application review shortlist, including moderate-risk candidates for
manual safety review, is in [Product leads report](product-leads-report.md).

## Reproduction

```bash
mfd evaluate-fusion-leaderboard \
  data/legacy/model.1000.json \
  data/legacy/features.1000.json \
  data/legacy/labels.tsv \
  --dense-source bacformer data/legacy/dense_model.bacformer.json /Users/miyuhoriuchi/microbe-foundation/data/bacformer_features_clean.npz \
  --dense-source hybrid_esm2_bacformer data/legacy/dense_model.hybrid_esm2_bacformer.json /Users/miyuhoriuchi/microbe-foundation/data/hybrid_esm2_bacformer.npz \
  --dense-source hybrid_v3 data/legacy/dense_model.hybrid_v3.json /Users/miyuhoriuchi/microbe-foundation/data/hybrid_v3.npz \
  --validation-split val \
  --test-split test \
  --k 10 --k 50 \
  --weight 0 --weight 0.25 --weight 0.5 --weight 0.75 --weight 1 \
  --select-k 10 \
  --out data/legacy/target_leaderboard.json

mfd rank-discovery-candidates \
  data/legacy/model.1000.json \
  data/legacy/features.1000.json \
  data/legacy/labels.tsv \
  data/legacy/target_leaderboard.json \
  --dense-source bacformer data/legacy/dense_model.bacformer.json /Users/miyuhoriuchi/microbe-foundation/data/bacformer_features_clean.npz \
  --dense-source hybrid_esm2_bacformer data/legacy/dense_model.hybrid_esm2_bacformer.json /Users/miyuhoriuchi/microbe-foundation/data/hybrid_esm2_bacformer.npz \
  --dense-source hybrid_v3 data/legacy/dense_model.hybrid_v3.json /Users/miyuhoriuchi/microbe-foundation/data/hybrid_v3.npz \
  --split test \
  --limit-per-target 10 \
  --min-precision 0.5 \
  --precision-k 10 \
  --max-targets 25 \
  --out data/legacy/discovery_candidates.top_targets.json

mfd export-safe-leads \
  data/legacy/discovery_candidates.annotated.json \
  --out data/legacy/safe_leads.tsv \
  --report-out docs/safe-leads-report.md \
  --format tsv \
  --allowed-risk low \
  --allowed-risk unknown \
  --min-precision 0.5 \
  --precision-k 10
```
