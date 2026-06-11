# Discovery Candidate Report

This report annotates ranked held-out candidates with taxonomy, genome accessions, evidence, and biosafety flags.

## Summary

- Candidate target sets: 25
- Split: test
- Minimum precision filter: 0.5
- Biosafety risk counts: high=14, low=7, moderate=174, unknown=13

## Top Candidates

| Target | Rank | Genome | Species | Source | Score | Safety | Evidence |
|---|---:|---|---|---|---:|---|---|
| `biofuels_industrial:carbon_utilization__trehalose` | 1 | 1268 | Mesobacillus foraminis | fusion:hybrid_esm2_bacformer | 1.000 | moderate | eggNOG:1RM9G |
| `biosafety:amr_phenotype__oxacillin` | 1 | 131334 | Hymenobacter latericoloratus | annotation | 1.000 | moderate | eggNOG:1RJMA |
| `biosafety:biosafety_level__bsl_1` | 1 | 10173 | Nakamurella multipartita | hybrid_esm2_bacformer | 1.000 | moderate | eggNOG:1RMV2 |
| `biosafety:gram_stain__negative` | 1 | 10423 | Nitratifractor salsuginis | fusion:hybrid_esm2_bacformer | 1.000 | moderate | eggNOG:1RMCS |
| `biosafety:gram_stain__positive` | 1 | 10173 | Nakamurella multipartita | fusion:hybrid_esm2_bacformer | 1.000 | moderate | eggNOG:1RMGZ |
| `environmental_terraforming:catalase_activity` | 1 | 10173 | Nakamurella multipartita | fusion:hybrid_esm2_bacformer | 1.000 | moderate | eggNOG:1RMD7 |
| `food_fermentation_agriculture:metabolite_production__acetate` | 1 | 132618 | Agathobaculum butyriciproducens | annotation | 1.000 | moderate | eggNOG:1RJ03 |
| `food_fermentation_agriculture:metabolite_production__butyrate` | 1 | 132618 | Agathobaculum butyriciproducens | annotation | 1.000 | moderate | eggNOG:1RJ5B |
| `food_fermentation_agriculture:metabolite_production__carbon_dioxide` | 1 | 14324 | Gracilinema caldarium | annotation | 1.000 | moderate | eggNOG:1RIZT |
| `food_fermentation_agriculture:metabolite_production__dihydrogen` | 1 | 14324 | Gracilinema caldarium | annotation | 1.000 | moderate | eggNOG:1RJ39 |
| `food_fermentation_agriculture:metabolite_production__formate` | 1 | 17805 | Alkaliflexus imshenetskii | annotation | 1.000 | moderate | eggNOG:1RJ03 |
| `food_fermentation_agriculture:metabolite_production__propionate` | 1 | 12660 | Propionicimonas paludicola | annotation | 1.000 | moderate | eggNOG:1RIZT |
| `biofuels_industrial:carbon_utilization__cellobiose` | 1 | 10173 | Nakamurella multipartita | fusion:hybrid_esm2_bacformer | 1.000 | moderate | eggNOG:1RMEC |
| `biofuels_industrial:carbon_utilization__d_mannose` | 1 | 130271 | Dyadobacter jejuensis | fusion:bacformer | 1.000 | moderate | eggNOG:1RNAW |
| `biofuels_industrial:temperature_class__mesophile` | 1 | 10173 | Nakamurella multipartita | fusion:hybrid_esm2_bacformer | 1.000 | moderate | eggNOG:1RP5X |
| `environmental_terraforming:cytochrome_oxidase_activity` | 1 | 10173 | Nakamurella multipartita | fusion:hybrid_esm2_bacformer | 1.000 | moderate | eggNOG:1RM8M |
| `environmental_terraforming:oxygen_tolerance__obligate_aerobe` | 1 | 10173 | Nakamurella multipartita | fusion:hybrid_esm2_bacformer | 1.000 | moderate | eggNOG:1RMCI |
| `environmental_terraforming:sporulation` | 1 | 10173 | Nakamurella multipartita | fusion:hybrid_v3 | 1.000 | moderate | eggNOG:1RMBA |
| `food_fermentation_agriculture:pigmentation` | 1 | 12659 | Propionibacterium acidifaciens | fusion:bacformer | 1.000 | high | eggNOG:1RJA2 |
| `biosafety:amr_phenotype__streptomycin` | 1 | 132999 | Pseudoalteromonas gelatinilytica | fusion:hybrid_esm2_bacformer | 1.000 | moderate | eggNOG:1RNC0 |
| `biofuels_industrial:ph_class__alkaliphile` | 1 | 10174 | Nakamurella lactea | fusion:hybrid_esm2_bacformer | 1.000 | moderate | eggNOG:1RJUK |
| `biosafety:amr_phenotype__gentamicin` | 1 | 140643 | Salinigranum salinum | bacformer | 1.000 | unknown | candidate_score |
| `biosafety:amr_phenotype__nalidixic_acid` | 1 | 1251 | Neobacillus fumarioli | hybrid_esm2_bacformer | 1.000 | low | candidate_score |
| `environmental_terraforming:halophily__halophile` | 1 | 11163 | Thermasporomyces composti | fusion:hybrid_esm2_bacformer | 1.000 | moderate | eggNOG:1RMD0 |
| `food_fermentation_agriculture:cell_shape__rod` | 1 | 10173 | Nakamurella multipartita | fusion:hybrid_esm2_bacformer | 1.000 | moderate | eggNOG:1RK33 |

## Notes

- These are ranked benchmark candidates, not wet-lab validated recommendations.
- Safety flags combine known BacDive fields with annotation-baseline biosafety scores when available.
- Evidence lists active annotation features with positive target weights; dense-only leads may have limited feature evidence.
