# Top 20 Microbial Validation Packets

These packets convert model-ranked microbial leads into review units for partner triage.
They are not release or deployment recommendations.

## Summary

- Packets: 20
- Risk mix: moderate=20
- Panels: biofuels_industrial, biosafety, environmental_terraforming, food_fermentation_agriculture

## Packet 1: Mesobacillus foraminis

- Genome ID: `1268`
- Accession: `GCA_004340465`
- Candidate function: carbon utilization: trehalose
- Target key: `biofuels_industrial:carbon_utilization__trehalose`
- Application panel: `biofuels_industrial`
- Why it matters: carbon utilization: trehalose can prioritize strains or enzymes for feedstock conversion and industrial biocatalysis.
- Model support: source `fusion:hybrid_esm2_bacformer`, target precision `1.0`, rank `1`, score `1.0`
- Evidence to review: `eggNOG:1RM9G`, `eggNOG:1RMCX`, `eggNOG:1RMCG`, `eggNOG:1RK5N`, `eggNOG:1RMGJ`
- Biosafety review: risk `moderate`, flags `predicted_human_pathogenicity,predicted_amr_signal`
- First validation step: Run a growth or activity assay with trehalose as the target carbon substrate.

## Packet 2: Hymenobacter latericoloratus

- Genome ID: `131334`
- Accession: `GCA_014199535`
- Candidate function: amr phenotype: oxacillin
- Target key: `biosafety:amr_phenotype__oxacillin`
- Application panel: `biosafety`
- Why it matters: amr phenotype: oxacillin is useful for safety characterization and go/no-go triage before application testing.
- Model support: source `annotation`, target precision `1.0`, rank `1`, score `1.0`
- Evidence to review: `eggNOG:1RJMA`, `eggNOG:1RK1X`, `eggNOG:1RM7T`, `eggNOG:1RMH2`, `eggNOG:1RMID`
- Biosafety review: risk `moderate`, flags `predicted_human_pathogenicity,predicted_amr_signal`
- First validation step: Confirm the predicted safety phenotype with strain-level metadata and standard biosafety screens.

## Packet 3: Nakamurella multipartita

- Genome ID: `10173`
- Accession: `GCA_000024365`
- Candidate function: biosafety level: bsl 1
- Target key: `biosafety:biosafety_level__bsl_1`
- Application panel: `biosafety`
- Why it matters: biosafety level: bsl 1 is useful for safety characterization and go/no-go triage before application testing.
- Model support: source `hybrid_esm2_bacformer`, target precision `1.0`, rank `1`, score `1.0`
- Evidence to review: `eggNOG:1RMV2`, `eggNOG:1RMDZ`, `eggNOG:1RN6V`, `eggNOG:1RJRA`, `eggNOG:1RJTG`
- Biosafety review: risk `moderate`, flags `predicted_human_pathogenicity,predicted_amr_signal`
- First validation step: Confirm the predicted safety phenotype with strain-level metadata and standard biosafety screens.

## Packet 4: Nitratifractor salsuginis

- Genome ID: `10423`
- Accession: `GCA_000186245`
- Candidate function: gram stain: negative
- Target key: `biosafety:gram_stain__negative`
- Application panel: `biosafety`
- Why it matters: gram stain: negative is useful for safety characterization and go/no-go triage before application testing.
- Model support: source `fusion:hybrid_esm2_bacformer`, target precision `1.0`, rank `1`, score `1.0`
- Evidence to review: `eggNOG:1RMCS`, `eggNOG:1RMHU`, `eggNOG:1RR0J`, `eggNOG:1RM8P`, `eggNOG:1RP3W`
- Biosafety review: risk `moderate`, flags `predicted_human_pathogenicity,predicted_amr_signal`
- First validation step: Confirm the predicted safety phenotype with strain-level metadata and standard biosafety screens.

## Packet 5: Nakamurella multipartita

- Genome ID: `10173`
- Accession: `GCA_000024365`
- Candidate function: gram stain: positive
- Target key: `biosafety:gram_stain__positive`
- Application panel: `biosafety`
- Why it matters: gram stain: positive is useful for safety characterization and go/no-go triage before application testing.
- Model support: source `fusion:hybrid_esm2_bacformer`, target precision `1.0`, rank `1`, score `1.0`
- Evidence to review: `eggNOG:1RMGZ`, `eggNOG:1RNAN`, `eggNOG:1RMI4`, `eggNOG:1RMHS`, `eggNOG:1RMFP`
- Biosafety review: risk `moderate`, flags `predicted_human_pathogenicity,predicted_amr_signal`
- First validation step: Confirm the predicted safety phenotype with strain-level metadata and standard biosafety screens.

## Packet 6: Nakamurella multipartita

- Genome ID: `10173`
- Accession: `GCA_000024365`
- Candidate function: catalase activity
- Target key: `environmental_terraforming:catalase_activity`
- Application panel: `environmental_terraforming`
- Why it matters: catalase activity can prioritize organisms for environmental stress response, remediation, or closed-system bioprocessing assays.
- Model support: source `fusion:hybrid_esm2_bacformer`, target precision `1.0`, rank `1`, score `1.0`
- Evidence to review: `eggNOG:1RMD7`, `eggNOG:1RM8M`, `eggNOG:1RM8P`, `eggNOG:1RMYV`, `eggNOG:1RJNQ`
- Biosafety review: risk `moderate`, flags `predicted_human_pathogenicity,predicted_amr_signal`
- First validation step: Run catalase and peroxide-stress assays, then compare activity against close relatives.

## Packet 7: Agathobaculum butyriciproducens

- Genome ID: `132618`
- Accession: `GCA_003096535`
- Candidate function: metabolite production: acetate
- Target key: `food_fermentation_agriculture:metabolite_production__acetate`
- Application panel: `food_fermentation_agriculture`
- Why it matters: metabolite production: acetate can prioritize microbes for fermentation, metabolite production, food, or agriculture screens.
- Model support: source `annotation`, target precision `1.0`, rank `1`, score `1.0`
- Evidence to review: `eggNOG:1RJ03`, `eggNOG:1RJ05`, `eggNOG:1RJ13`, `eggNOG:1RJ19`, `eggNOG:1RJ1Q`
- Biosafety review: risk `moderate`, flags `predicted_human_pathogenicity,predicted_amr_signal`
- First validation step: Run targeted acetate quantification under matched culture or enrichment conditions.

## Packet 8: Agathobaculum butyriciproducens

- Genome ID: `132618`
- Accession: `GCA_003096535`
- Candidate function: metabolite production: butyrate
- Target key: `food_fermentation_agriculture:metabolite_production__butyrate`
- Application panel: `food_fermentation_agriculture`
- Why it matters: metabolite production: butyrate can prioritize microbes for fermentation, metabolite production, food, or agriculture screens.
- Model support: source `annotation`, target precision `1.0`, rank `1`, score `1.0`
- Evidence to review: `eggNOG:1RJ5B`, `eggNOG:1RJ6A`, `eggNOG:1RJ7D`, `eggNOG:1RJ8K`, `eggNOG:1RJ8P`
- Biosafety review: risk `moderate`, flags `predicted_human_pathogenicity,predicted_amr_signal`
- First validation step: Run targeted butyrate quantification under matched culture or enrichment conditions.

## Packet 9: Gracilinema caldarium

- Genome ID: `14324`
- Accession: `GCA_000219725`
- Candidate function: metabolite production: carbon dioxide
- Target key: `food_fermentation_agriculture:metabolite_production__carbon_dioxide`
- Application panel: `food_fermentation_agriculture`
- Why it matters: metabolite production: carbon dioxide can prioritize microbes for fermentation, metabolite production, food, or agriculture screens.
- Model support: source `annotation`, target precision `1.0`, rank `1`, score `1.0`
- Evidence to review: `eggNOG:1RIZT`, `eggNOG:1RJ05`, `eggNOG:1RJ13`, `eggNOG:1RJ19`, `eggNOG:1RJ1K`
- Biosafety review: risk `moderate`, flags `predicted_human_pathogenicity,predicted_amr_signal`
- First validation step: Run targeted carbon dioxide quantification under matched culture or enrichment conditions.

## Packet 10: Gracilinema caldarium

- Genome ID: `14324`
- Accession: `GCA_000219725`
- Candidate function: metabolite production: dihydrogen
- Target key: `food_fermentation_agriculture:metabolite_production__dihydrogen`
- Application panel: `food_fermentation_agriculture`
- Why it matters: metabolite production: dihydrogen can prioritize microbes for fermentation, metabolite production, food, or agriculture screens.
- Model support: source `annotation`, target precision `1.0`, rank `1`, score `1.0`
- Evidence to review: `eggNOG:1RJ39`, `eggNOG:1RJ5B`, `eggNOG:1RJ7D`, `eggNOG:1RJ8K`, `eggNOG:1RJ8P`
- Biosafety review: risk `moderate`, flags `predicted_human_pathogenicity,predicted_amr_signal`
- First validation step: Run targeted dihydrogen quantification under matched culture or enrichment conditions.

## Packet 11: Alkaliflexus imshenetskii

- Genome ID: `17805`
- Accession: `GCA_000517065`
- Candidate function: metabolite production: formate
- Target key: `food_fermentation_agriculture:metabolite_production__formate`
- Application panel: `food_fermentation_agriculture`
- Why it matters: metabolite production: formate can prioritize microbes for fermentation, metabolite production, food, or agriculture screens.
- Model support: source `annotation`, target precision `1.0`, rank `1`, score `1.0`
- Evidence to review: `eggNOG:1RJ03`, `eggNOG:1RJ05`, `eggNOG:1RJ13`, `eggNOG:1RJ19`, `eggNOG:1RJ1Q`
- Biosafety review: risk `moderate`, flags `predicted_human_pathogenicity,predicted_amr_signal`
- First validation step: Run targeted formate quantification under matched culture or enrichment conditions.

## Packet 12: Propionicimonas paludicola

- Genome ID: `12660`
- Accession: `GCA_002563675`
- Candidate function: metabolite production: propionate
- Target key: `food_fermentation_agriculture:metabolite_production__propionate`
- Application panel: `food_fermentation_agriculture`
- Why it matters: metabolite production: propionate can prioritize microbes for fermentation, metabolite production, food, or agriculture screens.
- Model support: source `annotation`, target precision `1.0`, rank `1`, score `1.0`
- Evidence to review: `eggNOG:1RIZT`, `eggNOG:1RJ03`, `eggNOG:1RJ05`, `eggNOG:1RJ13`, `eggNOG:1RJ19`
- Biosafety review: risk `moderate`, flags `predicted_human_pathogenicity,predicted_amr_signal`
- First validation step: Run targeted propionate quantification under matched culture or enrichment conditions.

## Packet 13: Listeria monocytogenes

- Genome ID: `131740`
- Accession: `GCA_900637785`
- Candidate function: carbon utilization: trehalose
- Target key: `biofuels_industrial:carbon_utilization__trehalose`
- Application panel: `biofuels_industrial`
- Why it matters: carbon utilization: trehalose can prioritize strains or enzymes for feedstock conversion and industrial biocatalysis.
- Model support: source `fusion:hybrid_esm2_bacformer`, target precision `1.0`, rank `2`, score `1.0`
- Evidence to review: `eggNOG:1RM9G`, `eggNOG:1RMCX`, `eggNOG:1RMCG`, `eggNOG:1RK5N`, `eggNOG:1RMGJ`
- Biosafety review: risk `moderate`, flags `known_bsl-2,predicted_human_pathogenicity,predicted_amr_signal`
- First validation step: Run a growth or activity assay with trehalose as the target carbon substrate.

## Packet 14: Hymenobacter luteus

- Genome ID: `131335`
- Accession: `GCA_014202325`
- Candidate function: amr phenotype: oxacillin
- Target key: `biosafety:amr_phenotype__oxacillin`
- Application panel: `biosafety`
- Why it matters: amr phenotype: oxacillin is useful for safety characterization and go/no-go triage before application testing.
- Model support: source `annotation`, target precision `1.0`, rank `2`, score `1.0`
- Evidence to review: `eggNOG:1RJMA`, `eggNOG:1RK1X`, `eggNOG:1RM7T`, `eggNOG:1RMH2`, `eggNOG:1RMID`
- Biosafety review: risk `moderate`, flags `predicted_human_pathogenicity,predicted_amr_signal`
- First validation step: Confirm the predicted safety phenotype with strain-level metadata and standard biosafety screens.

## Packet 15: Nakamurella lactea

- Genome ID: `10174`
- Accession: `GCA_000426645`
- Candidate function: biosafety level: bsl 1
- Target key: `biosafety:biosafety_level__bsl_1`
- Application panel: `biosafety`
- Why it matters: biosafety level: bsl 1 is useful for safety characterization and go/no-go triage before application testing.
- Model support: source `hybrid_esm2_bacformer`, target precision `1.0`, rank `2`, score `1.0`
- Evidence to review: `eggNOG:1RMV2`, `eggNOG:1RMDZ`, `eggNOG:1RN6V`, `eggNOG:1RJRA`, `eggNOG:1RJTG`
- Biosafety review: risk `moderate`, flags `predicted_human_pathogenicity,predicted_amr_signal`
- First validation step: Confirm the predicted safety phenotype with strain-level metadata and standard biosafety screens.

## Packet 16: Caminibacter mediatlanticus

- Genome ID: `10428`
- Accession: `GCA_005843985`
- Candidate function: gram stain: negative
- Target key: `biosafety:gram_stain__negative`
- Application panel: `biosafety`
- Why it matters: gram stain: negative is useful for safety characterization and go/no-go triage before application testing.
- Model support: source `fusion:hybrid_esm2_bacformer`, target precision `1.0`, rank `2`, score `1.0`
- Evidence to review: `eggNOG:1RMCS`, `eggNOG:1RMHU`, `eggNOG:1RR0J`, `eggNOG:1RM8P`, `eggNOG:1RP3W`
- Biosafety review: risk `moderate`, flags `predicted_human_pathogenicity,predicted_amr_signal`
- First validation step: Confirm the predicted safety phenotype with strain-level metadata and standard biosafety screens.

## Packet 17: Nakamurella lactea

- Genome ID: `10174`
- Accession: `GCA_000426645`
- Candidate function: gram stain: positive
- Target key: `biosafety:gram_stain__positive`
- Application panel: `biosafety`
- Why it matters: gram stain: positive is useful for safety characterization and go/no-go triage before application testing.
- Model support: source `fusion:hybrid_esm2_bacformer`, target precision `1.0`, rank `2`, score `1.0`
- Evidence to review: `eggNOG:1RMGZ`, `eggNOG:1RNAN`, `eggNOG:1RMI4`, `eggNOG:1RMHS`, `eggNOG:1RMFP`
- Biosafety review: risk `moderate`, flags `predicted_human_pathogenicity,predicted_amr_signal`
- First validation step: Confirm the predicted safety phenotype with strain-level metadata and standard biosafety screens.

## Packet 18: Nakamurella lactea

- Genome ID: `10174`
- Accession: `GCA_000426645`
- Candidate function: catalase activity
- Target key: `environmental_terraforming:catalase_activity`
- Application panel: `environmental_terraforming`
- Why it matters: catalase activity can prioritize organisms for environmental stress response, remediation, or closed-system bioprocessing assays.
- Model support: source `fusion:hybrid_esm2_bacformer`, target precision `1.0`, rank `2`, score `1.0`
- Evidence to review: `eggNOG:1RMD7`, `eggNOG:1RM8M`, `eggNOG:1RM8P`, `eggNOG:1RMYV`, `eggNOG:1RJNQ`
- Biosafety review: risk `moderate`, flags `predicted_human_pathogenicity,predicted_amr_signal`
- First validation step: Run catalase and peroxide-stress assays, then compare activity against close relatives.

## Packet 19: Gracilinema caldarium

- Genome ID: `14324`
- Accession: `GCA_000219725`
- Candidate function: metabolite production: acetate
- Target key: `food_fermentation_agriculture:metabolite_production__acetate`
- Application panel: `food_fermentation_agriculture`
- Why it matters: metabolite production: acetate can prioritize microbes for fermentation, metabolite production, food, or agriculture screens.
- Model support: source `annotation`, target precision `1.0`, rank `2`, score `1.0`
- Evidence to review: `eggNOG:1RJ05`, `eggNOG:1RJ13`, `eggNOG:1RJ19`, `eggNOG:1RJ1Q`, `eggNOG:1RJ1T`
- Biosafety review: risk `moderate`, flags `predicted_human_pathogenicity,predicted_amr_signal`
- First validation step: Run targeted acetate quantification under matched culture or enrichment conditions.

## Packet 20: Anaerotruncus colihominis

- Genome ID: `17679`
- Accession: `GCA_025146135`
- Candidate function: metabolite production: butyrate
- Target key: `food_fermentation_agriculture:metabolite_production__butyrate`
- Application panel: `food_fermentation_agriculture`
- Why it matters: metabolite production: butyrate can prioritize microbes for fermentation, metabolite production, food, or agriculture screens.
- Model support: source `annotation`, target precision `1.0`, rank `2`, score `1.0`
- Evidence to review: `eggNOG:1RJ5B`, `eggNOG:1RJ6A`, `eggNOG:1RJ7D`, `eggNOG:1RJ8K`, `eggNOG:1RJ8P`
- Biosafety review: risk `moderate`, flags `predicted_human_pathogenicity,predicted_amr_signal`
- First validation step: Run targeted butyrate quantification under matched culture or enrichment conditions.

## Operating Notes

- Confirm taxonomy, accession metadata, and strain availability before experimental planning.
- Re-check biosafety and AMR/pathogenicity signals before culturing or partner handoff.
- Treat model evidence as prioritization support; require orthogonal assay confirmation.
