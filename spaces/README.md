---
title: Microbial Function Discovery
sdk: gradio
app_file: app.py
license: mit
---

# Microbial Function Discovery Demo

This is a no-GPU Gradio browser for exported useful-function lead TSVs. It is
designed for fast review of candidate microbes before wet-lab validation.

Run locally:

```bash
uv run --with gradio python spaces/app.py
```

The demo loads `examples/product_leads.review.tsv` by default and can also show
`examples/safe_leads.review.tsv`, plus smaller sample files for smoke tests. To
point it at a local export:

```bash
LEADS_TSV=data/legacy/product_leads_review.tsv uv run --with gradio python spaces/app.py
```

Expected TSV columns match the `mfd export-safe-leads` output:

```text
panel target_key label source target_precision rank genome_id species genus family accession score risk_level biosafety_flags evidence
```

Deploy this folder as a public Hugging Face Space:

```bash
hf repo create miyuiu/microbial-function-discovery --repo-type space --space_sdk gradio --exist-ok
hf upload miyuiu/microbial-function-discovery spaces . --repo-type space --commit-message "Deploy microbial lead review demo"
```

Current public app: https://miyuiu-microbial-function-discovery.hf.space
