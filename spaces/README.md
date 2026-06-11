# Microbial Function Discovery Demo

This is a no-GPU Gradio browser for exported useful-function lead TSVs. It is
designed for fast review of candidate microbes before wet-lab validation.

Run locally:

```bash
uv run --with 'gradio>=4.44,<6' python spaces/app.py
```

The demo loads `examples/product_leads.sample.tsv` by default and can also show
`examples/safe_leads.sample.tsv`. To point it at a local export:

```bash
LEADS_TSV=data/legacy/product_leads_review.tsv uv run --with 'gradio>=4.44,<6' python spaces/app.py
```

Expected TSV columns match the `mfd export-safe-leads` output:

```text
panel target_key label source target_precision rank genome_id species genus family accession score risk_level biosafety_flags evidence
```
