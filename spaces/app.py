"""Gradio demo for reviewing microbial useful-function discovery leads."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import gradio as gr

try:
    from spaces.lead_browser import (
        filter_leads,
        lead_detail,
        load_leads,
        option_values,
        summary_markdown,
        table_columns,
        table_rows,
    )
except ModuleNotFoundError:
    from lead_browser import (
        filter_leads,
        lead_detail,
        load_leads,
        option_values,
        summary_markdown,
        table_columns,
        table_rows,
    )


SPACE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = Path(__file__).resolve().parents[1]


def _example_path(filename: str) -> Path:
    for root in [SPACE_ROOT, REPO_ROOT]:
        path = root / "examples" / filename
        if path.exists():
            return path
    return REPO_ROOT / "examples" / filename


DEFAULT_DATASETS = {
    "Product review leads": _example_path("product_leads.review.tsv"),
    "Strict safe leads": _example_path("safe_leads.review.tsv"),
    "Product review sample": _example_path("product_leads.sample.tsv"),
    "Strict safe sample": _example_path("safe_leads.sample.tsv"),
}

CSS = """
:root {
  --mfd-ink: #17201b;
  --mfd-muted: #5b665f;
  --mfd-line: #d7ddd8;
  --mfd-surface: #f7f9f7;
  --mfd-accent: #1f6f54;
}
.gradio-container {
  color: var(--mfd-ink);
  font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
}
.mfd-title {
  border-bottom: 1px solid var(--mfd-line);
  margin-bottom: 14px;
  padding-bottom: 12px;
}
.mfd-title h1 {
  font-size: 28px;
  line-height: 1.15;
  margin: 0 0 6px;
}
.mfd-title p {
  color: var(--mfd-muted);
  margin: 0;
}
.mfd-summary {
  background: var(--mfd-surface);
  border: 1px solid var(--mfd-line);
  border-radius: 8px;
  padding: 12px 14px;
}
button.primary {
  background: var(--mfd-accent) !important;
  border-color: var(--mfd-accent) !important;
}
"""


def _dataset_map() -> dict[str, Path]:
    datasets = dict(DEFAULT_DATASETS)
    env_path = os.environ.get("LEADS_TSV")
    if env_path:
        datasets["Local LEADS_TSV"] = Path(env_path).expanduser()
    return datasets


def _load_dataset(dataset_name: str) -> list[dict[str, Any]]:
    datasets = _dataset_map()
    path = datasets.get(dataset_name) or next(iter(datasets.values()))
    return load_leads(path)


def _dropdown_update(values: list[str]) -> gr.Dropdown:
    return gr.Dropdown(choices=["All", *values], value="All")


def dataset_changed(dataset_name: str) -> tuple[Any, Any, Any, Any, list[list[Any]], str, str]:
    leads = _load_dataset(dataset_name)
    filtered = filter_leads(leads)
    return (
        _dropdown_update(option_values(leads, "panel")),
        _dropdown_update(option_values(leads, "risk_level")),
        _dropdown_update(option_values(leads, "source")),
        table_rows(filtered),
        summary_markdown(filtered, dataset_name),
        lead_detail(filtered[0] if filtered else None),
        "",
    )


def update_view(
    dataset_name: str,
    panel: str,
    risk: str,
    source: str,
    query: str,
    min_precision: float,
    genome_id: str,
) -> tuple[list[list[Any]], str, str]:
    leads = _load_dataset(dataset_name)
    filtered = filter_leads(
        leads,
        panel=panel,
        risk=risk,
        source=source,
        query=query,
        min_precision=min_precision,
    )
    selected = _find_detail(filtered, genome_id) or (filtered[0] if filtered else None)
    return table_rows(filtered), summary_markdown(filtered, dataset_name), lead_detail(selected)


def _find_detail(leads: list[dict[str, Any]], genome_id: str) -> dict[str, Any] | None:
    lookup = genome_id.strip()
    if not lookup:
        return None
    for lead in leads:
        if str(lead.get("genome_id", "")) == lookup or str(lead.get("accession", "")) == lookup:
            return lead
    return None


def create_app() -> gr.Blocks:
    datasets = _dataset_map()
    default_dataset = "Local LEADS_TSV" if "Local LEADS_TSV" in datasets else "Product review leads"
    default_leads = _load_dataset(default_dataset)
    default_filtered = filter_leads(default_leads)

    with gr.Blocks(css=CSS, title="Microbial Function Discovery") as demo:
        gr.HTML(
            """
            <div class="mfd-title">
              <h1>Microbial Function Discovery</h1>
              <p>Review useful-function leads from genome-scale microbial predictions.</p>
            </div>
            """
        )
        summary = gr.Markdown(
            summary_markdown(default_filtered, default_dataset),
            elem_classes=["mfd-summary"],
        )
        with gr.Row():
            dataset = gr.Dropdown(
                choices=list(datasets.keys()),
                value=default_dataset,
                label="Dataset",
            )
            panel = gr.Dropdown(
                choices=["All", *option_values(default_leads, "panel")],
                value="All",
                label="Panel",
            )
            risk = gr.Dropdown(
                choices=["All", *option_values(default_leads, "risk_level")],
                value="All",
                label="Risk",
            )
            source = gr.Dropdown(
                choices=["All", *option_values(default_leads, "source")],
                value="All",
                label="Source",
            )
        with gr.Row():
            query = gr.Textbox(label="Search", placeholder="species, accession, target, evidence")
            min_precision = gr.Slider(0.0, 1.0, value=0.0, step=0.05, label="Minimum target precision")
            genome_id = gr.Textbox(label="Detail genome/accession", placeholder="e.g. 1268 or GCA_004340465")
        table = gr.Dataframe(
            headers=table_columns(),
            value=table_rows(default_filtered),
            label="Ranked leads",
            wrap=True,
            interactive=False,
        )
        detail = gr.Markdown(lead_detail(default_filtered[0] if default_filtered else None))

        dataset.change(
            dataset_changed,
            inputs=[dataset],
            outputs=[panel, risk, source, table, summary, detail, genome_id],
        )
        for control in (panel, risk, source, query, min_precision, genome_id):
            control.change(
                update_view,
                inputs=[dataset, panel, risk, source, query, min_precision, genome_id],
                outputs=[table, summary, detail],
            )

    return demo


if __name__ == "__main__":
    create_app().launch()
