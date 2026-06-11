"""Application panels for microbial function discovery."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Panel:
    """A model output panel with targets and evidence databases."""

    key: str
    name: str
    description: str
    example_targets: tuple[str, ...]
    evidence_sources: tuple[str, ...]


PANELS: tuple[Panel, ...] = (
    Panel(
        key="environmental_terraforming",
        name="Environmental / Terraforming",
        description="Functions relevant to environmental engineering and extreme-condition survival.",
        example_targets=(
            "carbon fixation",
            "nitrogen cycling",
            "sulfur metabolism",
            "metal reduction",
            "plastic degradation",
            "desiccation tolerance",
            "radiation tolerance",
            "salinity tolerance",
        ),
        evidence_sources=("KEGG", "MetaCyc", "eggNOG", "COG", "GTDB", "MGnify"),
    ),
    Panel(
        key="therapeutics_antimicrobials",
        name="Therapeutics / Antimicrobials",
        description="Functions relevant to drug discovery, pathogen suppression, and microbiome effects.",
        example_targets=(
            "biosynthetic gene cluster potential",
            "antimicrobial peptide potential",
            "pathogen suppression",
            "microbiome-relevant functions",
            "toxin risk",
            "virulence risk",
        ),
        evidence_sources=("antiSMASH", "MIBiG", "BAGEL", "VFDB", "CARD", "AMRFinderPlus"),
    ),
    Panel(
        key="biofuels_industrial",
        name="Biofuels / Industrial Enzymes",
        description="Functions relevant to biomass conversion, fermentation, and robust enzymes.",
        example_targets=(
            "cellulose degradation",
            "xylan degradation",
            "lignin degradation",
            "lipid accumulation",
            "fermentation pathways",
            "thermostable enzymes",
            "acid-stable enzymes",
            "salt-stable enzymes",
        ),
        evidence_sources=("CAZy", "BRENDA", "KEGG", "MetaCyc", "Pfam", "UniProt"),
    ),
    Panel(
        key="food_fermentation_agriculture",
        name="Food / Fermentation / Agriculture",
        description="Functions relevant to fermentation, food safety, crops, and soil productivity.",
        example_targets=(
            "fermentation traits",
            "flavor metabolism",
            "aroma metabolism",
            "probiotic-relevant functions",
            "plant growth promotion",
            "nitrogen fixation",
            "phosphate solubilization",
            "spoilage risk",
        ),
        evidence_sources=("KEGG", "MetaCyc", "eggNOG", "BacDive", "MediaDive", "VFDB"),
    ),
    Panel(
        key="biosafety",
        name="Biosafety",
        description="Risk predictions used to filter discovery candidates before validation.",
        example_targets=(
            "pathogenicity",
            "virulence factors",
            "antimicrobial resistance",
            "toxin genes",
            "biosafety risk",
            "generalization caveats",
        ),
        evidence_sources=("VFDB", "CARD", "AMRFinderPlus", "BacDive", "NCBI"),
    ),
)

APPLICATION_AREAS: tuple[str, ...] = tuple(panel.key for panel in PANELS)
_PANEL_BY_KEY = {panel.key: panel for panel in PANELS}


def list_panels() -> tuple[Panel, ...]:
    """Return all application panels in display order."""

    return PANELS


def get_panel(key: str) -> Panel:
    """Return an application panel by key."""

    return _PANEL_BY_KEY[key]
