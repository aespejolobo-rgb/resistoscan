"""Access to the bundled reference panel and its annotation table."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

DATA_DIR = Path(__file__).parent / "data"
DEFAULT_FASTA = DATA_DIR / "amr_panel.fasta"
DEFAULT_METADATA = DATA_DIR / "amr_panel.tsv"
MARKER_GENE = "16S_rRNA"

REQUIRED_COLUMNS = {"gene", "category", "display_class"}


def load_metadata(path: str | Path = DEFAULT_METADATA) -> pd.DataFrame:
    """Load the gene annotation table (one row per reference sequence).

    Required columns: ``gene`` (FASTA identifier), ``category`` (ARG, MGE or
    marker) and ``display_class`` (antibiotic class used for grouping).
    """
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = REQUIRED_COLUMNS - set(table.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {sorted(missing)}")
    return table
