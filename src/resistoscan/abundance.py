"""Normalisation of per-gene read counts into comparable abundances.

Two complementary units are reported:

* **RPKM** - reads per kilobase of gene per million sequenced reads. Corrects
  for gene length and sequencing depth.
* **Copies per 16S rRNA gene** - the standard unit in environmental resistome
  studies (e.g. ARGs-OAP). It approximates how many copies of an ARG are
  carried per bacterial genome-equivalent, which corrects for differences in
  bacterial load and host DNA between samples::

      copies/16S = (reads_ARG / length_ARG) / (reads_16S / length_16S)
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

from .database import MARKER_GENE


def normalise(genes: pd.DataFrame, metadata: pd.DataFrame, marker: str = MARKER_GENE,
              min_breadth: float = 0.5, min_reads: float = 2) -> pd.DataFrame:
    """Add detection calls and normalised abundances to a long gene table.

    ``genes`` needs the columns produced by :func:`resistoscan.screen.screen_reads`:
    sample, gene, length, reads, breadth, total_reads. A gene is *detected* when
    at least ``min_breadth`` of its length is covered by ``min_reads`` or more
    reads; undetected genes get an abundance of zero, which removes most
    spurious signal from short shared k-mer stretches.
    """
    table = genes.copy()
    table["detected"] = (table["breadth"] >= min_breadth) & (table["reads"] >= min_reads)
    table["rpkm"] = table["reads"] / (table["length"] / 1e3) / (table["total_reads"] / 1e6)

    marker_rows = table[table["gene"] == marker].set_index("sample")
    if marker_rows.empty:
        warnings.warn(f"Marker gene '{marker}' not in the table; copies/16S will be NaN")
        table["copies_per_16S"] = np.nan
    else:
        marker_coverage = (marker_rows["reads"] / marker_rows["length"]).replace(0, np.nan)
        if marker_coverage.isna().any():
            warnings.warn("No 16S reads in samples: " + ", ".join(marker_coverage[marker_coverage.isna()].index))
        table["copies_per_16S"] = (table["reads"] / table["length"]) / table["sample"].map(marker_coverage)

    table.loc[~table["detected"], ["rpkm", "copies_per_16S"]] = 0.0
    annotation = metadata[["gene", "category", "display_class"]]
    return table.merge(annotation, on="gene", how="left")


def abundance_matrix(table: pd.DataFrame, value: str = "copies_per_16S",
                     categories: tuple[str, ...] = ("ARG",)) -> pd.DataFrame:
    """Pivot a normalised long table into a samples x genes matrix."""
    subset = table[table["category"].isin(categories)]
    matrix = subset.pivot_table(index="sample", columns="gene", values=value, aggfunc="sum", sort=False)
    return matrix.fillna(0.0)


def class_totals(table: pd.DataFrame, value: str = "copies_per_16S") -> pd.DataFrame:
    """Sum ARG abundance per antibiotic class (samples x classes)."""
    subset = table[table["category"] == "ARG"]
    return subset.pivot_table(index="sample", columns="display_class", values=value,
                              aggfunc="sum", sort=False).fillna(0.0)
