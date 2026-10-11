"""High-level workflows used by the command-line interface."""
from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from . import plots
from .abundance import abundance_matrix, normalise
from .diversity import alpha_diversity, bray_curtis, pcoa
from .index import KmerIndex
from .phenotype import interpret
from .screen import ReadScreenParams, screen_assembly, screen_reads

log = logging.getLogger("resistoscan")


@dataclass
class MetagenomeResults:
    table: pd.DataFrame
    alpha: pd.DataFrame
    distances: pd.DataFrame
    coords: pd.DataFrame
    explained: np.ndarray
    sheet: pd.DataFrame
    figures: dict[str, Path] = field(default_factory=dict)


@dataclass
class IsolateResults:
    hits: pd.DataFrame
    summary: pd.DataFrame
    sheet: pd.DataFrame


def run_metagenomes(sample_sheet: str | Path, outdir: str | Path, index: KmerIndex, metadata: pd.DataFrame,
                    params: ReadScreenParams | None = None, min_breadth: float = 0.5) -> MetagenomeResults:
    """Screen every sample in a sheet (columns: sample, group, reads[, description])."""
    sample_sheet, outdir = Path(sample_sheet), Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    sheet = pd.read_csv(sample_sheet, sep="\t", dtype=str)
    missing = {"sample", "reads"} - set(sheet.columns)
    if missing:
        raise ValueError(f"Sample sheet is missing columns: {sorted(missing)}")
    sheet["group"] = sheet.get("group", pd.Series("NA", index=sheet.index)).fillna("NA")
    sheet["description"] = sheet.get("description", sheet["sample"]).fillna(sheet["sample"])

    results = []
    for row in sheet.itertuples(index=False):
        paths = [sample_sheet.parent / p.strip() for p in row.reads.split(",")]
        start = time.perf_counter()
        result = screen_reads(paths, index, params, sample=row.sample)
        log.info("%-12s %9d reads screened in %.1fs", row.sample, result.total_reads, time.perf_counter() - start)
        results.append(result)

    table = normalise(pd.concat([r.genes for r in results], ignore_index=True), metadata, min_breadth=min_breadth)
    matrix = abundance_matrix(table).reindex(sheet["sample"]).fillna(0.0)
    alpha = alpha_diversity(matrix)
    distances = bray_curtis(matrix)
    coords, explained = pcoa(distances)

    table.to_csv(outdir / "gene_abundance.tsv", sep="\t", index=False, float_format="%.6g")
    matrix.to_csv(outdir / "abundance_matrix.tsv", sep="\t", float_format="%.6g")
    alpha.to_csv(outdir / "alpha_diversity.tsv", sep="\t", float_format="%.4f")
    distances.to_csv(outdir / "bray_curtis.tsv", sep="\t", float_format="%.4f")
    coords.to_csv(outdir / "pcoa.tsv", sep="\t", float_format="%.4f")

    samples = sheet["sample"].tolist()
    labels = dict(zip(sheet["sample"], sheet["description"]))
    groups = dict(zip(sheet["sample"], sheet["group"]))
    figdir = outdir / "figures"
    figures = {
        "heatmap": plots.resistome_heatmap(table, samples, labels, figdir / "resistome_heatmap.png"),
        "classes": plots.class_profile(table, samples, labels, figdir / "class_profile.png"),
    }
    if len(samples) >= 3:
        figures["pcoa"] = plots.pcoa_plot(coords, explained, groups, labels, figdir / "pcoa_bray_curtis.png")
    return MetagenomeResults(table, alpha, distances, coords, explained, sheet, figures)


def run_isolates(sample_sheet: str | Path, outdir: str | Path, index: KmerIndex, metadata: pd.DataFrame,
                 min_breadth: float = 0.9) -> IsolateResults:
    """Screen every assembly in a sheet (columns: isolate, assembly[, organism])."""
    sample_sheet, outdir = Path(sample_sheet), Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    sheet = pd.read_csv(sample_sheet, sep="\t", dtype=str)
    missing = {"isolate", "assembly"} - set(sheet.columns)
    if missing:
        raise ValueError(f"Isolate sheet is missing columns: {sorted(missing)}")

    hits = []
    for row in sheet.itertuples(index=False):
        found = screen_assembly(sample_sheet.parent / row.assembly, index, min_breadth=min_breadth, sample=row.isolate)
        log.info("%-12s %d genes found", row.isolate, (found["call"] == "present").sum())
        hits.append(found)
    hits = pd.concat(hits, ignore_index=True)
    summary = interpret(hits, metadata, samples=sheet["isolate"].tolist())
    if "organism" in sheet.columns:
        summary.insert(1, "organism", summary["isolate"].map(dict(zip(sheet["isolate"], sheet["organism"]))))

    hits.to_csv(outdir / "isolate_genes.tsv", sep="\t", index=False, float_format="%.4f")
    summary.to_csv(outdir / "isolate_summary.tsv", sep="\t", index=False)
    return IsolateResults(hits, summary, sheet)
