"""Markdown run report (renders directly on GitHub)."""
from __future__ import annotations

import datetime as dt
from pathlib import Path

import pandas as pd

from . import __version__
from .pipeline import IsolateResults, MetagenomeResults

CRITICAL_PATTERN = r"Carbapenemase|ESBL|MRSA|VRE|colistin"


def md_table(df: pd.DataFrame, floatfmt: str = "{:.3g}") -> str:
    """Render a DataFrame as a GitHub-flavoured Markdown table (no extra dependency)."""
    def fmt(value) -> str:
        if isinstance(value, float):
            return "–" if pd.isna(value) else floatfmt.format(value)
        return str(value).replace("|", "\\|")

    header = "| " + " | ".join(map(str, df.columns)) + " |"
    rule = "| " + " | ".join("---" for _ in df.columns) + " |"
    rows = ["| " + " | ".join(fmt(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join([header, rule, *rows])


def write_report(path: str | Path, metadata: pd.DataFrame, settings: dict,
                 metagenomes: MetagenomeResults | None = None, isolates: IsolateResults | None = None,
                 benchmark: dict[str, dict] | None = None) -> Path:
    path = Path(path)
    out = [
        "# resistoscan report",
        "",
        f"Generated {dt.datetime.now():%Y-%m-%d %H:%M} with resistoscan {__version__}. "
        + " · ".join(f"{k}: `{v}`" for k, v in settings.items()),
        "",
        "> Research use only. Gene detection does not replace phenotypic antimicrobial susceptibility testing.",
    ]

    if metagenomes is not None:
        table = metagenomes.table
        labels = dict(zip(metagenomes.sheet["sample"], metagenomes.sheet["description"]))
        out += ["", "## Environmental metagenomes", ""]
        for name, fig in metagenomes.figures.items():
            out.append(f"![{name}]({fig.relative_to(path.parent).as_posix()})")
            out.append("")

        alpha = metagenomes.alpha.reset_index()
        alpha.insert(1, "description", alpha["sample"].map(labels))
        alpha.columns = ["Sample", "Description", "ARG richness", "Shannon H'", "Gini-Simpson", "Pielou J",
                         "ARG copies / 16S"]
        out += ["### Resistome diversity", "", md_table(alpha), ""]

        critical_genes = metadata.loc[metadata["note"].str.contains(CRITICAL_PATTERN, regex=True), "gene"] \
            if "note" in metadata else pd.Series(dtype=str)
        critical = table[table["gene"].isin(critical_genes)].pivot_table(
            index="gene", columns="sample", values="copies_per_16S", aggfunc="sum", sort=False)
        critical = critical.reindex(columns=metagenomes.sheet["sample"])
        critical = critical[(critical > 0).any(axis=1)]
        if not critical.empty:
            critical = critical.where(critical > 0).reset_index().rename(columns={"gene": "Gene"})
            out += ["### Clinically critical genes in environmental samples (copies / 16S; – = not detected)", "",
                    md_table(critical), ""]

    if isolates is not None:
        summary = isolates.summary.copy()
        summary.columns = [c.replace("_", " ").capitalize() for c in summary.columns]
        hits = isolates.hits[isolates.hits["gene"] != "16S_rRNA"].copy()
        hits = hits[["sample", "gene", "call", "breadth", "identity", "est_snps", "contig", "start", "end", "strand"]]
        out += ["", "## Clinical isolates", "", "### Predicted resistance profile", "", md_table(summary), "",
                "### Gene hits", "", md_table(hits, floatfmt="{:.4g}"), ""]

    if benchmark:
        out += ["", "## Validation against simulated ground truth", ""]
        for title, metrics in benchmark.items():
            rows = pd.DataFrame({"Metric": list(metrics), "Value": list(metrics.values())})
            out += [f"### {title}", "", md_table(rows), ""]
        validation = path.parent / "figures" / "validation.png"
        if validation.exists():
            out += [f"![validation]({validation.relative_to(path.parent).as_posix()})", ""]

    path.write_text("\n".join(out) + "\n", encoding="utf-8")
    return path
