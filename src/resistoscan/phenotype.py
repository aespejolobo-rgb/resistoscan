"""Genotype-to-phenotype interpretation for clinical isolates.

Translates the resistance genes detected in an assembly into the antibiotic
classes they are expected to compromise, and raises alerts for the genotypes
behind WHO priority pathogens (carbapenemases, ESBLs, MRSA, VRE, mcr).

This is a research aid. Gene presence does not guarantee expression, and
chromosomal point mutations are not covered, so predictions must be
confirmed with phenotypic antimicrobial susceptibility testing.
"""
from __future__ import annotations

import pandas as pd

ALERTS = [
    # (rule, alert text)
    (lambda r: "CARBAPENEM" in r.ncbi_subclass, "Carbapenemase producer"),
    (lambda r: r.gene.startswith("blaCTX-M"), "ESBL producer"),
    (lambda r: r.gene == "mecA", "MRSA genotype (mecA)"),
    (lambda r: r.gene == "vanA", "VRE genotype (vanA)"),
    (lambda r: r.display_class == "Polymyxin", "Mobile colistin resistance (mcr)"),
]


def interpret(hits: pd.DataFrame, metadata: pd.DataFrame, samples: list[str] | None = None) -> pd.DataFrame:
    """One row per isolate: predicted resistant classes, key genes and alerts."""
    metadata = metadata.assign(ncbi_subclass=metadata.get("ncbi_subclass", ""))
    annotated = hits[hits["call"] == "present"].merge(metadata, on="gene", how="left")
    annotated = annotated[annotated["category"] == "ARG"]
    rows = []
    for sample in samples if samples is not None else hits["sample"].unique():
        group = hits[hits["sample"] == sample]
        genes = annotated[annotated["sample"] == sample]
        alerts = []
        for rule, text in ALERTS:
            matched = [r.gene for r in genes.itertuples() if rule(r)]
            if matched:
                alerts.append(f"{text}: {', '.join(matched)}")
        variants = [f"{r.gene} (~{int(r.est_snps)} SNP)" for r in genes.itertuples() if r.est_snps and r.est_snps > 0]
        partial = group.loc[group["call"] == "partial", "gene"].tolist()
        rows.append({
            "isolate": sample,
            "n_args": len(genes),
            "resistance_genes": ", ".join(genes["gene"]),
            "predicted_resistance": ", ".join(sorted(genes["display_class"].unique())),
            "alerts": "; ".join(alerts) or "none",
            "allelic_variants": ", ".join(variants) or "none",
            "partial_hits": ", ".join(partial) or "none",
        })
    return pd.DataFrame(rows)
