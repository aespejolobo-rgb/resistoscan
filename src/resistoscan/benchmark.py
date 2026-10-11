"""Validate pipeline output against the simulated ground truth."""
from __future__ import annotations

import numpy as np
import pandas as pd


def benchmark_metagenomes(table: pd.DataFrame, truth: pd.DataFrame,
                          detection_floor: float = 0.02) -> tuple[pd.DataFrame, dict]:
    """Detection and quantification accuracy for every (sample, gene) pair.

    ``detection_floor`` is the abundance above which a gene is expected to be
    detectable at the simulated sequencing depth (~3x coverage of the gene).
    """
    estimated = table[table["category"].isin(["ARG", "MGE"])][["sample", "gene", "copies_per_16S", "detected"]]
    targets = truth[truth["role"] == "target"][["sample", "gene", "true_copies_per_16S"]]
    merged = estimated.merge(targets, on=["sample", "gene"], how="left").fillna({"true_copies_per_16S": 0.0})

    present = merged["true_copies_per_16S"] > 0
    called = merged["detected"]
    tp, fn = int((present & called).sum()), int((present & ~called).sum())
    fp, tn = int((~present & called).sum()), int((~present & ~called).sum())
    above = merged["true_copies_per_16S"] >= detection_floor

    both = present & called
    log_true = np.log10(merged.loc[both, "true_copies_per_16S"])
    log_est = np.log10(merged.loc[both, "copies_per_16S"])

    decoys = truth[truth["role"] == "decoy"].merge(estimated, on=["sample", "gene"])
    summary = {
        "Sample-gene pairs evaluated": len(merged),
        "True positives / false negatives": f"{tp} / {fn}",
        "False positives / true negatives": f"{fp} / {tn}",
        "Sensitivity (all present genes)": tp / max(tp + fn, 1),
        f"Sensitivity (genes >= {detection_floor} copies/16S)": (above & called).sum() / max(above.sum(), 1),
        "Specificity": tn / max(tn + fp, 1),
        "Divergent decoy homologs wrongly called": f"{int(decoys['detected'].sum())} / {len(decoys)}",
        "Pearson r, log10 abundance (detected genes)": float(np.corrcoef(log_true, log_est)[0, 1]),
        "Median fold error (detected genes)": float(10 ** np.median(np.abs(log_est - log_true))),
    }
    return merged, summary


def benchmark_isolates(hits: pd.DataFrame, truth: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Gene calls, locations, SNP estimates and truncation flags for assemblies."""
    calls = hits.rename(columns={"sample": "isolate"})
    expected = truth[truth["expected_call"] == "present"]
    present_calls = calls[calls["call"] == "present"]
    merged = expected.merge(present_calls, on=["isolate", "gene"], how="outer",
                            suffixes=("_true", "_call"), indicator=True)
    merged = merged[~merged["gene"].isin(["16S_rRNA"])]
    tp = merged["_merge"] == "both"
    location_ok = tp & (merged["contig_true"] == merged["contig_call"]) & (merged["start_true"] == merged["start_call"]) \
        & (merged["end_true"] == merged["end_call"]) & (merged["strand_true"] == merged["strand_call"])
    snps_ok = tp & (merged["snps"] == merged["est_snps"])

    truncated = truth[truth["expected_call"] == "partial"].merge(calls, on=["isolate", "gene"], how="left")
    n_tp, n_fn = int(tp.sum()), int((merged["_merge"] == "left_only").sum())
    n_fp = int((merged["_merge"] == "right_only").sum())
    summary = {
        "Expected resistance genes": len(expected),
        "Precision": n_tp / max(n_tp + n_fp, 1),
        "Recall": n_tp / max(n_tp + n_fn, 1),
        "Exact location (contig, start, end, strand)": f"{int(location_ok.sum())} / {n_tp}",
        "SNP count estimated exactly": f"{int(snps_ok.sum())} / {n_tp}",
        "Truncated genes flagged as partial": f"{int((truncated['call'] == 'partial').sum())} / {len(truncated)}",
    }
    return merged, summary
