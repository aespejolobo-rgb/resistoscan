"""Screen sequencing reads or assembled contigs against a :class:`KmerIndex`.

Two modes are provided:

* :func:`screen_reads` - metagenomic / amplicon-free shotgun reads. Every read
  is assigned to the reference gene that shares the most k-mers with it, and
  per-gene read counts, depth and breadth of coverage are accumulated.
* :func:`screen_assembly` - assembled genomes (e.g. clinical isolates). Reports
  which genes are present, their breadth of coverage, an identity estimate and
  their location (contig, coordinates, strand).

Identity is estimated from k-mer containment *C* with the Mash model: if
every base mutates independently with probability *d*, a k-mer survives with
probability (1 - d)^k, so identity ~= C^(1/k)  (Ondov et al., 2016).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import islice
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .index import KmerIndex
from .kmers import canonical_kmers, encode
from .seqio import read_sequences


@dataclass
class ReadScreenParams:
    min_kmer_hits: int = 10        # minimum k-mers a read must share with a gene (~40 bp for k=31)
    min_containment: float = 0.25  # hits / k-mers in the matched span; 0.25 ~ 95.6% identity at k=31
    batch_size: int = 20_000       # reads processed per vectorised batch


@dataclass
class ReadScreenResult:
    sample: str
    total_reads: int
    total_bases: int
    genes: pd.DataFrame
    depth: dict[str, np.ndarray] = field(repr=False)


def _batched(iterable: Iterable, size: int):
    iterator = iter(iterable)
    while batch := list(islice(iterator, size)):
        yield batch


def _iter_reads(paths: Iterable[str | Path]):
    for path in paths:
        for _, seq in read_sequences(path):
            yield seq


def screen_reads(paths: str | Path | Iterable[str | Path], index: KmerIndex,
                 params: ReadScreenParams | None = None, sample: str | None = None) -> ReadScreenResult:
    """Assign reads to reference genes and accumulate per-gene coverage.

    Reads are concatenated into one long string separated by ``N`` so that a
    whole batch is k-merised and looked up with a handful of NumPy calls; the
    ``N`` separators make every k-mer spanning two reads invalid.
    """
    params = params or ReadScreenParams()
    if isinstance(paths, (str, Path)):
        paths = [paths]
    paths = [Path(p) for p in paths]
    sample = sample or paths[0].name.split(".")[0]

    k, n_genes = index.k, index.n_genes
    seg_start = np.concatenate(([0], np.cumsum(index.lengths + 1)))
    depth_diff = np.zeros(int(seg_start[-1]))
    reads_per_gene = np.zeros(n_genes)
    identity_sum = np.zeros(n_genes)
    total_reads = total_bases = 0

    for batch in _batched(_iter_reads(paths), params.batch_size):
        total_reads += len(batch)
        total_bases += sum(map(len, batch))
        starts = np.cumsum([0] + [len(s) + 1 for s in batch[:-1]])
        kmers, kpos, _ = canonical_kmers(encode("N".join(batch)), k)
        qi, gene, gpos, _ = index.lookup(kmers)
        if qi.size == 0:
            continue

        # One vote per (read k-mer, gene), even if the k-mer repeats inside the gene.
        _, first = np.unique(qi * n_genes + gene, return_index=True)
        qi, gene, gpos = qi[first], gene[first], gpos[first]
        read = np.searchsorted(starts, kpos[qi], side="right") - 1
        rpos = kpos[qi] - starts[read]

        # Aggregate hits per (read, gene) pair.
        key = read.astype(np.int64) * n_genes + gene
        order = np.argsort(key, kind="stable")
        key, rpos, gpos = key[order], rpos[order], gpos[order]
        group = np.flatnonzero(np.r_[True, key[1:] != key[:-1]])
        hits = np.diff(np.r_[group, key.size])
        span = np.maximum.reduceat(rpos, group) - np.minimum.reduceat(rpos, group) + 1
        gmin = np.minimum.reduceat(gpos, group)
        gmax = np.maximum.reduceat(gpos, group)
        pair_read, pair_gene = key[group] // n_genes, key[group] % n_genes
        containment = hits / span

        keep = (hits >= params.min_kmer_hits) & (containment >= params.min_containment)
        if not keep.any():
            continue
        pair_read, pair_gene, hits = pair_read[keep], pair_gene[keep], hits[keep]
        containment, gmin, gmax = containment[keep], gmin[keep], gmax[keep]

        # Best gene per read; exact ties split the read equally.
        best = np.zeros(len(batch))
        np.maximum.at(best, pair_read, hits)
        is_best = hits == best[pair_read]
        n_best = np.bincount(pair_read[is_best], minlength=len(batch))
        weight = 1.0 / n_best[pair_read[is_best]]
        g = pair_gene[is_best]

        reads_per_gene += np.bincount(g, weights=weight, minlength=n_genes)
        identity_sum += np.bincount(g, weights=weight * containment[is_best] ** (1 / k), minlength=n_genes)
        np.add.at(depth_diff, seg_start[g] + gmin[is_best], weight)
        np.add.at(depth_diff, seg_start[g] + gmax[is_best] + k, -weight)

    rows, depth = [], {}
    for gi, name in enumerate(index.names):
        length = int(index.lengths[gi])
        profile = np.cumsum(depth_diff[seg_start[gi]:seg_start[gi] + length])
        depth[name] = profile
        n = reads_per_gene[gi]
        rows.append({
            "sample": sample, "gene": name, "length": length, "reads": n,
            "mean_depth": profile.mean(), "breadth": float((profile > 1e-9).mean()),
            "identity": identity_sum[gi] / n if n > 0 else np.nan,
        })
    genes = pd.DataFrame(rows)
    genes["total_reads"] = total_reads
    return ReadScreenResult(sample, total_reads, total_bases, genes, depth)


def screen_assembly(path: str | Path, index: KmerIndex, min_breadth: float = 0.9,
                    min_partial: float = 0.2, sample: str | None = None) -> pd.DataFrame:
    """Find reference genes in an assembly (FASTA of contigs).

    Returns one row per gene with any k-mer support, with ``call`` set to
    ``"present"`` (breadth >= ``min_breadth``) or ``"partial"`` (truncated
    or highly divergent copy; breadth >= ``min_partial``).
    """
    path = Path(path)
    sample = sample or path.name.split(".")[0]
    k = index.k
    n_windows = np.maximum(index.lengths - k + 1, 1)
    window_start = np.concatenate(([0], np.cumsum(n_windows)))
    found = np.zeros(int(window_start[-1]), dtype=bool)
    best: dict[int, dict] = {}  # gene -> best supporting contig

    for contig, seq in read_sequences(path):
        kmers, kpos, kfwd = canonical_kmers(encode(seq), k)
        qi, gene, gpos, gfwd = index.lookup(kmers)
        if qi.size == 0:
            continue
        found[window_start[gene] + gpos] = True
        same_strand = kfwd[qi] == gfwd
        for g in np.unique(gene):
            m = gene == g
            n_hits = int(np.unique(gpos[m]).size)
            if n_hits > best.get(int(g), {}).get("hits", 0):
                cpos = kpos[qi[m]]
                best[int(g)] = {
                    "hits": n_hits, "contig": contig,
                    "start": int(cpos.min()) + 1, "end": int(cpos.max()) + k,
                    "strand": "+" if same_strand[m].mean() >= 0.5 else "-",
                }

    rows = []
    for g, loc in best.items():
        length = int(index.lengths[g])
        hit_windows = np.flatnonzero(found[window_start[g]:window_start[g + 1]])
        covered = np.zeros(length + 1)
        np.add.at(covered, hit_windows, 1)
        np.add.at(covered, np.minimum(hit_windows + k, length), -1)
        breadth = float((np.cumsum(covered[:length]) > 0).mean())
        containment = hit_windows.size / n_windows[g]
        missing = n_windows[g] - hit_windows.size
        if breadth < min_partial:
            continue
        rows.append({
            "sample": sample, "gene": index.names[g], "length": length,
            "breadth": breadth, "identity": containment ** (1 / k),
            # an isolated SNP removes k overlapping k-mers, so missing/k ~ number of SNPs
            "est_snps": int(round(missing / k)) if breadth >= min_breadth else np.nan,
            "contig": loc["contig"], "start": loc["start"], "end": loc["end"], "strand": loc["strand"],
            "call": "present" if breadth >= min_breadth else "partial",
        })
    columns = ["sample", "gene", "length", "breadth", "identity", "est_snps",
               "contig", "start", "end", "strand", "call"]
    return pd.DataFrame(rows, columns=columns).sort_values(["call", "gene"], ascending=[False, True],
                                                           ignore_index=True)
