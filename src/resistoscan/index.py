"""Sorted-array k-mer index of a reference gene database.

Instead of a Python dictionary (hundreds of bytes per k-mer) the index keeps
four flat NumPy arrays sorted by k-mer value. Lookups use binary search
(``np.searchsorted``) on whole batches of query k-mers at once, so memory is
~17 bytes per reference k-mer and queries run at C speed. K-mers shared by
several genes simply appear several times in the sorted array.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .kmers import canonical_kmers, encode
from .seqio import read_fasta


@dataclass
class KmerIndex:
    k: int
    names: list[str]
    lengths: np.ndarray   # (n_genes,) gene length in bp
    kmers: np.ndarray     # (n,) uint64, sorted
    gene: np.ndarray      # (n,) int32, gene index of each k-mer
    pos: np.ndarray       # (n,) int32, start position of the k-mer in its gene
    forward: np.ndarray   # (n,) bool, canonical form equals the gene's forward strand

    @classmethod
    def from_records(cls, records, k: int = 31) -> "KmerIndex":
        names, lengths, parts = [], [], []
        for gene_idx, (name, seq) in enumerate(records):
            kmers, pos, fwd = canonical_kmers(encode(seq), k)
            names.append(name)
            lengths.append(len(seq))
            parts.append((kmers, np.full(kmers.size, gene_idx, np.int32), pos.astype(np.int32), fwd))
        if not names:
            raise ValueError("Reference database is empty")
        kmers, gene, pos, fwd = (np.concatenate(cols) for cols in zip(*parts))
        order = np.argsort(kmers, kind="stable")
        return cls(k, names, np.asarray(lengths, np.int64), kmers[order], gene[order], pos[order], fwd[order])

    @classmethod
    def from_fasta(cls, path: str | Path, k: int = 31) -> "KmerIndex":
        return cls.from_records(read_fasta(path), k=k)

    @property
    def n_genes(self) -> int:
        return len(self.names)

    def gene_id(self, name: str) -> int:
        return self.names.index(name)

    def lookup(self, query: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Find every reference occurrence of each query k-mer.

        Returns aligned arrays ``(query_index, gene, gene_position, gene_forward)``
        with one row per (query k-mer, reference occurrence) pair.
        """
        left = np.searchsorted(self.kmers, query, side="left")
        right = np.searchsorted(self.kmers, query, side="right")
        counts = right - left
        total = int(counts.sum())
        if total == 0:
            empty = np.empty(0, np.int64)
            return empty, empty.astype(np.int32), empty.astype(np.int32), np.empty(0, bool)
        query_idx = np.repeat(np.arange(query.size), counts)
        within = np.arange(total) - np.repeat(np.cumsum(counts) - counts, counts)
        ref = np.repeat(left, counts) + within
        return query_idx, self.gene[ref], self.pos[ref], self.forward[ref]

    def gene_specificity(self) -> np.ndarray:
        """Fraction of each gene's k-mers that are unique to that gene in the database."""
        _, counts = np.unique(self.kmers, return_counts=True)
        shared = np.repeat(counts > 1, counts)  # self.kmers is sorted, so groups are contiguous
        total = np.bincount(self.gene, minlength=self.n_genes)
        unique = np.bincount(self.gene[~shared], minlength=self.n_genes)
        return unique / np.maximum(total, 1)

    def save(self, path: str | Path) -> None:
        np.savez_compressed(path, k=self.k, names=np.array(self.names), lengths=self.lengths,
                            kmers=self.kmers, gene=self.gene, pos=self.pos, forward=self.forward)

    @classmethod
    def load(cls, path: str | Path) -> "KmerIndex":
        data = np.load(path)
        return cls(int(data["k"]), data["names"].tolist(), data["lengths"], data["kmers"],
                   data["gene"], data["pos"], data["forward"])
