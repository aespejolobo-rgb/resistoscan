"""Minimal FASTA/FASTQ readers and writers (plain text or gzip).

Kept dependency-free on purpose: the screening engine only needs sequence
strings, and avoiding Biopython keeps installation light for teaching and CI.
"""
from __future__ import annotations

import gzip
from pathlib import Path
from typing import IO, Iterable, Iterator


def _open(path: str | Path, mode: str = "rt") -> IO:
    path = Path(path)
    if path.suffix == ".gz":
        # level 3 compresses sequence data almost as well as 9 at a fraction of the time
        return gzip.open(path, mode, compresslevel=3, encoding="utf-8")
    return open(path, mode, encoding="utf-8")


def read_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield ``(identifier, sequence)`` pairs. The identifier is the first header token."""
    name, chunks = None, []
    with _open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def read_fastq(path: str | Path) -> Iterator[tuple[str, str, str]]:
    """Yield ``(identifier, sequence, quality)`` triples from a 4-line FASTQ file."""
    with _open(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            seq = fh.readline().strip()
            fh.readline()
            qual = fh.readline().strip()
            if not header.startswith("@") or len(seq) != len(qual):
                raise ValueError(f"Malformed FASTQ record near '{header.strip()}' in {path}")
            yield header[1:].split()[0], seq.upper(), qual


def read_sequences(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield ``(identifier, sequence)`` from FASTA or FASTQ, detected from the first byte."""
    with _open(path) as fh:
        first = fh.read(1)
    if first == ">":
        yield from read_fasta(path)
    elif first == "@":
        for name, seq, _ in read_fastq(path):
            yield name, seq
    elif first == "":
        return
    else:
        raise ValueError(f"{path} does not look like FASTA or FASTQ")


def write_fasta(records: Iterable[tuple[str, str]], path: str | Path, width: int = 80) -> None:
    with _open(path, "wt") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def write_fastq(records: Iterable[tuple[str, str, str]], path: str | Path) -> None:
    with _open(path, "wt") as fh:
        fh.writelines(f"@{name}\n{seq}\n+\n{qual}\n" for name, seq, qual in records)


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]
