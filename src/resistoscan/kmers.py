"""Vectorised 2-bit k-mer encoding with NumPy.

Each nucleotide is packed into 2 bits (A=0, C=1, G=2, T=3), so any k-mer with
k <= 31 fits in a single unsigned 64-bit integer. K-mers are stored in their
*canonical* form, min(forward, reverse complement), which makes matching
independent of the strand a read or contig was sequenced from.
"""
from __future__ import annotations

import numpy as np

MAX_K = 31

_ENCODE = np.full(256, 4, dtype=np.uint8)
for _code, _base in enumerate("ACGT"):
    _ENCODE[ord(_base)] = _code
    _ENCODE[ord(_base.lower())] = _code


def encode(seq: str) -> np.ndarray:
    """Encode a nucleotide string as uint8 codes; anything but ACGT becomes 4 (invalid)."""
    return _ENCODE[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]


def canonical_kmers(codes: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute the canonical k-mers of an encoded sequence.

    Returns three aligned arrays describing every window that contains only
    A/C/G/T: the canonical k-mer value (uint64), its 0-based start position,
    and whether the canonical form is the forward strand.

    The forward and reverse-complement values are built with k vectorised
    shift/OR passes over the whole sequence, which is O(k * n) but runs in C,
    so even multi-megabase contigs take well under a second.
    """
    if not 1 <= k <= MAX_K:
        raise ValueError(f"k must be between 1 and {MAX_K}")
    n = codes.size - k + 1
    if n <= 0:
        return np.empty(0, np.uint64), np.empty(0, np.int64), np.empty(0, bool)

    invalid = codes > 3
    invalid_prefix = np.concatenate(([0], np.cumsum(invalid)))
    window_ok = (invalid_prefix[k:] - invalid_prefix[:-k]) == 0

    bases = np.where(invalid, 0, codes).astype(np.uint64)
    forward = np.zeros(n, dtype=np.uint64)
    reverse = np.zeros(n, dtype=np.uint64)
    two, three = np.uint64(2), np.uint64(3)
    for j in range(k):
        window = bases[j:j + n]
        forward = (forward << two) | window
        reverse |= (three - window) << np.uint64(2 * j)

    canonical = np.minimum(forward, reverse)
    positions = np.flatnonzero(window_ok)
    return canonical[window_ok], positions, (forward <= reverse)[window_ok]


def decode(value: int, k: int) -> str:
    """Turn a packed k-mer back into a string (handy for debugging and tests)."""
    return "".join("ACGT"[(int(value) >> (2 * (k - 1 - i))) & 3] for i in range(k))
