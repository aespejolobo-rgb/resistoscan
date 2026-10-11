import numpy as np
import pytest

from resistoscan.kmers import canonical_kmers, decode, encode
from resistoscan.seqio import reverse_complement


def test_encode_maps_bases_and_flags_ambiguous():
    assert encode("ACGTacgtN").tolist() == [0, 1, 2, 3, 0, 1, 2, 3, 4]


def test_forward_kmers_decode_back_to_sequence():
    seq = "ACGTTGCAAGGCT"
    kmers, pos, forward = canonical_kmers(encode(seq), 5)
    assert pos.tolist() == list(range(len(seq) - 4))
    for value, p, fwd in zip(kmers, pos, forward):
        window = seq[p:p + 5]
        expected = window if fwd else reverse_complement(window)
        assert decode(value, 5) == expected
        assert decode(value, 5) == min(window, reverse_complement(window))  # A<C<G<T matches 2-bit order


def test_canonical_kmers_are_strand_independent(rng):
    seq = "".join(rng.choice(list("ACGT"), 500))
    fwd, _, _ = canonical_kmers(encode(seq), 31)
    rev, _, _ = canonical_kmers(encode(reverse_complement(seq)), 31)
    assert np.array_equal(np.sort(fwd), np.sort(rev))


def test_windows_with_n_are_skipped():
    kmers, pos, _ = canonical_kmers(encode("AAAAANCCCCC"), 3)
    assert pos.tolist() == [0, 1, 2, 6, 7, 8]
    assert kmers.size == 6


def test_short_sequence_gives_no_kmers():
    kmers, pos, fwd = canonical_kmers(encode("ACGT"), 31)
    assert kmers.size == pos.size == fwd.size == 0


@pytest.mark.parametrize("k", [0, 32])
def test_invalid_k_raises(k):
    with pytest.raises(ValueError):
        canonical_kmers(encode("ACGT" * 20), k)
