import numpy as np

from resistoscan.index import KmerIndex
from resistoscan.kmers import canonical_kmers, encode


def test_lookup_returns_gene_and_position(references):
    index = KmerIndex.from_records([("sul1", references["sul1"]), ("tet(M)", references["tet(M)"])], k=21)
    query, _, _ = canonical_kmers(encode(references["tet(M)"][100:121]), 21)
    qi, gene, pos, _ = index.lookup(query)
    assert qi.tolist() == [0]
    assert index.names[gene[0]] == "tet(M)"
    assert pos[0] == 100


def test_shared_kmers_report_every_gene():
    shared = "ACGTACGGTCAGTTGACCATGCAAGT"
    index = KmerIndex.from_records([("a", "TTTTT" + shared), ("b", shared + "GGGGG")], k=21)
    query, _, _ = canonical_kmers(encode(shared[:21]), 21)
    _, gene, _, _ = index.lookup(query)
    assert sorted(index.names[g] for g in gene) == ["a", "b"]


def test_unknown_kmers_return_nothing(panel_index):
    qi, gene, pos, fwd = panel_index.lookup(np.array([0], dtype=np.uint64))  # poly-A 31-mer
    assert qi.size == gene.size == pos.size == fwd.size == 0


def test_panel_genes_are_fully_specific(panel_index):
    # the curated panel has no 31-mers shared between genes, so read assignment is unambiguous
    assert np.all(panel_index.gene_specificity() == 1.0)


def test_save_and_load_roundtrip(panel_index, tmp_path):
    path = tmp_path / "panel.npz"
    panel_index.save(path)
    loaded = KmerIndex.load(path)
    assert loaded.names == panel_index.names
    assert loaded.k == panel_index.k
    assert np.array_equal(loaded.kmers, panel_index.kmers)
