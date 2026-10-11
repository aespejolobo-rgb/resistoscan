import pandas as pd
import pytest

from resistoscan.abundance import abundance_matrix, normalise
from resistoscan.database import load_metadata


def test_copies_per_16s_and_detection():
    genes = pd.DataFrame({
        "sample": ["s"] * 3,
        "gene": ["16S_rRNA", "sul1", "blaKPC-2"],
        "length": [1500, 750, 900],
        "reads": [300, 30, 1],
        "breadth": [1.0, 0.9, 0.15],
        "total_reads": [1_000_000] * 3,
    })
    table = normalise(genes, load_metadata()).set_index("gene")
    # (30 / 750) / (300 / 1500) = 0.2 copies per 16S
    assert table.loc["sul1", "copies_per_16S"] == pytest.approx(0.2)
    assert table.loc["sul1", "rpkm"] == pytest.approx(40)
    assert not table.loc["blaKPC-2", "detected"]
    assert table.loc["blaKPC-2", "copies_per_16S"] == 0
    assert table.loc["sul1", "display_class"] == "Sulfonamide"

    matrix = abundance_matrix(table.reset_index())
    assert list(matrix.columns) == ["sul1", "blaKPC-2"]  # 16S marker is excluded
