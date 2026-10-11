import numpy as np
import pytest

from resistoscan.database import DEFAULT_FASTA
from resistoscan.index import KmerIndex
from resistoscan.seqio import read_fasta


@pytest.fixture(scope="session")
def references() -> dict[str, str]:
    return dict(read_fasta(DEFAULT_FASTA))


@pytest.fixture(scope="session")
def panel_index() -> KmerIndex:
    return KmerIndex.from_fasta(DEFAULT_FASTA)


@pytest.fixture
def rng() -> np.random.Generator:
    return np.random.default_rng(7)
