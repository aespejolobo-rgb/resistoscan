import numpy as np
import pandas as pd
import pytest

from resistoscan.diversity import alpha_diversity, bray_curtis, pcoa, shannon, simpson


def test_shannon_and_simpson_of_even_community():
    assert shannon(np.ones(4)) == pytest.approx(np.log(4))
    assert simpson(np.ones(4)) == pytest.approx(0.75)
    assert shannon(np.array([5, 0, 0])) == 0.0


def test_alpha_diversity_table():
    matrix = pd.DataFrame({"a": [1, 1], "b": [1, 0], "c": [2, 0]}, index=["s1", "s2"])
    alpha = alpha_diversity(matrix)
    assert alpha.loc["s1", "richness"] == 3
    assert alpha.loc["s2", "richness"] == 1
    assert np.isnan(alpha.loc["s2", "evenness"])


def test_bray_curtis_bounds():
    matrix = pd.DataFrame([[1, 2, 0], [1, 2, 0], [0, 0, 5]], index=["a", "b", "c"])
    dist = bray_curtis(matrix)
    assert dist.loc["a", "b"] == 0
    assert dist.loc["a", "c"] == 1
    assert np.allclose(dist, dist.T)


def test_pcoa_recovers_euclidean_configuration():
    points = np.array([[0, 0], [3, 0], [0, 4], [1, 1]], dtype=float)
    dist = pd.DataFrame(np.linalg.norm(points[:, None] - points[None], axis=2))
    coords, explained = pcoa(dist)
    recovered = np.linalg.norm(coords.to_numpy()[:, None] - coords.to_numpy()[None], axis=2)
    assert np.allclose(recovered, dist.to_numpy())
    assert explained.sum() == pytest.approx(1.0)
