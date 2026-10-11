"""Alpha and beta diversity of resistomes.

Implemented from first principles with NumPy (no scikit-bio dependency) so the
maths is visible: Shannon and Simpson indices, Pielou evenness, Bray-Curtis
dissimilarity and Principal Coordinates Analysis (classical MDS).
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def shannon(abundances: np.ndarray) -> float:
    """Shannon entropy H' = -sum(p_i ln p_i) over non-zero abundances."""
    x = np.asarray(abundances, float)
    x = x[x > 0]
    if x.size == 0:
        return 0.0
    p = x / x.sum()
    h = float(-(p * np.log(p)).sum())
    return h if h > 0 else 0.0  # a single-gene sample gives -0.0


def simpson(abundances: np.ndarray) -> float:
    """Gini-Simpson index 1 - sum(p_i^2): probability two draws are different genes."""
    x = np.asarray(abundances, float)
    if x.sum() == 0:
        return 0.0
    p = x / x.sum()
    return float(1 - (p ** 2).sum())


def alpha_diversity(matrix: pd.DataFrame) -> pd.DataFrame:
    """Per-sample richness, Shannon, Simpson and Pielou evenness."""
    rows = []
    for sample, values in matrix.iterrows():
        richness = int((values > 0).sum())
        h = shannon(values.to_numpy())
        rows.append({
            "sample": sample, "richness": richness, "shannon": h,
            "simpson": simpson(values.to_numpy()),
            "evenness": h / np.log(richness) if richness > 1 else np.nan,
            "total_abundance": float(values.sum()),
        })
    return pd.DataFrame(rows).set_index("sample")


def bray_curtis(matrix: pd.DataFrame) -> pd.DataFrame:
    """Bray-Curtis dissimilarity: BC(a, b) = sum|a_i - b_i| / sum(a_i + b_i)."""
    x = matrix.to_numpy(float)
    num = np.abs(x[:, None, :] - x[None, :, :]).sum(axis=2)
    den = (x[:, None, :] + x[None, :, :]).sum(axis=2)
    with np.errstate(invalid="ignore", divide="ignore"):
        dist = np.where(den > 0, num / den, 0.0)
    return pd.DataFrame(dist, index=matrix.index, columns=matrix.index)


def pcoa(distances: pd.DataFrame, n_axes: int = 2) -> tuple[pd.DataFrame, np.ndarray]:
    """Principal Coordinates Analysis of a distance matrix.

    Double-centres the squared distances (B = -1/2 J D^2 J), takes the
    eigen-decomposition of B and scales eigenvectors by sqrt(eigenvalue).
    Returns the sample coordinates and the fraction of variance explained by
    each returned axis (relative to the sum of positive eigenvalues).
    """
    d2 = distances.to_numpy(float) ** 2
    n = d2.shape[0]
    centering = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * centering @ d2 @ centering
    eigvals, eigvecs = np.linalg.eigh(b)
    order = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    positive = eigvals > 1e-12
    coords = eigvecs[:, :n_axes] * np.sqrt(np.clip(eigvals[:n_axes], 0, None))
    explained = eigvals[:n_axes] / eigvals[positive].sum()
    columns = [f"PCo{i + 1}" for i in range(n_axes)]
    return pd.DataFrame(coords, index=distances.index, columns=columns), explained
