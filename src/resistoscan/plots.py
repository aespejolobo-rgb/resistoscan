"""Publication-style figures (matplotlib, static PNG).

Colour choices follow a colour-vision-deficiency-validated palette:
categorical hues are assigned in a fixed order per entity, magnitudes use a
single-hue sequential ramp, and identity is never carried by colour alone
(legends, direct labels and the TSV tables back every figure).
"""
from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, LogNorm  # noqa: E402
from matplotlib.transforms import Bbox  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
EMPTY_CELL = "#f0efec"
CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
BLUE_RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]

CLASS_ORDER = ["Beta-lactam", "Sulfonamide", "Tetracycline", "MLS", "Aminoglycoside", "Quinolone", "Other"]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Segoe UI", "Helvetica Neue", "Arial", "DejaVu Sans"],
    "font.size": 9.5,
    "text.color": INK,
    "axes.labelcolor": INK_2,
    "axes.edgecolor": AXIS,
    "axes.linewidth": 0.8,
    "axes.titlesize": 11,
    "axes.titleweight": "semibold",
    "axes.titlelocation": "left",
    "axes.titlepad": 12,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "xtick.labelcolor": INK_2,
    "ytick.labelcolor": INK_2,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "savefig.dpi": 200,
    "legend.frameon": False,
})


def _clean(ax, grid_axis: str | None = None) -> None:
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=0)
    if grid_axis:
        ax.grid(axis=grid_axis, color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)


def _save(fig, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)
    return path


def resistome_heatmap(table: pd.DataFrame, samples: list[str], labels: dict[str, str], path: Path) -> Path:
    """Genes (grouped by antibiotic class) x samples, coloured by copies per 16S (log scale)."""
    genes = table[table["category"].isin(["ARG", "MGE"])]
    order = (genes[["gene", "display_class", "category"]].drop_duplicates()
             .assign(rank=lambda d: d["display_class"].map({c: i for i, c in enumerate(CLASS_ORDER)}).fillna(5.5))
             .assign(rank=lambda d: np.where(d["category"] == "MGE", 99, d["rank"]))
             .sort_values(["rank", "gene"], kind="stable"))
    matrix = genes.pivot_table(index="gene", columns="sample", values="copies_per_16S", aggfunc="sum")
    matrix = matrix.reindex(index=order["gene"], columns=samples)
    values = matrix.to_numpy()

    fig, ax = plt.subplots(figsize=(7.4, 7.6))
    cmap = LinearSegmentedColormap.from_list("blue", BLUE_RAMP).with_extremes(bad=EMPTY_CELL)
    shown = np.ma.masked_less_equal(values, 0)
    mesh = ax.pcolormesh(shown, cmap=cmap, norm=LogNorm(vmin=3e-3, vmax=1.0), edgecolors=SURFACE, linewidth=2)

    ax.set_xticks(np.arange(len(samples)) + 0.5, [textwrap.fill(labels.get(s, s), 11) for s in samples])
    ax.xaxis.tick_top()
    ax.set_yticks(np.arange(len(matrix)) + 0.5, matrix.index, fontstyle="italic")
    ax.invert_yaxis()
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)

    # class brackets on the right-hand side
    classes = order.assign(label=np.where(order["category"] == "MGE", "Integron (MGE)", order["display_class"]))["label"]
    pos = 0
    for label, size in classes.groupby(classes.ne(classes.shift()).cumsum()).agg(["first", "size"]).itertuples(index=False):
        ax.plot([len(samples) + 0.25] * 2, [pos + 0.15, pos + size - 0.15], color=AXIS, lw=1.2, clip_on=False)
        ax.text(len(samples) + 0.45, pos + size / 2, label, va="center", ha="left", color=INK_2, fontsize=9)
        pos += size

    cbar = fig.colorbar(mesh, ax=ax, orientation="horizontal", fraction=0.035, pad=0.03, aspect=40)
    cbar.set_ticks([3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1], labels=["0.003", "0.01", "0.03", "0.1", "0.3", "1"])
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(length=0, labelsize=8.5)
    cbar.set_label("Copies per 16S rRNA gene (log scale)   ·   grey = not detected", color=INK_2, fontsize=9)
    ax.set_title("Resistome profile across the One Health gradient", pad=40)
    return _save(fig, path)


def class_profile(table: pd.DataFrame, samples: list[str], labels: dict[str, str], path: Path) -> Path:
    """Left: total ARG load per sample. Right: composition by antibiotic class (100% stacked)."""
    args = table[table["category"] == "ARG"].copy()
    args["cls"] = np.where(args["display_class"].isin(CLASS_ORDER[:-1]), args["display_class"], "Other")
    totals = args.pivot_table(index="sample", columns="cls", values="copies_per_16S", aggfunc="sum")
    totals = totals.reindex(index=samples, columns=CLASS_ORDER).fillna(0)
    load = totals.sum(axis=1)
    share = totals.div(load.replace(0, np.nan), axis=0).fillna(0)

    fig, (ax_load, ax_mix) = plt.subplots(1, 2, figsize=(10, 3.9), gridspec_kw={"width_ratios": [1, 1.6], "wspace": 0.08})
    y = np.arange(len(samples))
    height = 0.42
    ax_load.barh(y, load.values, height=height, color=CATEGORICAL[0])
    for yi, v in zip(y, load.values):
        ax_load.text(v + load.max() * 0.02, yi, f"{v:.2f}", va="center", color=INK_2, fontsize=8.5)
    ax_load.set_yticks(y, [labels.get(s, s) for s in samples])
    ax_load.invert_yaxis()
    ax_load.set_xlim(0, load.max() * 1.18)
    ax_load.set_xlabel("ARG copies per 16S rRNA gene")
    ax_load.set_title("Total resistance gene load", pad=48)
    _clean(ax_load, "x")

    left = np.zeros(len(samples))
    for i, cls in enumerate(CLASS_ORDER):
        ax_mix.barh(y, share[cls].values, left=left, height=height, color=CATEGORICAL[i],
                    edgecolor=SURFACE, linewidth=1.5, label=cls)
        left += share[cls].values
    ax_mix.set_yticks(y, [""] * len(samples))
    ax_mix.invert_yaxis()
    ax_mix.set_xlim(0, 1)
    ax_mix.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
    ax_mix.set_xlabel("Share of ARG abundance")
    ax_mix.set_title("Composition by antibiotic class", pad=48)
    ax_mix.legend(ncol=4, loc="lower left", bbox_to_anchor=(0, 1.01), fontsize=8.5, handlelength=1, handleheight=1,
                  columnspacing=1.2, labelcolor=INK_2, borderaxespad=0)
    _clean(ax_mix, "x")
    return _save(fig, path)


def pcoa_plot(coords: pd.DataFrame, explained: np.ndarray, groups: dict[str, str], labels: dict[str, str],
              path: Path) -> Path:
    """Bray-Curtis PCoA, one colour per environmental compartment, every point labelled."""
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    group_names = list(dict.fromkeys(groups[s] for s in coords.index))
    for i, group in enumerate(group_names):
        members = [s for s in coords.index if groups[s] == group]
        ax.scatter(coords.loc[members, "PCo1"], coords.loc[members, "PCo2"], s=90, color=CATEGORICAL[i],
                   edgecolors=SURFACE, linewidths=2, label=group, zorder=3)
    ax.margins(x=0.3, y=0.25)
    ax.axhline(0, color=AXIS, lw=0.8, zorder=1)
    ax.axvline(0, color=AXIS, lw=0.8, zorder=1)
    ax.set_xlabel(f"PCo1 ({explained[0]:.0%} of variation)")
    ax.set_ylabel(f"PCo2 ({explained[1]:.0%} of variation)")
    ax.set_title("Resistome dissimilarity (Bray-Curtis PCoA)", pad=30)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=len(group_names), fontsize=8.5, labelcolor=INK_2,
              handletextpad=0.3, columnspacing=1.2, borderaxespad=0)
    _clean(ax, None)
    ax.grid(color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    _label_points(ax, coords["PCo1"].to_numpy(), coords["PCo2"].to_numpy(),
                  [labels.get(s, s) for s in coords.index])
    return _save(fig, path)


def _label_points(ax, xs, ys, texts, marker_radius_px: float = 9) -> None:
    """Direct-label scatter points without collisions.

    Tries candidate offsets around each point (near ones first, then farther
    ones with a thin leader line) and keeps the first label box that overlaps
    neither another label nor a marker; falls back to the least-overlapping.
    """
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    centers = ax.transData.transform(np.column_stack([xs, ys]))
    markers = [Bbox.from_bounds(cx - marker_radius_px, cy - marker_radius_px, 2 * marker_radius_px,
                                2 * marker_radius_px) for cx, cy in centers]
    near = [(8, 6), (8, -6), (-8, 6), (-8, -6), (0, 11), (0, -11)]
    far = [(sign * dx, dy) for sign in (1, -1) for dx, dy in ((22, 22), (22, -22), (34, 0), (26, 34), (26, -34))]
    placed = []
    for i, (x, y, text) in enumerate(zip(xs, ys, texts)):
        others = [m for j, m in enumerate(markers) if j != i]
        best, best_score = None, None
        for dx, dy in near + far:
            leader = (dx, dy) in far
            ann = ax.annotate(
                text, (x, y), xytext=(dx, dy), textcoords="offset points", color=INK_2, fontsize=8.5,
                ha="left" if dx > 0 else "right" if dx < 0 else "center", va="bottom" if dy >= 0 else "top",
                arrowprops=dict(arrowstyle="-", color=AXIS, lw=0.8, shrinkA=1, shrinkB=7) if leader else None)
            box = ann.get_window_extent(renderer).expanded(1.04, 1.2)
            inside = ax.bbox.contains(box.x0, box.y0) and ax.bbox.contains(box.x1, box.y1)
            # label collisions are worst, covering a marker next, spilling past the axes is tolerable
            score = (10 * sum(box.overlaps(b) for b in placed) + 5 * sum(box.overlaps(m) for m in others)
                     + (0 if inside else 1))
            if best_score is None or score < best_score:
                if best is not None:
                    best[0].remove()
                best, best_score = (ann, box), score
                if score == 0:
                    break
            else:
                ann.remove()
        placed.append(best[1])


def validation_plot(merged: pd.DataFrame, path: Path) -> Path:
    """Estimated vs simulated abundance for genes that are present and detected."""
    hits = merged[(merged["true_copies_per_16S"] > 0) & merged["detected"]]
    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    lims = (2e-3, 2)
    ax.plot(lims, lims, color=AXIS, lw=1, zorder=1)
    ax.scatter(hits["true_copies_per_16S"], hits["copies_per_16S"], s=46, color=CATEGORICAL[0],
               edgecolors=SURFACE, linewidths=1.5, zorder=3)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_aspect("equal")
    ticks = [0.01, 0.1, 1]
    ax.set_xticks(ticks, ["0.01", "0.1", "1"])
    ax.set_yticks(ticks, ["0.01", "0.1", "1"])
    ax.set_xlabel("Simulated copies per 16S (ground truth)")
    ax.set_ylabel("Estimated copies per 16S")
    ax.set_title("Quantification accuracy on simulated data", pad=24)
    ax.text(0, 1.02, "One dot per detected gene and sample · grey line marks perfect agreement",
            transform=ax.transAxes, color=INK_2, fontsize=8.5, va="bottom")
    _clean(ax, "both")
    return _save(fig, path)

