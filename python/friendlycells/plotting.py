"""Plots example neighborhoods, to check neighborhood settings."""

import numpy as np
import pandas as pd

from .similarity import _neighbor_pairs, _neighborhood_settings, _nn_distances

# Okabe-Ito, which suits color blindness, without black and grey
_COLORS = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]


def plot_neighborhoods(
    samples,
    samples_to_plot=None,
    neighborhood="auto",
    radius=None,
    k=6,
    n=7,
    seed=123,
):
    """Plot example neighborhoods.

    Draws a patch of each sample, highlighting a few spots/cells and their
    neighbors, to check that the neighborhood settings give the neighborhoods
    intended. Uses the same settings as spatial_similarity would given the same
    arguments. Needs matplotlib (pip install friendlycells[plot]).

    samples, neighborhood, radius, k: as for spatial_similarity.
    samples_to_plot: names of the samples to plot. By default, all samples.
    n: number of neighborhoods to highlight in each sample. Colors repeat after 7.
    seed: random seed, so the same patches and spots/cells are plotted each
        time. None for different ones each time.
    Returns the matplotlib Figure.
    """
    if samples_to_plot is None:
        samples_to_plot = list(samples)
    unknown = [name for name in samples_to_plot if name not in samples]
    if unknown:
        raise ValueError("Unknown samples_to_plot: " + ", ".join(map(str, unknown)))
    neighborhood, radius = _neighborhood_settings(samples, neighborhood, radius)
    fig, axes = _panels(len(samples_to_plot))
    rng = np.random.default_rng(seed)
    for ax, name in zip(axes.flat, samples_to_plot):
        xy = samples[name][["x", "y"]].to_numpy(dtype=float)
        _draw_patch(
            ax, xy, *_example_neighborhoods(xy, neighborhood, radius, k, n, rng)
        )
        ax.set_title(name)
    return fig


def _panels(count):
    """A figure with a grid of count panels, without axes."""
    try:
        # Optional dependency, so only imported when plotting
        # pylint: disable-next=import-outside-toplevel
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise ImportError(
            "plot_neighborhoods needs matplotlib: pip install friendlycells[plot]"
        ) from e
    cols = int(np.ceil(np.sqrt(count)))
    rows = int(np.ceil(count / cols))
    fig, axes = plt.subplots(
        rows, cols, squeeze=False, figsize=(4 * cols, 4 * rows), layout="constrained"
    )
    for ax in axes.flat:
        ax.set_axis_off()
    return fig, axes


def _draw_patch(ax, xy, patch, hoods):
    """Draws the patch's spots/cells, coloring each neighborhood."""
    # pylint: disable=import-outside-toplevel
    from matplotlib.colors import to_rgb
    from matplotlib.collections import EllipseCollection

    fill = np.full((len(xy), 3), 0.85)
    for i, hood in enumerate(hoods):
        color = np.array(to_rgb(_COLORS[i % len(_COLORS)]))
        # Halfway to white
        fill[hood[1:]] = (color + 1) / 2
        fill[hood[0]] = color
    spacing = np.median(
        _nn_distances({"patch": pd.DataFrame(xy[patch], columns=["x", "y"])})
    )
    if not np.isfinite(spacing):
        spacing = 1
    # Circles in data units, so spots/cells keep their spacing at any size
    ax.add_collection(
        EllipseCollection(
            widths=0.9 * spacing,
            heights=0.9 * spacing,
            angles=0,
            units="xy",
            offsets=xy[patch],
            offset_transform=ax.transData,
            facecolors=fill[patch],
            edgecolors="none",
        )
    )
    lo, hi = xy[patch].min(axis=0) - spacing, xy[patch].max(axis=0) + spacing
    ax.set_xlim(lo[0], hi[0])
    ax.set_ylim(lo[1], hi[1])
    ax.set_aspect("equal")


def _example_neighborhoods(xy, neighborhood, radius, k, n, rng, patch_size=1000):
    """A patch of the patch_size spots/cells nearest a random one, and up to n
    separate neighborhoods of spots/cells in its inner half, where they
    match the neighborhoods in the full sample."""
    start = xy[rng.integers(len(xy))]
    patch = np.argsort(((xy - start) ** 2).sum(axis=1), kind="stable")[:patch_size]
    pairs = _neighbor_pairs(xy[patch], neighborhood, radius, k)
    # pairs are sorted by their first column
    neighbors = np.split(
        pairs[:, 1], np.searchsorted(pairs[:, 0], np.arange(1, len(patch)))
    )
    used = np.zeros(len(patch), dtype=bool)
    hoods = []
    for i in rng.permutation((len(patch) + 1) // 2):
        if len(hoods) == n:
            break
        hood = np.concatenate([[i], neighbors[i]]).astype(int)
        if used[hood].any():
            continue
        # Their neighbors too, so neighborhoods don't touch
        used[np.concatenate([hood, *(neighbors[j] for j in hood)])] = True
        hoods.append(patch[hood])
    return patch, hoods
