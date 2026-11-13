"""Computes similarity between samples from counts of neighboring labels."""

import sys
import warnings

import numpy as np
import pandas as pd
from scipy.spatial import Delaunay, cKDTree


def spatial_similarity(
    samples,
    x,
    y,
    labels,
    normalize=True,
    pools=None,
    only_transitions=False,
    neighborhood="auto",
    radius=None,
    k=6,
):
    """Computes spatial similarity between spatial transcriptomics samples by
    counting connections between neighboring spots/cells.

    samples: dict of DataFrames, one per sample, with a row per spot/cell.
    x, y: names of the columns holding the positions of spots/cells, in any
        unit (e.g. pixels or micrometres) consistent between samples.
    labels: name of the column holding the labels (e.g. cell types).
    normalize: divide each sample's connection counts by their total to
        control for sample size.
    pools: optional dict mapping pool names to lists of sample names;
        connection counts are summed within each pool.
    only_transitions: if True, only count connections between spots with
        distinct labels (i.e. borders); does not count connections between
        spots with the same label.
    neighborhood: definition of the neighbourhood of each spot/cell: "radius"
        selects all within radius; "knn" selects the k nearest (symmetric:
        both are included in each other's neighborhood if either is among the
        other's nearest); or "delaunay" selects those surrounding it by
        Delaunay triangulation, dropping links longer than radius. The
        default, "auto", uses "radius" if the positions lie on a grid (at least
        75% of nearest-neighbor distances within 10% of their median) and
        "delaunay" otherwise, with a message saying which was used.
    radius: neighborhood radius for "radius", or maximum link length for
        "delaunay", in the units of the positions. By default, set from the
        spacing of the data: set to the median distance from each spot/cell to
        its nearest neighbor over all samples, multiplied by 1.2 for "radius"
        or 3 for "delaunay". Used to make neighborhoods wider or narrower.
    k: number of nearest neighbors for neighborhood="knn". Ties are broken
        arbitrarily.
    Returns a symmetric DataFrame of similarities in [0, 1], or NaN for
    samples with no neighboring spots when normalize=True.

    Under default settings, the neighbourhoods of each sample's spots/cells
    depend on technology:
    - Visium spots: six immediate neighbors.
    - Visium HD bins: four side neighbors. Set radius to 1.5 times the bin size
      to add the diagonals.
    - Segmented cells (Xenium, or Visium HD with cell segmentation): the cells
      surrounding them, ignoring links across gaps.

    Positions can be the pixel or micrometre coordinates from Space Ranger or
    Xenium outputs, as long as all samples use the same units. As the spacing
    is taken over all samples, don't mix samples in different units or from
    different technologies.
    """
    neighborhood, radius = _neighborhood_settings(
        samples, {"x": x, "y": y, "labels": labels}, neighborhood, radius
    )
    if any(s[labels].isna().any() for s in samples.values()):
        # e.g. pd.read_csv() reads a label literally called "NA" as missing
        warnings.warn(
            "Some labels are 'NA'. These will be treated as a label called "
            "'NA' and used in the comparison, rather than being treated as "
            "the absence of a label and dropped."
        )
    # Sorted as set order varies between runs, changing floating-point results
    all_labels = sorted(
        set().union(*(_na_as_label(s[labels]) for s in samples.values()))
    )
    counts = pd.DataFrame(
        {
            name: _connection_counts(
                s, x, y, labels, all_labels, only_transitions, neighborhood, radius, k
            )
            for name, s in samples.items()
        }
    ).T
    if pools is not None:
        counts = pd.DataFrame(
            {name: counts.loc[members].sum() for name, members in pools.items()}
        ).T
    if normalize:
        totals = counts.sum(axis=1)
        if (totals == 0).any():
            warnings.warn(
                "Similarities are NA for samples with no neighboring spots"
                + (
                    " with different labels (no transitions)"
                    if only_transitions
                    else ""
                )
                + ": "
                + ", ".join(map(str, counts.index[totals == 0]))
            )
        counts = counts.div(totals, axis=0)
    x = counts.to_numpy(dtype=float)
    d = np.sqrt(((x[:, None, :] - x[None, :, :]) ** 2).sum(axis=2))
    # Identical samples give a max distance of 0, so don't divide by it
    d_max = np.nanmax(d, initial=0)
    if d_max > 0:
        d = d / d_max
    return pd.DataFrame(1 - d, index=counts.index, columns=counts.index)


def _neighborhood_settings(samples, columns, neighborhood, radius):
    """Checks each sample has the named columns (a dict of x, y and maybe labels),
    then resolves "auto" and the default radius from all samples."""
    if neighborhood not in ("auto", "radius", "knn", "delaunay"):
        raise ValueError("neighborhood must be 'auto', 'radius', 'knn' or 'delaunay'.")
    for arg, column in columns.items():
        missing = [name for name, s in samples.items() if column not in s.columns]
        if missing:
            raise ValueError(
                f'Column "{column}" ({arg}) is missing from samples: '
                + ", ".join(map(str, missing))
            )
    if neighborhood == "auto" or (radius is None and neighborhood != "knn"):
        nn = _nn_distances(samples, columns["x"], columns["y"])
    if neighborhood == "auto":
        # pylint: disable-next=possibly-used-before-assignment
        neighborhood = "radius" if _on_grid(nn) else "delaunay"
        lie = "lie" if neighborhood == "radius" else "don't lie"
        print(
            f'Using neighborhood="{neighborhood}", as the positions {lie} on a grid.',
            file=sys.stderr,
        )
    if radius is None and neighborhood != "knn":
        radius = _default_radius(nn, neighborhood)
    return neighborhood, radius


def _nn_distances(samples, x, y):
    """Distance from each spot or cell to its nearest neighbor, over all samples."""
    xys = [s[[x, y]].to_numpy(dtype=float) for s in samples.values()]
    return np.concatenate(
        [cKDTree(xy).query(xy, k=2)[0][:, 1] for xy in xys if len(xy) >= 2]
        or [[np.nan]]
    )


def _on_grid(nn):
    """On a grid, nearly all spots are the same distance from their nearest neighbor,
    while segmented cells vary widely."""
    return bool(np.mean(np.abs(nn / np.median(nn) - 1) < 0.1) >= 0.75)


def _default_radius(nn, neighborhood):
    """The median nearest-neighbor distance, times 1.2 for "radius" (Visium's six neighbors,
    or HD bins' four, allowing for jitter) or 3 for "delaunay" (dropping only links across gaps).
    """
    return np.median(nn) * (1.2 if neighborhood == "radius" else 3)


def _connection_counts(
    sample,
    x,
    y,
    labels,
    all_labels,
    only_transitions=False,
    neighborhood="radius",
    radius=None,
    k=6,
):
    codes = pd.Categorical(_na_as_label(sample[labels]), categories=all_labels).codes
    pairs = _neighbor_pairs(
        sample[[x, y]].to_numpy(dtype=float), neighborhood, radius, k
    )
    counts = np.zeros((len(all_labels), len(all_labels)))
    np.add.at(counts, (codes[pairs[:, 0]], codes[pairs[:, 1]]), 1)
    if only_transitions:
        np.fill_diagonal(counts, 0)
    return counts.ravel()


def _na_as_label(labels):
    """Missing labels are counted as a label called "NA"."""
    return labels.astype(str).where(labels.notna(), "NA")


def _neighbor_pairs(xy, neighborhood, radius=None, k=6):
    """Pairs (i, j) of neighboring rows of xy, each pair listed in both directions."""
    n = len(xy)
    if n < 2:
        return np.empty((0, 2), dtype=int)
    if neighborhood == "radius":
        pairs = cKDTree(xy).query_pairs(radius, output_type="ndarray")
    elif neighborhood == "knn":
        # The nearest point to each is itself
        _, nbr = cKDTree(xy).query(xy, k=min(k, n - 1) + 1)
        pairs = np.column_stack(
            [np.repeat(np.arange(n), nbr.shape[1] - 1), nbr[:, 1:].ravel()]
        )
    else:
        pairs = _delaunay_edges(xy)
        if radius is not None:
            pairs = pairs[
                np.linalg.norm(xy[pairs[:, 0]] - xy[pairs[:, 1]], axis=1) <= radius
            ]
    # kNN links can be one-way, so add the reverse of every pair, then dedupe
    return np.unique(np.vstack([pairs, pairs[:, ::-1]]), axis=0)


def _delaunay_edges(xy):
    centred = xy - xy.mean(axis=0)
    if np.linalg.matrix_rank(centred) < 2:
        # Collinear points can't be triangulated, so link each to the next along the line
        order = np.argsort(centred @ np.linalg.svd(centred)[2][0])
        return np.column_stack([order[:-1], order[1:]])
    tri = Delaunay(xy).simplices
    return np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [0, 2]]])
