"""Tests for spatial_similarity and its helpers."""

import io
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from friendlycells import plot_neighborhoods, spatial_similarity
from friendlycells.plotting import _example_neighborhoods
from friendlycells.similarity import _default_radius, _neighbor_pairs, _nn_distances


def visium(array_row, array_col, label):
    """Visium spots, from array coordinates to micrometres."""
    array_row, array_col = np.asarray(array_row), np.asarray(array_col)
    return pd.DataFrame(
        {"x": 50 * array_col, "y": 50 * np.sqrt(3) * array_row, "label": label}
    )


SAMPLES = {
    "A": visium([0, 0, 1], [0, 2, 1], ["a", "b", "a"]),
    "B": visium([0, 0, 0, 2], [0, 2, 4, 0], ["a", "a", "b", "b"]),
    "C": visium([0, 0], [0, 2], ["b", "b"]),
}


def assert_sim(sim, ab, ac, bc):
    """Check sim is the 3x3 similarity matrix with off-diagonals ab, ac, bc."""
    expected = np.array([[1, ab, ac], [ab, 1, bc], [ac, bc, 1]])
    np.testing.assert_allclose(sim.to_numpy(), expected)
    assert list(sim.columns) == list(sim.index)


def test_neighbor_counts():
    """Hexagonal neighbor counts give similarity 1 - d / max(d)."""
    sim = spatial_similarity(SAMPLES, "x", "y", "label", sample_normalize=False)
    assert_sim(sim, 1 - np.sqrt(2) / 4, 0, 1 - np.sqrt(10) / 4)


def test_names_in_input_order():
    """Output row/column names respect input sample order."""
    sim = spatial_similarity(
        dict(reversed(SAMPLES.items())), "x", "y", "label", sample_normalize=False
    )
    assert list(sim.index) == ["C", "B", "A"]


def test_identical_samples():
    """Identical samples are seen as fully similar."""
    sim = spatial_similarity({"A": SAMPLES["A"], "B": SAMPLES["A"]}, "x", "y", "label")
    np.testing.assert_array_equal(sim.to_numpy(), np.ones((2, 2)))


def test_categorical_labels():
    """Categorical labels yield the same output as string labels."""
    sim = spatial_similarity(SAMPLES, "x", "y", "label", sample_normalize=False)
    categories = {
        n: s.assign(label=s["label"].astype("category")) for n, s in SAMPLES.items()
    }
    pd.testing.assert_frame_equal(
        spatial_similarity(categories, "x", "y", "label", sample_normalize=False), sim
    )


def test_missing_labels():
    """Labels read as missing (e.g. NA, nan, None) become one label called 'NA', with a warning."""
    sim = spatial_similarity(SAMPLES, "x", "y", "label", sample_normalize=False)
    # pd.read_csv() turns labels written as NA, nan, None, etc. into NaN
    csv = pd.read_csv(
        io.StringIO(
            "sample,x,y,label\n"
            "A,0,0,a\nA,100,0,NA\nA,50,86.60254037844386,a\n"
            "B,0,0,a\nB,100,0,a\nB,200,0,nan\nB,0,173.20508075688772,nan\n"
            "C,0,0,None\nC,100,0,None\n"
        )
    )
    from_file = {name: s.drop(columns="sample") for name, s in csv.groupby("sample")}
    with pytest.warns(UserWarning, match="treated as a label called 'NA'"):
        pd.testing.assert_frame_equal(
            spatial_similarity(from_file, "x", "y", "label", sample_normalize=False),
            sim,
        )


def test_normalize():
    """Normalize divides counts by their total."""
    dmax = np.sqrt(11 / 8)
    sim = spatial_similarity(SAMPLES, "x", "y", "label")
    assert_sim(sim, 1 - (np.sqrt(6) / 12) / dmax, 1 - np.sqrt(4 / 3) / dmax, 0)


def test_label_normalize():
    """label_normalize scales by self-connections, then by degree."""
    # A's counts (a-a 2, a-b 2, b-b 0), with a self-connection added (3, 2, 1)
    # and divided by sqrt(3 * 1) where a meets b, give each label 1 to itself
    # and 2 / sqrt(3) to the other; B's (2, 1, 0) give 1 / sqrt(3); C's are
    # all to itself. p and q are A's and B's shares of each label to itself
    p = 1 / (1 + 2 / np.sqrt(3))
    q = 1 / (1 + 1 / np.sqrt(3))
    sim = spatial_similarity(
        SAMPLES, "x", "y", "label", sample_normalize=False, label_normalize=True
    )
    assert_sim(sim, 1 - (q - p) / (1 - p), 0, 1 - (1 - q) / (1 - p))
    with pytest.raises(ValueError, match="not both"):
        spatial_similarity(SAMPLES, "x", "y", "label", label_normalize=True)


def test_no_connections_with_normalize():
    """Samples where no spot has any neighbor yield NA and a warning."""
    sim = spatial_similarity(SAMPLES, "x", "y", "label")
    # D's two spots don't touch, so it has no connections to normalize
    d = visium([0, 0], [0, 4], ["a", "a"])
    with pytest.warns(UserWarning, match="no neighboring spots: D"):
        with_d = spatial_similarity(
            {**SAMPLES, "D": d}, "x", "y", "label", neighborhood="radius"
        )
    pd.testing.assert_frame_equal(with_d.loc[["A", "B", "C"], ["A", "B", "C"]], sim)
    assert with_d["D"].isna().all() and with_d.loc["D"].isna().all()


def test_pools():
    """Pools sum counts of their samples."""
    pools = {"p1": ["A", "B"], "p2": ["C"], "p3": ["A"]}
    sim = spatial_similarity(
        SAMPLES, "x", "y", "label", sample_normalize=False, pools=pools
    )
    assert list(sim.index) == ["p1", "p2", "p3"]
    assert_sim(sim, 0, 1 - np.sqrt(6) / np.sqrt(38), 1 - 4 / np.sqrt(38))


def test_pig_esophagus_paper():
    """Reproduces the pig esophagus paper's similarities."""
    # Durkin et al. (2026)
    spots = pd.read_csv(Path(__file__).parent / "pig_esophagus_spots.csv.xz")
    samples = {
        name: visium(s["array_row"], s["array_col"], s["label"].to_numpy())
        for name, s in spots.groupby("sample_name")
    }
    pools = {
        "Pre-implant": ["Pre-implant"],
        "1 month": ["1 month (D)"],
        "3 months": ["3 months (E)", "3 months (F)"],
        "6 months": ["6 months (B)", "6 months (C)"],
        "Native": ["Native DP1", "Native DP2", "Native DP3", "Native (C)"],
    }
    sim = spatial_similarity(samples, "x", "y", "label", pools=pools)
    paper = np.array(
        [
            [1.00, 0.01, 0.03, 0.09, 0.00],
            [0.01, 1.00, 0.90, 0.68, 0.27],
            [0.03, 0.90, 1.00, 0.70, 0.28],
            [0.09, 0.68, 0.70, 1.00, 0.52],
            [0.00, 0.27, 0.28, 0.52, 1.00],
        ]
    )
    # The paper's 7.5 px neighbor cutoff missed a few neighbors in Native DP2
    paper[0, 1] = paper[1, 0] = 0.02
    assert list(sim.index) == list(pools)
    np.testing.assert_allclose(sim.round(2).to_numpy(), paper, atol=1e-12)


def test_only_transitions():
    """only_transitions drops same-label contacts, so only borders count."""
    sim = spatial_similarity(
        SAMPLES, "x", "y", "label", sample_normalize=False, only_transitions=True
    )
    assert_sim(sim, 0.5, 0, 0.5)


def test_only_transitions_single_label():
    """One-label samples have no transitions, giving NA and a warning."""
    with pytest.warns(
        UserWarning, match=r"with different labels \(no transitions\): C"
    ):
        sim = spatial_similarity(SAMPLES, "x", "y", "label", only_transitions=True)
    assert sim.loc["A", "B"] == 1
    assert sim["C"].isna().all() and sim.loc["C"].isna().all()


def test_only_transitions_default_unchanged():
    """only_transitions=False gives the default result."""
    pd.testing.assert_frame_equal(
        spatial_similarity(SAMPLES, "x", "y", "label", only_transitions=False),
        spatial_similarity(SAMPLES, "x", "y", "label"),
    )


def test_only_transitions_with_pools():
    """Pools sum the border counts of their samples."""
    pools = {"p1": ["A", "B"], "p2": ["C"], "p3": ["A"]}
    sim = spatial_similarity(
        SAMPLES,
        "x",
        "y",
        "label",
        sample_normalize=False,
        pools=pools,
        only_transitions=True,
    )
    assert_sim(sim, 0, 2 / 3, 1 / 3)


def test_radius():
    """Radius links cells up to and including the radius."""
    xy = np.array([[0, 0], [10, 0], [20.5, 0]])
    np.testing.assert_array_equal(
        _neighbor_pairs(xy, "radius", radius=10), [[0, 1], [1, 0]]
    )


def test_knn():
    """kNN links cells both ways if either is among the other's nearest."""
    # 2's nearest is 1, but 1's nearest is 0
    xy = np.array([[0, 0], [1, 0], [3, 0]])
    np.testing.assert_array_equal(
        _neighbor_pairs(xy, "knn", k=1), [[0, 1], [1, 0], [1, 2], [2, 1]]
    )


def test_delaunay():
    """Delaunay links surrounding cells, dropping links over radius."""
    # A triangle, and a far cell that Delaunay links to it
    xy = np.array([[0, 0], [10, 0], [5, 8], [100, 0]])
    np.testing.assert_array_equal(
        _neighbor_pairs(xy, "delaunay", radius=20),
        [[0, 1], [0, 2], [1, 0], [1, 2], [2, 0], [2, 1]],
    )


def test_delaunay_collinear():
    """Delaunay links collinear cells to the next along the line."""
    xy = np.array([[0, 0], [20, 0], [10, 0]])
    np.testing.assert_array_equal(
        _neighbor_pairs(xy, "delaunay"), [[0, 2], [1, 2], [2, 0], [2, 1]]
    )


def test_columns():
    """Columns are chosen by name, with an error naming missing ones."""
    renamed = {
        n: pd.DataFrame(
            {"col": s["x"], "row": s["y"], "cell_type": s["label"], "label": "z"}
        )
        for n, s in SAMPLES.items()
    }
    pd.testing.assert_frame_equal(
        spatial_similarity(renamed, "col", "row", "cell_type"),
        spatial_similarity(SAMPLES, "x", "y", "label"),
    )
    with pytest.raises(
        ValueError, match=r'"cell_type" \(labels\) is missing from samples: A, B, C'
    ):
        spatial_similarity(SAMPLES, "x", "y", "cell_type")


def test_default_radius():
    """Default radius scales the median spacing over all samples."""
    # Nearest-neighbor distances 10, 10 and 22, 22, 28, so the spacing is 22
    s = {
        "a": pd.DataFrame({"x": [0, 10], "y": 0}),
        "b": pd.DataFrame({"x": [0, 22, 50], "y": 0}),
    }
    nn = _nn_distances(s, "x", "y")
    assert _default_radius(nn, "radius") == pytest.approx(1.2 * 22)
    assert _default_radius(nn, "delaunay") == pytest.approx(3 * 22)


def test_auto(capsys):
    """Auto uses radius on a grid and delaunay otherwise, saying which."""
    spatial_similarity(SAMPLES, "x", "y", "label")
    assert 'neighborhood="radius"' in capsys.readouterr().err
    cells = {
        "A": pd.DataFrame(
            {"x": [0, 4, 11, 13, 30], "y": [0, 9, 2, 20, 6], "label": "a"}
        )
    }
    spatial_similarity(cells, "x", "y", "label")
    assert 'neighborhood="delaunay"' in capsys.readouterr().err


def test_default_radius_units():
    """With the default radius, results don't depend on units."""
    pixels = {n: s.assign(x=s["x"] / 0.3, y=s["y"] / 0.3) for n, s in SAMPLES.items()}
    pd.testing.assert_frame_equal(
        spatial_similarity(pixels, "x", "y", "label"),
        spatial_similarity(SAMPLES, "x", "y", "label"),
    )


def test_example_neighborhoods():
    """Example neighborhoods don't touch and match the full sample's."""
    rng = np.random.default_rng(1)
    xy = rng.uniform(0, 100, (3000, 2))
    _, hoods = _example_neighborhoods(xy, "delaunay", 6, 6, 10, rng)
    assert len(hoods) == 10
    flat = np.concatenate(hoods)
    assert len(np.unique(flat)) == len(flat)
    pairs = _neighbor_pairs(xy, "delaunay", radius=6)
    hood_of = np.zeros(len(xy), dtype=int)
    for i, hood in enumerate(hoods, start=1):
        hood_of[hood] = i
    linked = pairs[(hood_of[pairs[:, 0]] > 0) & (hood_of[pairs[:, 1]] > 0)]
    assert (hood_of[linked[:, 0]] == hood_of[linked[:, 1]]).all()
    for hood in hoods:
        assert set(hood[1:]) == set(pairs[pairs[:, 0] == hood[0], 1])


def test_unknown_samples_to_plot():
    """Unknown samples_to_plot give an error."""
    with pytest.raises(ValueError, match="Unknown"):
        plot_neighborhoods(SAMPLES, "x", "y", samples_to_plot=["D"])
