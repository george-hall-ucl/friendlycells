import io
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from friendlyspots import spatial_similarity

SAMPLES = {
    "A": pd.DataFrame({"array_row": [0, 0, 1], "array_col": [0, 2, 1], "label": ["a", "b", "a"]}),
    "B": pd.DataFrame({"array_row": [0, 0, 0, 2], "array_col": [0, 2, 4, 0], "label": ["a", "a", "b", "b"]}),
    "C": pd.DataFrame({"array_row": [0, 0], "array_col": [0, 2], "label": ["b", "b"]}),
}


def assert_sim(sim, ab, ac, bc):
    expected = np.array([[1, ab, ac], [ab, 1, bc], [ac, bc, 1]])
    np.testing.assert_allclose(sim.to_numpy(), expected)
    assert list(sim.columns) == list(sim.index)


def test_neighbor_counts():
    """Hexagonal neighbor counts give similarity 1 - d / max(d)."""
    sim = spatial_similarity(SAMPLES, normalize=False)
    assert_sim(sim, 1 - np.sqrt(2) / 4, 0, 1 - np.sqrt(10) / 4)


def test_names_in_input_order():
    """Output row/column names respect input sample order."""
    sim = spatial_similarity(dict(reversed(SAMPLES.items())), normalize=False)
    assert list(sim.index) == ["C", "B", "A"]


def test_identical_samples():
    """Identical samples are seen as fully similar."""
    sim = spatial_similarity({"A": SAMPLES["A"], "B": SAMPLES["A"]})
    np.testing.assert_array_equal(sim.to_numpy(), np.ones((2, 2)))


def test_categorical_labels():
    """Categorical labels yield the same output as string labels."""
    sim = spatial_similarity(SAMPLES, normalize=False)
    categories = {n: s.assign(label=s["label"].astype("category")) for n, s in SAMPLES.items()}
    pd.testing.assert_frame_equal(spatial_similarity(categories, normalize=False), sim)


def test_missing_labels():
    """Labels read as missing (e.g. NA, nan, None) become one label called 'NA', with a warning."""
    sim = spatial_similarity(SAMPLES, normalize=False)
    # pd.read_csv() turns labels written as NA, nan, None, etc. into NaN
    csv = pd.read_csv(io.StringIO(
        "sample,array_row,array_col,label\n"
        "A,0,0,a\nA,0,2,NA\nA,1,1,a\n"
        "B,0,0,a\nB,0,2,a\nB,0,4,nan\nB,2,0,nan\n"
        "C,0,0,None\nC,0,2,None\n"
    ))
    from_file = {name: s.drop(columns="sample") for name, s in csv.groupby("sample")}
    with pytest.warns(UserWarning, match="treated as a label called 'NA'"):
        pd.testing.assert_frame_equal(spatial_similarity(from_file, normalize=False), sim)


def test_normalize():
    """Normalize divides counts by their total."""
    dmax = np.sqrt(11 / 8)
    sim = spatial_similarity(SAMPLES)
    assert_sim(sim, 1 - (np.sqrt(6) / 12) / dmax, 1 - np.sqrt(4 / 3) / dmax, 0)


def test_no_connections_with_normalize():
    """Samples where no spot has any neighbor yield NA and a warning."""
    sim = spatial_similarity(SAMPLES)
    # D's two spots don't touch, so it has no connections to normalize
    d = pd.DataFrame({"array_row": [0, 0], "array_col": [0, 4], "label": ["a", "a"]})
    with pytest.warns(UserWarning, match="no neighboring spots: D"):
        with_d = spatial_similarity({**SAMPLES, "D": d})
    pd.testing.assert_frame_equal(with_d.loc[["A", "B", "C"], ["A", "B", "C"]], sim)
    assert with_d["D"].isna().all() and with_d.loc["D"].isna().all()


def test_pools():
    """Pools sum counts of their samples."""
    pools = {"p1": ["A", "B"], "p2": ["C"], "p3": ["A"]}
    sim = spatial_similarity(SAMPLES, normalize=False, pools=pools)
    assert list(sim.index) == ["p1", "p2", "p3"]
    assert_sim(sim, 0, 1 - np.sqrt(6) / np.sqrt(38), 1 - 4 / np.sqrt(38))


def test_pig_esophagus_paper():
    """Reproduces the pig esophagus paper's similarities."""
    # Durkin et al. (2026)
    spots = pd.read_csv(Path(__file__).parent / "pig_esophagus_spots.csv.xz")
    samples = {name: s for name, s in spots.groupby("sample_name")}
    pools = {
        "Pre-implant": ["Pre-implant"],
        "1 month": ["1 month (D)"],
        "3 months": ["3 months (E)", "3 months (F)"],
        "6 months": ["6 months (B)", "6 months (C)"],
        "Native": ["Native DP1", "Native DP2", "Native DP3", "Native (C)"],
    }
    sim = spatial_similarity(samples, pools=pools)
    paper = np.array([
        [1.00, 0.01, 0.03, 0.09, 0.00],
        [0.01, 1.00, 0.90, 0.68, 0.27],
        [0.03, 0.90, 1.00, 0.70, 0.28],
        [0.09, 0.68, 0.70, 1.00, 0.52],
        [0.00, 0.27, 0.28, 0.52, 1.00],
    ])
    # The paper's 7.5 px neighbor cutoff missed a few neighbors in Native DP2
    paper[0, 1] = paper[1, 0] = 0.02
    assert list(sim.index) == list(pools)
    np.testing.assert_allclose(sim.round(2).to_numpy(), paper, atol=1e-12)


def test_only_transitions():
    """only_transitions drops same-label contacts, so only borders count."""
    sim = spatial_similarity(SAMPLES, normalize=False, only_transitions=True)
    assert_sim(sim, 0.5, 0, 0.5)


def test_only_transitions_single_label():
    """One-label samples have no transitions, giving NA and a warning."""
    with pytest.warns(UserWarning, match=r"with different labels \(no transitions\): C"):
        sim = spatial_similarity(SAMPLES, only_transitions=True)
    assert sim.loc["A", "B"] == 1
    assert sim["C"].isna().all() and sim.loc["C"].isna().all()


def test_only_transitions_default_unchanged():
    """only_transitions=False gives the default result."""
    pd.testing.assert_frame_equal(spatial_similarity(SAMPLES, only_transitions=False), spatial_similarity(SAMPLES))


def test_only_transitions_with_pools():
    """Pools sum the border counts of their samples."""
    pools = {"p1": ["A", "B"], "p2": ["C"], "p3": ["A"]}
    sim = spatial_similarity(SAMPLES, normalize=False, pools=pools, only_transitions=True)
    assert_sim(sim, 0, 2 / 3, 1 / 3)
