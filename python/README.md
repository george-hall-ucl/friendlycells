<!-- README.md is generated from README.Rmd with knitr::knit("README.Rmd"). Please edit that file -->



<p align="center">
  <img src="https://raw.githubusercontent.com/george-hall-ucl/friendlycells/main/assets/logo.gif" width="210" alt="friendlycells logo">
</p>

Compute similarity between spatial transcriptomics samples based on which
labels (e.g. cell types) border one another.

## Installation

```bash
pip install "git+https://github.com/george-hall-ucl/friendlycells#subdirectory=python"
```

Requires `pandas`, `numpy` and `scipy`. Also available as an
[R package](https://github.com/george-hall-ucl/friendlycells/tree/main/R).

## Usage

Compatible with Visium, Visium HD (binned or segmented) and Xenium.

```python
from friendlycells import spatial_similarity

# One DataFrame per sample, with columns x, y (positions, e.g. in pixels or
# micrometres) and label (see Extracting positions below)
samples = {"sample_1": df1, "sample_2": df2, "sample_3": df3}

sim = spatial_similarity(samples)
```

Returns a sample × sample DataFrame with normalised similarity scores: 1 =
identical, 0 = least similar pair.

### Options

- `normalize=False`: compare raw connection counts instead of proportions. This
  will make samples with different numbers of spots be considered spatially
  different even if their underlying spatial features are similar.
- `pools={"group_1": ["sample_1", "sample_2"], "group_2": ["sample_3"]}`: sum
  connection counts within groups before comparing.
- `only_transitions=True`: only count connections between spots with distinct
  labels (i.e. borders), not between spots with the same label.
- `neighborhood`: how to find each spot's or cell's neighbors. The default,
  `"auto"`, uses `"radius"` (all within `radius`) if the positions lie on a
  grid, and `"delaunay"` (the cells surrounding each cell) otherwise, with a
  message saying which. `"knn"` (the `k` nearest, default 6) is also available.
- `radius`: by default it's set from the spacing of the data, the median
  distance from each spot or cell to its nearest neighbor over all samples,
  times 1.2 for `"radius"` or 3 for `"delaunay"` (where it's the maximum link
  length). Give it, in the units of `x` and `y`, to widen or narrow
  neighborhoods.

### Extracting positions

For Visium read with scanpy or squidpy, and for tables read with
[spatialdata-io](https://spatialdata.scverse.org/projects/io/) (Xenium's
`adata = sdata.tables["table"]`, or Visium HD bins' e.g.
`adata = sdata.tables["square_008um"]`), with cell types in
`adata.obs["cell_type"]`:

```python
x, y = adata.obsm["spatial"].T
df = pd.DataFrame({"x": x, "y": y, "label": adata.obs["cell_type"]})
```

For Visium HD segmented cells, `sdata.tables["cell_segmentations"]` has no
`obsm["spatial"]`, so take the centroids of the cells' outlines instead:

```python
adata = sdata.tables["cell_segmentations"]
cells = sdata.shapes[adata.obs["region"].iloc[0]].centroid.loc[adata.obs_names]
df = pd.DataFrame({"x": cells.x, "y": cells.y, "label": adata.obs["cell_type"]})
```

### Neighbors

With the defaults, each sample's spots or cells are linked to:

- Visium spots: their six immediate neighbors.
- Visium HD bins: their four side neighbors. Set `radius` to 1.5 times the bin
  size to add the diagonals.
- Segmented cells (Xenium, or Visium HD with cell segmentation): the cells
  surrounding them, ignoring links across gaps.

Positions can be the pixel or micrometre coordinates from Space Ranger or
Xenium outputs, as long as all samples use the same units. As the spacing is
taken over all samples, don't mix samples in different units, or from
different technologies.

To check them, `plot_neighborhoods(samples)` draws a patch of each sample with
seven spots or cells highlighted, each in its own color, and their neighbors in a
lighter shade. It uses the same settings as `spatial_similarity()`, so pass it
the same `neighborhood`, `radius` and `k`. Choose samples with
`samples_to_plot`, and change `seed` to draw different spots.
It needs matplotlib.

### Example

Uses four mouse brain sections from 10x Genomics (two anterior, two posterior),
with each spot labelled by its most enriched of five marker genes.

<details><summary>Setup: download and label the data, define plotting helpers (click to expand)</summary>


``` python
# To rebuild this README from README.Rmd, run knitr::knit("README.Rmd") in R
# from the python/ directory, with RETICULATE_PYTHON set to a Python that has
# friendlycells installed (rmarkdown::render breaks the <details> block).
# Requires squidpy, matplotlib and seaborn.
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import squidpy as sq

markers = ["Mbp", "Slc17a7", "Penk", "Hpca", "Gabra6"]
ids = {
    "anterior_1": "V1_Mouse_Brain_Sagittal_Anterior",
    "anterior_2": "V1_Mouse_Brain_Sagittal_Anterior_Section_2",
    "posterior_1": "V1_Mouse_Brain_Sagittal_Posterior",
    "posterior_2": "V1_Mouse_Brain_Sagittal_Posterior_Section_2",
}
cache = Path.home() / ".cache" / "friendlycells"

adatas = {}
labels = {}
for name, sample_id in ids.items():
    adata = sq.datasets.visium(sample_id, base_dir=cache)
    adata.var_names_make_unique()
    counts = adata[:, markers].X.toarray()
    top = (counts / counts.mean(axis=0)).argmax(axis=1)
    adatas[name] = adata
    labels[name] = np.array(markers)[top]


def plot_slices(samples):
    palette = ["#CC79A7", "#D55E00", "#0072B2", "#F0E442", "#009E73"]
    colors = dict(zip(sorted(markers), palette))
    fig, axes = plt.subplots(1, len(samples), layout="constrained")
    for ax, (name, df) in zip(axes, samples.items()):
        ax.scatter(
            df["x"],
            -df["y"],
            c=df["label"].map(colors),
            s=3,
            linewidths=0,
        )
        ax.set_title(name)
        ax.set_aspect("equal")
        ax.axis("off")
    handles = [
        plt.Line2D([], [], marker="o", linestyle="", color=color, label=label)
        for label, color in colors.items()
    ]
    fig.legend(
        handles=handles,
        loc="outside lower center",
        ncol=len(colors),
        frameon=False,
    )
    plt.show()


def plot_similarity(sim):
    style = {"font.size": 14, "font.weight": "bold", "figure.autolayout": True}
    with plt.rc_context(style):
        ax = sns.heatmap(
            sim,
            mask=np.tril(np.ones(sim.shape, dtype=bool), k=-1),
            annot=True,
            fmt=".2f",
            cmap="inferno",
            vmin=0,
            vmax=1,
            square=True,
            cbar_kws={"label": "Similarity"},
        )
        ax.invert_yaxis()
        plt.xticks(rotation=45, ha="right")
    plt.show()
```

</details>

Each sample is a DataFrame of spot positions (in full-resolution image
pixels) and labels:


``` python
from friendlycells import spatial_similarity

samples = {
    name: pd.DataFrame({
        "x": adata.obsm["spatial"][:, 0],
        "y": adata.obsm["spatial"][:, 1],
        "label": labels[name],
    }, index=adata.obs_names)
    for name, adata in adatas.items()
}
print(samples["anterior_1"].head())
#>                        x     y    label
#> AAACAAGTATCTCCCA-1  8500  7474      Mbp
#> AAACACCAATAACTGC-1  2788  8552  Slc17a7
#> AAACAGAGCGACTCCT-1  7950  3163  Slc17a7
#> AAACAGCTTTCAGAAG-1  2100  6636      Mbp
#> AAACAGGGTCTATATT-1  2375  7115     Penk
```


``` python
plot_slices(samples)
```

<img src="https://raw.githubusercontent.com/george-hall-ucl/friendlycells/main/python/figures/README-slices-1.png" alt="" width="75%" style="display: block; margin: auto;" />


``` python
sim = spatial_similarity(samples)
#> Using neighborhood="radius", as the positions lie on a grid.
print(sim.round(2))
#>              anterior_1  anterior_2  posterior_1  posterior_2
#> anterior_1         1.00        0.93         0.00         0.01
#> anterior_2         0.93        1.00         0.03         0.04
#> posterior_1        0.00        0.03         1.00         0.92
#> posterior_2        0.01        0.04         0.92         1.00
plot_similarity(sim)
```

<img src="https://raw.githubusercontent.com/george-hall-ucl/friendlycells/main/python/figures/README-similarity-3.png" alt="" width="50%" style="display: block; margin: auto;" />

## Citation

The algorithm and original implementation of `friendlycells` are based on the
methodology in:

- Durkin _et al_. "Functional integration of an autologous engineered esophagus
  in a large-animal model." _Nature Biotechnology_ (2026).
- Budinger _et al_. "An in vivo and in vitro spatiotemporal profile of human
  midbrain development." _Nature Communications_ (2026).

## Licence

Copyright (C) 2026 University College London. Licensed under GPL-3.0-or-later.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see http://www.gnu.org/licenses/.
