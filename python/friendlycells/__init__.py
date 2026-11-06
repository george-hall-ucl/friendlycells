"""
friendlycells computes similarity between spatial transcriptomics samples based
on which labels (e.g. cell types) border one another.
"""

from .plotting import plot_neighborhoods
from .similarity import spatial_similarity

__all__ = ["plot_neighborhoods", "spatial_similarity"]
