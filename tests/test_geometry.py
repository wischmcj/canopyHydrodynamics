from __future__ import annotations

import pytest
from shapely.geometry import Polygon

from canhydro.geometry import get_projected_overlap
from tests.expected_results_shapes import squares_projection_overlap


def test_2D_overlap():
    polys_bottom = [
        Polygon([(0, 0), (2, 0), (2, 2), (0, 2)]),
        Polygon([(1, 1), (3, 1), (3, 3), (1, 3)]),
    ]
    polys_mid = [
        Polygon([(2, 2), (4, 2), (4, 4), (2, 4)]),
        Polygon([(3, 3), (5, 3), (5, 5), (3, 5)]),
    ]
    polys_top = [Polygon([(4, 4), (6, 4), (6, 6), (4, 6)])]

    poly_matrix = [polys_bottom, polys_mid, polys_top]
    labels = ["bottom", "mid", "top"]
    actual_overlap_dict = get_projected_overlap(poly_matrix, labels)
    assert squares_projection_overlap == actual_overlap_dict


@pytest.mark.parametrize("flexible_collection", ["5_SmallTree.csv"], indirect=True)
def test_collection_overlap(flexible_collection):
    flexible_collection.project_cylinders("XY")
    flexible_collection.draw(plane="XY", filter_lambda=lambda: branch_order > 0)
    overlap_dict = flexible_collection.find_overlap_by_percentile()
    breakpoint()
