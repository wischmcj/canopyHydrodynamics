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

def test_project_cyl():
    """
           Tests projection of cylinders parallel with:
            the XY plane, Z axis the line x=y (45 deg)
    """

def test_lu_factor():
    """
        tests circumcenter func with some mock examples
        (one, two, many points, size zero array, cooincindent coords)
    """

def test_lapack_v_lu_factor():
    """
        Tests thse answers are the same so that we only really need
        to test one
    """


def test_proj_overlap():
    """
        Tests for diff number of  (outer list) and diff number of polys (outer list)
    """

def test_curcumradius():
    """
        test with input center, and without
        corrd list of one two and many lengths (div by three and non div by three)
        test with 3+ dims
    """

def simplicies():
    """

        probs would need a mock to do this one? maybe skip for now...
    """

def maximal_alpha():
    """
        maybe find a diff package to test against
        not very important atm
    """

def test_concave_hull():
    """
        see test_alpha_shape.py
    """

def
