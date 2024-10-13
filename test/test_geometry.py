from __future__ import annotations

from test.expected_results_shapes import (
    example_coords,
    funky_squares_overlap_areas,
    squares_projection_overlap,
    star_poly,
    start_poly_coords,
)
from test.utils import within_range

import pytest
from shapely.geometry import LineString, Point, Polygon

import src.canopyhydro.geometry as geometry


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
    actual_overlap_dict = geometry.get_projected_overlap(poly_matrix, labels)
    assert squares_projection_overlap == actual_overlap_dict


def test_proj_overlap():
    """
    Tests for diff number of  (outer list) and diff number of polys (outer list)
    """
    square1 = Polygon([(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)])
    square2 = Polygon([(0.25, 0.25), (0.25, 0.5), (0.5, 0.5), (0.5, 0.25)])
    square3 = Polygon([(1.0, 1.0), (1.0, 2.0), (2.0, 2.0), (2.0, 1.0)])
    square4 = Polygon([(1.5, 1.5), (1.5, 2.5), (2.5, 2.5), (2.5, 1.5)])
    bottom = [square1, square2]
    mid = [square1, square2, square3]
    top = [square1, square2, square3, square4]
    poly_matrix = [bottom, mid, top]
    actual = geometry.get_projected_overlap(poly_matrix, ["bottom", "mid", "top"])
    assert actual == funky_squares_overlap_areas


@pytest.mark.parametrize(
    "vector,magnitude,radius,e_angle,e_area",
    [
        ([(2, -1), (4, -1), (5, -1)], [-3, -5, -6], 1, -0.785398, 13.91),
        ([(1, -0.5), (1, -1.5), (1, -3)], [1.5, 2.5, 3], 1, 0.785398, 8.08),
        ([(-2, 1), (-4, 1), (-5, 1)], [3, 5, 6], 1, 0.785398, 13.91),
        ([(1, 4), (1, 6), (1, 7)], [3, 5, 6], 1, 0.785398, 13.91),
        ([(-1, -4), (-1, -6), (-1, -7)], [3, 5, 6], 1, 0.785398, 13.91),
        ([(1, -2), (1, -4), (1, -5)], [3, 5, 6], 1, 0.785398, 13.91),
    ],
)
def test_project_cyl(vector, magnitude, radius, e_angle, e_area):
    """
    Tests projection of cylinders parallel with:
    the XY plane, Z axis the line x=y (45 deg)
    """
    actual = geometry.get_projection(vector, magnitude, radius)
    assert within_range(e_angle, actual["angle"], 0.02)
    assert within_range(e_area, actual["area"], 0.02)


@pytest.mark.parametrize(
    "points,expected",
    [
        ([(1, 1), (-1, 1), (1, -1)], [0.0, 0.5, 0.5]),
        ([(0.5, 0.5), (0, 0), (0.5, 0)], [0.5, 0.5, 0.0]),
        ([(1, 0), (0, 0), (0, 1)], [0.5, 0.0, 0.5]),
    ],
)
def test_lu_factor(points, expected):
    """
    tests circumcenter func with some mock examples
    (one, two, many points, size zero array, cooincindent coords)
    """
    actual = geometry.circumcenter_lapack(points)
    for idx, exp in enumerate(expected):
        assert within_range(exp, actual[idx], 0.02)


@pytest.mark.parametrize(
    "points",
    [
        ([(1, 1), (-1, 1), (1, -1)]),
        ([(0.5, 0.5), (0, 0), (0.5, 0)]),
        ([(1, 0), (0, 0), (0, 1)]),
    ],
)
def test_lapack_v_lu_factor(points):
    """
    Tests thse answers are the same so that we only really need
    to test one
    """
    actual = geometry.circumcenter_lu_factor(points)
    expected = geometry.circumcenter_lapack(points)
    for idx, exp in enumerate(expected):
        assert within_range(exp, actual[idx], 0.02)


@pytest.mark.parametrize(
    "points,expected",
    [
        ([(1, 1), (-1, 1), (1, -1)], 1.41),
        ([(5, 6), (4, 4), (5, 4)], 1.11),
        ([(0.5, 0.5), (0, 0), (0.5, 0)], 0.35),
    ],
)
def test_curcumradius(points, expected):
    """
    test with input center, and without
    corrd list of one two and many lengths (div by three and non div by three)
    test with 3+ dims
    """
    center = geometry.circumcenter_lapack(points)
    actual = geometry.circumradius(points, center)
    assert within_range(actual, expected, 0.02)


# Alpha Shape tests below
# Credit for the coordinates utilized in these tests goes to
# user 'bellokk' on git hub (https://github.com/bellockk/alphashape/tree/master)


def test_given_a_point_return_a_point():
    """
    Given a point, the alphashape function should return the same point
    """
    point_a, _ = geometry.concave_hull([(0.0, 0.0)], 0)
    point_c, _ = geometry.concave_hull([(1.0, 0.0)], 0)
    point_b, _ = geometry.concave_hull([(0.0, 1.0)], 0)
    point_d, _ = geometry.concave_hull([(0.0, 0.0)], 99)
    point_e, _ = geometry.concave_hull([(1.0, 0.0)], 99)
    point_f, _ = geometry.concave_hull([(0.0, 1.0)], 99)
    assert Point([0.0, 0.0]) == point_a
    assert Point([1.0, 0.0]) == point_c
    assert Point([0.0, 1.0]) == point_b
    assert Point([0.0, 0.0]) == point_d
    assert Point([1.0, 0.0]) == point_e
    assert Point([0.0, 1.0]) == point_f


def test_given_a_line_with_dupicate_points_return_a_point():
    """
    Given a line with duplicate points, the alphashape function should
    return a point
    """
    actual, _ = geometry.concave_hull([(0.0, 1.0), (0.0, 1.0)], 0)
    assert Point([0.0, 1.0]) == actual


def test_given_a_line_with_unique_points_return_a_line():
    """
    Given a line with unique points, the alphashape function should return
    the same line
    """
    actual_horz, _ = geometry.concave_hull([(0.0, 0.0), (0.0, 1.0)], 0)
    actual_diag, _ = geometry.concave_hull([(1.0, 0.0), (0.0, 1.0)], 0)
    assert LineString([(0.0, 0.0), (0.0, 1.0)]) == actual_horz
    assert LineString([(1.0, 0.0), (0.0, 1.0)]) == actual_diag


def test_given_a_triangle_with_duplicate_points_returns_a_point():
    """
    Given a triangle with two unique points, the alphashape function should
    return a point
    """
    dupe_tri, _ = geometry.concave_hull([(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)], 0)
    assert Point((0.0, 1.0)) == dupe_tri


def test_given_a_triangle_with_two_duplicate_points_returns_a_line():
    """
    Given a line with two unique points, the alphashape function should
    return a line with the unique points
    """
    tri_two_dupe, _ = geometry.concave_hull([(1.0, 0.0), (0.0, 1.0), (0.0, 1.0)], 0)
    assert LineString([(1.0, 0.0), (0.0, 1.0)]) == tri_two_dupe


def test_star_poly():
    """
    Given a 3-dimensional data set, return an expected set of edges.
    """
    coords_2d = sorted(
        [
            (0.0, 0.0),
            (0.0, 0.0),
            (0.0, 1.0),
            (1.0, 0.0),
            (1.0, 1.0),
            (1.0, 0.0),
            (0.0, 1.0),
            (1.0, 1.0),
            (0.25, 0.5),
            (0.5, 0.25),
            (0.5, 0.5),
            (0.75, 0.5),
            (0.5, 0.75),
            (0.5, 0.5),
        ]
    )
    points_2d = [Point((x, y)) for (x, y) in coords_2d]
    # expected_vertices = [
    #     [0.0, 0.0],
    #     [0.0, 0.0],
    #     [0.0, 1.0],
    #     [0.0, 1.0],
    #     [1.0, 0.0],
    #     [1.0, 0.0],
    #     [1.0, 1.0],
    #     [1.0, 1.0],
    #     [0.25, 0.5],
    #     [0.5, 0.25],
    #     [0.5, 0.5],
    #     [0.5, 0.5],
    #     [0.5, 0.75],
    #     [0.75, 0.5],
    # ]
    expected = star_poly
    results, _, _ = geometry.concave_hull(points_2d, 2.1)
    assert results.contains(expected)
    assert expected.contains(results)
    # breakpoint()


@pytest.mark.parametrize(
    "point,points,num, expected",
    [((2, 2), example_coords, 1, (1, 2)), ((9, 7), start_poly_coords, 1, (1, 1))],
)
def test_closest_point(point, points, num, expected):
    actual = geometry.closest_points(point, points, num)
    assert actual == expected
