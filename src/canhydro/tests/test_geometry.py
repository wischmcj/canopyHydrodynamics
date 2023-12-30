from __future__ import annotations

import pytest
from shapely.geometry import Polygon

import canhydro.geometry as geometry
from tests.expected_results_shapes import (funky_squares_overlap_areas,
                                           squares_projection_overlap)
from tests.utils import within_range


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
    actual = geometry.numba_get_projection(vector, magnitude, radius)
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

# def test_given_a_point_return_a_point():
#     """
#     Given a point, the alphashape function should return the same point
#     """
#     point_a, _ = geometry.concave_hull([(0.0, 0.0)], 0)
#     point_c, _ = geometry.concave_hull([(1.0, 0.0)], 0)
#     point_b, _ = geometry.concave_hull([(0.0, 1.0)], 0)
#     point_d, _ = geometry.concave_hull([(0.0, 0.0)], 99)
#     point_e, _ = geometry.concave_hull([(1.0, 0.0)], 99)
#     point_f, _ = geometry.concave_hull([(0.0, 1.0)], 99)
#     assert Point([0.0, 0.0]) == point_a
#     assert Point([1.0, 0.0]) == point_c
#     assert Point([0.0, 1.0]) == point_b
#     assert Point([0.0, 0.0]) == point_d
#     assert Point([1.0, 0.0]) == point_e
#     assert Point([0.0, 1.0]) == point_f

# def test_given_a_line_with_dupicate_points_return_a_point():
#     """
#     Given a line with duplicate points, the alphashape function should
#     return a point
#     """
#     actual, _ = geometry.concave_hull([(0.0, 1.0), (0.0, 1.0)], 0)
#     assert Point([0.0, 1.0]) == actual
# def test_given_a_line_with_unique_points_return_a_line():
#     """
#     Given a line with unique points, the alphashape function should return
#     the same line
#     """
#     actual_horz, _ = geometry.concave_hull([(0.0, 0.0), (0.0, 1.0)], 0)
#     actual_diag, _ = geometry.concave_hull([(1.0, 0.0), (0.0, 1.0)], 0)
#     assert LineString([(0.0, 0.0), (0.0, 1.0)]) == actual_horz
#     assert LineString([(1.0, 0.0), (0.0, 1.0)]) == actual_diag
# def test_given_a_triangle_with_duplicate_points_returns_a_point():
#     """
#     Given a triangle with two unique points, the alphashape function should
#     return a point
#     """
#     dupe_tri, _ = geometry.concave_hull([(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)], 0)
#     assert Point((0.0, 1.0)) == dupe_tri
# def test_given_a_triangle_with_two_duplicate_points_returns_a_line():
#     """
#     Given a line with two unique points, the alphashape function should
#     return a line with the unique points
#     """
#     tri_two_dupe, _ = geometry.concave_hull([(1.0, 0.0), (0.0, 1.0), (0.0, 1.0)], 0)
#     assert LineString([(1.0, 0.0), (0.0, 1.0)]) == tri_two_dupe
# # def test_given_a_four_point_polygon_with_small_alpha_return_input():
# #     """
# #     Given a polygon with four points, and an alpha value of zero, return
# #     the input as a polygon.
# #     """
# #     pts = [Point((x,y)) for (x,y) in  [(0., 0.), (0., 1.), (1., 1.), (1., 0.)]]
# #     quad_small, _= concave_hull(pts, .0001)
# #     assert shapely.geometry.Polygon([(0., 0.), (0., 1.), (1., 1.), (1., 0.), (0., 0.)]).equals(quad_small)
# @pytest.mark.parametrize("flexible_collection", ["5_SmallTree.csv"], indirect=True)
# def test_small_tree(flexible_collection):
#     flexible_collection.project_cylinders(plane="XY")
#     flexible_collection.initialize_graph()
#     flexible_collection.watershed_boundary(flexible_collection.graph, plane="XY")
#     expected = small_tree_wateshed_poly
#     actual = flexible_collection.hull
#     assert str(actual) == expected

# def test_star_poly():
#     """
#     Given a 3-dimensional data set, return an expected set of edges.
#     """
#     coords_2d = sorted(
#         [
#             (0.0, 0.0),
#             (0.0, 0.0),
#             (0.0, 1.0),
#             (1.0, 0.0),
#             (1.0, 1.0),
#             (1.0, 0.0),
#             (0.0, 1.0),
#             (1.0, 1.0),
#             (0.25, 0.5),
#             (0.5, 0.25),
#             (0.5, 0.5),
#             (0.75, 0.5),
#             (0.5, 0.75),
#             (0.5, 0.5),
#         ]
#     )
#     points_2d = [Point((x, y)) for (x, y) in coords_2d]
#     expected_vertices = [
#         [0.0, 0.0],
#         [0.0, 0.0],
#         [0.0, 1.0],
#         [0.0, 1.0],
#         [1.0, 0.0],
#         [1.0, 0.0],
#         [1.0, 1.0],
#         [1.0, 1.0],
#         [0.25, 0.5],
#         [0.5, 0.25],
#         [0.5, 0.5],
#         [0.5, 0.5],
#         [0.5, 0.75],
#         [0.75, 0.5],
#     ]
#     expected = star_poly
#     results, _, _ = geometry.concave_hull(points_2d, 2.1)
#     assert results.contains(expected)
#     assert expected.contains(results)
#     # breakpoint()
# def test_star_poly():
#     """
#     Given a 3-dimensional data set, return an expected set of edges.
#     """
#     coords_2d = sorted(
#         [
#             (0.0, 0.0),
#             (0.0, 0.0),
#             (0.0, 1.0),
#             (1.0, 0.0),
#             (1.0, 1.0),
#             (1.0, 0.0),
#             (0.0, 1.0),
#             (1.0, 1.0),
#             (0.25, 0.5),
#             (0.5, 0.25),
#             (0.5, 0.5),
#             (0.75, 0.5),
#             (0.5, 0.75),
#             (0.5, 0.5),
#         ]
#     )
#     points_2d = [Point((x, y)) for (x, y) in coords_2d]
#     expected_vertices = [
#         [0.0, 0.0],
#         [0.0, 0.0],
#         [0.0, 1.0],
#         [0.0, 1.0],
#         [1.0, 0.0],
#         [1.0, 0.0],
#         [1.0, 1.0],
#         [1.0, 1.0],
#         [0.25, 0.5],
#         [0.5, 0.25],
#         [0.5, 0.5],
#         [0.5, 0.5],
#         [0.5, 0.75],
#         [0.75, 0.5],
#     ]
#     expected = star_poly
#     results, _, _ = geometry.concave_hull(points_2d, 2.1)
#     assert results.contains(expected)
#     assert expected.contains(results)
#     # breakpoint()
# def test_zero_curve_alpha():
#     points_2d = sorted(
#         [
#             (0.0, 0.0),
#             (0.0, 1.0),
#             (1.0, 1.0),
#             (1.0, 0.0),
#             (0.5, 0.25),
#             (0.5, 0.75),
#             (0.25, 0.5),
#             (0.75, 0.5),
#         ]
#     )
#     alpha_shape = alphashape.alphashape(points_2d, 2.0)
#     boundary = [Point(x, y) for (x, y) in points_2d]
#     hull, _, _ = geometry.concave_hull(boundary, 2.0)
#     # fig, ax = plt.subplots()
#     # geopoly = geo.GeoSeries(hull)
#     # geopoly.plot(ax =ax)
#     # plt.show()
#     assert alpha_shape == hull
# # def test_3_dimensional_regression():
# #     """
# #     Given a 3-dimensional data set, return an expected set of edges.
# #     """
# #     points_3d = [
# #         (0.0, 0.0, 0.0),
# #         (0.0, 0.0, 1.0),
# #         (0.0, 1.0, 0.0),
# #         (1.0, 0.0, 0.0),
# #         (1.0, 1.0, 0.0),
# #         (1.0, 0.0, 1.0),
# #         (0.0, 1.0, 1.0),
# #         (1.0, 1.0, 1.0),
# #         (0.25, 0.5, 0.5),
# #         (0.5, 0.25, 0.5),
# #         (0.5, 0.5, 0.25),
# #         (0.75, 0.5, 0.5),
# #         (0.5, 0.75, 0.5),
# #         (0.5, 0.5, 0.75),
# #     ]
# #     expected = {
# #     }
# #     expected_vertices = [
# #         [0., 0., 0.], [0., 0., 1.], [0., 1., 0.],
# #         [0., 1., 1.], [1., 0., 0.], [1., 0., 1.],
# #         [1., 1., 0.], [1., 1., 1.], [0.25, 0.5, 0.5],
# #         [0.5, 0.25, 0.5], [0.5, 0.5, 0.25], [0.5, 0.5, 0.75],
# #         [0.5, 0.75, 0.5], [0.75, 0.5, 0.5]]
# #     expected_faces = [
# #         (13, 10, 6), (13, 9, 4), (6, 12, 13),
# #         (13, 12, 7), (5, 11, 9), (8, 10, 0),
# #         (3, 12, 8), (0, 10, 9), (5, 9, 13),
# #         (12, 11, 7), (9, 10, 4), (8, 9, 1),
# #         (12, 10, 2), (13, 11, 5), (1, 11, 8),
# #         (4, 10, 13), (9, 11, 1), (2, 10, 8),
# #         (8, 12, 2), (3, 11, 12), (0, 9, 8),
# #         (7, 11, 13), (6, 10, 12), (8, 11, 3)]
# #     results = alphashape(points_3d, 2.1)
# #     assertTrue(len(results.vertices) == len(expected_vertices))
# #     assertTrue(len(points_3d) == len(expected_vertices))
# #     assertTrue(len(results.faces) == len(expected_faces))
# #     vertex_map = {i: expected_vertices.index(
# #         list(vertex)) for i,vertex in enumerate(results.vertices)}
# #     for edge in list(results.faces):
# #        assertTrue(any([(
# #            vertex_map[e[0]],
# #            vertex_map[e[1]],
# #            vertex_map[e[2]]) in expected_faces \
# #                 for e in itertools.combinations(edge, r=len(edge))]))
# # def test_4_dimensional_regression():
# #     """
# #     Given a 4-dimensional data set, return an expected set of edges.
# #     """
# #     points_4d = [
# #        (0., 0., 0., 0.), (0., 0., 0., 1.), (0., 0., 1., 0.),
# #        (0., 1., 0., 0.), (0., 1., 1., 0.), (0., 1., 0., 1.),
# #        (0., 0., 1., 1.), (0., 1., 1., 1.), (1., 0., 0., 0.),
# #        (1., 0., 0., 1.), (1., 0., 1., 0.), (1., 1., 0., 0.),
# #        (1., 1., 1., 0.), (1., 1., 0., 1.), (1., 0., 1., 1.),
# #        (1., 1., 1., 1.), (.25, .5, .5, .5), (.5, .25, .5, .5),
# #        (.5, .5, .25, .5), (.5, .5, .5, .25), (.75, .5, .5, .5),
# #        (.5, .75, .5, .5), (.5, .5, .75, .5), (.5, .5, .5, .75)
# #     ]
# #     expected = {
# #         (16, 1, 2, 0), (16, 1, 3, 0), (16, 2, 3, 0),
# #         (16, 4, 2, 3), (16, 4, 7, 2), (16, 4, 7, 3),
# #         (16, 5, 1, 3), (16, 5, 7, 1), (16, 5, 7, 3),
# #         (16, 6, 1, 2), (16, 6, 7, 1), (16, 6, 7, 2),
# #         (17, 1, 2, 0), (17, 1, 8, 0), (17, 2, 8, 0),
# #         (17, 6, 1, 2), (17, 6, 14, 1), (17, 6, 14, 2),
# #         (17, 9, 1, 8), (17, 9, 14, 1), (17, 9, 14, 8),
# #         (17, 10, 2, 8), (17, 10, 14, 2), (17, 10, 14, 8),
# #         (18, 1, 3, 0), (18, 1, 8, 0), (18, 3, 8, 0),
# #         (18, 5, 1, 3), (18, 5, 13, 1), (18, 5, 13, 3),
# #         (18, 9, 1, 8), (18, 9, 13, 1), (18, 9, 13, 8),
# #         (18, 11, 3, 8), (18, 11, 13, 3), (18, 11, 13, 8),
# #         (19, 2, 3, 0), (19, 2, 8, 0), (19, 3, 8, 0),
# #         (19, 4, 2, 3), (19, 4, 12, 2), (19, 4, 12, 3),
# #         (19, 10, 2, 8), (19, 10, 12, 2), (19, 10, 12, 8),
# #         (19, 11, 3, 8), (19, 11, 12, 3), (19, 11, 12, 8),
# #         (20, 9, 13, 8), (20, 9, 14, 8), (20, 9, 14, 13),
# #         (20, 10, 12, 8), (20, 10, 14, 8), (20, 10, 14, 12),
# #         (20, 11, 12, 8), (20, 11, 13, 8), (20, 11, 13, 12),
# #         (20, 13, 12, 15), (20, 14, 12, 15), (20, 14, 13, 15),
# #         (21, 4, 7, 3), (21, 4, 7, 12), (21, 4, 12, 3),
# #         (21, 5, 7, 3), (21, 5, 7, 13), (21, 5, 13, 3),
# #         (21, 7, 12, 15), (21, 7, 13, 15), (21, 11, 12, 3),
# #         (21, 11, 13, 3), (21, 11, 13, 12), (21, 13, 12, 15),
# #         (22, 4, 7, 2), (22, 4, 7, 12), (22, 4, 12, 2),
# #         (22, 6, 7, 2), (22, 6, 7, 14), (22, 6, 14, 2),
# #         (22, 7, 12, 15), (22, 7, 14, 15), (22, 10, 12, 2),
# #         (22, 10, 14, 2), (22, 10, 14, 12), (22, 14, 12, 15),
# #         (23, 5, 7, 1), (23, 5, 7, 13), (23, 5, 13, 1),
# #         (23, 6, 7, 1), (23, 6, 7, 14), (23, 6, 14, 1),
# #         (23, 7, 13, 15), (23, 7, 14, 15), (23, 9, 13, 1),
# #         (23, 9, 14, 1), (23, 9, 14, 13), (23, 14, 13, 15)}
# #     results = alphashape(points_4d, 1.0)
# #     assertTrue(len(results) == len(expected))
# #     for edge in list(results):
# #        assertTrue(any([e in expected for e in itertools.combinations(
# #             edge, r=len(edge))]))

# def test_voronoi():
#     """
#         Placeholder to test this functionality when it is built
#     """

# @pytest.mark.parametrize(
#     "point,points,num, expected",
#     [((2,2), example_coords,1,(1,2)),
#      ((9,7),start_poly_coords, 1,(1,1))
#      ]
# )
# def test_closest_point(point,points, num,expected):
#     actual = geometry.closest_points(point,points, num)
#     assert actual ==expected


# def simplicies():
#     """
#         probs would need a mock to do this one?
#           maybe better for integtation? maybe skip for now...
#     """

# def maximal_alpha():
#     """
#         maybe find a diff package to test against
#         not very important atm
#     """
