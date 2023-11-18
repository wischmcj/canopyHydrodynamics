"""This is a sample python file for testing functions from the source code."""

from __future__ import annotations

import os
import shutil
import stat
from pathlib import Path

import pytest

import sys
import pandas as pd
import geopandas as geo
import numpy as np
from descartes import PolygonPatch
from shapely.geometry import Point
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.getcwd()))
import alphashape



from canhydro.global_vars import DIR, test_input_dir
from canhydro.utils import lam_filter, read_file_names
from canhydro.geometry import concave_hull, furthest_point, union, draw_cyls
from tests.expected_results import (ez_projection_xy_angle, hp_edges,
                                    ten_cyls_bo_and_rad, ten_cyls_bo_one,
                                    ten_cyls_edges, ten_cyls_id_one,
                                    ten_cyls_is_stem, ten_cyls_rows,
                                    small_tree_wateshed_poly)

import alphashape

# DIR = DIR
# test_input_dir = test_input_dir


# def on_rm_error(func, path, exc_info):
#     # path contains the path of the file that couldn't be removed
#     # let's just assume that it's read-only and unlink it.
#     os.chmod(path, stat.S_IWRITE)
#     os.unlink(path)


# def create_dir_and_file(filename) -> None:
#     print(type(filename))
#     os.makedirs(filename, exist_ok=True)
#     f = open("test\\demofile2.csv", "w")
#     f.write("Now the file has more content!")
#     f.close()


# def test_file_names():
#     file_path = "".join([DIR, "test\\"])
#     file = os.path.dirname(file_path)
#     create_dir_and_file(file)
#     file_names = read_file_names(Path("".join([DIR, "test"])))
#     assert file_names == ["demofile2.csv"]
#     shutil.rmtree(file, onerror=on_rm_error)
#     print("File Names Successful")


# def within_range(expected, actual, err):
#     return actual > expected - expected * err and actual < expected + expected * err


# # Needs tested for various filters, as well as for XZ, YZ
# @pytest.mark.parametrize("flexible_collection", ["2_EZ_projection.csv"], indirect=True)
# def test_project_cylinders(flexible_collection, accepted_err=0.03):
#     flexible_collection.project_cylinders(plane="XZ")
#     actual = flexible_collection.cylinders[0].projected_data["XZ"]["angle"]
#     actual = flexible_collection.cylinders[0].angle
#     expected = ez_projection_xy_angle
#     assert within_range(expected, actual, accepted_err)


# def test_split(self):
#     s = "hello world"
#     """This is an example of a type error test"""
#     pytest.assertEqual(s.split(), ["hello", "world"])
#     # check that s.split fails when the separator is not a string
#     with self.assertRaises(TypeError):
#         s.split(2)

def test_zero_curve_alpha():
    points_2d = sorted([(0., 0.), (0., 1.), (1., 1.), (1., 0.),
        (0.5, 0.25), (0.5, 0.75), (0.25, 0.5), (0.75, 0.5)])
    alpha_shape =  alphashape.alphashape(points_2d, 2.0)
    boundary = [Point((x,y)) for (x,y) in points_2d]
    
    hull, _ = concave_hull(boundary,2.0)
    # fig, ax = plt.subplots()
    # geopoly = geo.GeoSeries(hull)
    # geopoly.plot(ax =ax)
    # plt.show()
    assert alpha_shape == hull
    breakpoint()

# @pytest.mark.parametrize(
#     "flexible_collection", ["5_SmallTree.csv"], indirect=True
# )
# def test_small_tree(flexible_collection):
#     flexible_collection.project_cylinders(plane="XY")
#     flexible_collection.initialize_graph()
#     flexible_collection.watershed_boundary(flexible_collection.graph, plane = 'XY')
#     expected = small_tree_wateshed_poly
#     actual = flexible_collection.hull
#     assert str(actual) == expected

# @pytest.mark.parametrize(
#     "flexible_collection", ["3_HappyPathWTrunk.csv"], indirect=True
# )
# def test_find_flows(flexible_collection):
#     flexible_collection.project_cylinders("XZ")
#     flexible_collection.initialize_graph()
#     flexible_collection.find_flow_components()
#     _, stem_bool = lam_filter(
#         flexible_collection.cylinders, lambda: is_stem, return_all=True
#     )
#     breakpoint()
#     assert stem_bool == ten_cyls_is_stem


# @pytest.mark.parametrize("flexible_collection", ["1_TenCyls.csv"], indirect=True)
# def test_lam_filter(flexible_collection):
#     bo_one = lam_filter(flexible_collection.cylinders, lambda: branch_order == 1)
#     bo_zero = lam_filter(
#         flexible_collection.cylinders, lambda: branch_order == 0 or length <= 0.22447
#     )
#     id_one = lam_filter(flexible_collection.cylinders, lambda: cyl_id == 1)
#     assert ten_cyls_bo_one == str(bo_one)
#     assert ten_cyls_bo_and_rad == str(bo_zero)
#     assert ten_cyls_id_one == str(id_one)


# @pytest.mark.parametrize("flexible_collection", ["1_TenCyls.csv"], indirect=True)
# def test_create_cylinders(flexible_collection):
#     actual = flexible_collection.get_collection_data()
#     expected = ten_cyls_rows
#     assert expected == actual


# @pytest.mark.parametrize("flexible_collection", ["1_TenCyls.csv"], indirect=True)
# def test_highlight_filt_draw(flexible_collection, accepted_err=0.03):
#     # flexible_collection.project_cylinders(plane="XZ")
#     # flexible_collection.draw('XZ')
#     # flexible_collection.draw('XZ', a_lambda = lambda: branch_order ==1)
#     # flexible_collection.draw('XZ', a_lambda = lambda: branch_order ==1, highlight = True)

#     assert 1 == 1


# @pytest.mark.parametrize("flexible_collection", ["1_TenCyls.csv"], indirect=True)
# def test_create_line_graph(flexible_collection):
#     flexible_collection.initialize_graph()
#     actual_edges = [edge for edge in flexible_collection.graph.edges]
#     assert actual_edges == ten_cyls_edges


# def test_create_happy_path_graph(happy_path_projection):
#     happy_path_projection.initialize_graph()
#     actual_edges = [edge for edge in happy_path_projection.graph.edges]
#     expected_edges = hp_edges
#     assert actual_edges == expected_edges


# def test_dbh():
#     assert 1 == 1


    
# @pytest.mark.parametrize("flexible_collection", ["5_SmallTree.csv"], indirect=True)
# def test_watershed(flexible_collection):
#     flexible_collection.initialize_graph()
#     breakpoint()
#     flexible_collection.watershed_boundary()
#     assert actual_edges == expected_edges
