"""This is a sample python file for testing functions from the source code."""

from __future__ import annotations

import os
import shutil
import stat
from pathlib import Path

import pytest

import tests.expected_results
from canhydro.Forester import Forester
from canhydro.global_vars import DIR, test_input_dir
from canhydro.utils import concave_hull, read_file_names, save_file
from tests.expected_results import ten_cyls_rows, ten_cyls_edges, ez_projection_xy_angle, hp_edges

DIR = DIR
test_input_dir = test_input_dir


def on_rm_error(func, path, exc_info):
    # path contains the path of the file that couldn't be removed
    # let's just assume that it's read-only and unlink it.
    os.chmod(path, stat.S_IWRITE)
    os.unlink(path)


def create_dir_and_file(filename) -> None:
    print(type(filename))
    os.makedirs(filename, exist_ok=True)
    f = open("test\\demofile2.csv", "w")
    f.write("Now the file has more content!")
    f.close()


def test_file_names():
    file_path = "".join([DIR, "test\\"])
    file = os.path.dirname(file_path)
    create_dir_and_file(file)
    file_names = read_file_names(Path("".join([DIR, "test"])))
    assert file_names == ["demofile2.csv"]
    shutil.rmtree(file, onerror=on_rm_error)
    print("File Names Successfull")


# def test_split(self):
#     s = "hello world"
#     """This is an example of a type error test"""
#     pytest.assertEqual(s.split(), ["hello", "world"])
#     # check that s.split fails when the separator is not a string
#     with self.assertRaises(TypeError):
#         s.split(2)

# expected_result = {}


# def test_create_cyliders(basic_forest):
#     ten_cyls = basic_forest.cylinder_collections[0]
#     actual = ten_cyls.get_collection_data()
#     expected = ten_cyls_rows
#     assert expected == actual

def within_range(expected,actual,err):
    return actual > expected - expected*err and actual < expected + expected*err


# Needs tested for various filters, as well as for XZ, YZ
def test_project_cyliders(ez_projection, accepted_err = 0.03):
    ez_projection.project_cylinders(plane="XY")
    actual = ez_projection.cylinders[0].projected_data["XY"]["angle"]
    expected = ez_projection_xy_angle
    assert within_range(expected,actual,accepted_err) 

def test_create_line_graph(ten_cyls_col):
    ten_cyls_col.initialize_graph()
    actual_edges = [edge for edge in ten_cyls_col.graph.edges]
    assert actual_edges == ten_cyls_edges

    
def test_create_happy_path_graph(happy_path_projection):
    happy_path_projection.initialize_graph()
    actual_edges = [edge for edge in happy_path_projection.graph.edges]
    expected_edges = hp_edges
    breakpoint()
    assert actual_edges == expected_edges