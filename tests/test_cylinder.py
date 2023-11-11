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
from tests.expected_results import ten_cyls_rows

DIR = DIR
test_input_dir = test_input_dir

# def test_create_cyliders(basic_forest):
#     actual = basic_forest.get_collection_data("1_TenCyls.csv")
#     expected = ten_cyls_rows
#     breakpoint()
#     assert expected == actual


def test_project_cyliders(ten_cyls_col):
    ten_cyls_col.project_cylinders("XZ")
    # basic_forest.cylinder_collections[0].project_cylinders("XZ")
    # test = basic_forest.get_collection_data("1_TenCyls.csv")
    ten_cyls_col.draw_cyls(plane="XZ")
    assert 1 == 0
