"""This is a sample python file for testing functions from the source code."""

from __future__ import annotations

import os
import shutil
import stat
from pathlib import Path

import pytest
import tests.expected_results
from canhydro.Forester import Forester
from canhydro.utils import concave_hull, read_file_names, save_file
from tests.expected_results import ten_cyls_rows
from canhydro.global_vars import (
    DIR, 
    test_input_dir
)

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

expected_result = {}


def test_create_cyliders(basic_forest):
    actual = basic_forest.get_collection_data("1_TenCyls.csv")
    expected = ten_cyls_rows
    breakpoint()
    assert expected == actual

def test_project_cyliders(basic_forest):
    collection = basic_forest.cylinder_collections[0]
    collection.project_cylinders('XY')
    breakpoint()
    assert 1==0
