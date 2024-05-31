from __future__ import annotations

import shutil
import toml
import os
from pathlib import Path

import pytest 
from numpy import array as arr, all
from src.canhydro.utils import create_dir_and_file, read_file_names, on_rm_error, stack, njit_stack


with open("src/canhydro/user_def_config.toml") as f:
    config = toml.load(f)
    DIR = config["directories"]['root_dir']

njit_stack_test_cases = [
    pytest.param(arr([arr([1,2,3]),arr([4,5,6])]),False ,arr([[1., 2., 3.],[4., 5., 6.]]), id = "3x3 row stack"),
    pytest.param(arr([arr([1,2,3]),arr([4,5,6]),arr([7,8,9])]),False ,arr([[1., 2., 3.],[4., 5., 6.],[7., 8., 9.]]), id = "3x3x3 row stack"),
    pytest.param(arr([arr([1,2]),arr([4,5])]),False ,arr([[1., 2.],[4., 5.]]), id = "2x2 row stack"),
    pytest.param(arr([arr([1,2,3]),arr([4,5,6])]),True ,arr([[1., 4.],[2., 5.],[3., 6.]]), id = "3x3 col stack"),
    pytest.param(arr([arr([1,2]),arr([4,5])]),True, arr([[1., 4.],[2., 5.]]), id = "2x2 col stack"),
]

stack_test_cases = [
    pytest.param([arr([1]),arr([4,5]),arr([4,5])], False, ValueError, id = "1x2x2 ValueError"),
    pytest.param([arr([1]),arr([4,5])], False, ValueError, id = "1x2 ValueError"),
]

def test_file_funcs():
    file_path = "".join([DIR, "test_data/"])
    file = os.path.dirname(file_path)
    create_dir_and_file(file)
    f = open("test_data/demofile2.csv", "w")
    f.write("Now the file has more content!")
    f.close()
    file_names = read_file_names(Path("".join([DIR, "test_data"])))
    assert file_names == ["demofile2.csv"]
    shutil.rmtree(file, onerror=on_rm_error)
    print("File Names Successful")

@pytest.mark.parametrize(
    "arr_list,is_col_stack,expected",
    njit_stack_test_cases
)
def test_njit_stack(arr_list,is_col_stack,expected):
    """
    Test vstack and col stack capabilities v. numpy
    thow in a couple exenmplar custom tests to show functionality
    """
    actual = njit_stack(arr_list,is_col_stack)
    assert all(actual == expected)


@pytest.mark.parametrize(
    "arr_list,is_col_stack,exception_type",
    stack_test_cases
)
def test_stack_err(arr_list,is_col_stack,exception_type):
    """
    Test vstack and col stack capabilities v. numpy
    thow in a couple exenmplar custom tests to show functionality
    """
    with pytest.raises(exception_type):
        actual = stack(arr_list,is_col_stack)
