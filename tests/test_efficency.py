"""This is a sample python file for testing functions from the source code."""

from __future__ import annotations

import os
import shutil
import stat
from pathlib import Path

import pytest

from canhydro.Forester import Forester
from canhydro.CylinderCollection import CylinderCollection
from canhydro.global_vars import test_input_dir

test_input_dir = test_input_dir

# initialize_graph
# initialize_object_graph
# initialize_minimal_graph
# sum_over_graph
# sum_over_object_graph
# sum_over_min_graph
from memory_profiler import profile

# @pytest.mark.report_uss
# @pytest.mark.report_tracemalloc
# @pytest.mark.report_duration

# def test_base_efficency(large_collection):  
#     large_collection.initialize_graph()
#     proj_area = large_collection.sum_over_graph()
#     print(str(proj_area))
#     assert 1==1

# # @pytest.mark.report_uss
# # @pytest.mark.report_tracemalloc
# # @pytest.mark.report_duration
# def test_object_efficency(large_collection):  
#     large_collection.initialize_object_graph()
#     proj_area = large_collection.sum_over_object_graph()
#     print(str(proj_area))
#     assert 1==1

# # @pytest.mark.report_uss
# # @pytest.mark.report_tracemalloc
# # @pytest.mark.report_duration
# def test_min_efficency(large_collection):  
#     large_collection.initialize_minimal_graph()
#     proj_area = large_collection.sum_over_min_graph()
#     print(str(proj_area))
#     assert 1==1