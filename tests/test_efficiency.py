"""This is a sample python file for testing functions from the source code."""

from __future__ import annotations

import pytest
from canhydro.global_vars import test_input_dir

test_input_dir = test_input_dir

# initialize_graph
# initialize_object_graph
# initialize_minimal_graph
# sum_over_graph
# sum_over_object_graph
# sum_over_min_graph

@pytest.mark.parametrize(
    "flexible_collection", ["4_LargeCollection.csv"], indirect=True
)
def test_base_graph_test(flexible_collection):
    flexible_collection.initialize_graph()
    proj_area = flexible_collection.sum_over_graph()
    flexible_collection.find_flow_components()

@pytest.mark.parametrize(
    "flexible_collection", ["4_LargeCollection.csv"], indirect=True
)
def test_min_graph_test(flexible_collection):
    flexible_collection.initialize_minimal_graph()
    proj_area = flexible_collection.sum_over_min_graph()
    flexible_collection.find_flow_components_minimal()

@pytest.mark.parametrize(
    "flexible_collection", ["4_LargeCollection.csv"], indirect=True
)   
def test_obj_graph_test(flexible_collection):
    flexible_collection.initialize_object_graph()
    proj_area = flexible_collection.sum_over_object_graph()
    flexible_collection.find_flow_components_object()

# @pytest.mark.report_uss
# @pytest.mark.report_tracemalloc
# @pytest.mark.report_duration

# def test_base_efficiency(large_collection):
#     large_collection.initialize_graph()
#     proj_area = large_collection.sum_over_graph()
#     print(str(proj_area))
#     assert 1==1

# # @pytest.mark.report_uss
# # @pytest.mark.report_tracemalloc
# # @pytest.mark.report_duration
# def test_object_efficiency(large_collection):
#     large_collection.initialize_object_graph()
#     proj_area = large_collection.sum_over_object_graph()
#     print(str(proj_area))
#     assert 1==1

# # @pytest.mark.report_uss
# # @pytest.mark.report_tracemalloc
# # @pytest.mark.report_duration
# def test_min_efficiency(large_collection):
#     large_collection.initialize_minimal_graph()
#     proj_area = large_collection.sum_over_min_graph()
#     print(str(proj_area))
#     assert 1==1
