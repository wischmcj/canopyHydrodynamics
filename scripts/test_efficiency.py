"""This is a sample python file for testing functions from the source code."""

from __future__ import annotations
import toml

from line_profiler import profile  # noqa

with open("src/canhydro/user_def_config.toml") as f:
    config = toml.load(f)
    test_input_dir = config["directories"]['test_input_dir']


# @pytest.mark.parametrize(
#     "flexible_collection", ["4_LargeCollection.csv"], indirect=True
# )
# def test_base_graph_test(flexible_collection):
#     flexible_collection.initialize_graph_from()
#     proj_area = flexible_collection.sum_over_graph()
#     flexible_collection.find_flow_components()
#     flexible_collection.calculate_flows()


# @pytest.mark.parametrize(
#     "flexible_collection", ["4_LargeCollection.csv"], indirect=True
# )
# def test_min_graph_test(flexible_collection):
#     flexible_collection.initialize_minimal_graph_from()
#     proj_area = flexible_collection.sum_over_min_graph()
#     flexible_collection.find_flow_components_minimal()
#     flexible_collection.calculate_flows_min()


# @pytest.mark.parametrize(
#     "flexible_collection", ["4_LargeCollection.csv"], indirect=True
# )
# # @pytest.mark.parametrize(
# #     "flexible_collection", ["3_HappyPathWTrunk.csv"], indirect=True
# # )
# def test_obj_graph_test(flexible_collection):
#     flexible_collection.initialize_object_graph_from()
#     proj_area = flexible_collection.sum_over_object_graph()
#     flexible_collection.find_flow_components_object()
#     flexible_collection.calculate_flows_obj()

# #@profile
# def test_small_tree_proj(small_tree, ez_projection):
#     ez_projection.numba_project_cylinders("XY")
#     small_tree.project_cylinders("XY")
#     assert 1 == 1


# #@profile
# def test_numba_small_tree_proj(small_tree, ez_projection):
#     ez_projection.numba_project_cylinders("XY")
#     small_tree.numba_project_cylinders("XY")
#     assert 1 == 1


# #@profile
# def test_happy_proj(happy_path_projection, ez_projection):
#     ez_projection.numba_project_cylinders("XY")
#     happy_path_projection.project_cylinders("XZ")
#     assert 1 == 1


# #@profile
# def test_numba_happy_proj(happy_path_projection, ez_projection):
#     ez_projection.numba_project_cylinders("XY")
#     happy_path_projection.numba_project_cylinders("XZ")
#     assert 1 == 1


# #@profile
# def test_large_proj(large_collection, ez_projection):
#     ez_projection.numba_project_cylinders("XY")
#     large_collection.project_cylinders("XZ")
#     assert 1 == 1


# #@profile
# def test_numba_large_proj(large_collection, ez_projection):
#     ez_projection.numba_project_cylinders("XY")
#     large_collection.numba_project_cylinders("XZ")
#     assert 1 == 1

# @pytest.mark.parametrize(
#     "flexible_collection", ["4_LargeCollection.csv"], indirect=True
# )
# # @pytest.mark.parametrize(
# #     "flexible_collection", ["3_HappyPathWTrunk.csv"], indirect=True
# # )
# def test_obj_graph_test(flexible_collection):
#     flexible_collection.initialize_object_graph_from()
#     proj_area = flexible_collection.sum_over_object_graph()
#     flexible_collection.find_flow_components_object()
#     flexible_collection.calculate_flows_obj()

# @pytest.mark.report_uss
# @pytest.mark.report_tracemalloc
# @pytest.mark.report_duration

# def test_base_efficiency(large_collection):
#     large_collection.initialize_graph_from()
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
