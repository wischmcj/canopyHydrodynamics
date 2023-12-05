"""This is a sample python file for testing functions from the source code."""

from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.getcwd()))


from canhydro.global_vars import DIR

DIR = DIR


# @pytest.mark.parametrize("flexible_collection", ["5_SmallTree.csv"], indirect=True)
# def test_sandbox(flexible_collection, happy_path_projection):
#     assert 1==1

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


# @pytest.mark.parametrize("flexible_collection", ["2_EZ_projection.csv"], indirect=True)
# def test_project_cylinders(flexible_collection, accepted_err=0.03):
#     flexible_collection.project_cylinders(plane="XZ")
#     actual = flexible_collection.cylinders[0].projected_data["XZ"]["angle"]
#     actual = flexible_collection.cylinders[0].angle
#     expected = ez_projection_xy_angle
#     assert within_range(expected, actual, accepted_err)


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
#     assert stem_bool == ten_cyls_is_stem


# @pytest.mark.parametrize("flexible_collection", ["1_TenCyls.csv"], indirect=True)
# def test_single_lam_filter(flexible_collection):
#     bo_one, _ = lam_filter(flexible_collection.cylinders, lambda: branch_order == 1)
#     bo_zero_and_rad, _ = lam_filter(
#         flexible_collection.cylinders, lambda: branch_order == 0 or length <= 0.22447
#     )
#     id_one, _ = lam_filter(flexible_collection.cylinders, lambda: cyl_id == 1)

#     bo_one_ids = [cyl.cyl_id for cyl in bo_one]
#     bo_zero_and_rad_ids = [cyl.cyl_id for cyl in bo_zero_and_rad]
#     id_one_ids = [cyl.cyl_id for cyl in id_one]
#     assert ten_cyls_bo_one == str(bo_one_ids)
#     assert ten_cyls_bo_and_rad == str(bo_zero_and_rad_ids)
#     assert ten_cyls_id_one == str(id_one_ids)


# @pytest.mark.parametrize("flexible_collection", ["1_TenCyls.csv"], indirect=True)
# def test_create_cylinders(flexible_collection):
#     actual = flexible_collection.get_collection_data()
#     expected = ten_cyls_rows
#     breakpoint()
#     assert expected == actual


# @pytest.mark.parametrize("flexible_collection", ["1_TenCyls.csv"], indirect=True)
# def test_highlight_filt_draw(flexible_collection, accepted_err=0.03):
#     flexible_collection.project_cylinders(plane="XZ")
#     flexible_collection.draw("XZ")
#     flexible_collection.draw("XZ", highlight_lambda=lambda: branch_order == 1)
#     flexible_collection.draw(
#         "XZ", highlight_lambda=lambda: branch_order == 1,
#         filter_lambda=lambda: cyl_id > 3
#     )
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


# @pytest.mark.parametrize("flexible_collection", ["5_SmallTree.csv"], indirect=True)
# def test_dbh(flexible_collection):
#     flexible_collection.get_dbh()
#     actual = flexible_collection.treeQualities["dbh"]
#     expected = small_tree_dbh
#     assert expected == actual


@pytest.mark.parametrize("flexible_collection", ["5_SmallTree.csv"], indirect=True)
def test_drip_map(flexible_collection):
    flexible_collection.find_flow_components()
    flexible_collection.calculate_flows()
    flexible_collection.drip_map()

    # assert expected == str(actual)


# @pytest.mark.parametrize("flexible_collection", ["5_SmallTree.csv"], indirect=True)
# def test_watershed(flexible_collection):
#     flexible_collection.initialize_graph()
#     flexible_collection.project_cylinders("XY")
#     flexible_collection.watershed_boundary(flexible_collection.graph)
#     actual_poly = flexible_collection.hull
#     expected_poly = small_tree_wateshed_poly
#     assert actual_poly == expected_poly


# @pytest.mark.parametrize("flexible_collection", ["1_TenCyls.csv"], indirect=True)
# def test_find_flows_ten_cyls(flexible_collection):
#     flexible_collection.project_cylinders("XY")
#     flexible_collection.initialize_graph()
#     flexible_collection.find_flow_components()
#     flexible_collection.calculate_flows()
#     actual = flexible_collection.flows
#     expected = ten_cyls_flows
#     breakpoint()
#     assert expected == str(actual)


# # needs verification
# @pytest.mark.parametrize("flexible_collection", ["5_SmallTree.csv"], indirect=True)
# def test_find_flows_small_tree(flexible_collection):
#     flexible_collection.find_flow_components()
#     flexible_collection.calculate_flows()
#     actual = flexible_collection.flows
#     flexible_collection.draw(
#         highlight_lambda=lambda: is_stem, filter_lambda=lambda: cyl_id > 100
#     )
#     expected = small_tree_flows
#     breakpoint()
#     assert expected == str(actual)
