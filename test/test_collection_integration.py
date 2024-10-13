"""This is a sample python file for testing functions from the source code."""

from __future__ import annotations

import os
import sys
from test.expected_results import (drip_adj_flows, drip_adj_stem_map,
                                   drip_mid_flows, drip_mid_stem_map,
                                   drip_on_trunk_flows, drip_on_trunk_stem_map,
                                   ez_projection_cyls, ez_projection_xy,
                                   ez_projection_xz, ez_projection_yz,
                                   happy_path_cyls, happy_path_dbh,
                                   happy_path_edges, happy_path_flows,
                                   happy_path_is_stem, small_tree_dbh,
                                   small_tree_edges, small_tree_flows,
                                   small_tree_is_stem, ten_cyls_bo_and_len,
                                   ten_cyls_bo_one, ten_cyls_cyls,
                                   ten_cyls_dbh, ten_cyls_edges,
                                   ten_cyls_flows, ten_cyls_id_one,
                                   ten_cyls_is_stem)
from test.expected_results_shapes import (small_tree_overlap,
                                          small_tree_wateshed_poly)
from test.utils import within_range

import pytest

from canopyhydro.configuration import test_input_dir
from canopyhydro.CylinderCollection import (CylinderCollection,
                                            pickle_collection,
                                            unpickle_collection)
from canopyhydro.utils import lam_filter

sys.path.insert(0, os.path.dirname(os.getcwd()))

create_cylinders_cases = [
    # (file, expected_cylinders )
    pytest.param("1_TenCyls.csv", ten_cyls_cyls, id="Ten Cyls"),
    pytest.param("2_EZ_projection.csv", ez_projection_cyls, id="EZ Projection"),
    pytest.param("3_HappyPathWTrunk.csv", happy_path_cyls, id="Happy Path"),
]

ez_projection_cases = [
    # (file, angles, projection )
    pytest.param("2_EZ_projection.csv", ez_projection_xy, "XY", id="XY projection"),
    pytest.param("2_EZ_projection.csv", ez_projection_xz, "XZ", id="XZ projection"),
    pytest.param("2_EZ_projection.csv", ez_projection_yz, "YZ", id="YZ projection"),
]

lam_filter_cases = [
    # (file, angles, projection )
    pytest.param(
        "1_TenCyls.csv",
        ten_cyls_bo_one,
        lambda: branch_order == 1,  # noqa
        id="XY projection",
    ),
    pytest.param(
        "1_TenCyls.csv",
        ten_cyls_bo_and_len,
        lambda: branch_order == 0 or length <= 0.22447,  # noqa
        id="XZ projection",
    ),
    pytest.param(
        "1_TenCyls.csv",
        ten_cyls_id_one,
        lambda: cyl_id == 1,  # noqa
        id="YZ projection",
    ),
]

find_flows_cases = [
    # (file, expected_stem_map, expected_flows )
    pytest.param("1_TenCyls.csv", ten_cyls_is_stem, ten_cyls_flows, id="Ten Cyls"),
    pytest.param(
        "7_DripPathAdjToTrunk.csv",
        drip_adj_stem_map,
        drip_adj_flows,
        id="Drip Adjacent Trunk",
    ),
    pytest.param(
        "8_DripPathMidBranch.csv",
        drip_mid_stem_map,
        drip_mid_flows,
        id="Drip Mid Branch",
    ),
    pytest.param(
        "9_DripOnTrunk.csv",
        drip_on_trunk_stem_map,
        drip_on_trunk_flows,
        id="Drip On Trunk",
    ),
    pytest.param(
        "3_HappyPathWTrunk.csv", happy_path_is_stem, happy_path_flows, id="Happy Path"
    ),
    pytest.param(
        "5_SmallTree.csv", small_tree_is_stem, small_tree_flows, id="Small Tree"
    ),
]


pickle_cases = [
    # (file, expected_stem_map, expected_flows )
    pytest.param("5_SmallTree.csv", id="Small Tree"),
    pytest.param("7_DripPathAdjToTrunk.csv", id="Drip Adjacent Trunk"),
    pytest.param("8_DripPathMidBranch.csv", id="Drip Mid Branch"),
    pytest.param("9_DripOnTrunk.csv", id="Drip On Trunk"),
    pytest.param("3_HappyPathWTrunk.csv", id="Happy Path"),
]

create_graph_cases = [
    # (file, expected_edges )
    pytest.param("1_TenCyls.csv", ten_cyls_edges, id="Ten Cyls"),
    pytest.param("3_HappyPathWTrunk.csv", happy_path_edges, id="Happy Path"),
    pytest.param("5_SmallTree.csv", small_tree_edges, id="Small Tree"),
]

dbh_cases = [
    # (file, expected_edges )
    pytest.param("1_TenCyls.csv", ten_cyls_dbh, id="Ten Cyls (Too Short)"),
    pytest.param("3_HappyPathWTrunk.csv", happy_path_dbh, id="Happy Path"),
    pytest.param("5_SmallTree.csv", small_tree_dbh, id="Small Tree"),
]


@pytest.mark.parametrize(
    "basic_collection, expected_cylinders",
    create_cylinders_cases,
    indirect=["basic_collection"],
)
def test_create_cylinders(basic_collection, expected_cylinders):
    actual = basic_collection.get_collection_data()
    expected = expected_cylinders
    assert expected == actual


@pytest.mark.parametrize("file_name, expected_cylinders", create_cylinders_cases)
def test_create_cylinders_from_csv(file_name, expected_cylinders):
    csv_collection = CylinderCollection()
    file_path = "/".join([str(test_input_dir), file_name])
    file_obj = open(file_path)
    csv_collection.from_csv(file_obj)
    actual = csv_collection.get_collection_data()
    assert actual == expected_cylinders


@pytest.mark.parametrize(
    "flexible_collection, angles, projection_axis",
    ez_projection_cases,
    indirect=["flexible_collection"],
)
def test_project_cylinders(
    flexible_collection, angles, projection_axis, accepted_err=0.03
):
    for idx, cyl in enumerate(flexible_collection.cylinders):
        expected_angle = angles[idx]
        actual_angle = cyl.projected_data[projection_axis]["angle"]
        if projection_axis == "XY":
            actual_setup_angle = cyl.angle
            assert within_range(expected_angle, actual_setup_angle, accepted_err)
        assert within_range(expected_angle, actual_angle, accepted_err)


@pytest.mark.parametrize(
    "basic_collection, expected_result, lam_func",
    lam_filter_cases,
    indirect=["basic_collection"],
)
def test_lam_filter(basic_collection, expected_result, lam_func):
    cyls_returned, _ = lam_filter(basic_collection.cylinders, lam_func)
    actual_result = [cyl.cyl_id for cyl in cyls_returned]
    assert actual_result == expected_result


@pytest.mark.parametrize(
    "basic_collection, expected_stem_map, expected_flows",
    find_flows_cases,
    indirect=["basic_collection"],
)
def test_find_flows(basic_collection, expected_stem_map, expected_flows):
    basic_collection.project_cylinders("XY")
    basic_collection.initialize_digraph_from()
    basic_collection.find_flow_components()
    basic_collection.calculate_flows()
    actual_flows = basic_collection.flows
    _, actual_stem_map = lam_filter(
        basic_collection.cylinders,
        lambda: is_stem,  # noqa
        return_all=True,
    )
    print(actual_flows)
    print(expected_flows)
    try:
        assert actual_flows == expected_flows
        assert actual_stem_map == expected_stem_map
    except AssertionError as e:
        # breakpoint()
        raise e


@pytest.mark.parametrize(
    "basic_collection, expected_stem_map, expected_flows",
    find_flows_cases,
    indirect=["basic_collection"],
)
def test_pickle(basic_collection, expected_stem_map, expected_flows):
    basic_collection.project_cylinders("XY")
    basic_collection.initialize_digraph_from()
    basic_collection.find_flow_components()
    basic_collection.calculate_flows()
    pickle_file = pickle_collection(basic_collection)

    unpickled_collection = unpickle_collection(pickle_file)

    actual_flows = unpickled_collection.flows
    _, actual_stem_map = lam_filter(
        unpickled_collection.cylinders,
        lambda: is_stem,  # noqa
        return_all=True,
    )
    print(actual_flows)
    print(expected_flows)
    assert actual_flows == expected_flows
    assert actual_stem_map == expected_stem_map


@pytest.mark.parametrize(
    "basic_collection, expected_dbh", dbh_cases, indirect=["basic_collection"]
)
def test_dbh(basic_collection, expected_dbh):
    basic_collection.get_dbh()
    actual_dbh = basic_collection.treeQualities["dbh"]
    assert actual_dbh == expected_dbh


@pytest.mark.parametrize("flexible_collection", ["5_SmallTree.csv"], indirect=True)
def test_collection_overlap(flexible_collection):
    flexible_collection.project_cylinders("XY")
    actual = flexible_collection.find_overlap_by_percentile()
    assert actual == small_tree_overlap


@pytest.mark.parametrize("flexible_collection", ["5_SmallTree.csv"], indirect=True)
def test_watershed(flexible_collection):
    flexible_collection.initialize_digraph_from()
    flexible_collection.project_cylinders("XY")
    flexible_collection.watershed_boundary("XY")
    actual_poly = flexible_collection.hulls["XY"]
    expected_poly = small_tree_wateshed_poly
    assert actual_poly == expected_poly
