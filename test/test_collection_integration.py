"""This is a sample python file for testing functions from the source code."""

from __future__ import annotations
import logging
import os
import sys
import toml
import pytest
import numpy as np
sys.path.insert(0, os.path.dirname(os.getcwd()))

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
                                   ten_cyls_is_stem,
                                   drip_adj_flows_rust)

from test.expected_results_shapes import (small_tree_overlap,
                                          small_tree_wateshed_poly)
from test.utils import within_range

from src.canhydro.utils import lam_filter
from src.canhydro.DataClasses import Flow
from src.canhydro.CylinderCollection import CylinderCollection, pickle_collection, unpickle_collection
from scripts.basic_recipies import initialize_collection
from scripts.compare_utils import compare_flows

log = logging.getLogger("script")

with open("src/canhydro/user_def_config.toml") as f:
    config = toml.load(f)
    test_input_dir = config["directories"]['test_input_dir']
    DIR = config["directories"]['root_dir']


with open("src/canhydro/user_def_config.toml") as f:
    config = toml.load(f)
    test_input_dir = config["directories"]['test_input_dir']
    DIR = config["directories"]['root_dir']


with open("src/canhydro/user_def_config.toml") as f:
    config = toml.load(f)
    test_input_dir = config["directories"]['test_input_dir']
    DIR = config["directories"]['root_dir']


with open("src/canhydro/user_def_config.toml") as f:
    config = toml.load(f)
    test_input_dir = config["directories"]['test_input_dir']
    DIR = config["directories"]['root_dir']

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
        "1_TenCyls.csv", ten_cyls_bo_one, lambda: branch_order == 1, id="XY projection"
    ),
    pytest.param(
        "1_TenCyls.csv",
        ten_cyls_bo_and_len,
        lambda: branch_order == 0 or length <= 0.22447,
        id="XZ projection",
    ),
    pytest.param(
        "1_TenCyls.csv", ten_cyls_id_one, lambda: cyl_id == 1, id="YZ projection"
    ),
]
find_flows_cases = [
    # (file, expected_stem_map, expected_flows )
    # pytest.param("1_TenCyls.csv", ten_cyls_is_stem, ten_cyls_flows, id="Ten Cyls"),
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
        "5_SmallTree.csv", small_tree_is_stem, small_tree_flows, id="Small Tree"
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
]

find_flows_cases_rust = [
    # (file, expected_stem_map, expected_flows )
    # pytest.param("1_TenCyls.csv", ten_cyls_is_stem, ten_cyls_flows, id="Ten Cyls"),
    pytest.param(
        "7_DripPathAdjToTrunk.csv",
        drip_adj_stem_map,
        drip_adj_flows_rust,
        id="Drip Adjacent Trunk",
    ),
    pytest.param(
        "8_DripPathMidBranch.csv",
        drip_mid_stem_map,
        drip_mid_flows,
        id="Drip Mid Branch",
    ),
    pytest.param(
        "5_SmallTree.csv", small_tree_is_stem, small_tree_flows, id="Small Tree"
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
]

pickle_cases = [
    # (file, expected_stem_map, expected_flows )
    # pytest.param("1_TenCyls.csv", ten_cyls_is_stem, ten_cyls_flows, id="Ten Cyls"),
    pytest.param(
        "5_SmallTree.csv", 
        id="Small Tree"
    ),
    pytest.param(
        "7_DripPathAdjToTrunk.csv",
        id="Drip Adjacent Trunk"
    ),
    pytest.param(
        "8_DripPathMidBranch.csv",
        id="Drip Mid Branch"
    ),
    pytest.param(
        "9_DripOnTrunk.csv",
        id="Drip On Trunk"
    ),
    pytest.param(
        "3_HappyPathWTrunk.csv",
        id="Happy Path"
    ),
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


# @pytest.mark.parametrize("basic_collection", ["7_DripPathAdjTrunk.csv"], indirect=True)
# @pytest.mark.parametrize("basic_collection", ["5_SmallTree.csv"], indirect=True)
# def test_sandbox(basic_collection):

#     basic_collection.initialize_graph_from()
#     basic_collection.project_cylinders("XY")
#     basic_collection.draw()
#     breakpoint()
#     basic_collection.draw(plane = 'XZ', filter_lambda = lambda: cyl_id>100, highlight_lambda = lambda:is_stem)


@pytest.mark.parametrize(
    "basic_collection, expected_cylinders",
    create_cylinders_cases,
    indirect=["basic_collection"],
)
def test_create_cylinders(basic_collection, expected_cylinders):
    actual = basic_collection.get_collection_data()
    expected = expected_cylinders
    assert expected == actual
    
@pytest.mark.parametrize(
    "file_name, expected_cylinders",
    create_cylinders_cases
)
def test_create_cylinders_from_csv(file_name, expected_cylinders):
    csv_collection = CylinderCollection()
    file_path = "/".join([str(test_input_dir), file_name])
    file_obj = open(file_path,'r')
    csv_collection.from_csv(file_obj, DIR)
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

def compare_flows(c1,c2):
    f1 = c1.flows
    f2 = c2.flows
    # _, map1 = lam_filter(
    #     c1.cylinders, lambda: is_stem, return_all=True
    # )
    # _, map2 = lam_filter(
    #     c2.cylinders, lambda: is_stem, return_all=True
    # )
    diffs= []
    unique_drip_nodes= [[],[]]
    
    for idf, flow in f1:
        drip_node = flow.drip_node_id
        compare_to = [flow2 for flow2 in f2 if flow2.drip_node_id ==drip_node ]
        if len(compare_to) == 0:
            diffs[idf] = flow - compare_to
        else:
            unique_drip_nodes[0].append(drip_node)
    f1_drip_nodes = [flow.drip_node_id for flow in f1]
    unique_drip_nodes[1].extend([flow.drip_node_id for flow in f2 if flow.drip_node_id not in f1_drip_nodes])
    return diffs, unique_drip_nodes

@pytest.mark.parametrize(
    "basic_collection, expected_stem_map, expected_flows",
    find_flows_cases,
    indirect=["basic_collection"],
)
def test_find_flows(basic_collection, expected_stem_map, expected_flows):
    basic_collection.project_cylinders("XY")
    basic_collection.initialize_digraph_from()
    basic_collection.find_flow_components_new()
    basic_collection.calculate_flows()
    actual_flows = basic_collection.flows
    _, actual_stem_map = lam_filter(
        basic_collection.cylinders, lambda: is_stem, return_all=True
    )
    print(actual_flows)
    print(expected_flows)
    try:
        assert actual_flows == expected_flows
        assert actual_stem_map == expected_stem_map
    except AssertionError as e:
        breakpoint()
        raise e

@pytest.mark.parametrize(
    "basic_collection, expected_stem_map, expected_flows",
    find_flows_cases_rust,
    indirect=["basic_collection"],
)
def test_find_flows_rust(basic_collection, expected_stem_map, expected_flows):
    basic_collection.project_cylinders("XY")
    basic_collection.initialize_digraph_from_rust()
    basic_collection.find_flow_components_rust()
    basic_collection.calculate_flows_rust()
    actual_flows = basic_collection.flows
    _, actual_stem_map = lam_filter(
        basic_collection.cylinders, lambda: is_stem, return_all=True
    )
    print(actual_flows)
    print(expected_flows)
    try:
        assert actual_flows == expected_flows
        assert actual_stem_map == expected_stem_map
    except AssertionError as e:
        collection = initialize_collection(basic_collection.file_name,from_pickle = False)
        collection.initialize_digraph_from()
        collection.find_flow_components_new()
        collection.calculate_flows()
        # diffs = compare_flows(collection,basic_collection)
        f1 = collection.flows
        f2 = basic_collection.flows
        # _, map1 = lam_filter(
        #     c1.cylinders, lambda: is_stem, return_all=True
        # )
        # _, map2 = lam_filter(
        #     c2.cylinders, lambda: is_stem, return_all=True
        # )
        diffs= []
        unique_drip_nodes= [[],[]]
        for idf, flow in enumerate(f1):
            drip_node = flow.drip_node_id
            compare_to = [flow2 for flow2 in f2 if flow2.drip_node_id ==drip_node ]
            if len(compare_to) == 0:
                unique_drip_nodes[0].append(drip_node)
                compare_to = f2[idf]
            diffs.append(flow.compare(compare_to[0]))
            
        f1_drip_nodes = [flow.drip_node_id for flow in f1]
        unique_drip_nodes[1].extend([flow.drip_node_id for flow in f2 if flow.drip_node_id not in f1_drip_nodes])
        
        missed_cyls = [[cyl for cyl in collection.cylinders if cyl.cyl_id == flow.drip_node_id][0] for flow in f1]
        correction = [ Flow(    1,    np.float64(cyl.projected_data['XY']["area"]),    cyl.surface_area,cyl.angle,cyl.volume,cyl.sa_to_vol,cyl.cyl_id,[0,0,0])   for cyl in missed_cyls ]
        def get_cyl_as_flow(node):
            cyl = [cyl for cyl in collection.cylinders if cyl.cyl_id ==node][0]
            correction =  Flow(    1,    np.float64(cyl.projected_data['XY']["area"]),    cyl.surface_area,cyl.angle,cyl.volume,cyl.sa_to_vol,cyl.cyl_id,[0,0,0])
        log.info(f'{missed_cyls=}')
        log.info(f'{diffs=}')

        breakpoint()    
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
    try:
        pickle_file = pickle_collection(basic_collection)
        unpickled_collection = unpickle_collection(pickle_file)

        actual_flows = unpickled_collection.flows
        _, actual_stem_map = lam_filter(
            unpickled_collection.cylinders, lambda: is_stem, return_all=True
        )
        print(actual_flows)
        print(expected_flows)
        assert actual_flows == expected_flows
        assert actual_stem_map == expected_stem_map
    except Exception as e:
        breakpoint()
        raise e


@pytest.mark.parametrize("flexible_collection", ["1_TenCyls.csv"], indirect=True)
def test_highlight_filt_draw(flexible_collection, accepted_err=0.03):
    flexible_collection.project_cylinders(plane="XZ")
    flexible_collection.draw("XZ")
    # flexible_collection.draw("XZ", highlight_lambda=lambda: branch_order == 1)
    # flexible_collection.draw(
    #     "XZ",
    #     highlight_lambda=lambda: branch_order == 1,
    #     filter_lambda=lambda: cyl_id > 3,
    # )
    assert 1 == 1


@pytest.mark.parametrize(
    "basic_collection, expected_dbh", dbh_cases, indirect=["basic_collection"]
)
def test_dbh(basic_collection, expected_dbh):
    basic_collection.get_dbh()
    actual_dbh = basic_collection.dbh
    assert actual_dbh == expected_dbh


# @pytest.mark.parametrize("flexible_collection", ["5_SmallTree.csv"], indirect=True)
# # @pytest.mark.parametrize("flexible_collection", ["4_LargeCollection.csv"], indirect=True)
# def test_drip_map(flexible_collection):
#     flexible_collection.find_flow_components()
#     flexible_collection.calculate_flows()
#     flexible_collection.drip_map(lambda: cyl_id > 400)
#     breakpoint()
#     assert expected == str(actual)


@pytest.mark.parametrize("flexible_collection", ["5_SmallTree.csv"], indirect=True)
def test_collection_overlap(flexible_collection):
    flexible_collection.project_cylinders("XY")
    # flexible_collection.draw(plane="XY", filter_lambda=lambda: branch_order > 0)
    actual = flexible_collection.find_overlap_by_percentile()
    assert actual == small_tree_overlap


@pytest.mark.parametrize("flexible_collection", ["5_SmallTree.csv"], indirect=True)
def test_watershed(flexible_collection):
    flexible_collection.initialize_digraph_from()
    flexible_collection.project_cylinders("XY")
    flexible_collection.watershed_boundary(flexible_collection.graph)
    actual_poly = flexible_collection.hull
    expected_poly = small_tree_wateshed_poly
    assert actual_poly == expected_poly
