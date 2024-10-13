from __future__ import annotations

# from test.expected_results import ez_projection_vectors, ez_projection_xy_angle
from test.utils import within_range

import numpy as np
import pytest

from canopyhydro.Cylinder import Cylinder

# Test Config
accepted_err = 0.02

# general test_cyl
cyl_list = np.array(
    [
        None,
        0.000,
        -1.000,
        -0.299115,
        2.537844,
        -0.598273,
        -0.299115,
        2.537844,
        -0.552,
        0.557981,
        0.045261,
        7.269781,
        0.046273,
        7955.52145,
        None,
        0.000,
        -1.000,
        0.0529898,
        5.983393,
        7954.0735,
        0.000,
        67.0,
        1.35908057e02,
        18471,
        0.000,
        0.000,
        None,
        0.000,
        0.000,
        0.000,
        0.01,
        0.0401606,
        0.000,
        0.01,
        0.0401606,
    ]
)

cyl_obj_param = Cylinder(
    cyl_id=0.0,
    x=np.array([-0.299115, -0.299115], dtype=object),
    y=np.array([2.537844, 2.537844], dtype=object),
    z=np.array([-0.598273, -0.552], dtype=object),
    radius=0.557981,
    length=0.046273,
    branch_order=0.0,
    branch_id=0.0,
    volume=0.045261,
    parent_id=-1.0,
    reverse_branch_order=67.0,
    segment_id=0.0,
    projected_data={},
    flow_id=None,
    flow_type=None,
    begins_at_drip_point=None,
    begins_at_divide_point=None,
    dx=0.0,
    dy=0.0,
    dz=0.04627300000000001,
    surface_area=0.16222841912042885,
    sa_to_vol=3.584287115185896,
    slope=0.0,
    is_stem=False,
)

cyl_sa = 0.162

# projection_test_cyl
proj_cyl_list = np.array(
    [
        None,
        0.00000000e00,
        -1.00000000e00,
        1.00000000e00,
        1.00000000e00,
        1.00000000e00,
        4.00000000e00,
        6.00000000e00,
        7.00000000e00,
        1.00000000e00,
        1.00210000e-02,
        3.09237800e00,
        6.44330000e-02,
        2.13541289e03,
        None,
        0.00000000e00,
        -1.00000000e00,
        1.97783000e-01,
        2.94657100e00,
        2.13431216e03,
        0.00000000e00,
        3.20000000e01,
        5.67450440e01,
        3.22000000e03,
        0.00000000e00,
        0.00000000e00,
        None,
        0.00000000e00,
        0.00000000e00,
        0.00000000e00,
        1.00000000e-02,
        4.01606000e-01,
        0.00000000e00,
        1.00000000e-02,
        4.01606000e-01,
    ]
)

proj_cyl = Cylinder(
    cyl_id=0.0,
    x=np.array([1.0, 4.0]),
    y=np.array([1.0, 6.0]),
    z=np.array([1.0, 7.0]),
    radius=1.0,
    length=0.064433,
    branch_order=0.0,
    branch_id=0.0,
    volume=0.010021,
    parent_id=-1.0,
    reverse_branch_order=32.0,
    segment_id=0.0,
    projected_data={},
    flow_id=None,
    flow_type=None,
    begins_at_drip_point=None,
    begins_at_divide_point=None,
    dx=3.0,
    dy=5.0,
    dz=6.0,
    surface_area=0.40484447889750186,
    sa_to_vol=40.399608711456125,
    slope=0.0,
    is_stem=False,
)

ez_projection_vectors = {
    "XY": [np.array([1.0, 1.0, 1.0]), np.array([4.0, 6.0, 7.0])],
    "XZ": [np.array([1.0, 1.0, 1.0]), np.array([4.0, 7.0, 6.0])],
    "YZ": [np.array([1.0, 1.0, 1.0]), np.array([6.0, 7.0, 4.0])],
}
ez_projection_xy_angle = 0.7996
ez_projection_xz_angle = 0.6405
ez_projection_yz_angle = 0.3667


@pytest.mark.parametrize(
    "test_cyl, cyl_obj", [(cyl_list, cyl_obj_param)], indirect=["test_cyl"]
)
def test_create_test_cyl(test_cyl, cyl_obj):
    """
    Tests the creation of the test cylinder fixture
    """
    assert test_cyl == cyl_obj


@pytest.mark.parametrize(
    "test_cyl, ex_surface_area", [(cyl_list, cyl_sa)], indirect=["test_cyl"]
)
def test_calc_surface_area(test_cyl, ex_surface_area):
    surface_area = test_cyl.surface_area
    assert within_range(ex_surface_area, surface_area, accepted_err)


@pytest.mark.parametrize("proj_cyl", [proj_cyl])
def test_project_cylinder(proj_cyl):
    proj_cyl.get_projection("XY")
    actual = proj_cyl.projected_data["XY"]["angle"]
    expected = ez_projection_xy_angle
    assert within_range(expected, actual, accepted_err)


# @pytest.mark.parametrize("proj_cyl", [proj_cyl])
# def test_project_cylinder_numba(proj_cyl):
#     proj_cyl.numba_get_projection("XY")
#     actual = proj_cyl.projected_data["XY"]["angle"]
#     expected = ez_projection_xy_angle
#     assert within_range(expected, actual, accepted_err)


@pytest.mark.parametrize("test_cyl", [(proj_cyl_list)], indirect=["test_cyl"])
def test_create_cyl_angle(test_cyl):
    """
    Tests the creation of the test cylinder fixture
    """
    actual = test_cyl.angle
    expected = ez_projection_xy_angle
    assert within_range(expected, actual, accepted_err)


@pytest.mark.parametrize("test_cyl", [(proj_cyl_list)], indirect=["test_cyl"])
def test_create_cyl_vectors(test_cyl):
    """
    Tests the creation of the test cylinder fixture
    """
    actual = {
        "XY": [
            np.array([test_cyl.x[0], test_cyl.y[0], test_cyl.z[0]]),
            np.array([test_cyl.x[1], test_cyl.y[1], test_cyl.z[1]]),
        ],
        "XZ": [
            np.array([test_cyl.x[0], test_cyl.z[0], test_cyl.y[0]]),
            np.array([test_cyl.x[1], test_cyl.z[1], test_cyl.y[1]]),
        ],
        "YZ": [
            np.array([test_cyl.y[0], test_cyl.z[0], test_cyl.x[0]]),
            np.array([test_cyl.y[1], test_cyl.z[1], test_cyl.x[1]]),
        ],
    }
    expected = ez_projection_vectors
    # Probably should mae a  data class for these cylinder vectors... oh well
    for k, v in actual.items():
        for idx, vector in enumerate(v):
            assert np.all(expected[k][idx] == vector)
