import pytest
from numpy import float16  # assuming you are using numpy.float16
from canopyhydro.DataClasses import Flow  # Adjust import based on your module


@pytest.fixture
def flow1():
    return Flow(
        num_cylinders=4,
        projected_area=float16(20.0),
        surface_area=float16(30.0),
        angle_sum=float16(15.0),
        volume=float16(100.0),
        sa_to_vol=float16(0.3),
        drip_node_id=101,
        drip_node_loc=(1.0, 1.0),
        cyl_list={1, 2, 3, 4},
    )


@pytest.fixture
def flow2():
    return Flow(
        num_cylinders=2,
        projected_area=float16(10.0),
        surface_area=float16(20.0),
        angle_sum=float16(10.0),
        volume=float16(50.0),
        sa_to_vol=float16(0.4),
        drip_node_id=102,
        drip_node_loc=(2.0, 2.0),
        cyl_list={5, 6},
    )


@pytest.fixture
def int_flow():
    return Flow(
        num_cylinders=3,
        projected_area=30,
        surface_area=40,
        angle_sum=20,
        volume=100,
        sa_to_vol=0.3,
        drip_node_id=101,
        drip_node_loc=(2.0, 2.0),
        cyl_list={5, 6, 7},
    )


def test_post_init(int_flow):
    assert isinstance(int_flow.projected_area, float16)
    assert isinstance(int_flow.surface_area, float16)
    assert isinstance(int_flow.angle_sum, float16)
    assert isinstance(int_flow.volume, float16)
    assert isinstance(int_flow.sa_to_vol, float16)


def test_getitem(flow1):
    assert flow1["num_cylinders"] == 4
    assert flow1["projected_area"] == float16(20.0)


def test_eq(flow1):
    assert flow1 == flow1  # Expecting flow1 and flow2 to be different


def test_add(flow1, flow2):
    flow_sum = flow1 + flow2
    assert flow_sum.num_cylinders == 6
    assert flow_sum.projected_area == float16(30.0)
    assert flow_sum.surface_area == float16(50.0)
    assert flow_sum.angle_sum == float16(25.0)
    assert flow_sum.volume == float16(150.0)
    assert flow_sum.cyl_list == {1, 2, 3, 4, 5, 6}


def test_sub(int_flow, flow2):
    flow_diff = int_flow - flow2
    assert flow_diff.num_cylinders == 1
    assert flow_diff.projected_area == float16(20.0)
    assert flow_diff.surface_area == float16(20.0)
    assert flow_diff.angle_sum == float16(10.0)
    assert flow_diff.volume == float16(50.0)
    assert flow_diff.cyl_list == {7}


def test_pct_diff(flow1, flow2):
    pct_diff = flow1.pct_diff(flow2)
    assert abs(pct_diff["projected_area"] - (-0.5)) < 1e-2  # 50% difference
    assert abs(pct_diff["surface_area"] - (-0.333)) < 1e-2  # 33.3% difference
    assert abs(pct_diff["volume"] - (-0.5)) < 1e-2  # 50% difference


def test_within_range(flow1, flow2):
    assert not flow1.within_range(
        flow2, 0.02
    )  # Should return False since they're > 2% different
    assert flow1.within_range(flow2, 0.5)  # Should return True if within 50%


def test_add_cyl(flow1):
    class MockCylinder:
        def __init__(self, cyl_id, projected_area, surface_area, angle_sum, volume):
            self.cyl_id = cyl_id
            self.projected_area = projected_area
            self.surface_area = surface_area
            self.angle_sum = angle_sum
            self.volume = volume

    cyl = MockCylinder(6, float16(5.0), float16(7.0), float16(3.0), float16(10.0))
    flow1.add_cyl(cyl)
    assert flow1.num_cylinders == 5
    assert flow1.projected_area == float16(25.0)
    assert flow1.surface_area == float16(37.0)
    assert flow1.angle_sum == float16(18.0)
    assert flow1.volume == float16(110.0)
    assert flow1.cyl_list == {1, 2, 3, 4, 6}


def test_add_cyls(flow1):
    class MockCylinder:
        def __init__(self, cyl_id, projected_area, surface_area, angle_sum, volume):
            self.cyl_id = cyl_id
            self.projected_area = projected_area
            self.surface_area = surface_area
            self.angle_sum = angle_sum
            self.volume = volume

    cylinders = [
        MockCylinder(7, float16(2.0), float16(3.0), float16(1.0), float16(5.0)),
        MockCylinder(8, float16(3.0), float16(4.0), float16(2.0), float16(6.0)),
    ]

    flow1.add_cyls(cylinders)
    assert flow1.num_cylinders == 4  # Cylinders count does not change
    assert flow1.projected_area == float16(25.0)  # 20.0 + 2.0 + 3.0
    assert flow1.surface_area == float16(37.0)  # 30.0 + 3.0 + 4.0
    assert flow1.angle_sum == float16(18.0)  # 15.0 + 1.0 + 2.0
    assert flow1.volume == float16(111.0)  # 100.0 + 5.0 + 6.0
    assert flow1.cyl_list == {1, 2, 3, 4, 7, 8}
