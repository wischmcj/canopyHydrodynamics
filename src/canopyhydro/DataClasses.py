from __future__ import annotations

from dataclasses import dataclass, field
from math import dist
from typing import Union

from numpy import float16, ndarray, sum
from shapely.geometry import Polygon

coord_list = Union[list[tuple[float16]], ndarray]


@dataclass
class Projection:
    plane: str
    polygon: Polygon
    base_vector: list[int]
    anti_vector: list[int]
    angle: int


@dataclass
class Flow:
    num_cylinders: int
    projected_area: float16
    surface_area: float16
    angle_sum: float16
    volume: float16
    sa_to_vol: float16 | None = None
    drip_node_id: int | None = None
    drip_node_loc: tuple | None = ()
    cyl_list: set | None = field(default_factory=set)

    def __post_init__(self):
        self.projected_area = float16(self.projected_area)
        self.surface_area = float16(self.surface_area)
        self.angle_sum = float16(self.angle_sum)
        self.volume = float16(self.volume)
        if self.sa_to_vol:
            self.sa_to_vol = float16(self.sa_to_vol)
        elif self.volume != 0:
            self.sa_to_vol = self.surface_area / self.volume or 0

    def __getitem__(self, attr):
        return getattr(self, attr)

    def __eq__(self, flow: Flow):
        """Two flows are considered equal if all of their properties are <2% different"""
        diff = self.pct_diff(flow)
        for value in diff.values():
            if value > 0.02:
                return False
        return True

    def __add__(self, flow: Flow):
        return Flow(
            self.num_cylinders + flow.num_cylinders,
            self.projected_area + flow.projected_area,
            self.surface_area + flow.surface_area,
            self.angle_sum + flow.angle_sum,
            self.volume + flow.volume,
            self.sa_to_vol + flow.sa_to_vol,
            self.drip_node_id,
            self.drip_node_loc,
            self.cyl_list.union(flow.cyl_list),
        )

    def __sub__(self, flow: Flow):
        cyls = self.cyl_list.difference(flow.cyl_list)
        diff = Flow(
            self.num_cylinders - flow.num_cylinders,
            self.projected_area - flow.projected_area,
            self.surface_area - flow.surface_area,
            self.angle_sum - flow.angle_sum,
            self.volume - flow.volume,
            self.sa_to_vol - flow.sa_to_vol,
            self.drip_node_id,
            dist(self.drip_node_loc, flow.drip_node_loc),  # Distance between points
            cyls,
        )
        return diff

    def __repr__(self) -> str:
        return f"""Flow(num_cylinders={self.num_cylinders},
                    projected_area={self.projected_area:2f},
                    surface_area={self.surface_area:2f},
                    volume={self.volume:2f},
                    drip_node_id={self.drip_node_id},
                    drip_node_loc={self.drip_node_loc}"""

    def pct_diff(self, flow):
        """
        Calculates the % difference between each attribute
        of the flows
        """
        diff_attrs = ["num_cylinders", "projected_area", "surface_area", "volume"]
        diff = {}
        for attr in diff_attrs:
            if self[attr] == 0:
                # in this case we consider self.attr to be 0.00001 instead
                diff[attr] = (
                    0
                    if flow[attr] == 0
                    else (flow[attr] - self[attr]) / float16(0.00001)
                )
            else:
                diff[attr] = (flow[attr] - self[attr]) / self[attr]
        return diff

    def within_range(self, flow, rng):
        diff = self.pct_diff(flow)
        rng = abs(rng)
        for key, value in diff.items():
            if key != "num_cylinders":
                if abs(value) > rng:
                    return False
        return True

    def add_cyls(self, cylinders: list):
        cyls = {cyl.cyl_id for cyl in cylinders}
        cyls_sum = sum(
            (
                (
                    cylinder.projected_area,
                    cylinder.surface_area,
                    cylinder.angle_sum,
                    cylinder.volume,
                )
                for cylinder in cylinders
            ),
            axis=0,
        )
        self.projected_area += cyls_sum[0]
        self.surface_area += cyls_sum[1]
        self.angle_sum += cyls_sum[2]
        self.volume += cyls_sum[3]
        self.sa_to_vol = self.surface_area / self.volume
        self.cyl_list = self.cyl_list.union(cyls)

    def add_cyl(self, cylinder):
        self.num_cylinders += 1
        self.projected_area += cylinder.projected_area
        self.surface_area += cylinder.surface_area
        self.angle_sum += cylinder.angle_sum
        self.volume += cylinder.volume
        self.sa_to_vol = self.surface_area / self.volume
        self.cyl_list.add(cylinder.cyl_id)
