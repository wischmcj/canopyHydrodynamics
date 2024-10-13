from __future__ import annotations

from dataclasses import dataclass
from typing import Union

import numpy as np
from shapely.geometry import Polygon


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
    projected_area: np.float64
    surface_area: np.float64
    angle_sum: np.float64
    volume: np.float64
    sa_to_vol: np.float64
    drip_node_id: int
    drip_node_loc: tuple

    def __eq__(self, flow: Flow):
        return (
            np.round(self.num_cylinders, decimals=1)
            == np.round(flow.num_cylinders, decimals=1)
            and np.round(self.projected_area, decimals=3)
            == np.round(flow.projected_area, decimals=3)
            and np.round(self.surface_area, decimals=3)
            == np.round(flow.surface_area, decimals=3)
            and np.round(self.angle_sum, decimals=3)
            == np.round(flow.angle_sum, decimals=3)
            and np.round(self.volume, decimals=3) == np.round(flow.volume, decimals=3)
            and np.round(self.sa_to_vol, decimals=3)
            == np.round(flow.sa_to_vol, decimals=3)
        )


coord_list = Union[list[tuple[np.float64]], np.ndarray]
