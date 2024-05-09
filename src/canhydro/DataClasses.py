from __future__ import annotations

from dataclasses import dataclass
from typing import Union, Optional

import numpy as np
from shapely.geometry import Polygon


@dataclass
class Projection:
    plane: str
    polygon: Polygon()
    base_vector: list[int]
    anti_vector: list[int]
    angle: int()


@dataclass
class Flow:
    num_cylinders: int()
    projected_area: np.float64()
    surface_area: np.float64()
    angle_sum: np.float64()
    volume: np.float64()
    sa_to_vol: np.float64()
    drip_node_id: int()
    drip_node_loc: tuple()


coord_list = Union[list[tuple[np.float64]], np.ndarray]

