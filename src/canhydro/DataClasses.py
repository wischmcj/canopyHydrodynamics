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

    def __eq__(self,flow:Flow):
        return (np.round(self.num_cylinders, decimals=1) == np.round(flow.num_cylinders, decimals=1)
                and np.round(self.projected_area, decimals=3)   == np.round(flow.projected_area, decimals=3)
                and np.round(self.surface_area  , decimals=3)       == np.round(flow.surface_area, decimals=3)
               and np.round(self.angle_sum     , decimals=3)       == np.round(flow.angle_sum, decimals=3)
                and np.round(self.volume        , decimals=3)       == np.round(flow.volume, decimals=3)
                and np.round(self.sa_to_vol     , decimals=3)       == np.round(flow.sa_to_vol, decimals=3)
                and np.round(self.drip_node_id   , decimals=0)       == np.round(flow.drip_node_id, decimals=0))

    def compare(self,flow):
        same_drip= self.drip_node_id == flow.drip_node_id
        ret = {  'num_cylinders':np.round(self.num_cylinders, decimals=1) - np.round(flow.num_cylinders, decimals=1)
                ,'projected_area':np.round(self.projected_area, decimals=3)  - np.round(flow.projected_area, decimals=3)
                ,'surface_area' :np.round(self.surface_area  , decimals=3)  - np.round(flow.surface_area, decimals=3)
                ,'angle_sum'    :np.round(self.angle_sum     , decimals=3)  - np.round(flow.angle_sum, decimals=3)
                ,'volume'       :np.round(self.volume        , decimals=3)  - np.round(flow.volume, decimals=3)
                ,'sa_to_vol'    :np.round(self.sa_to_vol     , decimals=3)  - np.round(flow.sa_to_vol, decimals=3)
                ,'drip_node_id' :self.drip_node_id
                ,'drip_node_loc':self.drip_node_loc 
                }
        comp_flow = Flow(**ret)
        return comp_flow


    def add_cyls(self,cylinders: list[Cylinder]):
        for cylinder in cylinders:
            add_cyl(cylinder)

    def add_cyl(self,cylinder):
        self.num_cylinders += 1
        self.projected_area += cylinder.projected_area
        self.surface_area += cylinder.surface_area
        self.angle_sum += cylinder.angle_sum
        self.volume += cylinder.volume
        self.sa_to_vol = self.surface_area/self.volume

coord_list = Union[list[tuple[np.float64]], np.ndarray]

