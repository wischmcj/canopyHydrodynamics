from __future__ import annotations

from dataclasses import dataclass
from typing import Union
import math

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
    projected_area: np.float16
    surface_area: np.float16
    angle_sum: np.float16
    volume: np.float16
    sa_to_vol: np.float16
    drip_node_id: int
    drip_node_loc: tuple
    
    def __post_init__(self):
        self.projected_area= np.float16(self.projected_area)
        self.surface_area  = np.float16(self.surface_area)  
        self.angle_sum     = np.float16(self.angle_sum)     
        self.volume        = np.float16(self.volume)       
        self.sa_to_vol     = np.float16(self.sa_to_vol)     

    def __eq__(self,flow:Flow):
        diff = self.pct_diff(flow)
        for value in diff.values():
            if value > 0.01:
                return False
        return True
    
    def __add__(self,flow:Flow):
        return Flow( self.num_cylinders  + flow.num_cylinders     
                    ,self.projected_area + flow.projected_area    
                    ,self.surface_area   + flow.surface_area      
                    ,self.angle_sum      + flow.angle_sum         
                    ,self.volume         + flow.volume            
                    ,self.sa_to_vol      + flow.sa_to_vol         
                    ,self.drip_node_id
                    ,self.drip_node_loc)
    
    def __sub__(self,flow:Flow):
        return self.compare(flow)
    
    def __repr__(self) -> str:
        return f'''Flow(num_cylinders={self.num_cylinders}, 
                    projected_area={self.projected_area:2f}, 
                    surface_area={self.surface_area:2f}, 
                    volume={self.volume:2f}, 
                    drip_node_id={self.drip_node_id},
                    drip_node_loc={self.drip_node_loc}'''
    
    def compare(self,flow):
        return Flow( self.num_cylinders  - flow.num_cylinders   
             ,self.projected_area - flow.projected_area  
             ,self.surface_area   - flow.surface_area    
             ,self.angle_sum      - flow.angle_sum       
             ,self.volume         - flow.volume          
             ,self.sa_to_vol      - flow.sa_to_vol       
             ,self.drip_node_id
             ,math.dist(self.drip_node_loc, flow.drip_node_loc))
    
    def pct_diff(self,flow):
        diff_attrs = ['num_cylinders','projected_area','surface_area','volume']
        diff ={}
        for attr in diff_attrs:
            if getattr(self,attr)==0:
                diff[attr] = 0 if getattr(flow,attr)==0 else (getattr(flow,attr)-getattr(self,attr))/np.float16(0.00001)
            else:
                diff[attr] = (getattr(flow,attr)-getattr(self,attr))/getattr(self,attr)
        return diff

    def within_range(self,flow,rng):
        diff = self.pct_diff(flow)
        rng = abs(rng)
        for key,value in diff.items():
            if key != 'num_cylinders':
                if abs(value)> rng:
                    return False
        return True

    def add_cyls(self,cylinders: list[Cylinder]):
        for cylinder in cylinders:
            self.add_cyl(cylinder)

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

