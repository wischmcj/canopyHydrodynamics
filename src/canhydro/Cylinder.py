"""Defines the component parts of the ingested QSM"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from src.canhydro.DataClasses import Projection
from src.canhydro.geometry import (draw_cyls, get_projection)
                                #    numba_get_projection)
from src.canhydro.global_vars import qsm_cols

# from descartes import PolygonPatch
# from mpl_toolkits import mplot3d


# import time
# import copy
# import math
# import openpyxl
# import geopandas as geo
NAME = "Cylinder"


def create_cyl(arr: np.array):
    cols = qsm_cols
    attrs = {k: arr[v] for (k, v) in cols.items()}
    cyl = Cylinder(**attrs)
    cyl.create_from_list(arr, cols)
    return cyl


@dataclass
class Cylinder:  # (defaultdict):
    cyl_id: int
    x: np.ndarray[np.float32]  # len 2 array
    y: np.ndarray[np.float32]  # len 2 array
    z: np.ndarray[np.float32]  # len 2 array
    radius: np.float32
    length: np.float32
    branch_order: int
    branch_id: int
    volume: np.float32
    parent_id: int
    reverse_branch_order: int
    segment_id: int

    projected_data: dict(Projection) = field(default_factory=dict)
    flow_id: int() = None
    flow_type: str = None
    drip_node: int() = None
    begins_at_drip_point: bool = None
    begins_at_divide_point: bool = None

    stem_path_id = int

    dx: np.float32 = 0
    dy: np.float32 = 0
    dz: np.float32 = 0

    surface_area: np.float32 = 0.0
    sa_to_vol: np.float32 = 0.0
    slope: np.float32 = 0.0

    is_stem: bool = False
    is_divide: bool = False

    # #blood for the blood god, software eng for the filter func
    # class_attrs = self.__get_class_attributes(type(self))
    # self.__init_instance(class_attrs, kwargs)
    def __repr__(self):
        return f"Cylinder( cyl_id={self.cyl_id}, x={self.x}, y={self.y}, z={self.z}, radius={self.radius}, length={self.length}, branch_order={self.branch_order}, branch_id={self.branch_id}, volume={self.volume}, parent_id={self.parent_id}, reverse_branch_order={self.reverse_branch_order}, segment_id={self.segment_id}"

    def __eq__(self, other):
        return type(self) == type(other) and self.__repr__() == other.__repr__()

    def calc_surface_area(self):
        radius = self.radius
        length = self.length
        sa = 2 * np.pi * radius * (radius + length) - 2 * np.pi * radius * radius
        return sa

    def create_from_list(self, attrs: list, columns=qsm_cols):
        """creates a cylinder corrosponding to that defined by a given row of the qsm (attrs)"""

        extract = (
            lambda attr: attrs[columns[attr]]
        )  # pulls a column from the qsm row (attrs) corrosponding to the input attribute

        self.dx = self.x[1] - self.x[0]
        self.dy = self.y[1] - self.y[0]
        self.dz = self.z[1] - self.z[0]
        self.vectors = {
            "XY": [
                np.array([self.x[0], self.y[0], self.z[0]]),
                np.array([self.x[1], self.y[1], self.z[1]]),
            ],
            "XZ": [
                np.array([self.x[0], self.z[0], self.y[0]]),
                np.array([self.x[1], self.z[1], self.y[1]]),
            ],
            "YZ": [
                np.array([self.y[0], self.z[0], self.x[0]]),
                np.array([self.y[1], self.z[1], self.x[1]]),
            ],
        }
        self.surface_area = self.calc_surface_area()
        self.sa_to_vol = 0 if self.volume == 0 else self.surface_area / self.volume
        run = np.sqrt(self.dx**2 + self.dy**2)
        self.angle = (
            np.arctan(0)
            if run == 0
            else np.arctan(self.dz / np.sqrt(self.dx**2 + self.dy**2))
        )
        self.xy_area = 0
        # log.info(str(self.__repr__()))

    def get_projection(self, plane="XY"):
        if plane == "XY":
            magnitude = [self.dx, self.dy, self.dz]
            vector = [np.transpose(self.x), np.transpose(self.y), np.transpose(self.z)]
        elif plane == "XZ":
            magnitude = [self.dx, self.dz, self.dy]
            vector = [np.transpose(self.x), np.transpose(self.z), np.transpose(self.y)]
        else:
            magnitude = [self.dy, self.dz, self.dx]
            vector = [np.transpose(self.y), np.transpose(self.z), np.transpose(self.x)]

        projection = get_projection(vector, magnitude, self.radius)
        self.projected_data[plane] = projection
        if plane == "XY":
            self.xy_area = self.projected_data["XY"]["area"]
        return projection["polygon"]

    # def numba_get_projection(self, plane="XY"):
    #     if plane == "XY":
    #         magnitude = [self.dx, self.dy, self.dz]
    #         vector = [np.transpose(self.x), np.transpose(self.y), np.transpose(self.z)]
    #     elif plane == "XZ":
    #         magnitude = [self.dx, self.dz, self.dy]
    #         vector = [np.transpose(self.x), np.transpose(self.z), np.transpose(self.y)]
    #     else:
    #         magnitude = [self.dy, self.dz, self.dx]
    #         vector = [np.transpose(self.y), np.transpose(self.z), np.transpose(self.x)]

    #     projection = numba_get_projection(vector, magnitude, self.radius)
    #     self.projected_data[plane] = projection
    #     if plane == "XY":
    #         self.xy_area = self.projected_data["XY"]["area"]
    #     return projection["polygon"]
    

    def draw(self, plane: str = "XY"):
        poly = self.projected_data[plane]["polygon"]
        draw_cyls([poly])

    def get_flow_data():
        """Returns the flow ID and flow characteristics of the flow the cyl is contained in"""
        print("Get flow data not written")
