"""Defines the component parts of the ingested QSM"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from canhydro.DataClasses import Projection
from canhydro.global_vars import qsm_cols
from canhydro.Plotter import draw_cyls
from canhydro.utils import get_projection

# from descartes import PolygonPatch
# from mpl_toolkits import mplot3d


# import time
# import copy
# import math
# import openpyxl
# import geopandas as geo
NAME = "Cylinder"


@dataclass
class Cylinder:  # (defaultdict):
    cyl_id: int
    x: np.ndarray  # len 2 array
    y: np.ndarray  # len 2 array
    z: np.ndarray  # len 2 array
    radius: float
    length: float
    branch_order: int
    branch_id: int
    volume: float
    parent_id: int
    reverse_branch_order: int
    segment_id: int

    projected_data: dict(Projection) = field(default_factory=dict)
    flow_id: int() = None
    flow_type: str = None
    begins_at_drip_point: bool = None
    begins_at_divide_point: bool = None

    stem_path_id = int

    dx: float = 0
    dy: float = 0
    dz: float = 0

    surface_area: float = 0.0
    sa_to_vol: float = 0.0
    slope: float = 0.0

    is_stem: bool = False
    is_divide: bool = False

    # #blood for the blood god, software eng for the filter func
    # class_attrs = self.__get_class_attributes(type(self))
    # self.__init_instance(class_attrs, kwargs)

    def calc_surface_area(self):
        radius = self.radius
        length = self.length
        sa = 2 * np.pi * radius * (radius + length) - 2 * np.pi * radius * radius
        return sa

    def weight_dict(self, plane: str = "XZ"):
        attr_dict = {
            "cyl_id": self.cyl_id,
            "parent_id": self.parent_id,
            "x": self.x,
            "y": self.y,
            "z": self.z,
            "radius": self.radius,
            "length": self.length,
            "xz_area": self.xz_area,
            "polygon": self.projected_data[plane]["polygon"],
            "angle": self.angle,
            "surface_area": self.surface_area,
            "volume": self.volume,
            "sa_to_vol": self.sa_to_vol,
            "branch_order": self.branch_order,
        }
        return attr_dict

    def create_from_list(self, attrs: list, columns=qsm_cols):
        """creates a cylinder corrosponding to that defined by a given row of the qsm (attrs)"""

        extract = (
            lambda attr: attrs[columns[attr]]
        )  # pulls a column from the qsm row (attrs) corrosponding to the input attribute

        self.dx = self.x[1] - self.x[0]
        self.dy = self.y[1] - self.y[0]
        self.dz = self.z[1] - self.z[0]
        self.surface_area = self.calc_surface_area()
        self.sa_to_vol = self.surface_area / self.volume
        run = np.sqrt(self.dx**2 + self.dy**2)
        self.angle = (
            np.arctan(self.dz / np.sqrt(self.dx**2 + self.dy**2))
            if run > 0
            else np.arctan(0)
        )
        self.xz_area = None
        # log.info(str(self.__repr__()))

    def get_projection(self, plane="XZ"):
        noCirPoints = 360
        tCir = np.linspace(
            0, 2 * np.pi, noCirPoints
        )  # 360 evenly spaced points between 0 - 2pi (radian degrees)
        a_ortho = np.cos(tCir)  # x coordinates of the points on a circle
        b_ortho = np.sin(tCir)  # y coordinates of the points on a circle

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
        if plane == "XZ":
            self.xz_area = (
                self.projected_data["XZ"]["area"] 
            )
        return projection["polygon"]

    def draw(self, plane: str = "XZ"):
        poly = self.projected_data[plane]["polygon"]
        draw_cyls([poly])

    def get_flow_data():
        """Returns the flow ID and flow characteristics of the flow the cyl is contained in"""
        print("Get flow data not written")
