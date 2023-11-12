"""Defines the component parts of the ingested QSM"""

from __future__ import annotations

import calendar
import os
from collections import defaultdict
from dataclasses import dataclass, field
from math import isnan, pi, sqrt
from multiprocessing import Pool
from pathlib import Path
from random import random
from time import sleep

import geopandas as geo
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.pyplot import cm
from shapely.geometry import Point, Polygon
from shapely.ops import unary_union

from canhydro.DataClasses import Projection, Model
from canhydro.global_vars import input_dir, log, output_dir, qsm_cols
from canhydro.Plotter import draw_cyls

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
    cyl_id: int()
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

    # #blood for the blood god, software eng for the filter func
    # class_attrs = self.__get_class_attributes(type(self))
    # self.__init_instance(class_attrs, kwargs)

    def calc_surface_area(self):
        radius = self.radius
        length = self.length
        sa = 2 * np.pi * radius * (radius + length) - 2 * np.pi * radius * radius
        return sa

    def weight_dict(self, cols: str = ""):
        if cols == "":
            attr_dict = {
                "x": self.x,
                "y": self.y,
                "z": self.z,
                "radius": self.radius,
                "unit_vector": self.projected_data[""],
                "length": self.length,
                "surface_area": self.surface_area,
                "volume": self.volume,
                "sa_to_vol": self.sa_to_vol,
                "branch_order": self.branch_order,
                "branch_id": self.branch_id,
                "parent_id": self.parent_id,
                "reverse_branch_order": self.reverse_branch_order,
                "segment_id": self.segment_id,
                "angle": self.angle,
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
        self.angle = np.arctan(self.dz / np.sqrt(self.dx**2 + self.dy**2))
        # log.info(str(self.__repr__()))

    def meets_criteria(
        self,
        plane: str = "",
        highlight: bool = False,
        color: str = "",
        branch_order: int = -1,
        min_radius: int = -1,
        max_radius: int = -1,
        branch_id: int = -1,
        segment_id: int = -1,
    ) -> str:
        order_match = (
            True if branch_order == -1 else (self.branch_order == branch_order)
        )
        radius_match = (
            self.radius >= min_radius
            if max_radius == -1
            else (self.radius >= min_radius and self.radius <= max_radius)
        )
        branch_id_match = True if branch_id == -1 else (self.branch_id == branch_id)
        seg_id_match = True if segment_id == -1 else self.segment_id == segment_id
        proj_match = self.projected_data.get(plane, "") != ""
        match = (
            order_match
            and radius_match
            and branch_id_match
            and seg_id_match
            and proj_match
        )
        return_color = color if color else "blue"
        to_return = return_color if match else "grey" if highlight else ""
        return to_return

    def get_projection(self, plane="XY"):
        noCirPoints = 360
        tCir = np.linspace(
            0, 2 * np.pi, noCirPoints
        )  # 360 evenly spaced points between 0 - 2pi (radian degrees)
        a_ortho = np.cos(tCir)  # x coordinates of the points on a circle
        b_ortho = np.sin(tCir)  # y coordinates of the points on a circle

        if plane == "XY":
            delt_a = self.dx
            delt_b = self.dy
            delt_c = self.dz
            dim_a = np.transpose(self.x)
            dim_b = np.transpose(self.y)
            dim_c = np.transpose(self.z)
        elif plane == "XZ":
            delt_a = self.dx
            delt_b = self.dz
            delt_c = self.dy
            dim_a = np.transpose(self.x)
            dim_b = np.transpose(self.z)
            dim_c = np.transpose(self.y)
        else:
            delt_a = self.dy
            delt_b = self.dz
            delt_c = self.dx
            dim_a = np.transpose(self.y)
            dim_b = np.transpose(self.z)
            dim_c = np.transpose(self.x)

        # unit vector at base of cylinder, pointing up cylinder axis
        vNorm = np.sqrt(delt_a**2 + delt_b**2 + delt_c**2)
        aV = np.hstack((delt_a, delt_b, delt_c)) / vNorm
        bV = -aV  # unit vector looking down from top circle (but not translated)

        # function to find orthgonal vectors
        oVz = lambda v, a, b: ((-v[0] * a - v[1] * b) / v[2])

        # initializing min max arrays+
        min_c = np.zeros_like(delt_c)
        max_c = np.zeros_like(delt_c)

        pSV = []
        try:
            # for each cylinder
            if not np.isnan(dim_a[0]):
                if np.logical_and(delt_a == 0, delt_b == 0):
                    pX = dim_a[0] + self.radius * a_ortho
                    pY = dim_b[0] + self.radius * b_ortho
                    cPS = Polygon(list(zip(pX, pY)))
                    min_c = np.min(dim_c[:])
                    max_c = np.max(dim_c[:])
                else:
                    # find orthogonal vectors @ endpoints
                    # Identifies corners of projected rectangle
                    aVp1 = np.hstack((aV[1], -aV[0]))
                    aVp2 = np.hstack((-aV[1], aV[0]))
                    bVp1 = np.hstack((bV[1], -bV[0]))
                    bVp2 = np.hstack((-bV[1], bV[0]))

                    aVp1 = aVp1 / np.linalg.norm(aVp1)
                    aVp2 = aVp2 / np.linalg.norm(aVp2)
                    bVp1 = bVp1 / np.linalg.norm(bVp1)
                    bVp2 = bVp2 / np.linalg.norm(bVp2)

                    # from each endpoint, use radius to find vertices of the rectangle
                    x1 = dim_a[0] + self.radius * aVp1[0]
                    y1 = dim_b[0] + self.radius * aVp1[1]
                    x2 = dim_a[0] + self.radius * aVp2[0]
                    y2 = dim_b[0] + self.radius * aVp2[1]
                    x3 = dim_a[1] + self.radius * bVp1[0]
                    y3 = dim_b[1] + self.radius * bVp1[1]
                    x4 = dim_a[1] + self.radius * bVp2[0]
                    y4 = dim_b[1] + self.radius * bVp2[1]

                    # calculate set of orthgonal vectors using lambda function
                    # That is 360 orthogonal vectors ending at eqidistant points along
                    # a circle of radius self.radius with the starting point of our cylinder
                    # at is center
                    ZOrtho = oVz(aV[:], a_ortho, b_ortho)

                    # unit-ify the orthgonal vectors
                    uovd = np.sqrt(a_ortho**2 + b_ortho**2 + ZOrtho**2)
                    # breakpoint()
                    # Confounded - why does removing the first three [:,None]'s below lead to non-circular projections
                    # for XZ?
                    uov = (
                        np.hstack((a_ortho[:, None], b_ortho[:, None], ZOrtho[:, None]))
                        / uovd[:, None]
                    )

                    # donot re unit-fy, you only want the horizontal component, not the
                    # renormalized horizontal component

                    # using only the X and Y components, find circle coods in plane of
                    # interest
                    xaC = dim_a[0] + uov[:, 0] * self.radius
                    yaC = dim_b[0] + uov[:, 1] * self.radius
                    zaC = dim_c[0] + uov[:, 2] * self.radius

                    xbC = dim_a[1] + uov[:, 0] * self.radius
                    ybC = dim_b[1] + uov[:, 1] * self.radius
                    zbC = dim_c[1] + uov[:, 2] * self.radius

                    min_c = np.min(np.vstack((zaC, zbC)))
                    max_c = np.max(np.vstack((zaC, zbC)))

                    # assymble total package
                    rX = np.vstack((x1, x2, x3, x4))
                    rY = np.vstack((y1, y2, y3, y4))

                    # test for circle parts in polygon
                    try:
                        c1 = Polygon(
                            list(
                                zip(
                                    [0 if isnan(x) else x for x in xaC],
                                    [0 if isnan(y) else y for y in yaC],
                                )
                            )
                        )
                        bBox = Polygon(
                            list(
                                zip(
                                    [0 if isnan(x) else x for x in rX],
                                    [0 if isnan(y) else y for y in rY],
                                )
                            )
                        )
                        c2 = Polygon(
                            list(
                                zip(
                                    [0 if isnan(x) else x for x in xbC],
                                    [0 if isnan(y) else y for y in ybC],
                                )
                            )
                        )

                        partsPS = [c1, bBox, c2]
                    except:
                        log.info(f"Error creating projection polygons")

                    try:
                        cPS = unary_union(partsPS)
                    except:
                        print(np.any(np.isnan(xaC)))
                        log.info(f"Error unioning projection polygons")
                        print(yaC)
                        print(rX)
                        print(rY)
                        print(xbC)
                        print(ybC)
                    # get angle away from plane projected on to
                    run = sqrt(delt_b**2 + delt_a**2)
                    rise = delt_c
                    if run == 0:
                        slope = 1  # straightDown e.g. is in flow
                    else:
                        slope = rise / run
                    ang = np.arctan(slope)
                    projection = {
                        "polygon": cPS,
                        "base_vector": aV,
                        "anti_vector": bV,
                        "angle": ang,
                    }
                    self.projected_data[plane] = projection
            else:
                log.info(str(self.__dict__))
                breakpoint()
                log.info(f"dim_a[0] is null, unable to project cylinder {self.cyl_id}")
            return cPS
        except UnboundLocalError:
            log.info(f"UnboundLocalError for {self.__dict__}")

    def draw(self, plane: str = "XY"):
        poly = self.projected_data[plane]["poly"]
        draw_cyls([poly])

    def get_flow_data():
        """Returns the flow ID and flow characteristics of the flow the cyl is contained in"""
        print("Get flow data not written")

    
@dataclass
class CylinderList(Cylinder,Model):

    def __init__(self, **kwargs):
        super(CylinderList,self).__init__(**kwargs)
