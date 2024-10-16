"""Defines the component parts of the ingested QSM"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
import toml
from src.canhydro.DataClasses import Projection
from src.canhydro.geometry import (draw_cyls, get_projection,
                                   get_projection_scikit)

log = logging.getLogger("model")

with open("src/canhydro/user_def_config.toml") as f:
    config = toml.load(f)
# QSM column order
qsm_cols = {}
for column in config["qsm"]:
    qsm_cols[column] = config["qsm"][column]


NAME = "Cylinder"


def create_cyl(arr: np.array):
    cols = qsm_cols
    attrs = {k: arr[v] for (k, v) in cols.items()}
    cyl = Cylinder(**attrs)
    cyl.create_from_list(arr, cols)
    return cyl


def create_empty_cyl():
    cols = qsm_cols
    attrs = {
        "cyl_id": -1,
        "x": 0,
        "y": 0,
        "z": 0,
        "radius": 0,
        "length": 0,
        "branch_order": 0,
        "branch_id": 0,
        "volume": 0,
        "parent_id": -1,
        "reverse_branch_order": 0,
        "segment_id": 0,
    }
    cyl = Cylinder(**attrs)
    cyl.create_from_list(attrs.values(), cols)
    return cyl


@dataclass
class Cylinder:
    cyl_id: int
    x: np.ndarray[np.float16]
    y: np.ndarray[np.float16]
    z: np.ndarray[np.float16]
    radius: np.float16
    length: np.float16
    branch_order: int
    branch_id: int
    volume: np.float16
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

    dx: np.float16 = 0
    dy: np.float16 = 0
    dz: np.float16 = 0

    surface_area: np.float16 = 0.0
    sa_to_vol: np.float16 = 0.0
    slope: np.float16 = 0.0

    is_stem: bool = False
    is_divide: bool = False

    def __repr__(self):
        return f"Cylinder( cyl_id={self.cyl_id}, x={self.x}, y={self.y}, z={self.z}, radius={self.radius}, length={self.length}, branch_order={self.branch_order}, branch_id={self.branch_id}, volume={self.volume}, parent_id={self.parent_id}, reverse_branch_order={self.reverse_branch_order}, segment_id={self.segment_id}"

    def __eq__(self, other):
        return type(self) == type(other) and self.__repr__() == other.__repr__()

    def __post_init__(self):
        self.cyl_id = int(self.cyl_id)

    def calc_surface_area(self):
        radius = self.radius
        length = self.length
        sa = 2 * np.pi * radius * (radius + length) - 2 * np.pi * radius * radius
        return sa

    def create_from_list(self, attrs: list, columns=qsm_cols):
        """creates a cylinder corrosponding to that defined by a given row of the qsm (attrs)"""

        extract = lambda attr: attrs[columns[attr]]
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
        log.debug(str(self.__repr__()))

    def get_flow_arr(self, plane="XY"):
        return np.array(
            [
                1,
                np.float16(self.projected_data["XY"]["area"]),
                self.surface_area,
                self.angle,
                self.volume,
                self.sa_to_vol,
            ]
        )

    def get_loc(self):
        return (self.x[1], self.y[1], self.z[1])

    def get_projection(self, plane="XY", scikit=False):
        noCirPoints = 360
        tCir = np.linspace(
            0, 2 * np.pi, noCirPoints
        )  # 360 evenly spaced points between 0 - 2pi (radian degrees)

        if plane == "XY":
            magnitude = [self.dx, self.dy, self.dz]
            vector = [np.transpose(self.x), np.transpose(self.y), np.transpose(self.z)]
        elif plane == "XZ":
            magnitude = [self.dx, self.dz, self.dy]
            vector = [np.transpose(self.x), np.transpose(self.z), np.transpose(self.y)]
        else:
            magnitude = [self.dy, self.dz, self.dx]
            vector = [np.transpose(self.y), np.transpose(self.z), np.transpose(self.x)]

        if scikit:
            from sklearn.preprocessing import StandardScaler

            scaler = StandardScaler()
            vector = scaler.fit_transform(vector)
            magnitude = scaler.fit_transform(magnitude)
        if scikit:
            projection = get_projection_scikit(vector, magnitude, self.radius)
        else:
            projection = get_projection(vector, magnitude, self.radius)
        self.projected_data[plane] = projection
        if plane == "XY":
            self.xy_area = self.projected_data["XY"]["area"]
        return projection["polygon"]

    def draw(self, plane: str = "XY"):
        poly = self.projected_data[plane]["polygon"]
        draw_cyls([poly])
