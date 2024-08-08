"""Defines the component parts of the ingested QSM"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
import toml

from src.canhydro.DataClasses import Projection
from src.canhydro.geometry import draw_cyls, get_projection

log = logging.getLogger("model")

with open("src/canhydro/user_def_config.toml") as f:
    config = toml.load(f)

# QSM column order
qsm_cols = {}
for column in config["qsm"]:
    qsm_cols[column] = config["qsm"][column]


NAME = "Cylinder"


def create_cyl(arr: np.array):
    """Creates a cylinder based off of an input array.
        Enables user defined file structure via configuration yml file 

    Args:
        arr (np.array): 
            Array containing attribute data for a single cylinder.
            Attributes must be ordered as per qsm_cols in user_def_config.yml
    """
    cols = qsm_cols
    attrs = {k: arr[v] for (k, v) in cols.items()}
    cyl = Cylinder(**attrs)
    cyl.create_from_list(arr, cols)
    return cyl


@dataclass
class Cylinder:
    cyl_id: int
    x: np.ndarray[np.float32]
    y: np.ndarray[np.float32]
    z: np.ndarray[np.float32]
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

    def __repr__(self):
        """Defines a readable string to represent a cylinder object
            *Utilized in tests to compare expected and actual results
        """
        return f"Cylinder( cyl_id={self.cyl_id}, x={self.x}, y={self.y}, z={self.z}, radius={self.radius}, length={self.length}, branch_order={self.branch_order}, branch_id={self.branch_id}, volume={self.volume}, parent_id={self.parent_id}, reverse_branch_order={self.reverse_branch_order}, segment_id={self.segment_id}"

    def __eq__(self, other):
        """Defines the minimum requirements for equality between two cylinders. 
            Allows use of '=' for comparing cylinders
        """
        return type(self) == type(other) and self.__repr__() == other.__repr__()

    def create_from_list(self, attrs: list, columns=qsm_cols):
        """Initializes remaining attributes based off of the 
                attributes provided at object creation
        """

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
        radius = self.radius
        length = self.length
        self.surface_area = 2 * np.pi * radius * (radius + length) - 2 * np.pi * radius * radius

        self.sa_to_vol = 0 if self.volume == 0 else self.surface_area / self.volume
        run = np.sqrt(self.dx**2 + self.dy**2)
        self.angle = (
            np.arctan(0)
            if run == 0
            else np.arctan(self.dz / np.sqrt(self.dx**2 + self.dy**2))
        )
        self.xy_area = 0

    def get_projection(self, plane="XY"):
        """Calculates the projection of the cylinder on the requested plane 

        Args:
            self (Cylinder): 
                Cylinder to be projected 
            plane (str, optional): 
                Plane on which to project the cylinder. Defaults to "XY".
        
        Sets Attributes:
                - projected_data (dict)
                    - key->value:
                        - polygon (shapely.Polygon):
                            represents the shape formed by projecting 
                            the cylinder (self) onto the given plane
                        - base_vector (list[float]): 
                            cylinder vector *after projection 
                            (oriented away from cylinder base)
                        - anti_vector (list[float]):
                            cylinder vector *after projection 
                            (oriented towards cylinder base)
                        - angle (float):
                            angle (in radians) of cylinder vector with 
                            given axis. Between -pi/2 and pi/2.
                        - area (float) - area of polygon (see above)
                - xy_area: float 
                    - polygon formed by projecting the cylinder
                        onto the XY plane

        Returns:
            (shapely.Polygon): projected_data[polygon]

        """
        if plane == "XY":
            magnitude = [self.dx, self.dy, self.dz]
            ranges = [self.x, self.y, self.z]
        elif plane == "XZ":
            magnitude = [self.dx, self.dz, self.dy]
            ranges = [self.x, self.z, self.y]
        else:
            magnitude = [self.dy, self.dz, self.dx]
            ranges = [self.y, self.z, self.x]

        vector = list(map(np.transpose, ranges))

        projection = get_projection(vector, magnitude, self.radius)
        self.projected_data[plane] = projection
        if plane == "XY":
            self.xy_area = self.projected_data["XY"]["area"]
        return projection["polygon"]

    def draw(self, plane: str = "XY"):
        """A wrapper around the draw_cyl function allowing 
            more readable code for drawing a single cylinder 
                - e.g. for some Cylider 'cyl'
                    cyl.draw() rather than geometry.draw_cyls([cyl]) """
        poly = self.projected_data[plane]["polygon"]
        draw_cyls([poly])
