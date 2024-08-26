"""Defines the component parts of the ingested QSM"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np

from canopyhydro.configuration import qsm_cols
from canopyhydro.DataClasses import Projection
from canopyhydro.geometry import draw_cylinders_3D, draw_cyls, get_projection

log = logging.getLogger("model")

NAME = "Cylinder"


def create_cyl(arr: np.array):
    """Creates a cylinder based off of an input array.
        Enables user defined file structure via configuration yml file

    Args:
        arr (np.array):
            Array containing attribute data for a single cylinder.
            Attributes must be ordered as per qsm_cols in canopyhydro_config.yml
    """
    cols = qsm_cols
    attrs = {k: arr[v] for (k, v) in cols.items()}
    cyl = Cylinder(**attrs)
    # cyl.initialize(arr, cols)
    return cyl


@dataclass
class Cylinder:
    """
    The Cylinder class is used to represent the 3-D cylinders that make up a QSM.
    Contains several wrappers for functions in 'geometry'.

    Attributes:
        cyl_id (int): Desc.
        x (np.ndarray[np.float32]): Desc.
        y (np.ndarray[np.float32]): Desc.
        z (np.ndarray[np.float32]): Desc.
        radius (np.float32): Desc.
        length (np.float32): Desc.
        branch_order (int): Desc.
        branch_id (int): Desc.
        volume (np.float32): Desc.
        parent_id (int): Desc.
        reverse_branch_order (int): Desc.
        segment_id (int): Desc.
        projected_data (dict(Projection)): Desc. default_factory=dict
        flow_id (int): Desc. Defaults to None
        flow_type (str): Desc. Defaults to None
        drip_node (int): Desc. Defaults to None
        begins_at_drip_point (bool): Desc. Defaults to None
        begins_at_divide_point (bool): Desc. Defaults to None
        stem_path_id (int): Desc. Defaults to None
        dx (np.float32): Desc. Defaults to 0
        dy (np.float32): Desc. Defaults to 0
        dz (np.float32): Desc. Defaults to 0
        surface_area (np.float32): Desc. Defaults to 0.0
        sa_to_vol (np.float32): Desc. Defaults to 0.0
        slope (np.float32): Desc. Defaults to 0.0
        is_stem (bool): Desc. Defaults to False
        is_divide (bool): Desc. Defaults to False
    """

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
    flow_id: int = None
    flow_type: str = None
    drip_node: int = None
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

    def __post_init__(self):
        """Initializes the Cylinder object after all attributes are set"""
        self.initialize()

    def initialize(self):
        """Initializes remaining attributes based off of the
        attributes provided at object creation

        References:
            - self.x, self.y, self.z
            - self.dx, self.dy, self.dz
            - self.radius, self.length
            - self.volume
            - self.surface_area
            - self.sa_to_vol
            - self.angle
            - self.xy_area

        """

        self.dx = self.x[1] - self.x[0]
        self.dy = self.y[1] - self.y[0]
        self.dz = self.z[1] - self.z[0]
        radius = self.radius
        length = self.length
        self.surface_area = (
            2 * np.pi * radius * (radius + length) - 2 * np.pi * radius * radius
        )

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

    def draw(self, plane: str = "XY", **kwargs):
        """
        A wrapper around the draw_cyl function allowing
        more readable code for drawing a single cylinder
        e.g. For some Cylinder 'cyl' cyl.draw() replaces geometry.draw_cyls([cyl])

        Args:
            - plane:
                The projection of the cylinder to draw:
                'XY, 'XZ', or 'YZ'. Defaults to "XY".
        """
        poly = self.projected_data[plane]["polygon"]
        draw_cyls([poly], **kwargs)

    def draw_3D(self, **kwargs):
        """Draws the cylinder in 3D space"""
        vector_start_end = np.array(
            [
                np.array([self.x[0], self.y[0], self.z[0]]),
                np.array([self.x[1], self.y[1], self.z[1]]),
            ]
        )
        fig = draw_cylinders_3D([self.radius], [vector_start_end], **kwargs)
        return fig
