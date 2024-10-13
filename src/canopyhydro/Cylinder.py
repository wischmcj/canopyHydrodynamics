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
    Examples:
        >>> import numpy as np
        >>> from canopyhydro.configuration import qsm_cols
        >>> print(qsm_cols)
        {'cyl_id': 1, 'parent_id': 2, 'x': [3, 6], 'y': [4, 7], 'z': [5, 8],
        'radius': 9, 'volume': 10, 'length': 12, 'segment_id': 15, 'branch_order': 20,
        'reverse_branch_order': 21, 'branch_id': 24}
        >>> myCyl = create_cyl(
        >>>    np.array([np.nan, 1, 0, -2.372914, 2.875943, -0.55, -2.382034, 2.887813,
        >>>    -0.452896, 0.277545, 0.023777, 4.755711, 0.098251, 3880.839667,
        >>>    np.nan, 0, -1, 0.240459, 4.604282, 3880.04742, 0, 41, 84.320816,
        >>>    7110, 0, 0, np.nan, 0, 0, 0, 0.01, 0.401606, 0, 0.01, 0.401606])
        >>> )
        >>> print(myCyl)
        Cylinder( cyl_id=1.0, x=[-2.372914 -2.382034], y=[2.875943 2.887813], z=[-0.55     -0.452896], radius=0.277545, length=0.098251, branch_order=0.0, branch_id=0.0, volume=0.023777, parent_id=0.0, reverse_branch_order=41.0, segment_id=0.0
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
        cyl_id (int):
            The ID of the cylinder.
        x (np.ndarray[np.float32]):
            The x-coordinates of the cylinder vertices. Length 2, ordered (start, end)
        y (np.ndarray[np.float32]):
            The y-coordinates of the cylinder vertices. Length 2, ordered (start, end)
        z (np.ndarray[np.float32]):
            The z-coordinates of the cylinder vertices. Length 2, ordered (start, end)
        radius (np.float32):
            The radius of the cylinder.
        length (np.float32):
            The length of the cylinder.
        branch_order (int):
            The branch order of the cylinder. See QSM documentation for more information.
        branch_id (int):
            The ID of the branch the cylinder belongs to. See QSM documentation for more information.
        volume (np.float32):
            The volume of the cylinder.
        parent_id (int):
            The ID of the parent cylinder. -1 if cylinder is a root cylinder
        reverse_branch_order (int):
            The reverse branch order of the cylinder. See QSM documentation for more information.
        segment_id (int):
            The ID of the segment the cylinder belongs to. See QSM documentation for more information.
        projected_data (dict(Projection)):
            A dictionary containing projected data of the cylinder
            on the 'XY', 'XZ' and 'YZ' planes. Defaults to an empty dictionary.
        flow_id (int):
            The ID of the flow to which the cylinder contributes intercepted precipitation.
            Defaults to None.
        flow_type (str):
            The type of flow associated with the cylinder (drip or stem). Defaults to None.
        drip_node (int):
            The ID of the drip node associated with the cylinder.
            Defined by an ID correlating to the node_id of the node in the
            containing CylinderCollection's graph . Defaults to None.
        begins_at_drip_point (bool):
            Indicates if (x[0],y[0],z[0]) at a drip point. Defaults to None.
        begins_at_divide_point (bool):
            Indicates if (x[0],y[0],z[0]) at a dividepoint. Defaults to None.
        dx (np.float32):
            The change in x-coordinate: x[1] - x[0]. Defaults to 0.
        dy (np.float32):
            The change in y-coordinate: y[1] - y[0]. Defaults to 0.
        dz (np.float32):
            The change in z-coordinate: z[1] - z[0]. Defaults to 0.
        surface_area (np.float32):
            The surface area of the cylinder. Defaults to 0.0.
        sa_to_vol (np.float32):
            The ratio of surface area to volume of the cylinder. Defaults to 0.0.
        slope (np.float32):
            The slope of the cylinder, the rise over run in 3 dimensions. Defaults to 0.0.
        is_stem (bool):
            Indicates if the cylinder contributes to stem flow. Defaults to False.
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

    projected_data: dict[str, Projection] = field(default_factory=dict)
    flow_id: int = None
    flow_type: str = None
    drip_node: int = None
    begins_at_drip_point: bool = None
    begins_at_divide_point: bool = None

    dx: np.float32 = 0
    dy: np.float32 = 0
    dz: np.float32 = 0

    surface_area: np.float32 = 0.0
    sa_to_vol: np.float32 = 0.0
    slope: np.float32 = 0.0

    is_stem: bool = False

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

        References:
            projected_data.polygon (shapely.Polygon):
              represents the shape formed by projecting
              the cylinder (self) onto the given plane
            projected_data.base_vector (list[float]):
              cylinder vector *after projection
              (oriented away from cylinder base)
            projected_data.anti_vector (list[float]):
              cylinder vector *after projection
              (oriented towards cylinder base)
            projected_data.angle (float):
              angle (in radians) of cylinder vector with
              given axis. Between -pi/2 and pi/2.
            projected_data.area (float) - area of polygon (see above)
            projected_data.xy_area: float
              - polygon formed by projecting the cylinder
                  onto the XY plane

        Returns:
            (shapely.Polygon): projected_data[polygon]

        Examples:
            >>> import numpy as np
            >>> from canopyhydro.Cylinder import Cylinder
            >>> cyl = Cylinder(1, np.array([0, 1]), np.array([0, 1]), np.array([0, 1]), 1, 1, 0, 0, 1, 0, 0, 0)
            >>> cyl.get_projection("XY")
            >>> print(cyl.projected_data['XY']['polygon'])
            POLYGON ((0.7055547485197222 -0.,....)
            >>> print(cyl.projected_data['XY']['base_vector'])
            [0.57735027 0.57735027 0.57735027]
            >>> print(cyl.projected_data['XY']['anti_vector'])
            [-0.57735027 -0.57735027 -0.57735027]
            >>> print(cyl.projected_data['XY']['angle'])
            0.6154797086703873
            >>> print(cyl.projected_data['XY']['area'])
            4.642087600836025
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
        """A wrapper around the draw_cyl function allowing
        more readable code for drawing a single cylinder.
        e.g. For some Cylinder 'cyl' cyl.draw() replaces geometry.draw_cyls([cyl])

        Args:
            plane:
                The projection of the cylinder to draw.
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
