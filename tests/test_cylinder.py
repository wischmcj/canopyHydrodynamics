from __future__ import annotations

import pytest


@pytest.mark.parametrize(
    "points,e_angle,e_area",
    [([(1, 4), (1, 6), (1, 7)], 0.785398), ([(-1, -4), (-1, -6), (-1, -7)], 1)],
)
def test_project_cyl():
    """
    Tests projection of cylinders parallel with:
     the XY plane, Z axis the line x=y (45 deg)
    """

    # numba_get_projection
    ang = 0.785398
