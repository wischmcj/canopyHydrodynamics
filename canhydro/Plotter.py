from __future__ import annotations

import geopandas as geo  # only import what we need
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
from shapely.ops import transform, unary_union

from canhydro.global_vars import log, qsm_cols
from canhydro.utils import intermitent_log


def draw_cyls(collection: list[Polygon], colors: list[str] = None):
    log.info(f"Plotting cylinder collection")
    fig, ax = plt.subplots()
    geoPolys = geo.GeoSeries(collection)
    breakpoint()
    geoPolys.plot(ax=ax, color=colors)
    plt.show()
