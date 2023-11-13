from __future__ import annotations

import geopandas as geo  # only import what we need
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from shapely.ops import unary_union as union

from canhydro.global_vars import log


def draw_cyls(collection: list[Polygon], colors: list[bool] = None):
    log.info("Plotting cylinder collection")
    fig, ax = plt.subplots()
    geoPolys = geo.GeoSeries(collection)
    colors = ["Blue" if col else "Grey" for col in colors]
    geoPolys.plot(ax=ax, color=colors)
    plt.show()


def unary_union(polys):
    return union(polys)
