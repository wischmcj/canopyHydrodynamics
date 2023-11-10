from __future__ import annotations

import matplotlib
import geopandas as geo  # only import what we need
from shapely.geometry import Point, Polygon
import matplotlib.pyplot as plt
from shapely.ops import transform, unary_union



from canhydro.Cylinder import Cylinder
from canhydro.DataClasses import Model
from canhydro.global_vars import log, qsm_cols
from canhydro.utils import intermitent_log


def draw_cyls(
    collection: list[Cylinder]
):
    log.info(f"Plotting cylinder collection")
    fig, ax = plt.subplots()
    geoPolys = geo.GeoSeries(collection)
    breakpoint()
    geoPolys.plot(ax=ax)
    plt.show()