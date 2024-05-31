from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.getcwd()))
sys.path.insert(0, os.getcwd())

if __name__ == "__main__":
    points_2d = [
        (0.0, 0.0),
        (0.0, 1.0),
        (1.0, 1.0),
        (1.0, 0.0),
        (0.5, 0.25),
        (0.5, 0.75),
        (0.25, 0.5),
        (0.75, 0.5),
    ]
    print("test")
    # alpha_shape = alphashape.alphashape(
    #     points_2d, 2.0
    # )  # alphashape.alphashape(points_2d, 0.)
    # breakpoint()
    # fig, ax = plt.subplots()
    # ax.scatter(*zip(*points_2d))
    # #     ax.add_patch(PolygonPatch(alpha_shape))
    # geopoly = geo.GeoSeries(alpha_shape)
    # geopoly.plot(ax=ax) 
    # plt.show()
    # breakpoint()
