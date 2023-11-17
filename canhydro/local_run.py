# from __future__ import annotations

import calendar
import os
import shutil
import stat
import time
import sys
# import numpy as np

# from canhydro.Forester import Forester
# from canhydro.global_vars import log, test_input_dir
# from memory_profiler import profile


# import pandas as pd
import numpy as np
import geopandas as geo
from descartes import PolygonPatch
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.getcwd()))
import alphashape

if __name__ == "__main__":
      points_2d = [(0., 0.), (0., 1.), (1., 1.), (1., 0.),
          (0.5, 0.25), (0.5, 0.75), (0.25, 0.5), (0.75, 0.5)]
      alpha_shape =  alphashape.alphashape(points_2d, 2.0)# alphashape.alphashape(points_2d, 0.)
      [Point]
      fig, ax = plt.subplots()
      ax.scatter(*zip(*points_2d))
#     ax.add_patch(PolygonPatch(alpha_shape))
      geopoly = geo.GeoSeries(alpha_shape)
      geopoly.plot(ax =ax)
      plt.show()
      breakpoint()