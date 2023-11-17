from __future__ import annotations

import calendar
import os
import shutil
import stat
import time
import numpy as np

from canhydro.Forester import Forester
from canhydro.global_vars import log, test_input_dir
from memory_profiler import profile


if __name__ == "__main__":
    nodes = [(0,0,0), (1,1,1), (2,2,2),(3,3,3)]
    node = (1.5,1.5,1.6)
    nodes = np.asarray(nodes)
    dist_2 = np.sum((nodes - node) ** 2, axis=1)
    breakpoint()
    # return np.argmax(dist_2)