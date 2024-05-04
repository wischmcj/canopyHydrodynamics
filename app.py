from __future__ import annotations

import os
import sys

import debugpy
debugpy.listen(("0.0.0.0", 5678))
# import geopandas as geo
# import matplotlib.pyplot as plt

# import numpy as np

from src.canhydro.Forester import Forester
from src.canhydro.global_vars import log, test_input_dir


sys.path.insert(0, os.path.dirname(os.getcwd()))

def initialize_forester(dir, file = None):
    forester = Forester()
    forester.get_file_names(dir=test_input_dir)
    forester.qsm_from_file_names(file_name=request.param)
    return forester

if __name__ == "__main__":
    forest = initialize_forester(test_input_dir,"5_SmallTree.csv")
    tree = forest.cylinder_collections[0]
    print("Waiting for client to attach...")
    debugpy.wait_for_client()
    # breakpoint()