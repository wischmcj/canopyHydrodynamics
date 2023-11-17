from __future__ import annotations

from canhydro.global_vars import DIR, test_input_dir

import calendar
import os
import shutil
import stat
import time
import numpy as np

    
        # return np.argmax(dist_2)
# def test_create_cylinders(basic_forest):
#     actual = basic_forest.get_collection_data("1_TenCyls.csv")
#     expected = ten_cyls_rows
#     assert expected == actual


def test_sandpit():
        nodes = [(0,0,0), (1,1,1), (2,2,2),(3,3,3)]
        node = (1.5,1.5,1.6)
        nodes = np.asarray(nodes)
        dist_2 = np.sum((nodes - node) ** 2, axis=1)
        breakpoint()