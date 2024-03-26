from __future__ import annotations

import os
import sys

print(sys.path)
import pytest

sys.path.insert(0, os.path.dirname(os.getcwd()))
sys.path.insert(0, os.getcwd())
print(sys.path)

# from test.expected_results_shapes import (small_tree_overlap,
#                                           small_tree_wateshed_poly)
# from test.utils import within_range

from src.canhydro.global_vars import DIR, test_input_dir
from src.canhydro.utils import lam_filter
from src.canhydro.CylinderCollection import CylinderCollection
from src.canhydro.Forester import Forester


if __name__ == "__main__":
    forest = Forester()
    forest.get_file_names(dir=test_input_dir)
    forest.qsm_from_file_names(file_name="5_SmallTree")
    basic_collection = forest.cylinder_collections[0]

    # forest_old = Forester()
    # forest_old.get_file_names(dir=test_input_dir)
    # forest_old.qsm_from_file_names(file_name="5_SmallTree.csv")
    # basic_collection_old = forest_old.cylinder_collections[0]


    basic_collection.project_cylinders("XY")
    # basic_collection_old.project_cylinders("XY")

    basic_collection.initialize_digraph_from()
    # basic_collection_old.old_initialize_digraph_from()

    basic_collection.find_flow_components()
    # basic_collection_old.old_find_flow_components()
    
    print(basic_collection.drip_summary())
    # print(basic_collection_old.drip_summary())


    basic_collection.calculate_flows()
    # actual_flows = basic_collection.flows
    # _, actual_stem_map = lam_filter(
    #     basic_collection.cylinders, lambda: is_stem, return_all=True
    # )