from __future__ import annotations

import os
import shutil
import stat

import global_vars
from Forester import Forester
from memory_profiler import profile

DIR = global_vars.DIR
test_input_dir = global_vars.test_input_dir


@profile
class Townie:
    def on_rm_error(func, path, exc_info):
        # path contains the path of the file that couldn't be removed
        # let's just assume that it's read-only and unlink it.
        os.chmod(path, stat.S_IWRITE)
        os.unlink(path)

    def create_dir_and_file(filename) -> None:
        print(type(filename))
        os.makedirs(filename, exist_ok=True)
        f = open(r"test\demofile2.csv", "w")
        f.write("Now the file has more content!")
        f.close()

    def del_dir(filename) -> None:
        shutil.rmtree(filename, onerror=on_rm_error)

    @profile
    def min_graph_test(flexible_collection):
        flexible_collection.initialize_minimal_graph()
        proj_area = flexible_collection.sum_over_min_graph()
        flexible_collection.find_flow_components_minimal()
        print(proj_area)

    @profile
    def base_graph_test(flexible_collection):
        flexible_collection.initialize_graph()
        proj_area = flexible_collection.sum_over_graph()
        flexible_collection.find_flow_components()
        print(proj_area)

    @profile
    def obj_graph_test(flexible_collection):
        flexible_collection.initialize_object_graph()
        proj_area = flexible_collection.sum_over_object_graph()
        flexible_collection.find_flow_components_object()
        print(proj_area)

    if __name__ == "__main__":
        forest = Forester()
        forest.get_file_names(dir=test_input_dir)
        forest.qsm_from_file_names(file_name="4_LargeCollection.csv")
        # forest.qsm_from_file_names(file_name="3_HappyPathWTrunk.csv")
        forest.cylinder_collections[0].project_cylinders("XZ")
        # base_collection = forest.cylinder_collections[0]
        # min_collection = forest.cylinder_collections[0]
        obj_collection = forest.cylinder_collections[0]
        # base_graph_test(base_collection)
        # min_graph_test(min_collection)
        obj_graph_test(obj_collection)
