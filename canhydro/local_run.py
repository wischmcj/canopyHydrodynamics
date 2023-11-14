from __future__ import annotations

import calendar
import os
import shutil
import stat
import time

from Forester import Forester
from global_vars import log, test_input_dir
from memory_profiler import profile


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
    def min_graph_test():
        forest = Forester()
        forest.get_file_names(dir=test_input_dir)
        forest.qsm_from_file_names(file_name="4_LargeCollection.csv")
        flexible_collection = forest.cylinder_collections[0]
        flexible_collection.project_cylinders("XZ")
        flexible_collection.initialize_minimal_graph_from()
        proj_area = flexible_collection.sum_over_min_graph()
        flexible_collection.find_flow_components_minimal()
        print(proj_area)

    @profile
    def base_graph_test():
        forest = Forester()
        forest.get_file_names(dir=test_input_dir)
        forest.qsm_from_file_names(file_name="4_LargeCollection.csv")
        flexible_collection = forest.cylinder_collections[0]
        flexible_collection.project_cylinders("XZ")
        flexible_collection.initialize_graph_from()
        proj_area = flexible_collection.sum_over_graph()
        flexible_collection.find_flow_components()
        print(proj_area)

    @profile
    def obj_graph_test():
        forest = Forester()
        forest.get_file_names(dir=test_input_dir)
        forest.qsm_from_file_names(file_name="4_LargeCollection.csv")
        flexible_collection = forest.cylinder_collections[0]
        flexible_collection.project_cylinders("XZ")
        flexible_collection.initialize_object_graph_from()
        proj_area = flexible_collection.sum_over_object_graph()
        flexible_collection.find_flow_components_object()
        print(proj_area)

    if __name__ == "__main__":
        current_GMT = time.gmtime()
        time_stamp = str(calendar.timegm(current_GMT))
        log.info(f"base_graph_test started at {time_stamp}")

        base_graph_test()

        current_GMT = time.gmtime()
        time_stamp = str(calendar.timegm(current_GMT))
        log.info(f"min_graph_test started at {time_stamp}")

        min_graph_test()

        current_GMT = time.gmtime()
        time_stamp = str(calendar.timegm(current_GMT))
        log.info(f"obj_graph_test started at {time_stamp}")

        obj_graph_test()

        current_GMT = time.gmtime()
        time_stamp = str(calendar.timegm(current_GMT))
        log.info(f"finished at {time_stamp}")
