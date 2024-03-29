from __future__ import annotations

import os
import sys
import multiprocessing as mp
from time import time 
from itertools import product
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
from src.canhydro.CylinderCollection import CylinderCollection, pickle_collection, unpickle_collection
from src.canhydro.Forester import Forester
from test.utils import within_range

from src.canhydro.global_vars import config_vars, log, output_dir


def pickle(collection, designation = ""):
    log.info(f"attempting to pickle for file {collection.file_name}, des {designation}")
    # pickle_file = pickle_collection(collection, designation)
    try: 
        pickle_file = pickle_collection(collection, designation)
        log.info(f"successfully created pickle for file {collection.file_name}")
    except Exception as e:
        print(f"Error pickling file {collection.file_name}: {e}")
        return

def initialize_collection(file = "5_SmallTree"):
    log.info(f"Initializing collection...{file}")
    # forest = Forester()
    # forest.get_file_names(dir=test_input_dir)
    # forest.qsm_from_file_names(file_name=file)
    # basic_collection = forest.cylinder_collections[0]
    # basic_collection.project_cylinders("XY")
    try:
        forest = Forester()
        forest.get_file_names(dir=test_input_dir)
        forest.qsm_from_file_names(file_name=file)
        basic_collection = forest.cylinder_collections[0]
        basic_collection.project_cylinders("XY")
    except Exception as e:
        print(f"Error initializing collection for file {file}: {e}")
        return None
    log.info(f"Successfully initialized  collection {file}")
    return basic_collection

def prep_for_stats(collection, case_angle, case_name):
    log.info(f"attempting to prep for stats for case {case_name}")
    try: 
        collection.initialize_digraph_from(in_flow_grade_lim=case_angle)
        collection.find_flow_components()
        collection.calculate_flows()
    except Exception as e:
        print(f"Error preping for  stats for case {case_name}: {e}")
        return None
    log.info(f"successfully prepped for stats {case_name}")
    return True
    
def generate_statistics(collection, case_name):
    log.info(f"attempting to generate stats for file {collection.file_name}, case_name {case_name}")
    # statistics = collection.statistics(file_ext = case_name)
    try: 
        statistics = collection.statistics(file_ext = case_name)
    except Exception as e:
        print(f"Error gernerating stats for case {case_name} : {e}")
        return None
    log.info(f"attempting to generate flow file for {case_name}")
    try: 
        collection.generate_flow_file(file_ext = case_name)
    except Exception as e:
        print(f"Error gernerating flow file for case {case_name}: {e}")
        return None
    log.info("successfully created flow and stats files")
    return True

def run_test_cases(cases_to_run):

    start = time()
    collection = None
    case_name = f"inital_case_name"
    for case in cases_to_run:
        file_name, angle = case
        if not collection:
            collection = initialize_collection(file_name)
            if not collection:  
                dur = time() - start
                return None, f'{file_name}_{case_name}', dur
        log.info(f" Running case {file_name}_{angle}")
        case_name = f"{angle}"
        log.info(f"running case {file_name}_{case_name}")
        preped = prep_for_stats(collection, angle, case_name)
        if preped:
            generate_statistics(collection, case_name)
        else:
            pickle(collection,f'_prep_{case_name}')
            dur = time() - start
            return None, f'{file_name}_{case_name}',dur
        log.info(f"successfully ran case {case_name}")
        pickle(collection,f'_stats_{case_name}')
        collection=None
        dur = time() - start
        return True, f'{file_name}_{case_name}',dur
    

from data.output.run_so_far import already_run
run_cases = already_run

angles = set([tup[1] for tup in already_run])

def get_cases(file_names, already_run, angles_to_tests):
    cases = product(file_names,angles_to_tests)
    return [case for case in cases if case not in already_run]

def sensitivity_analysis():
    files_to_test = ["5_SmallTree"]
    files_to_test = [#"Secrest27-05_000000","Secrest32-06_000000"]
                         "Secrest02-26_000000"
                        ,"Secrest02-30_000000"
                        ,"Secrest03-12_000000"
                        ,"Secrest07-32_000000"
                        ,"Secrest08-24c_000000"
                        ,"Secrest10-02_000000"
                        ,"Secrest10-08_000000"
                        ,"Secrest11-27_000000"
                        ,"Secrest14-09_000000"
                        ,"Secrest16-3TI-CO_000000"
                        ,"Secrest16-14LI-ST_000000"
                        ,"Secrest18-13_000000"
                        ,"Secrest23-23_000000"
                        ,"Secrest24-03_000000"
                        ,"Secrest24-07_000000"
                        ,"Secrest26-03_000000"
                        ,"Secrest28-31_000000"
                        ,"Secrest29-20_000000"
                        ,"Secrest29-25_000000"
                        ,"Secrest31-05_000000"
                        ,"Secrest32-01_000000"
                        ,"Secrest32-03_000000"
                        ,"Secrest32-06_000000"
                        ,"Secrest32-14_000000"]
    cases_to_run = get_cases(files_to_test,run_cases,angles)
    log.info(f'Will run {len(cases_to_run)} cases : {cases_to_run}')
    start = time()
    # for file in files_to_test:
    #     success = run_test_cases(file)
    #     if not success:
    #         log.info(f"Failed run cases")
    #     else:
    #         log.info(f"suceeded running cases")
    # dur = time() - start

    # log.info(f"total time old method - {dur}")


    with mp.Pool(2) as p:
        task_pool = [p.apply_async(run_test_cases, args=(cases_to_run,)) for file in files_to_test]
        results = [task.get() for task in task_pool]

    for success, case, dur in results:
        log.info(f"total time running {case} - {dur}")
        if not success:
            log.info(f"Failed run case {case}")
        else:
            log.info(f"suceeded running cases {case}")


if __name__ == "__main__":
    sensitivity_analysis()
    
    
    # forest = Forester()
    # forest.get_file_names(dir=test_input_dir)
    # forest.qsm_from_file_names(file_name="5_SmallTree")
    # basic_collection = forest.cylinder_collections[0]

    # forest_old = Forester()
    # forest_old.get_file_names(dir=test_input_dir)
    # forest_old.qsm_from_file_names(file_name="5_SmallTree.csv")
    # basic_collection_old = forest_old.cylinder_collections[0]

    # basic_collection_old.project_cylinders("XY")

    # angles = [ (cyl.cyl_id,cyl.angle) for cyl in basic_collection.cylinders]
    # new_angles = [ (cyl.cyl_id,cyl.angle) for cyl in basic_collection.cylinders]
    # for idx, angle in angles:
    #     new_angle = [angle for cyl, angle in new_angles if cyl == idx][0]
    #     try:
    #         assert within_range(angle, new_angle, .03)
    #     except Exception as e:    
    #         print("Failure in projection {e}")
    #         breakpoint()
    #         #here 
    # # breakpoint()

    # # basic_collection.initialize_digraph_from(in_flow_grade_lim=-0.3)
    # basic_collection_old.initialize_digraph_from()

    # # nodes = [n for n in basic_collection_old.digraph.nodes()]
    # # # edges_new = [ edge for edge in basic_collection.digraph.nodes()]
    # # for node in nodes:
    # #     neighbors = basic_collection_old.digraph.edges(node)
    # #     new_neighbors = basic_collection.digraph.edges(node)
    # #     try:
    # #        assert [x for x in new_neighbors] == [y for y in neighbors]
    # #     except Exception as e:    
    # #         print("edges no equal {e}")
    # #         breakpoint()
    #         #here 

    # # breakpoint()
    # print("edges equal")


    # basic_collection_old.find_flow_components()
    # # breakpoint()
    # basic_collection_old.calculate_flows()
    # breakpoint()
    # breakpoint()
    # accepted_err = .01
   
    # #there


    # basic_collection.initialize_digraph_from_new()
    # basic_collection_old.initialize_digraph_from()

    # breakpoint()

    # basic_collection.find_flow_components_new()
    # basic_collection_old.old_find_flow_components()
    # breakpoint()
    
    # print(basic_collection.drip_summary())
    # print(basic_collection_old.drip_summary())


    # basic_collection.calculate_flows()
    # actual_flows = basic_collection.flows
    # _, actual_stem_map = lam_filter(
    #     basic_collection.cylinders, lambda: is_stem, return_all=True
    # )