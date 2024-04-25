from __future__ import annotations

import os
import sys

from time import time 
from itertools import product

print(sys.path)
sys.path.insert(0, os.path.dirname(os.getcwd()))
print(sys.path)
sys.path.insert(0, os.getcwd())
print(sys.path)

from shapely.ops import unary_union
from data.output.run_so_far import already_run
from src.canhydro.global_vars import DIR, test_input_dir
from src.canhydro.CylinderCollection import (
    pickle_collection, unpickle_collection
)
from src.canhydro.Forester import Forester

from src.canhydro.global_vars import log
from scripts.sftp_utils import sftp


def try_pickle_collection(collection, designation = ""):
    log.info(f"attempting to pickle for file {collection.file_name}, des {designation}")
    # pickle_file = pickle_collection(collection, designation)
    try: 
        pickle_file = pickle_collection(collection, designation)
        log.info(f"successfully created pickle for file {collection.file_name}")
    except Exception as e:
        log.info(f"Error pickling file {collection.file_name}: {e}")
        return

def initialize_collection(file = "5_SmallTree", from_pickle = False, **kwargs):
    if from_pickle:
        collection = load_from_pickle(**kwargs)
    else:
        log.info(f"initializing collection...")
        try:
            forest = Forester()
            forest.get_file_names(dir=test_input_dir)
            forest.qsm_from_file_names(file_name=file)
            basic_collection = forest.cylinder_collections[0]
            basic_collection.project_cylinders("XY")
        except Exception as e:
            log.info(f"Error initializing collection for file {file}: {e}")
            return None
        log.info(f"successfully initialized collection")
        # try_pickle_collection(basic_collection,file)
    return basic_collection

def prep_for_stats(collection, case_angle, case_name, calculate:bool = True):
    log.info(f"attempting to prep for stats {case_name}")
    try: 
        collection.initialize_digraph_from(in_flow_grade_lim=case_angle)
        collection.find_flow_components()
        if calculate: collection.calculate_flows()
    except Exception as e:
        log.info(f"Error preping for  stats for case {case_name}: {e}")
        return None
    log.info(f"successfully prepped for stats {case_name}")
    return True
    
def generate_statistics(collection, case_name):
    log.info(f"attempting to generate stats for file {collection.file_name}, case_name {case_name}")
    # statistics = collection.statistics(file_ext = case_name)
    try: 
        stat_file = collection.statistics(file_ext = case_name)
        send_file_to_pputer(stat_file)
    except Exception as e:
        log.info(f"Error gernerating stats for case {case_name} : {e}")
        return None
    log.info(f"attempting to generate flow file for {case_name}")
    try: 
        flow_file = collection.generate_flow_file(file_ext = case_name)
        send_file_to_pputer(flow_file)
    except Exception as e:
        log.info(f"Error gernerating flow file for case {case_name}: {e}")
        return None
    log.info("successfully created flow and stats files")
    return True

def run_test_cases(cases_to_run, stats :bool = True, fig:bool = False, from_pickle:bool =False):

    start = time()
    collection = None
    case_name = f"inital_case_name"
    ret = []
    for case in cases_to_run:
        file_name, angle = case
        case_name = f"{angle}"
        if not collection:
            collection = initialize_collection(file_name)
            if not collection:  
                dur = time() - start
                ret.append((None, f'{file_name}_{case_name}', dur))
                continue
        log.info(f"running case {file_name}_{case_name}")

        preped = prep_for_stats(collection, angle, case_name, calculate=stats)
        if stats:
            if preped:
                stat_file = generate_statistics(collection, f'{case_name}')
            else:
                log.info(f'Error prepping, pickling and ending ')
                # try_pickle_collection(collection,f'_prep_{case_name}')
                dur = time() - start
                ret.append((None, f'{file_name}_{case_name}', dur))
                continue
        if fig:
            draw_case(collection, angle = angle, file = file_name)
            drip_map(collection, angle = angle, file = file_name)
        log.info(f"successfully ran case {case_name}")
        try_pickle_collection(collection,f'_stats_{case_name}')
        collection=None
        dur = time() - start
        ret.append((True, f'{file_name}_{case_name}', dur))
    return ret

def load_from_pickle(file, pickle_point, angle):
    pickle_file = f'{file}_pickle__{pickle_point}_{angle}'
    collection = None
    try:
        log.info(f'Loading collection for {file}, case {angle} from pickle file {pickle_file}')
        collection = unpickle_collection(pickle_file)  
        if collection: 
            log.info(f'Collection for {file}, case {angle} successfully loaded from pickle file {pickle_file}')
    except Exception as e:
        log.info(f'failed to load pickle for {file}, case {angle}  :{e}')
    return collection


def draw_case(collection = None, file:str = '', pickle_point:str = '', angle = ''):
    try:
        if not collection:
            collection = load_from_pickle(file, pickle_point,angle)

        stem_flow_fig = collection.draw(plane = 'XZ', highlight_lambda = lambda:is_stem,
                                        save = True, file_ext = f'{file}_{angle}_XZ.png', show = False)
        stemflow_and_trunk_fig = collection.draw(plane = 'XZ', filter_lambda = lambda: is_stem, highlight_lambda = lambda:branch_order==0,
                                        save = True, file_ext = f'{file}_{angle}_XZ.png', show= False)
        log.info(f'Draw commented out {file}, case {angle} ')
    except Exception as e:  
        log.info(f'Failed to draw and save pickle for {file}, case {angle}  :{e}')

def drip_map(collection = None, file:str = '', pickle_point:str = '', angle = ''):
    try:
        if not collection:
            collection = load_from_pickle(file, pickle_point,angle)

        drip_map = collection.drip_map(file_ext = f'{file}_{angle}_XZ.png')
    except Exception as e:  
        log.info(f'Failed to draw and save pickle for {file}, case {angle}  :{e}')

def alpha_shape(collection = None, file:str = '', pickle_point:str = '', angle = ''):
    try:
        if not collection:
            collection = load_from_pickle(file, pickle_point,angle)

        drip_map = collection.drip_map(file_ext = f'{file}_{angle}_XZ.png')
    except Exception as e:  
        log.info(f'Failed to draw and save pickle for {file}, case {angle}  :{e}')



run_cases = already_run

angles = [ 0.96, 0.88,
         0.8, 0.72, 0.64, 0.56, 0.48, 0.4, 0.32, 0.24, 0.16,
           0.08, 0, -0.08, -0.16, -0.24, -0.32, -0.4, -0.48,
           -0.56, -0.64, -0.72, -0.8, -0.88, -0.96, -1.02,

            2.04,1.96,1.86, 1.8,1.72,1.66, 1.64, 1.58, 1.5, 1.42, 1.34,
           1.26, 1.18, 1.1, 1.02,
           -1.1, -1.18, -1.26, -1.34, -1.42, -1.5, -1.58
           ,-1.66,-1.64,-1.72,-1.8,-1.88,-1.96,-2.04]

def get_cases(file_names, already_run, angles_to_tests):
    already_run = [(x,float(y)) for x,y in already_run]
    cases = product(file_names,angles_to_tests)
    return [case for case in cases if case not in already_run]


def send_file_to_pputer(file,parent_folder= ''):
    log.info(f' file w/ name {file} sent for sftp')
    qualified_name = f'{parent_folder}/{file}' if parent_folder else file
    try:
        sftp(qualified_name, get = False)
    except Exception as e:
        log.info(f'difficulty sharing file {qualified_name}')


def sensitivity_analysis():
    files_to_test = ["5_SmallTree"] 
    # files_to_test  = ["Secrest02-26_000000"
    #                 ,"Secrest18-13_000000"
    #                 ,"Secrest10-08_000000"

    #                 ,"Secrest16-14LI-ST_000000"
    #                 ,"Secrest11-27_000000"
    #                 ,"Secrest14-09_000000"
    #                 ,"Secrest23-23_000000"
    #                 ,"Secrest24-03_000000"
    #                 ,"Secrest11-27_000000"]
                        

            #run or running ontowr 
                    # ,"Secrest27-05_000000"
                        # ,"Secrest26-03_000000"
                        # ,"Secrest28-31_000000"
                        # ,"Secrest29-20_000000"
                        # ,"Secrest29-25_000000"
                        # ,"Secrest31-05_000000"
                        # ,"Secrest32-01_000000"
                        # ,"Secrest32-03_000000"
                        # ,"Secrest32-06_000000"
                        # ,"Secrest32-14_000000"]
                        # ,"Secrest24-07_000000"
                    #       ["Secrest02-30_000000"
                    #     ,"Secrest03-12_000000"
                    #     ,"Secrest07-32_000000"
                    #     ,"Secrest08-24c_000000"
                    #     ,"Secrest10-02_000000"]
    #                     ,"Secrest16-3TI-CO_000000"


    cases_to_run = get_cases(files_to_test,run_cases,angles)
    log.info(f'Will run {len(cases_to_run)} cases : {cases_to_run}')
    success = run_test_cases(cases_to_run, fig = False)
    return success
    # for file, angle in cases_to_run:
    #     success = run_test_cases(cases_to_run)
    #     if not success:
    #         log.info(f"Failed run cases")
    #     else:
    #         log.info(f"suceeded running cases")

    # log.info(f"total time old method - {dur}")


    # with mp.Pool(1) as p:
    #     task_pool = [p.apply_async(run_test_cases, args=(cases_to_run,)) for file in files_to_test]
    #     results = [task.get() for task in task_pool]

    # for success, case, dur in results:
    #     log.info(f"total time running {case} - {dur}")
    #     if not success:
    #         log.info(f"Failed run case {case}")
    #     else:
    #         log.info(f"suceeded running cases {case}")

if __name__ == "__main__":
    # data/output/pickle/5_SmallTree_pickle__prep_-0.1

    # load_from_pickle('Secrest32-06_000000', 'stats', 0.3666)
    # load_from_pickle('Secrest32-06_000000', 'stats', 0.3666)
    # load_from_pickle('Secrest32-06_000000_1', 'stats', -1.5)
    # load_from_pickle('5_SmallTree_1', 'stats', 0.36666)
    sensitivity_analysis()       
    # draw_for_paper()
    # run_test_cases([('5_SmallTree',-0.16666),
    #                 ('5_SmallTree',-0.16666),
    #                 ('5_SmallTree',-0.16666)],
    #                 stats = True, fig = False, from_pickle = False)
    # run_test_cases([('Secrest32-06_000000',1.2),
    #                 ('Secrest32-06_000000',1.4),
    #                 ('Secrest32-06_000000',1.5),c
    #                 ('Secrest27-05_000000',1.2),
    #                 ('Secrest27-05_000000',1.4),
    #                 ('Secrest27-05_000000',1.5)],
    #                 stats = True, from_pickle = False)