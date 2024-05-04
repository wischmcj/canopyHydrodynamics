from __future__ import annotations

import os
import sys
import multiprocessing as mp
from time import time 
from itertools import product
import pytest


sys.path.insert(0, os.path.dirname(os.getcwd()))
sys.path.insert(0, os.getcwd())
# print(sys.path)


# from data.output.run_so_far import already_ru
from src.canhydro.global_vars import DIR, test_input_dir
from src.canhydro.utils import lam_filter
from src.canhydro.CylinderCollection import CylinderCollection, pickle_collection, unpickle_collection
from src.canhydro.Forester import Forester
from test.utils import within_range
from data.output.run_so_far import already_run
from src.canhydro.global_vars import config_vars, log, output_dir


def pickle(collection, designation = ""):
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
        # pickle(basic_collection,file)
    return basic_collection

def prep_for_stats(collection, case_angle, case_name, calculate:bool = True):
    log.info(f"attempting to prep for stats for case {case_name}")
    collection.initialize_digraph_from(in_flow_grade_lim=case_angle)
    collection.find_flow_components()
    if calculate: collection.calculate_flows()
    # try: 
    #     collection.initialize_digraph_from(in_flow_grade_lim=case_angle)
    #     collection.find_flow_components()
    #     if calculate: collection.calculate_flows()
    # except Exception as e:
    #     log.info(f"Error preping for  stats for case {case_name}: {e}")
    #     return None
    log.info(f"successfully prepped for stats {case_name}")
    return True
    
def generate_statistics(collection, case_name):
    log.info(f"attempting to generate stats for file {collection.file_name}, case_name {case_name}")
    # statistics = collection.statistics(file_ext = case_name)
    try: 
        statistics = collection.statistics(file_ext = case_name)
    except Exception as e:
        log.info(f"Error gernerating stats for case {case_name} : {e}")
        return None
    log.info(f"attempting to generate flow file for {case_name}")
    try: 
        collection.generate_flow_file(file_ext = case_name)
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
            # pickle_args = { "from_pickle": from_pickle,
            #                 'file':file_name,
            #                 'pickle_point':'stats', 
            #                 'angle': angle}
            collection = initialize_collection(file_name)
            # collection = load_from_pickle('Secrest27-05_000000', 'stats', -0.86666)
            if not collection:  
                dur = time() - start
                ret.append((None, f'{file_name}_{case_name}', dur))
                continue
        log.info(f"running case {file_name}_{case_name}")

        preped = prep_for_stats(collection, angle, case_name, calculate=stats)
        if stats:
            if preped:
                generate_statistics(collection, case_name)
            else:
                log.info(f'Error prepping, pickling and ending ')
                pickle(collection,f'_prep_{case_name}')
                dur = time() - start
                ret.append((None, f'{file_name}_{case_name}', dur))
                continue
        if fig:
            draw_case(collection, angle = angle, file = file_name)
        log.info(f"successfully ran case {case_name}")
        # pickle(collection,f'_stats_{case_name}')
        collection=None
        dur = time() - start
        ret.append((True, f'{file_name}_{case_name}', dur))
    return ret

def run_test_case( case, stats :bool = True):
    start = time()
    file_name, angle = case
    case_name = f"{angle}"
    # collection = load_from_pickle('Secrest27-05_000000",1', 'prep', -0.5)
    collection = initialize_collection(file_name)
    log.info(f"running case {file_name}_{case_name}")

    preped = prep_for_stats(collection, angle, case_name, calculate=stats)
    if stats:
        if preped:
            generate_statistics(collection, case_name)
        else:
            log.info(f'Error prepping, pickling and ending ')
            # pickle(collection,f'_prep_{case_name}')
            dur = time() - start
            return (None, f'{file_name}_{case_name}', dur)
    log.info(f"successfully ran case {case_name}")
    # pickle(collection,f'_stats_{case_name}')
    collection=None
    dur = time() - start
    log.info(f"running case {file_name}_{case_name} took {dur} seconds")
    return True, f'{file_name}_{case_name}', dur

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
                                        save = True, file_ext = f'{file}_{angle}.png', show = False)
        trunk_fig = collection.draw(plane = 'XZ', filter_lambda = lambda: is_stem, highlight_lambda = lambda:branch_order==0,
                                        save = True, file_ext = f'{file}_{angle}', show= False)
    except Exception as e:  
        log.info(f'Failed to draw and save pickle for {file}, case {angle}  :{e}')



run_cases = already_run

# angles = set([tup[1] for tup in already_run])
# angles = [-.8,-.88,-.86,-.84,-.82,-.80,-.78,-.76,-.74,-.72,-.7,-.68,-.66,-.64,-.62,-.6,-.58,-.56,-.54,-.52,-.5,-.48,-.46,-.44
#             ,-.42,-.44,-.46,-.48,-.4,-.38,-.36,-.34,-.32,-.30,-.28,-.26,-.24,-.22,-.2,-.18,-.16,-.14,-.12,-.1,-.08,-.06,-.04,-.02,
# angles = [  0.42,0.44,0.46,0.48,  0.8,0.88,0.86,0.84,0.82,0.80,0.78,0.76,0.74,0.72,0.7,0.68,0.66,
#            .4,.38,.36,.34,.32,.30,.28,.26,.24,.22,.2,.18,.16,.14,.12,.1,.08,.06,.04,.02,0.64,0.,
#            0.6,0.58,0.56,0.54,0.52,0.5,0.48,0.46,0.44]
# angles = [  0.42, 0.44, 0.46, 0.48, .4,.38,.36,.34,.32,.30,.28,.26,.24,.22,.2,.18,.16,.14,.12,.1,.08, 0.58,0.56,0.54,0.52,0.5,0.48,0.46,0.44]
# angles = [0.96, 0.64, 0.32, 0.24, -0.08,-0.24,-0.64,-0.72,-0.88]
# angles = [ -0.1, -0.24, -0.02,0.3]
# angles = [-.58, -.56, -.54, -.52, -0.24]

# angles = [-1.52,-1.48,-1.44,-1.4,-1.36,-1.32,-1.28,-1.24,-1.2,-1.16,-1.12,-1.08,-1.04,-1,
# -0.96,-0.92,-0.88,-0.84,-0.8,-0.76,-0.72,-0.68,-0.64,-0.6,-0.56,-0.52,-0.48,-0.44,-0.4,
# -0.36,-0.32,-0.28,-0.24,-0.2,-0.16,-0.12,-0.08,-0.04,0,0.04,0.08,0.12,0.16,0.2
# ,0.24,0.28,0.32,0.36,0.4,0.44,0.48,0.52,0.56,0.6,0.64,0.68,0.72,0.76,0.8,0.84,0.88
# ,0.92,0.96,1,1.04,1.08,1.12,1.16,1.2,1.24,1.28,1.32,1.36,1.4,1.44,1.48,1.52]

angles = [-.56, -.66, -.74, -0.82,-0.96, 0.96,0.16, -1.02,-1.1,-1.18,-1.26,-1.34,-1.42,-1.5, 1.02,1.1,1.18,1.26,1.34,1.42,1.5]

def get_cases(file_names, already_run, angles_to_tests):
    already_run = [(x,float(y)) for x,y in already_run]
    cases = product(file_names,angles_to_tests)
    return [case for case in cases if case not in already_run]

def sensitivity_analysis():
    files_to_test = ["Secrest27-05_000000","Secrest32-06_000000"]
    # files_to_test =["Secrest03-12_000000"
    #                     ,"Secrest07-32_000000"
    #                     ,"Secrest08-24c_000000"]
    # files_to_test = ["Secrest32-06_000000", "Secrest27-05_000000","Secrest03-12_000000"
    #                     ,"Secrest07-32_000000"]
                        #  "Secrest02-26_000000"1,-.08,-.06,-.1,-.08,-.06,-.04,-.02,
#           .4,.38,.36,.34,.32,.30,.28,.26,.24,.22,.2,.18,.16,.14,.12,.1,.08,.06,.04,.0204,-.02,
#           .4,.38,.36,.34,.32,.30,.28,.26,.24,.22,.2,.18,.16,.14,.12,.1,.08,.06,.04,.021,-.08,-.06,-.04,-.02,
#           .4,.38,.36,.34,.32,.30,.28,.26,.24,.22,.2,.18,.16,.14,.12,.1,.08,.06,.04,.02
    #  files_to_test =    ["Secrest02-30_000000"
    #                     ,"Secrest03-12_000000"
                        # ,"Secrest07-32_000000"
                        # ,"Secrest08-24c_000000"
                        # ,"Secrest10-02_000000"
                        # ,"Secrest10-08_000000"
                        # ,"Secrest11-27_000000"
                        # ,"Secrest14-09_000000"
                        # ,"Secrest16-3TI-CO_000000"
                        # ,"Secrest16-14LI-ST_000000"
                        # ,"Secrest18-13_000000"
                        # ,"Secrest23-23_000000"
                        # ,"Secrest24-03_000000"
                        # ,"Secrest24-07_000000"
                        # ,"Secrest26-03_000000"
                        # ,"Secrest28-31_000000"
                        # ,"Secrest29-20_000000"
                        # ,"Secrest29-25_000000"
                        # ,"Secrest31-05_000000"
                        # ,"Secrest32-01_000000"
                        # ,"Secrest32-03_000000"
                        # ,"Secrest32-06_000000"
                        # ,"Secrest32-14_000000"]
    cases_to_run = get_cases(files_to_test,run_cases,angles)
    breakpoint()
    # log.info(f'Will run {len(cases_to_run)} cases : {cases_to_run}')
    # success = run_test_cases(cases_to_run, fig = True)
    # for file, angle in cases_to_run:
    #     success = run_test_cases(cases_to_run)
    #     if not success:
    #         log.info(f"Failed run cases")
    #     else:
    #         log.info(f"suceeded running cases")


    with mp.Pool(4) as p:
        task_pool = [p.apply_async(run_test_case, args=(case,)) for case in cases_to_run]
        results = [task.get() for task in task_pool]

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
    # load_from_pickle('Secrest32-06_000000",1', 'stats', -1.5)
    # load_from_pickle('5_SmallTree_1', 'stats', 0.36666)
    sensitivity_analysis()       

    # forest = Forester()
    # forest.get_file_names(dir=test_input_dir)
    # forest.qsm_from_file_names(file_name="5_SmallTree")
    # basic_collection = forest.cylinder_collections[0]or cyl in collection.cylinders: rint(cyl.angle>=-1.5) if cyl.cyl_id in collection.trunk_nodes print(cyl.cyl_id

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