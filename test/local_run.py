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


# from data.output.run_so_far import already_run
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
            collection = initialize_collection(file_name)
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
            drip_map(collection, angle = angle, file = file_name)
        log.info(f"successfully ran case {case_name}")
        pickle(collection,f'_stats_{case_name}')
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
    except Exception as e:  
        log.info(f'Failed to draw and save pickle for {file}, case {angle}  :{e}')

def drip_map(collection = None, file:str = '', pickle_point:str = '', angle = ''):
    try:
        if not collection:
            collection = load_from_pickle(file, pickle_point,angle)

        drip_map = collection.drip_map(file_ext = f'{file}_{angle}_XZ.png')
    except Exception as e:  
        log.info(f'Failed to draw and save pickle for {file}, case {angle}  :{e}')



run_cases = [('Secrest02-30_000000', -0.2666), ('Secrest08-24c_000000', -1.5), ('Secrest02-30_000000', -1.0), ('Secrest32-06_000000', -0.82), ('Secrest03-12_000000', -0.3), ('Secrest10-08_000000', 0.1666), ('Secrest32-06_000000', -0.78), ('Secrest07-32_000000', 0.1666), ('Secrest24-03_000000', -0.95), ('Secrest32-06_000000', -0.33333), ('Secrest27-05_000000', -0.32), ('Secrest02-30_000000', -0.3), ('Secrest29-25_000000', -1.5), ('Secrest32-06_000000', 0.26), ('Secrest32-06_000000', 0.28), ('Secrest07-32_000000', -1.0), ('Secrest23-23_000000', 0.0), ('Secrest02-30_000000', -0.16666), ('Secrest32-06_000000', -0.8), ('Secrest02-26_000000', 0.0), ('Secrest32-03_000000', -0.0666), ('Secrest32-03_000000', -0.16666), ('Secrest31-05_000000', -2.0), ('Secrest24-03_000000', -0.0666), ('Secrest32-06_000000', 0.3666), ('Secrest10-02_000000', 0.0), ('Secrest10-08_000000', -0.3), ('Secrest27-05_000000', -1.0), ('Secrest18-13_000000', -0.16666), ('Secrest32-14_000000', -0.3), ('Secrest32-14_000000', -1.5), ('Secrest02-30_000000', 0.1666), ('Secrest10-02_000000', -0.16666), ('Secrest11-27_000000', 0.56666), ('Secrest24-03_000000', 0.4), ('Secrest02-30_000000', 0.0), ('Secrest24-03_000000', 0.3), ('Secrest11-27_000000', -2.0), ('Secrest14-09_000000', -1.5), ('Secrest11-27_000000', 0.4), ('Secrest14-09_000000', -0.0666), ('Secrest32-14_000000', 0.1666), ('Secrest29-20_000000', 0.1666), ('Secrest16-14LI-ST_000000', -1.5), ('Secrest11-27_000000', 0.3), ('Secrest27-05_000000', 0.56666), ('Secrest24-03_000000', 0.36666), ('Secrest32-03_000000', 0.1666), ('Secrest02-26_000000', -1.5), ('Secrest29-20_000000', 0.0), ('Secrest32-06_000000', 0.04), ('Secrest32-06_000000', 0.12), ('Secrest31-05_000000', -0.16666), ('Secrest32-06_000000', -0.5666), ('Secrest02-26_000000', -0.16666), ('Secrest16-3TI-CO_000000', -1.5), ('Secrest24-07_000000', -1.5), ('Secrest32-06_000000', -0.34), ('Secrest03-12_000000', 0.0), ('Secrest32-06_000000', -0.66666), ('Secrest27-05_000000', -0.7), ('Secrest31-05_000000', -1.0), ('Secrest32-14_000000', 0.0), ('Secrest32-06_000000', 0.08), ('Secrest32-03_000000', -1.5), ('Secrest16-14LI-ST_000000', 0.0), ('Secrest14-09_000000', -1.0), ('Secrest32-06_000000', -0.12), ('Secrest29-20_000000', -1.5), ('Secrest26-03_000000', -1.5), ('Secrest03-12_000000', -0.16666), ('Secrest03-12_000000', 0.1666), ('Secrest31-05_000000', -1.5), ('Secrest08-24c_000000', -0.16666), ('Secrest24-03_000000', 0.56666), ('Secrest10-02_000000', 0.46666), ('Secrest31-05_000000', 0.1666), ('Secrest18-13_000000', -2.5), ('Secrest27-05_000000', 0.36666), ('Secrest23-23_000000', -0.3), ('Secrest32-06_000000', 0.5), ('Secrest16-3TI-CO_000000', -0.3), ('Secrest32-03_000000', -0.2666), ('Secrest24-03_000000', -2.5), ('Secrest24-03_000000', -0.3), ('Secrest07-32_000000', 0.0), ('Secrest32-06_000000', -1.0), ('Secrest32-06_000000', -0.48), ('Secrest32-06_000000', -4.0), ('Secrest32-06_000000', 0.4), ('Secrest27-05_000000', 0.2), ('Secrest32-06_000000', 0.14), ('Secrest02-26_000000', -0.3), ('Secrest27-05_000000', -0.38), ('Secrest27-05_000000', -0.5), ('Secrest27-05_000000', 0.0), ('Secrest32-06_000000', -0.42), ('Secrest16-14LI-ST_000000', -0.16666), ('Secrest10-02_000000', -0.2666), ('Secrest08-24c_000000', -0.3), ('Secrest23-23_000000', -0.16666), ('Secrest16-3TI-CO_000000', -1.0), ('Secrest32-06_000000', -0.28), ('Secrest32-06_000000', -1.5), ('Secrest32-06_000000', -0.08), ('Secrest32-06_000000', -0.44), ('Secrest10-02_000000', -2.0), ('Secrest32-06_000000', -0.8), ('Secrest18-13_000000', -0.95), ('Secrest29-20_000000', -0.0666), ('Secrest32-06_000000', -0.06), ('Secrest32-06_000000', -0.22), ('Secrest32-06_000000', -0.16), ('Secrest32-06_000000', -0.1), ('Secrest10-02_000000', -2.5), ('Secrest32-06_000000', -0.18), ('Secrest27-05_000000', -0.3), ('Secrest11-27_000000', -0.16666), ('Secrest02-26_000000', -1.0), ('Secrest10-08_000000', -0.0666), ('Secrest14-09_000000', 0.1666), ('Secrest27-05_000000', -0.33333), ('Secrest32-06_000000', 0.0), ('Secrest27-05_000000', -0.28), ('Secrest32-06_000000', 0.56666), ('Secrest18-13_000000', -0.2666), ('Secrest27-05_000000', 0.1), ('Secrest18-13_000000', -2.0), ('Secrest32-06_000000', -0.36666), ('Secrest27-05_000000', -0.44), ('Secrest16-14LI-ST_000000', -0.3), ('Secrest27-05_000000', -0.95), ('Secrest32-06_000000', 0.16), ('Secrest16-3TI-CO_000000', 0.0), ('Secrest24-03_000000', 0.1666), ('Secrest32-03_000000', -2.5), ('Secrest32-03_000000', -2.0), ('Secrest32-06_000000', -0.2), ('Secrest11-27_000000', 0.36666), ('Secrest32-06_000000', -0.14), ('Secrest27-05_000000', -0.9), ('Secrest32-06_000000', 0.3), ('Secrest11-27_000000', -0.95), ('Secrest32-06_000000', -0.76666), ('Secrest11-27_000000', -0.0666), ('Secrest32-06_000000', 0.24), ('Secrest18-13_000000', 0.36666), ('Secrest10-02_000000', -1.5), ('Secrest32-06_000000', 0.26666), ('Secrest14-09_000000', 0.0), ('Secrest31-05_000000', -0.3), ('Secrest27-05_000000', -0.6), ('Secrest27-05_000000', -0.46), ('Secrest29-25_000000', -0.3), ('Secrest16-14LI-ST_000000', -1.0), ('Secrest07-32_000000', -0.16666), ('Secrest28-31_000000', -1.5), ('Secrest23-23_000000', 0.1666), ('Secrest24-03_000000', 0.46666), ('Secrest32-06_000000', 0.1666), ('Secrest32-06_000000', -0.24), ('Secrest32-01_000000', 0.1666), ('Secrest24-07_000000', -1.0), ('Secrest32-06_000000', -0.72), ('Secrest24-03_000000', -0.26666), ('Secrest02-30_000000', -1.5), ('Secrest32-06_000000', -2.0), ('Secrest08-24c_000000', 0.1666), ('Secrest10-08_000000', -1.0), ('Secrest32-06_000000', -0.46), ('Secrest24-03_000000', -0.9), ('Secrest16-3TI-CO_000000', 0.1666), ('Secrest10-02_000000', 0.1666), ('Secrest10-08_000000', -1.5), ('Secrest14-09_000000', -0.16666), ('Secrest27-05_000000', -0.42), ('Secrest16-14LI-ST_000000', 0.1666), ('Secrest32-06_000000', -0.9), ('Secrest08-24c_000000', -1.0), ('Secrest18-13_000000', -1.0), ('Secrest27-05_000000', -0.36), ('Secrest32-06_000000', 0.36666), ('Secrest24-03_000000', -1.0), ('Secrest32-06_000000', -0.26), ('Secrest32-06_000000', -0.95), ('Secrest32-06_000000', -2.5), ('Secrest32-06_000000', -0.3), ('Secrest32-06_000000', -0.68), ('Secrest32-06_000000', 0.02), ('Secrest07-32_000000', -1.5), ('Secrest10-08_000000', -0.16666), ('Secrest27-05_000000', -0.48), ('Secrest27-05_000000', 0.3), ('Secrest27-05_000000', -0.4), ('Secrest27-05_000000', -0.8), ('Secrest32-06_000000', -0.5), ('Secrest14-09_000000', -0.3), ('Secrest32-06_000000', -0.04), ('Secrest27-05_000000', -0.5666), ('Secrest11-27_000000', -1.0), ('Secrest11-27_000000', 0.5), ('Secrest18-13_000000', 0.4), ('Secrest10-02_000000', -0.95), ('Secrest32-06_000000', -0.02), ('Secrest11-27_000000', 0.1666), ('Secrest23-23_000000', -1.5), ('Secrest32-06_000000', -0.6), ('Secrest24-07_000000', 0.1666), ('Secrest03-12_000000', -1.5), ('Secrest27-05_000000', 0.46666), ('Secrest27-05_000000', 0.1666), ('Secrest32-06_000000', -0.32), ('Secrest32-06_000000', 0.36), ('Secrest29-25_000000', -1.0), ('Secrest32-06_000000', -0.4), ('Secrest08-24c_000000', 0.0), ('Secrest10-02_000000', -1.0), ('Secrest18-13_000000', -0.9), ('Secrest18-13_000000', -1.5), ('Secrest27-05_000000', -0.16666), ('Secrest32-06_000000', -0.84), ('Secrest32-06_000000', -0.76), ('Secrest27-05_000000', -0.26), ('Secrest32-06_000000', 0.06), ('Secrest32-06_000000', 0.38), ('Secrest32-03_000000', -1.0), ('Secrest24-03_000000', -1.5), ('Secrest32-06_000000', 0.1), ('Secrest32-06_000000', -0.74), ('Secrest27-05_000000', -0.76666), ('Secrest29-20_000000', -0.16666), ('Secrest32-06_000000', -1.5), ('Secrest26-03_000000', -1.0), ('Secrest32-14_000000', -1.0), ('Secrest18-13_000000', 0.46666), ('Secrest10-02_000000', 0.5), ('Secrest32-14_000000', -0.16666), ('Secrest07-32_000000', -0.3), ('Secrest11-27_000000', -0.9), ('Secrest23-23_000000', -0.0666), ('Secrest32-03_000000', 0.0), ('Secrest11-27_000000', 0.46666), ('Secrest10-02_000000', 0.36666), ('Secrest24-03_000000', 0.0), ('Secrest11-27_000000', -2.5), ('Secrest29-25_000000', 0.0), ('Secrest29-20_000000', -0.3), ('Secrest27-05_000000', -1.5), ('Secrest24-03_000000', -2.0), ('Secrest24-03_000000', -0.16666), ('Secrest24-03_000000', -0.2666), ('Secrest32-06_000000', 0.18), ('Secrest29-25_000000', 0.1666), ('Secrest10-02_000000', 0.3), ('Secrest16-3TI-CO_000000', -0.16666), ('Secrest32-03_000000', -0.3), ('Secrest27-05_000000', -2.5), ('Secrest29-20_000000', -1.0), ('Secrest27-05_000000', -0.34), ('Secrest27-05_000000', -2.0), ('Secrest18-13_000000', 0.1666), ('Secrest32-06_000000', 0.22), ('Secrest32-06_000000', 0.32), ('Secrest29-25_000000', -0.16666), ('Secrest27-05_000000', -0.66666), ('Secrest31-05_000000', 0.0), ('Secrest18-13_000000', -0.0666), ('Secrest27-05_000000', -0.24), ('Secrest32-06_000000', 0.2), ('Secrest24-03_000000', 0.5), ('Secrest02-30_000000', -0.0666), ('Secrest32-06_000000', -0.86), ('Secrest10-08_000000', 0.0), ('Secrest02-26_000000', 0.1666), ('Secrest28-31_000000', -0.3), ('Secrest18-13_000000', 0.0), ('Secrest10-02_000000', -0.0666), ('Secrest32-06_000000', -0.16666), ('Secrest10-02_000000', -0.3), ('Secrest32-06_000000', -0.36), ('Secrest32-06_000000', -0.26666), ('Secrest27-05_000000', -0.86666), ('Secrest10-02_000000', 0.4), ('Secrest23-23_000000', -1.0), ('Secrest32-06_000000', -0.7), ('Secrest11-27_000000', -1.5), ('Secrest29-25_000000', -2.0), ('Secrest18-13_000000', -0.3), ('Secrest18-13_000000', 0.3), ('Secrest27-05_000000', 0.4), ('Secrest23-23_000000', -0.26666), ('Secrest27-05_000000', 0.26666), ('Secrest11-27_000000', -0.2666), ('Secrest27-05_000000', -1.5), ('Secrest10-02_000000', -0.9), ('Secrest32-06_000000', -0.38), ('Secrest10-02_000000', 0.56666), ('Secrest11-27_000000', 0.0), ('Secrest11-27_000000', -0.3), ('Secrest32-06_000000', 0.34), ('Secrest03-12_000000', -1.0), ('Secrest29-20_000000', -2.0), ('Secrest32-06_000000', -0.88), ('Secrest27-05_000000', 0.5), ('Secrest32-06_000000', 0.46666)]


# angles = set([tup[1] for tup in already_run])
# angles = [-.8,-.88,-.86,-.84,-.82,-.80,-.78,-.76,-.74,-.72,-.7,-.68,-.66,-.64,-.62,-.6,-.58,-.56,-.54,-.52,-.5,-.48,-.46,-.44
#             ,-.42,-.44,-.46,-.48,-.4,-.38,-.36,-.34,-.32,-.30,-.28,-.26,-.24,-.22,-.2,-.18,-.16,-.14,-.12,-.1,-.08,-.06,-.04,-.02,
#   0.42,0.44,0.46,0.48, 0.8,0.88,0.86,0.84,0.82,0.80,0.78,0.76,0.74,0.72,0.7,0.68,0.66,
#            .4,.38,.36,.34,.32,.30,.28,.26,.24,.22,.2,.18,.16,.14,.12,.1,.08,.06,.04,.02,0.64,0.,
#            0.6,0.58,0.56,0.54,0.52,0.5,0.48,0.46,0.44]

# angles = [-.64,-.62,-.6,-.58,-.56,-.54,-.52,-.5,-.48,-.46,-.44
#           ,-.42,-.44,-.46,-.48,-.4,-.38,-.36,-.34,-.32,-.30,-.28,-.26
#           ,-.24,-.22,-.2,-.18,-.16,-.14,-.12,-.1,-.08,-.06,-.04,-.02
#           ,0.66,0.64,0.62, 0.6,0.58,0.56,0.54,0.52,0.5,0.48,0.46,0.44,0.42
#            ,.4,.38,.36,.34,.32,.30,.28,.26,.24,.22,.2,.18,.16,.14,.12,.1,.08,.06,.04,.02 ]

# angles = [0.8,0.88,0.86,0.84,0.82,0.80,0.78,0.76,0.74,0.72,0.7,0.68,0.66,0.64,0.,0.6,0.58,0.56,0.54,0.52,0.5]
# angles = [.08,-1.02,1.02,-1.18,-1.26,-1.32,1.18,1.26,1.32]
# angles = [0.94, 0.62, 0.30, 0.22, -0.1, -0.24, -0.64, -0.72, -0.88]
angles = [0.96, 0.88, 0.8, 0.64, 0.72, 0.48, 0.4, 0.32, 0.24, 0,
  1.1, -1.32, 1.5, -0.24, -0.32, -0.4, -1.4, -0.08, -0.56,
    -0.48, -0.64, -1.3, -0.72, -0.8, 1.2, -0.88, 0.16, -0.16, 
    -0.96, 1.3, -1.1, 0.08, 1.4, -1.18, -1.26, -1.5, -1.02]

def get_cases(file_names, already_run, angles_to_tests):
    already_run = [(x,float(y)) for x,y in already_run]
    cases = product(file_names,angles_to_tests)
    return [case for case in cases if case not in already_run]

def sensitivity_analysis():
    # files_to_test = ["5_SmallTree"]
    # files_to_test = [ "Secrest27-05_000000",
    #                  "Secrest32-06_000000"]
    # , "Secrest27-05_000000","Secrest03-12_000000"
                    #     ,"Secrest07-32_000000"]
                        #  "Secrest02-26_000000"1,-.08,-.06,-.1,-.08,-.06,-.04,-.02,
#           .4,.38,.36,.34,.32,.30,.28,.26,.24,.22,.2,.18,.16,.14,.12,.1,.08,.06,.04,.0204,-.02,
#           .4,.38,.36,.34,.32,.30,.28,.26,.24,.22,.2,.18,.16,.14,.12,.1,.08,.06,.04,.021,-.08,-.06,-.04,-.02,
#           .4,.38,.36,.34,.32,.30,.28,.26,.24,.22,.2,.18,.16,.14,.12,.1,.08,.06,.04,.02
    files_to_test =    ["Secrest02-30_000000"
                        ,"Secrest03-12_000000"
                        ,"Secrest07-32_000000"
                        ,"Secrest08-24c_000000"
                        ,"Secrest10-02_000000"]
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
    # log.info(f'Will run {len(cases_to_run)} cases : {cases_to_run}')
    start = time()
    # # breakpoint()
    success = run_test_cases(cases_to_run, fig = True)
    # for file, angle in cases_to_run:
    #     success = run_test_cases(cases_to_run)
    #     if not success:
    #         log.info(f"Failed run cases")
    #     else:
    #         log.info(f"suceeded running cases")
    dur = time() - start

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

    # run_test_cases([('Secrest32-06_000000',-0.16666),
    #                 ('Secrest32-06_000000',-0.25),
    #                 ('Secrest27-05_000000',-0.16666),
    #                 ('Secrest27-05_000000',-0.25)],
    #                 stats = True, fig = True, from_pickle = False)
    # run_test_cases([('Secrest32-06_000000',1.2),
    #                 ('Secrest32-06_000000',1.4),
    #                 ('Secrest32-06_000000',1.5),
    #                 ('Secrest27-05_000000',1.2),
    #                 ('Secrest27-05_000000',1.4),
    #                 ('Secrest27-05_000000',1.5)],
    #                 stats = True, from_pickle = False)


    # forest = Forester()
    # forest.get_file_names(dir=test_input_dir)
    # forest.qsm_from_file_names(file_name="5_SmallTree")
    # basic_collection = forest.cylinder_collections[0]


    # forest_old = Forester()
    # forest_old.get_file_names(dir=test_input_dir)
    # forest_old.qsm_from_file_names(file_name="5_SmallTree.csv")
    # basic_collection_old = forest_old.cylinder_collections[0]

    # basic_collection.project_cylinders("XY")

    # angles = [ (cyl.cyl_id,cyl.angle) for cyl in basic_collection.cylinders]
    # new_angles = [ (cyl.cyl_id,cyl.angle) for cyl in basic_collection.cylinders]
    # for idx, angle in angles:
    #     new_angle = [angle for cyl, angle in new_angles if cyl == idx][0]
    #     try:
    #         assert within_range(angle, new_angle, .03)
    #     except Exception as e:    
    #         print("Failure in projection {e}")
    #         breakpoint()
            #here 
    # breakpoint()

    # basic_collection.initialize_digraph_from(in_flow_grade_lim=-0.3)
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


    # basic_collection.find_flow_components()
    # # breakpoint()
    # basic_collection.calculate_flows()
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