from timeit import Timer 
import functools
from copy import deepcopy
import sys
import os
import yaml
import logging
from memory_profiler import profile
sys.path.insert(0, os.path.dirname(os.getcwd()))
sys.path.insert(0, os.getcwd())

mem_logs =open('memory_profiler.log','w')
mem_logs =open('memory_profiler_rust.log','w')

from scripts.basic_recipies import initialize_collection,try_pickle_collection
from src.canhydro.CylinderCollection import unpickle_collection
from src.canhydro.utils import lam_filter
from src.canhydro.DataClasses import Flow

with open('logging_config.yml', 'rt') as f:
    config = yaml.safe_load(f.read())
    logging.config.dictConfig(config)

script_logger = logging.getLogger('script')

log = script_logger

@profile(stream=mem_logs)
def find_flows(collection):
    collection.find_flow_components_pre_ancestors()
    collection.calculate_flows()
    collection.digraph=None
    collection.flows=None

def find_flows_new(collection):
    collection.initialize_digraph_from()
    collection.find_flow_components()
    collection.calculate_flows()

@profile(stream=mem_logs)
def find_flows_rust(collection):
    collection.initialize_digraph_from_rust()
    collection.find_flow_components_rust()
    collection.digraph=None
    collection.flows=None
    # collection.calculate_flows_rust()

def compare_flows(c1,c2):
    breakpoint()
    f1 = c1.flows
    f2 = c2.flows
    diffs= []
    corrections= []
    corrected= []
    unique_drip_nodes= [[],[]]
    for idf, flow in enumerate(f1):
        drip_node = flow.drip_node_id
        compare_to = [flow2 for flow2 in f2 if flow2.drip_node_id ==drip_node ]
        if len(compare_to) == 0:
            unique_drip_nodes[0].append(drip_node)
            compare_to = f2[idf]
        else: 
            compare_to = compare_to[0]
        if flow != compare_to:
            cyl = [drp_cyl for drp_cyl in c1.cylinders if drp_cyl.cyl_id == flow.drip_node_id][0]
            correction = Flow(    1,    np.float16(cyl.projected_data['XY']["area"]),    cyl.surface_area,cyl.angle,cyl.volume,cyl.sa_to_vol,cyl.cyl_id,[0,0,0])
            diff = flow.compare(compare_to)
            corrections.append(correction)
            if flow - diff == compare_to:
                corrected.append(diff)
            else:
                diffs.append(diff)

    breakpoint()


def set_up_and_run(file):
    collection = initialize_collection(file,from_pickle = False)
    other_collection = initialize_collection(file,from_pickle = False)
    # collection.find_flow_components_pre_ancestors()
    # other_collection.find_flow_components_rust()
    # breakpoint()
    # try_pickle_collection(collection)

    # collection = unpickle_collection(f'data/output/pickle/{file}.pickle')
    # other_collection = deepcopy(collection)
    # other_collection.initialize_digraph_from_rust()
    log.info(f'Starting Testing for find_flows')#function {func.__name__}')


    t = Timer(functools.partial(find_flows_rust,other_collection))
    rust_func_sec_to_complete =t.timeit(1)
    log.info('Completed testing new function')
    
    t = Timer(functools.partial(find_flows,collection))
    original_func_time_to_complete =t.timeit(1)
    log.info('Completed testing base line function')





    log.info('Completed 30 trials for both functions')
    log.info(f'{rust_func_sec_to_complete=}')
    log.info(f'{original_func_time_to_complete=}')
    breakpoint()
    # log.info(f'Itr {i} for {func.__name__} complete')

if __name__ == "__main__":
    # 5_SmallTree
    # set_up_and_run('10_MediumCollection')
    # set_up_and_run('5_SmallTree')
    # set_up_and_run('Secrest04-19_000000')
    set_up_and_run('Secrest27-05_000000')
