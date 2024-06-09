from timeit import Timer 
import functools
from copy import deepcopy
import sys
import os
import yaml
import logging

sys.path.insert(0, os.path.dirname(os.getcwd()))
sys.path.insert(0, os.getcwd())


from scripts.basic_recipies import initialize_collection,try_pickle_collection
from src.canhydro.CylinderCollection import unpickle_collection
from src.canhydro.utils import lam_filter

with open('logging_config.yml', 'rt') as f:
    config = yaml.safe_load(f.read())
    logging.config.dictConfig(config)

script_logger = logging.getLogger('script')

log = script_logger


def find_flows(collection):
    collection.find_flow_components_pre_ancestors()

def find_flows_new(collection):
    collection.find_flow_components()

def find_flows_rust(collection):
    collection.find_flow_components_rust()

def compare_flows(c1,c2):
    breakpoint()
    f1 = c1.flows
    f2 = c2.flows
    # _, map1 = lam_filter(
    #     c1.cylinders, lambda: is_stem, return_all=True
    # )
    # _, map2 = lam_filter(
    #     c2.cylinders, lambda: is_stem, return_all=True
    # )
    diffs= []
    unique_drip_nodes= [[],[]]
    
    for idf, flow in enumerate(f1):
        drip_node = flow.drip_node_id
        compare_to = [flow2 for flow2 in f2 if flow2.drip_node_id ==drip_node ]
        if len(compare_to) == 0:
            diffs[idf] = flow.compare(compare_to)
        else:
            unique_drip_nodes[0].append(drip_node)
    f1_drip_nodes = [flow.drip_node_id for flow in f1]
    unique_drip_nodes[1].extend([flow.drip_node_id for flow in f2 if flow.drip_node_id not in f1_drip_nodes])

    breakpoint()


def set_up_and_run(file):
    collection = initialize_collection(file,from_pickle = False)
    other_collection = initialize_collection(file,from_pickle = False)
    collection.initialize_digraph_from()
    other_collection.initialize_digraph_from_rust()
    # breakpoint()
    # try_pickle_collection(collection)

    # collection = unpickle_collection(f'data/output/pickle/{file}.pickle')
    # other_collection = deepcopy(collection)
    # other_collection.initialize_digraph_from_rust()
    log.info(f'Starting Testing for find_flows')#function {func.__name__}')



    t = Timer(functools.partial(find_flows_rust,other_collection))
    rust_func_sec_to_complete =t.timeit(10)
    log.info('Completed testing new function')

    t = Timer(functools.partial(find_flows,collection))
    original_func_time_to_complete =t.timeit(10)
    log.info('Completed testing base line function')



    log.info('Completed 30 trials for both functions')
    log.info(f'{rust_func_sec_to_complete=}')
    log.info(f'{original_func_time_to_complete=}')
    breakpoint()
    # log.info(f'Itr {i} for {func.__name__} complete')


# if __name__ == "__main__":
    # 5_SmallTree
    # set_up_and_run('10_MediumCollection')
    # set_up_and_run('5_SmallTree')
    # set_up_and_run('Secrest04-19_000000')
    # set_up_and_run('Secrest27-05_000000')