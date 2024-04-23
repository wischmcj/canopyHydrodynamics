
from time import time 

from .Forester import Forester
from .Cylinder import create_cyl
from .global_vars import test_input_dir , log
import numpy as np

from timeit import timeit

def initialize_forester(dir, file = None):
    forester = Forester()
    forester.get_file_names(dir=test_input_dir)
    print(forester.file_names)
    forester.qsm_from_file_names(file_name=file)
    return forester

def classic_init(file):
    forester = Forester()
    forester.get_file_names(dir=test_input_dir)
    arr = np.genfromtxt(file, delimiter=",", skip_header=True)[0:, :-1]
    # cylinders = [create_cyl(row) for row in arr]
    # cylinders = cylinders
    # cyl_index = dict((cyl.cyl_id , cyl ) for cyl in cylinders)
    # if 1==1:
    #     min_x = np.min([cyl.x[0] for cyl in cylinders])
    #     min_y = np.min([cyl.y[0] for cyl in cylinders])
    #     min_z = np.min([cyl.z[0] for cyl in cylinders])
    #     max_x = np.max([cyl.x[1] for cyl in cylinders])
    #     max_y = np.max([cyl.y[1] for cyl in cylinders])
    #     max_z = np.max([cyl.z[1] for cyl in cylinders])
    #     # Aggregate values from file
    #     no_cylinders = len(cylinders)
    #     surface_area = np.sum([cyl.surface_area for cyl in cylinders])
    #     volume = np.sum([cyl.volume for cyl in cylinders])
    #     max_branch_order = np.max([cyl.branch_order for cyl in cylinders])
    #     max_reverse_branch_order = np.max(
    #         [cyl.reverse_branch_order for cyl in cylinders]
    #     )
    #     avg_sa_to_vol = (
    #         np.sum([cyl.sa_to_vol for cyl in cylinders]) / no_cylinders
    #     )


def map_init(file):
    forester = Forester()
    forester.get_file_names(dir=test_input_dir)
    arr = np.genfromtxt(file, delimiter=",", skip_header=True)[0:, :-1]
    # start = time()
    # cylinders = [create_cyl(row) for row in arr]
    # cylinders = cylinders
    # cyl_index = dict((cyl.cyl_id , cyl ) for cyl in cylinders)
    # for
    # if 1==1:
    #     min_x = np.min([cyl.x[0] for cyl in cylinders])
    #     min_y = np.min([cyl.y[0] for cyl in cylinders])
    #     min_z = np.min([cyl.z[0] for cyl in cylinders])
    #     max_x = np.max([cyl.x[1] for cyl in cylinders])
    #     max_y = np.max([cyl.y[1] for cyl in cylinders])
    #     max_z = np.max([cyl.z[1] for cyl in cylinders])
    #     # Aggregate values from file
    #     no_cylinders = len(cylinders)
    #     surface_area = np.sum([cyl.surface_area for cyl in cylinders])
    #     volume = np.sum([cyl.volume for cyl in cylinders])
    #     max_branch_order = np.max([cyl.branch_order for cyl in cylinders])
    #     max_reverse_branch_order = np.max(
    #         [cyl.reverse_branch_order for cyl in cylinders]
    #     )
    #     avg_sa_to_vol = (
    #         np.sum([cyl.sa_to_vol for cyl in cylinders]) / no_cylinders
    #     )
        
    
    # tree = forest.cylinder_collections[0]

def old_flows(input):
    # 4_LargeCollection
    # 5_SmallTree
    start = time()
    input.find_flow_components()
    input.old_calculate_flows()
    dur = time() - start
    return dur

def new_flows(input):
    start = time()
    input.find_flow_components()
    input.calculate_flows()
    dur = time() - start
    return dur


def test(file, func):
    forest = initialize_forester(test_input_dir,file)
    log.info(f'Initializing for function {func.__name__}')
    tree = forest.cylinder_collections[0]
    tree.initialize_digraph_from()
    # tree.find_flow_components()
    tree.project_cylinders()
    log.info(f'Starting Testing for function {func.__name__}')
    func(tree)
    log.info(f'Itr {i} for {func.__name__} complete')


def compare(file = '5_Small_Collection', old = old_flows, new = new_flows):
    iterations = 1
    new_time = timeit('test(file, new)',number =1)
    old_time = timeit('test(file, old)',number =1)

    new_line ='\n'
    return  f"""Old average time is {old_time / iterations:.2f} seconds 

                {new_line.join(['',''])} New average time is {new_time / iterations:.2f} seconds

                {new_line.join(['',''])} Thats a {old_time / new_time:.2f} times speedup"""

if __name__ == "__main__":
    print('run')