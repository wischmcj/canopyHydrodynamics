from __future__ import annotations

import os
import sys
import multiprocessing as mp
from time import time 
from itertools import product
import gc
import pytest

sys.path.insert(0, os.path.dirname(os.getcwd()))
sys.path.insert(0, os.getcwd())
# print(sys.path)


from shapely.ops import unary_union
from data.output.run_so_far import already_run
from src.canhydro.global_vars import DIR, test_input_dir
from src.canhydro.utils import lam_filter
from src.canhydro.CylinderCollection import CylinderCollection, pickle_collection, unpickle_collection
from src.canhydro.Forester import Forester
from test.utils import within_range

from src.canhydro.global_vars import config_vars, log, output_dir


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
    # try: 
    statistics = collection.statistics(file_ext = case_name)
    # except Exception as e:
    #     log.info(f"Error gernerating stats for case {case_name} : {e}")
    #     return None
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
                generate_statistics(collection, f'{case_name}')
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

def sensitivity_analysis():
    # files_to_test = ["5_SmallTree"] 
    # angles = [.96,.16,.9]
    # files_to_test = [ "Secrest27-05_000000",
    #                  "Secrest32-06_000000"]
                    # ["Secrest02-30_000000"
                    #     ,"Secrest03-12_000000"
                    #     ,"Secrest07-32_000000"
                    #     ,"Secrest08-24c_000000"
                    #     ,"Secrest10-02_000000"]
    # , "Secrest27-05_000000","Secrest03-12_000000"
                    #     ,"Secrest07-32_000000"]
                        #  "Secrest02-26_000000"
    files_to_test =   [
                        "Secrest10-08_000000"
                        ,"Secrest11-27_000000"
                        ,"Secrest14-09_000000"
                        ,"Secrest16-3TI-CO_000000"
                        ,"Secrest16-14LI-ST_000000"]
                        # ,"Secrest18-13_000000"
                        # ,"Secrest23-23_000000"
                        # ,"Secrest24-03_000000"
                        # ,"Secrest24-07_000000"
                        # ,"Secrest26-03_000000"
                        

            #run or running ontowr 
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
    log.info(f'Will run {len(cases_to_run)} cases : {cases_to_run}')
    start = time()
    breakpoint()
    success = run_test_cases(cases_to_run, fig = False)
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
import pickle
import math
import matplotlib.pyplot as plt
import numpy as np
from geopandas import GeoSeries
import pandas as pd
import matplotlib.colors as colors

def retain_quantile(df, field, percentile):
    percentile_val = df[field].quantile(percentile)
    # print(f'percentile_val = {percentile_val} found for {field}, percentile {percentile}')
    return df[df[field] >= percentile_val]

def return_quantile(df, field, percentile):
    percentile_val = df[field].quantile(percentile)
    # print(f'percentile_val = {percentile_val} found for {field}, percentile {percentile}')
    return df[df[field] >= percentile_val][field]

def plot_drip_points(collection, percentile):
    scale = 1
    flows = pd.DataFrame([flow.__dict__ for flow in collection.flows])
    flows = flows[flows['drip_node_id']!=0]
    drip_points  = retain_quantile(flows, 'projected_area', percentile)
    drip_point_locs_x = [pt[0] * scale for pt in drip_points['drip_node_loc']]
    drip_point_locs_y = [pt[1] * scale for pt in drip_points['drip_node_loc']]
    drip_point_size = [pt for pt in drip_points['projected_area']]
    drip_node = [pt for pt in drip_points['drip_node_id']]
    return drip_point_locs_x,drip_point_locs_y,drip_point_size

def draw_for_paper():
    import matplotlib.pyplot as plt
    class nlcmap(object):
        def __init__(self, cmap, levels):
            self.cmap = cmapr
            self.N = cmap.N
            self.monochrome = self.cmap.monochrome
            self.levels = np.asarray(levels, dtype='float64')
            self._x = self.levels
            self.levmax = self.levels.max()
            self.transformed_levels = np.linspace(0.0, self.levmax,
                len(self.levels))

        def __call__(self, xi, alpha=1.0, **kw):
            yi = np.interp(xi, self._x, self.transformed_levels)
            return self.cmap(yi / self.levmax, alpha)
    
    
    # import numpy as np
    # import matplotlib.pyplot as plt

    # x = y = np.linspace(1, 10, 10)

    # t1mean, t2mean = 2, 9
    # sigma1, sigma2 = .3, .01
    # t1 = np.random.normal(t1mean, sigma1, 10)
    # t2 = np.random.normal(t2mean, sigma2, 10)

    class nlcmap(object):
        def __init__(self, cmap, levels):
            self.cmap = cmap
            self.N = cmap.N
            self.monochrome = self.cmap.monochrome
            self.levels = np.asarray(levels, dtype='float64')
            self._x = self.levels
            self.levmax = self.levels.max()
            self.transformed_levels = np.linspace(0.0, self.levmax,
                len(self.levels))

        def __call__(self, xi, alpha=1.0, **kw):
            yi = np.interp(xi, self._x, self.transformed_levels)
            return self.cmap(yi / self.levmax, alpha)

    # tmax = max(t1.max(), t2.max())
    # #the choice of the levels depends on the data:
    # levels = np.concatenate((
    #     [0, tmax],
    #     np.linspace(t1mean - 4 * sigma1, t1mean + 4 * sigma1, 5),
    #     np.linspace(t2mean - 4 * sigma2, t2mean + 4 * sigma2, 5),
    #     ))

    # levels = levels[levels <= tmax]
    # levels.sort()
    # # breakpoint()

    # cmap_nonlin = nlcmap(plt.cm.jet, levels)

    # fig, (ax1, ax2) = plt.subplots(1, 2)

    # ax1.scatter(x, y, edgecolors=cmap_nonlin(t1), s=15, linewidths=4)
    # ax2.scatter(x, y, edgecolors=cmap_nonlin(t2), s=15, linewidths=4)

    # fig.subplots_adjust(left=.25)
    # cbar_ax = fig.add_axes([0.10, 0.15, 0.05, 0.7])

    # #for the colorbar we map the original colormap, not the nonlinear one:
    # sm = plt.cm.ScalarMappable(cmap=plt.cm.jet, 
    #                 norm=plt.Normalize(vmin=0, vmax=tmax))
    # sm._A = []
    # plt.show()


    # breakpoint()
    files = [ 
        # ('/code/code/canopyHydrodynamics/data/output/pickle/Secrest32-06_000000_pickle__stats_-0.36','32-06_low_drip_points_36')
               # #  ('/code/code/canopyHydrodynamics/data/output/pickle/Secrest32-06_000000_pickle__stats_0.04','32-06_high_drip_points_04')
               # #  ('/code/code/canopyHydrodynamics/data/output/pickle/Secrest32-06_000000_pickle__stats_-0.14','32-06_mid_drip_points_14')
               ('/code/code/canopyHydrodynamics/data/output/pickle/Secrest27-05_000000_pickle__stats_-0.34','27-05_low_end_drip_points_34'),
            #    ('/code/code/canopyHydrodynamics/data/output/pickle/Secrest27-05_000000_pickle__stats_0.04.pickle','27-05_high_end_drip_points_04'),
            #    ('/code/code/canopyHydrodynamics/data/output/pickle/Secrest27-05_000000_pickle__stats_-0.14.pickle','27-05_mi_drip_points_0014')
                ]
    # db_files = [(open(file_i[0],'rb'),file_i[0]) for file_i in files]
    for file_i in files:
        col = open(file_i[0],'rb')
        file = file_i[1]
    # #         open('/code/code/canopyHydrodynamics/data/output/pickle/Secrest32-06_000000_pickle__stats_-0.36', 'rb')]
    # # file = '32-06_low_drip_points_36'
    # #         open('/code/code/canopyHydrodynamics/data/output/pickle/Secrest32-06_000000_pickle__stats_0.04', 'rb')]
    # # file = '32-06_high_drip_points_04'
    # #         open('/code/code/canopyHydrodynamics/data/output/pickle/Secrest32-06_000000_pickle__stats_-0.14', 'rb')]
    # # file = '32-06_mid_drip_points_14

    #     open('/code/code/canopyHydrodynamics/data/output/pickle/Secrest27-05_000000_pickle__stats_-0.34', 'rb')]
    # file = '27-05_low_end_drip_points_34'
    # # pick_file = 'Secrest27-05_000000'
    # # pik_designation = 'stats_-0.34'
    # #   open('/code/code/canopyHydrodynamics/data/output/pickle/Secrest27-05_000000_pickle__stats_0.04.pickle', 'rb')]
    # # file = '27-05_high_end_drip_points_04'
    # #     open('/code/code/canopyHydrodynamics/data/output/pickle/Secrest27-05_000000_pickle__stats_-0.14.pickle', 'rb')]
    # # file = '27-05_mi_drip_points_14'
    

        # collections = [pickle.load(db_file[0]) for db_file in db_files]
        collections = [pickle.load(col)]
    # collections[0].project_cylinders('XZ')
    # pickle_collection(collections[0],designation=pik_designation)
    # breakpoint()

        flow_info = [plot_drip_points(col, .98) for col in collections]
        drip_info = [x for x in flow_info]
        # breakpoint()
        # max_drip = 0
        # for x,y,drip_point_size in drip_info:
        #     potential_max_drip = np.max(drip_point_size)
        #     if potential_max_drip>max_drip: 
        #         max_drip = potential_max_drip
    #32-06
    #10.613
    #19.25299317980633
    #28.45471360468393

    #27-05
    #6.603
    #6.603
    #6.603
        ghost_tree = True
        only_hulls = False
        figs_to_draw = any([ghost_tree, only_hulls])

        if figs_to_draw:
            for idx, graph_data in enumerate(drip_info):
                x,y,drip_point_size = graph_data
                if '32' in file:
                    x_trans = 9
                    y_trans = 4
                    x_lim = [0,14]
                    y_lim = [0,13]

                    # x.append(-20)
                    # y.append(-20)
                    # drip_point_size.append(28.4547)

                if '27' in file:
                    x_trans = 2
                    y_trans = 8
                    x_lim = [0,12]
                    y_lim = [0,11]
                    # x.append(-20)
                    # y.append(-20)
                    
                    # drip_point_size.append(6.603)
                    
                col = collections[idx]
                polys = unary_union([cyl.projected_data['XY']['polygon'] for cyl in col.cylinders])
                stem_hull = col.stem_hull
                geopolys = GeoSeries(polys)
                geohull = GeoSeries(stem_hull)
                tot_geoHull = GeoSeries(col.hull)

                # fig, ax = plt.subplots()
                # ext = 'stem_highlight_xz'
                # try:
                #     polys = [cyl.projected_data['XZ']['polygon'] for cyl in col.cylinders]
                #     geopolys = GeoSeries(polys)
                #     colors = ["Blue" if cyl.is_stem else "Grey" for cyl in col.cylinders]
                #     geopolys_trans = geopolys.translate( xoff=x_trans, yoff=y_trans, zoff=0.0)
                #     geopolys_trans.plot(ax=ax, color=colors)
                #     ax.set_xlim(x_lim)
                #     ax.set_ylim(y_lim)
                #     plt.show()
                #     # plt.savefig(f'/code/code/canopyHydrodynamics/data/output/draw/{file}_{ext}.svg')
                #     plt.savefig(f'/media/penguaman/Healthy/BranchHighlight/branchHighlight/{file}_{ext}.svg')
                # except KeyError as e:
                #     log.error(f'error getting XZ projected data {e}')
                #             #   trying to run project_cylinders {e}')
                #     col.project_cylinders('XZ')
                #     polys = [cyl.projected_data['XZ']['polygon'] for cyl in col.cylinders]
                

                
                # fig, ax = plt.subplots()
                # ext = 'stem_highlight_xy'
                # polys = [cyl.projected_data['XY']['polygon'] for cyl in col.cylinders]
                # geopolys = GeoSeries(polys)
                # colors = ["Blue" if cyl.is_stem else "Grey" for cyl in col.cylinders]
                # geopolys_trans = geopolys.translate( xoff=x_trans, yoff=y_trans, zoff=0.0)
                # geopolys_trans.plot(ax=ax, color=colors)
                # ax.set_xlim(x_lim)
                # ax.set_ylim(y_lim)
                # plt.savefig(f'/media/penguaman/Healthy/BranchHighlight/branchHighlight/{file}_{ext}.svg')

                # plt.show()


                
                # t1mean, t2mean = np.mean(drip_info), 9
                # sigma1, sigma2 = .3, .01
                # tmax = np.max(drip_info)
                # #the choice of the levels depends on the data:
                # levels = np.concatenate((
                #     [0, tmax],
                #     np.linspace(t1mean - 4 * sigma1, t1mean + 4 * sigma1, 5),
                #     np.linspace(t2mean - 4 * sigma2, t2mean + 4 * sigma2, 5),
                #     ))

                # levels = levels[levels <= tmax]
                # levels.sort()

                # cmap_nonlin = nlcmap(plt.cm.jet, levels)


                t1mean = np.mean(drip_info)
                sigma1 = np.std(drip_info)
                t1 = np.random.normal(t1mean, sigma1, 10)

                tmax = np.max(drip_info)
                #the choice of the levels depends on the data:
                levels = np.concatenate((
                    [0, tmax],
                    np.linspace(t1mean - 4 * sigma1, t1mean + 4 * sigma1, 5),
                    ))

                levels = levels[levels <= tmax]
                levels.sort()

                cmap_nonlin = nlcmap(plt.cm.jet, levels)

                # if ghost_tree:
                fig, ax = plt.subplots()
                ext = 'ghost_tree_cbar'        
                # geopolys_trans = geopolys.translate( xoff=x_trans, yoff=y_trans, zoff=0.0)
                geohull_trans = geohull.translate( xoff=x_trans, yoff=y_trans, zoff=0.0)

               
                geohull_trans.plot(ax=ax, color='darkgrey', alpha = 0.5)

                # geopolys_trans.plot(ax=ax, color='lightgrey', alpha = 0.7)
                CS = ax.scatter([x_val+x_trans for x_val in x], 
                                [y_val+y_trans for y_val in y], 
                                edgecolors=cmap_nonlin(drip_info),
                                 s=15, linewidths=4
                                # c = [x if x > .1  else .1 for x in drip_point_size], 
                                # cmap = 'Greys', edgecolors='slategray'
                                )
                cbar = plt.colorbar(CS)
                ax.set_xlim(x_lim)
                ax.set_ylim(y_lim)
                plt.show()
                breakpoint()
                # norm=colors.LogNorm(vmin=np.min(drip_point_size), vmax=np.min(drip_point_size)),
                # plt.savefig(f'/media/penguaman/Healthy/BranchHighlight/dripsWithTrees/{file}_{ext}.svg')

                # fig, ax = plt.subplots()
                # ext = 'two_hulls'
                # tot_geoHull = GeoSeries(col.hull)
                # tot_geoHull_trans = tot_geoHull.translate( xoff=x_trans, yoff=y_trans, zoff=0.0)
                # geohull_trans = geohull.translate( xoff=x_trans, yoff=y_trans, zoff=0.0)
                # tot_geoHull_trans.plot(ax=ax, color='darkgrey', alpha = 0.3)
                # geohull_trans.plot(ax=ax, color='darkgrey', alpha = 0.3)
                # ax.scatter([x_val+x_trans for x_val in x], [y_val+y_trans for y_val in y], facecolors='none', edgecolors='slategray', s=[x*5 for x in drip_point_size] )
                # ax.set_xlim(x_lim)
                # ax.set_ylim(y_lim)
                # plt.savefig(f'/media/penguaman/Healthy/BranchHighlight/dripsWithHulls/{file}_{ext}.svg')
                # plt.show()

            
                # fig, ax = plt.subplots()
                # ext = 'stem_hull'
                # geohull_trans = geohull.translate( xoff=x_trans, yoff=y_trans, zoff=0.0)
                # geohull_trans.plot(ax=ax, color='darkgrey', alpha = 0.3)
                # ax.scatter([x_val+x_trans for x_val in x], [y_val+y_trans for y_val in y], facecolors='none', edgecolors='slategray', s=[x*5 for x in drip_point_size] )
                # ax.set_xlim(x_lim)
                # ax.set_ylim(y_lim)
                # plt.savefig(f'/media/penguaman/Healthy/BranchHighlight/dripsWithHulls/{file}_{ext}.svg')
                # plt.show()

                # if only_hulls:
                #     fig, ax = plt.subplots()
                #     ext = 'only_hulls'
                #     tot_geoHull = GeoSeries(col.hull)
                #     tot_geoHull_trans = tot_geoHull.translate( xoff=x_trans, yoff=y_trans, zoff=0.0)
                #     geohull_trans = geohull.translate( xoff=x_trans, yoff=y_trans, zoff=0.0)
                #     tot_geoHull_trans.plot(ax=ax, color='darkgrey', alpha = 0.3)
                #     geohull_trans.plot(ax=ax, color='darkgrey', alpha = 0.3)
                #     ax.set_xlim(x_lim)
                #     ax.set_ylim(y_lim)
                #     plt.savefig(f'/media/penguaman/Healthy/BranchHighlight/hulls/{file}_{ext}.svg')
                #     plt.show()

        # except Exception as e:
        #     print(e)
            # plt.show()

    # filename='32_med_hull_dark'
    # plt.savefig('./data/output/draw/{filename}.png')
    # plt.savefig(f'/home/penguaman/Desktop/sensitivityAnalysis/Hulls/{filename}.png')
    # breakpoint()


    # fig, ax = plt.subplots()

    # geopolys.plot(ax=ax,color='lightgray',alpha=0.5)
    # geohull.plot(ax=ax, color='darkgray',alpha=0.6)




    # forest = Forester()
    # forest.get_file_names(dir=test_input_dir)
    # forest.qsm_from_file_names(file_name="Secrest32-06_000000")
    # basic_collection = forest.cylinder_collections[0]
    # basic_collection.project_cylinders("XY")
    # basic_collection.initialize_digraph_from(in_flow_grade_lim=-0.3)
    # basic_collection.find_flow_components()
    # basic_collection.calculate_flows()
    # breakpoint()
    # pickle_collection(basic_collection,basic_collection.file_name)

    
    # forest_old = Forester()
    # forest_old.get_file_names(dir=test_input_dir)
    # forest_old.qsm_from_file_names(file_name="5_SmallTree.csv")
    # basic_collection_old = forest_old.cylinder_collections[0]

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


    # # breakpoint()
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
