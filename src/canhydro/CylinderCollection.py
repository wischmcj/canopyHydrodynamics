"""These constitute our way of 'copying' the QSMs"""


from __future__ import annotations

import multiprocessing as mp
import copy
import math
import sys
import pickle
import os

from itertools import chain

import networkx as nx
# import rustworkx as rx
import numpy as np
# from memory_profiler import LogFile
# from scipy.spatial import distance
from shapely.geometry import Point
from shapely.ops import unary_union

from src.canhydro.Cylinder import create_cyl
from src.canhydro.DataClasses import Flow
from src.canhydro.geometry import (concave_hull, draw_cyls, furthest_point,
                                   get_projected_overlap,
                                    pool_get_projection)
# vectorized_get_projection
from src.canhydro.global_vars import config_vars, log, output_dir, output_dir
from src.canhydro.utils import intermitent_log, lam_filter, save_file, create_dir_and_file, create_dir_and_file
# sys.stdout = LogFile()

# import matplotlib.pyplot as plt

NAME = "CylinderCollection"


def pickle_collection(collection, designation: str = ""):   
    # file_path = "".join([output_dir, "pickle\\", f'{collection.file_name}_pickle'])
    if designation == "": designation = collection.file_name
    out_dir = output_dir.__str__()
    file_path ="".join([f'{out_dir}/pickle/', f'{collection.file_name.replace(".csv","")}_pickle_{designation}'])
    directory = os.path.dirname(file_path)
    create_dir_and_file(directory)
    pickle_file = open(file_path, 'ab')
    pickle.dump(collection, pickle_file) 
    pickle_file.close()
    return pickle_file

def unpickle_collection(file_name:str):
    out_dir = output_dir.__str__()
    file_path ="".join([f'{out_dir}/pickle/', file_name])
    dbfile = open(file_path, 'rb')    
    db = pickle.load(dbfile)
    dbfile.close()
    return db

# By inheriting the Model class, lambda cyl : cyl.branch_order = br CC gains managed functionality- lidirke lambda searching
class CylinderCollection:
    cylinders = {}

    # initialize our object level variables for cylinder objects
    def __init__(self) -> None:
        self.file = ""
        # self.collection = CylinderList()
        # Aggregate values from file
        self.surface_area = np.nan
        self.cyl_index = {}
        self.file_name = ""
        self.volume = np.nan
        self.avg_sa_to_vol = np.nan
        self.max_branch_order = np.nan
        self.max_reverse_branch_order = np.nan
        self.canopy_scope = np.nan  # desc of canopy
        self.extent = {
            "min": [np.nan, np.nan, np.nan],
            "max": [np.nan, np.nan, np.nan],
        }
        # to populate with x,y,z mins and maxs
        self.aggregate_angle = np.nan
        self.descriptive_vectors = np.nan  # Average, median, mode vectors
        self.treeQualities = {
            "total_psa": None,
            "tot_hull_area": None,
            "stem_flow_hull_area": None,
            "stem_psa": None,
            "flowStats": None,
            "dbh": None,
            "tot_surface_area": None,
            "stem_surface_area": None,
        }

        # Projection Attrs
        self.projections = {"XY": False, "XZ": False, "YZ": False}
        self.stem_path_lengths = []
        self.hull = None
        self.stem_hull = None

        # Special case tree attributes
        self.stem_paths = [[]]  # Cyl collection?
        self.trunk = []  # Collection of cylinders? id list?

        # Graph and Attributes
        self.graph = None
        self.digraph = None
        self.min_graph = None
        self.flows = None
        self.divide_nodes = None
        self.drip_nodes = None
        self.drip_node_to_cyl = None
        self.drip_points = {"x": np.nan, "y": np.nan, "z": np.nan, "flow_id": np.nan}
        self.flow_to_drip = {
            0: 1
        }  # A dictionary of flow ids with values equal to their drip node ids
        self.trunk_nodes = []
        self.drip_point_loc = None
        self.stem_flow_component = None
        self.drip_graph = None
        self.drip_flow_components = None
        # Calculations using graph results
        self.flow_chars = {}
        self.divide_points = []
        self.stemPolys = []
        self.compGraphs = []
        self.trunk_lean = None
        self.cyl_to_drip_node = {}

    def drip_summary(self)->str: 
        ret = ''
        ret += f'cyl_to_drip_node {self.cyl_to_drip_node} '
        ret += f'drip_flow_components {self.drip_flow_components} '
        ret += f'drip_points {self.drip_points} '
        ret += f'drip_nodes {self.drip_nodes} '
        return ret

    def get_trunk_nodes(self) -> list[int]:
        if self.trunk_nodes:
            return self.trunk_nodes
        elif len(self.cylinders) > 0:
            trunk_cyls, _ = lam_filter(self.cylinders, lambda: branch_order == 0)
            trunk_nodes = [cyl.cyl_id for cyl in trunk_cyls]
            self.trunk_nodes = trunk_nodes
            return trunk_nodes
        else:
            log.warning(
                "Graph not initialized, run <CylinderCollection>.initialize_graph(**args)"
            )
            return list(None)

    def from_csv(self, file, aggregate_cyls=True):
        """Initializes a new Cyl Collection based on the data in a QSM
        with the configured column locations"""
        self.file = file
        self.file_name = file.name
        log.info(f"Processing {str(file)}")
        # self.arr = pd.read_csv(file, header=0)
        self.arr = np.genfromtxt(file, delimiter=",", skip_header=True)[0:, :-1]
        cylinders = [create_cyl(row) for row in self.arr]
        self.cylinders = cylinders
        # self.cyl_index = dict((cyl.cyl_id , cyl ) for cyl in cylinders)
        if aggregate_cyls:
            min_x = np.min([cyl.x[0] for cyl in cylinders])
            min_y = np.min([cyl.y[0] for cyl in cylinders])
            min_z = np.min([cyl.z[0] for cyl in cylinders])
            max_x = np.max([cyl.x[1] for cyl in cylinders])
            max_y = np.max([cyl.y[1] for cyl in cylinders])
            max_z = np.max([cyl.z[1] for cyl in cylinders])
            # Aggregate values from file
            self.no_cylinders = len(cylinders)
            self.surface_area = np.sum([cyl.surface_area for cyl in cylinders])
            self.volume = np.sum([cyl.volume for cyl in cylinders])
            self.max_branch_order = np.max([cyl.branch_order for cyl in cylinders])
            self.max_reverse_branch_order = np.max(
                [cyl.reverse_branch_order for cyl in cylinders]
            )
            self.avg_sa_to_vol = (
                np.sum([cyl.sa_to_vol for cyl in cylinders]) / self.no_cylinders
            )
            self.extent = {
                "min": [min_x, min_y, min_z],
                "max": [max_x, max_y, max_z],
            }

        self.descriptive_vectors = np.nan  # Average, median, mode vectors

        self.theta = np.nan
        log.info(f"{file.name} initialized with {self.no_cylinders} cylinders")

    def project_cylinders(self, plane: str = "XY", force_rerun: bool = False):
        """Projects cylinders onto the specified plane"""
        if plane not in ("XY", "XZ", "YZ"):
            log.info(f"{plane}: invalid value for plane")
        elif not force_rerun and self.projections[plane]:
            log.info(
                "cached projections exist, pass 'force_rerun=True to calculate new projections "
            )
        else:
            polys = [None] * self.no_cylinders
            for idx, cyl in enumerate(self.cylinders):
                poly = cyl.get_projection(plane)
                polys[idx] = poly

            self.projections[plane] = True
            self.pSV = polys

    def pool_project_cylinders(self, plane: str = "XY", force_rerun: bool = False):
        """Projects cylinders onto the specified plane"""
        if plane not in ("XY", "XZ", "YZ"):
            log.info(f"{plane}: invalid value for plane")
        elif not force_rerun and self.projections[plane]:
            log.info(
                "cached projections exist, pass 'force_rerun=True to calculate new projections "
            )
        else:
            log.info(f"Projection into {plane} axis begun for file {self.file_name}")
            # log.info(f'Container has {mp.cpu_count()} cores')

            with mp.Pool(5) as p:
                task_pool = [p.apply_async(pool_get_projection, args=(cyl, plane,)) for cyl in self.cylinders]
                polys = [task.get() for task in task_pool]
            
            self.projections[plane] = True
            self.pSV = polys

    # def numba_project_cylinders(self, plane: str = "XY", force_rerun: bool = False):
    #     """Projects cylinders onto the specified plane"""
    #     if plane not in ("XY", "XZ", "YZ"):
    #         log.info(f"{plane}: invalid value for plane")
    #     elif not force_rerun and self.projections[plane]:
    #         log.info(
    #             "cached projections exist, pass 'force_rerun=True to calculate new projections "
    #         )
    #     else:
    #         polys = []
    #         log.info(f"Projection into {plane} axis begun for file {self.file_name}")
    #         for idx, cyl in enumerate(self.cylinders):
    #             poly = cyl.numba_get_projection(plane)
    #             polys.append(poly)
    #             # print a progress update once every 10 thousand or so cylinders
    #             intermitent_log(idx, self.no_cylinders, "Cylinder projection: ")
    #         # Used by other functions to know what projections have been run
    #         self.projections[plane] = True
    #         self.pSV = polys

    # def vectorized_project_cylinders(
    #     self, plane: str = "XY", force_rerun: bool = False
    # ):
    #     """Projects cylinders onto the specified plane"""
    #     if plane not in ("XY", "XZ", "YZ"):
    #         log.info(f"{plane}: invalid value for plane")
    #     # elif not force_rerun and self.projections[plane]:
    #     #     log.info(
    #     #         "cached projections exist, pass 'force_rerun=True to calculate new projections "
    #     #     )
    #     else:
    #         polys = []
    #         log.info(f"Projection into {plane} axis begun for file {self.file_name}")
    #         starts = np.array([cyl.vectors[plane][0] for cyl in self.cylinders])
    #         ends = np.array([cyl.vectors[plane][1] for cyl in self.cylinders])
    #         radii = np.array([cyl.radius for cyl in self.cylinders])
    #         vectorized_get_projection(starts, ends, radii)

    #         self.projections[plane] = True
    #         self.pSV = polys

    def get_collection_data(self):
        cyl_desc = [cyl.__repr__() for cyl in self.cylinders]
        return cyl_desc

    def __eq__(self, other):
        if len(self.cylinders) == 0 or len(other.cylinders) == 0:
            raise AttributeError(
                "One or both Cylinder Collections contain no Cylinders (Did you forget to initialize?)"
            )
        if type(self) == type(other):
            raise TypeError(
                "CylinderCollections may only be compared to other cylinder collections"
            )
        # order matters here
        return np.all([cyl == other[idx] for idx, cyl in enumerate(self.cylinders)])

    def draw(
        self,
        plane: str = "XY",
        highlight_lambda: function = lambda: True,
        filter_lambda: function = lambda: True,
        include_drips: bool = False,
        include_contour: bool = False,
        include_alpha_shape:bool = False ,
        stem = False,
        **args,
    ):
        """Draws cylinders meeting given characteristics onto the specified plane"""
        if plane not in ("XY", "XZ", "YZ"):
            log.info(f"{plane}: invalid value for plane")
        if not self.projections[plane]:
            self.project_cylinders(plane)
        cylinders, _ = lam_filter(self.cylinders, filter_lambda)
        filtered_cyls, matches = lam_filter(
            cylinders, highlight_lambda, return_all=True
        )
        to_draw = [cyl.projected_data[plane]["polygon"] for cyl in filtered_cyls]
        log.info(f"{len(to_draw)} cylinders matched criteria")
        if include_drips:
            self.drip_map()
        if include_contour:
            self.drip_map()
        # if include_alpha_shape:
        #     self.drip_map()
        # if stem:
        #     self.drip_map()
        fig =draw_cyls(collection=to_draw, colors=matches, **args)
        return fig

    def get_dbh(self):
        """a real trainwreck of a function to find dbh"""
        if self.treeQualities["dbh"]:
            return self.treeQualities["dbh"]
        else:
            desired_height = 1.3
            approx_height_range_start = desired_height - 0.1
            start_z = self.extent["min"][2]
            # trunk = lam_filter(
            #     self.cylinders,
            #     a_lambda=lambda: branch_order == 0,
            # )
            diff_from_desired_height = [
                (abs(cyl.z[0] - desired_height), cyl.cyl_id)
                for cyl in self.cylinders
                if cyl.z[0] >= start_z + approx_height_range_start
                and cyl.branch_order == 0
            ]
            closest_to_target_id = np.argmin(np.array(diff_from_desired_height[:][0]))
            rbh = self.cylinders[
                int(diff_from_desired_height[closest_to_target_id][1])
            ].radius
            self.treeQualities["dbh"] = 2 * rbh

    def get_end_nodes(self) -> list[int]:
        g = self.digraph
        if self.end_nodes:
            return self.end_nodes
        elif len(g.nodes) > 0:
            end_nodes = [n for n in g.nodes if g.degree(n) == 1 and n != -1]
            self.end_nodes = end_nodes
            return end_nodes
        else:
            log.warning(
                "Graph not initialized, run <CylinderCollection>.initialize_graph(**args)"
            )
            return list(None)


    def watershed_boundary(
        self,
        component=None,
        plane: str = "XY",
        curvature_alpha: np.float64 = 2,
        stem: bool = False,
        draw: bool = False,
        save: bool = False,
        file_ext: str = ''
    ) -> None:
        """Generates tightly fit concave_hull (alpha shape) around the passed component"""
        """Alpha determines the tightness of the fit of the shape. Too low an alpha results in multiple"""
        if not self.projections[plane]:
            self.project_cylinders(plane)
        # want to get a ring of points around the component
        # however, the root node also has degree one, so we exclude it with n!= -1
        component = component if component else self.digraph
        end_nodes = [n for n in component.nodes()
                        if component.in_degree(n) + component.out_degree(n)  == 1 
                            and n != 0]

        endCyls = [
            cyl
            for cyl in self.cylinders
            if cyl.cyl_id in end_nodes or cyl.branch_order == 0
        ]
        boundary_points = [Point(cyl.x[1], cyl.y[1]) for cyl in endCyls]

        hull, _, _ = concave_hull(boundary_points, curvature_alpha)
        if draw:
            draw_cyls([hull], save = save, file_ext = file_ext)
        if stem:
            self.stem_hull = hull
        else:
            self.hull = hull

    def initialize_graph_from(self):
        """This function creates an undirected_graph and initializes edge attributes as cylinder objects"""
        gr = nx.Graph()
        edges = (
            (int(cyl.cyl_id), int(cyl.parent_id), {"cylinder": cyl})
            for cyl in self.cylinders
        )
        gr.add_edges_from(edges)
        self.graph = gr

    def initialize_digraph_from(
        self, in_flow_grade_lim=config_vars["in_flow_grade_lim"]
    ):
        """This function creates a directed graph and its undirected counterpart.
        Initializes edge attributes as cylinder objects"""
        gr = nx.DiGraph()
        trunk_nodes = self.get_trunk_nodes()
        edges = (
            (
                (int(cyl.cyl_id), int(cyl.parent_id), {"cylinder": cyl})
                if (cyl.angle >= in_flow_grade_lim or cyl.cyl_id in trunk_nodes)
                else (int(cyl.parent_id), int(cyl.cyl_id), {"cylinder": cyl})
            )
            for cyl in self.cylinders
        )
        gr.add_edges_from(edges)
        self.digraph = gr
        self.graph = gr.to_directed()

    def contracted_nodes(self, node_list, next_node, G):
        """Returns a graph made from removing nodes from G,
        and connecting the removed nodes' neighbors via an edge"""
        # finding and removing the listed nodes to find their edges
        contracted = copy.deepcopy(G)
        contracted.remove_nodes_from(node_list)

        neighbors = [
            node
            for node in contracted.nodes()
            if (G.degree(node) != contracted.degree(node))
        ]

        contracted.add_node(next_node)
        edges_to_trunk = [(u, next_node) for u in neighbors]
        contracted.add_edges_from(edges_to_trunk)

        return contracted, neighbors

    # def has_path(self, end, start):
    #     g = self.digraph
    #     ancestors = nx.ancestors(g, end )
    #     if start in ancestors:
    #         return ancestors
    #     else:
    #         return None

    # def get_path(self, end, start):
    #     g = self.digraph
    #     ancestors = self.has_path(end, start)
    #     if not ancestors:
    #         log.error(f'No path exists between {end} and {start}') 
    #         raise ValueError(f'No path exists between {end} and {start}')
    #     else:
    #         g.edges(end)
    #     return start in ancestors

    def find_trunk_distances(self):
        """
        Finds the distance in the graph (in number of nodes) between each node and the closest trunk node
        """
        trunk_nodes = self.get_trunk_nodes()

        trunk_contraction, titans = self.contracted_nodes(trunk_nodes, 0, self.graph)
        trunk_paths = nx.shortest_path(trunk_contraction, target=0)
        dists = {node: len(path) - 1 for node, path in trunk_paths.items()}
        return dists

    # def initialize_digraph_from_new(
    #         self, in_flow_grade_lim=config_vars["in_flow_grade_lim"]
    #     ):
    #     """This function creates a directed graph and its undirected counterpart.
    #     Initializes edge attributes as cylinder objects"""
    #     log.info(
    #             "Begining initializing graphs"
    #     )
    #     gr = rx.PyDiGraph()
       
    #     edges = [
    #         (
    #             (int(cyl.cyl_id +1), int(cyl.parent_id+1), cyl.cyl_id)
    #             if cyl.angle >= in_flow_grade_lim
    #             else (int(cyl.parent_id+1), int(cyl.cyl_id+1), cyl.cyl_id)
    #         )
    #         for cyl in self.cylinders
    #     ]
    #     gr.add_nodes_from([int(cyl.cyl_id + 1) for cyl in self.cylinders])
    #     gr.add_nodes_from([0])
    #     gr.extend_from_weighted_edge_list(edges)
    #     self.digraph = gr

    # def get_end_nodes_new(self) -> list[int]:
    #     g = self.digraph
    #     if self.end_nodes:
    #         return self.end_nodes
    #     elif len(g.nodes) > 0:
    #         end_nodes = [n for n in g.nodes if (g.in_degree(n) + g.out_degree(n)) == 1 and n != -1]
    #         self.end_nodes = end_nodes
    #         return end_nodes
    #     else:
    #         log.warning(
    #             "Graph not initialized, run <CylinderCollection>.initialize_graph(**args)"
    #         )
    #         return list(None)

    # def get_trunk_nodes_new(self) -> list[int]:
    #     g = self.digraph
    #     if self.trunk_nodes:a
    #         return self.trunk_nodes
    #     elif len(g.nodes) > 0:
    #         trunk_cyls, _ = lam_filter(self.cylinders, lambda: branch_order == 0)
    #         trunk_nodes = [cyl.cyl_id for cyl in trunk_cyls]
    #         self.trunk_nodes = trunk_nodes
    #         return trunk_nodes
    #     else:
    #         log.warning(
    #             "Graph not initialized, run <CylinderCollection>.initialize_graph(**args)"
    #         )
    #         return list(None)
        


    # def find_drip_component_new(self, idx, pair):
    #     drip_node, source_divides = pair
    #     # import pdb; pdb.set_trace()
    #     paths = [
    #         rx.digraph_dijkstra_shortest_paths(self.drip_graph, source =div_node, target = drip_node)[drip_node]
    #         for div_node in source_divides
    #     ]
    #     component_nodes = [node for node in set(chain.from_iterable(paths))]
    #     component_graph = self.drip_graph.subgraph(component_nodes,preserve_attrs=True).copy()
    #     component_cyls = [
    #         attr
    #         for _, _, attr in component_graph.weighted_edge_list()
    #     ]

    #     for cyl_id in component_cyls:
    #        self.cyl_to_drip_node[int(cyl_id)].append(drip_node)
    #     log.info(
    #         f"component cyls  {component_cyls}"
    #     )
    #     if idx%50 == 0:
    #         log.info(
    #             f"completed drip component {idx}"
    #         )

    #     return (drip_node,component_cyls)

    # def find_flow_components_new(self, inFlowGradeLim=-1 / 6):
    #     log.info(
    #             "Finding flow components"
    #     )
    #     self.cyl_to_drip_node = {int(cyl.cyl_id): [] for cyl in self.cylinders}
    #     g = self.digraph
    #     root_node = 0


    #     g_drip = copy.deepcopy(g)
    #     divide_nodes = [
    #         node
    #         for node in g.nodes()
    #         if node and ( g.out_degree(node) > 1
    #         or (g.out_degree(node) == 1 and g.in_degree(node) == 0))
    #     ]
    #     drip_nodes = [
    #         node
    #         for node in g.nodes()
    #         if node and g.out_degree(node) == 0
    #     ]
    #     stem_comp_nodes = list(rx.ancestors(g, root_node))
    #     stem_component_edges = [(u,v,cyl_id) for (u,v,cyl_id)  
    #                                 in g.weighted_edge_list() 
    #                                     if u in stem_comp_nodes and v in stem_comp_nodes]
    #     stem_flow_component = rx.PyDiGraph()

    #     stem_flow_component.add_nodes_from(stem_comp_nodes)
    #     stem_flow_component.add_nodes_from([0])
    #     stem_flow_component.extend_from_weighted_edge_list(stem_component_edges)

    #     # [tup for tup in self.digraph.edge_list() if tup not in [tup for tup in stem_flow_component.edge_list()]]
    #     stem_cylinders = [ cyl_id for _,_,cyl_id in stem_component_edges ]
    #     log.info(
    #             "Setting is stem"
    #     )
    #     for cyl in self.cylinders:
    #         if cyl.cyl_id in stem_cylinders:
    #             cyl.is_stem = True
    #     log.info(
    #                     f"creating g drip: divde nodes {len(divide_nodes)} drip nodes {len(drip_nodes)}"
    #             )
    #     g_drip = copy.deepcopy(g)

    #     g_drip.remove_edges_from([(u,v) for (u,v) in stem_flow_component.edge_list()])


    #     drip_divide_pairings = [
    #         (
    #             drip_node,
    #             [node for node in divide_nodes if rx.has_path(g_drip, node, drip_node)],
    #         )
    #         for drip_node in drip_nodes
    #     ]
    #     log.info(
    #             "starting proecessing of drip components "
    #     )
    #     drip_components = []

    #     log.info(
    #             "Assessing drip points"
    #     )

    #     self.drip_graph = g_drip

    #     num_pairings = len(drip_divide_pairings)

    #     for d, source_divides in drip_divide_pairings:
    #     # import pdb; pdb.set_trace()
    #         breakpoint()
    #         paths = [
    #             rx.digraph_dijkstra_shortest_paths(self.drip_graph, source =div_node, target = drip_node)[drip_node]
    #             for div_node in source_divides
    #         ]
    #         component_nodes = [node for node in set(chain.from_iterable(paths))]
    #         component_graph = self.drip_graph.subgraph(component_nodes,preserve_attrs=True).copy()
    #         component_cyls = [    attr    for _, _, attr in component_graph.weighted_edge_list()]
    #         for cyl_id in component_cyls:
    #             self.cyl_to_drip_node[int(cyl_id)].append(drip_node)
    #         log.info(
    #             f"component cyls  {component_cyls}"
    #         )
    #         if idx%50 == 0:
    #             log.info(
    #                 f"completed drip component {idx}"
    #             )
        
    #     # with mp.Pool(5) as p:
    #     #         task_pool = [p.apply_async(self.find_drip_component, args=(idx,pair)) 
    #     #                         for idx, pair in enumerate(drip_divide_pairings)]
    #     #         component_cyl_tuples = [task.get() for task in task_pool]

    #     self.drip_node_to_cyl = { k:v for k, v in component_cyl_tuples}


    #     log.info(
    #         f"{self.file_name} found to have {len(drip_components)} drip components"
    #     )

    #     log.info(self.cyl_to_drip_node)

    #     self.stem_flow_component = stem_flow_component
    #     # self.drip_flow_components = drip_components
    #     self.divide_nodes = divide_nodes
    #     self.drip_nodes = drip_nodes


    def find_flow_components(self):
        g = self.digraph
        if type(g) == nx.Graph:
            msg = "Find Flow Digraph invoked for undirected graph"
            log.error(msg)
            raise TypeError(msg)
        root_node = -1

        divide_nodes = [
            node
            for node in g.nodes()
            if g.out_degree(node) > 1
            or (g.out_degree(node) == 1 and g.in_degree(node) == 0)
        ]
        trunk_nodes = self.get_trunk_nodes()
        drip_nodes = [
            node
            for node, out_degree in g.out_degree()
            if out_degree == 0 and node != -1 and node not in trunk_nodes
        ]
        stem_flow_component = g.subgraph(nx.ancestors(g, root_node) | {0}).copy()

        stem_cylinders = [
            node
            for node in stem_flow_component.nodes()
            if stem_flow_component.degree(node) > 0
        ]
        for cyl in self.cylinders:
            if cyl.cyl_id in stem_cylinders:
                cyl.is_stem = True

        g_drip = copy.deepcopy(g)
        g_drip.remove_edges_from(stem_flow_component.edges())
        # g_drip.remove_nodes_from(stem_flow_component.nodes())

        drip_divide_pairings = [
            (
                drip_node,
                [node for node in divide_nodes if nx.has_path(g_drip, node, drip_node)],
            )
            for drip_node in drip_nodes
            if drip_node != -1
        ]

        drip_components = []
        cyl_to_drip_node = {cyl.cyl_id: [] for cyl in self.cylinders}
        for pair in drip_divide_pairings:
            drip_node, source_divides = pair
            paths = [
                nx.shortest_path(g_drip, div_node, drip_node)
                for div_node in source_divides
            ]
            component_nodes = chain.from_iterable(paths)
            component_graph = g_drip.subgraph(component_nodes).copy()
            component_cyls = [
                attr["cylinder"].cyl_id
                for _, _, attr in component_graph.edges(data=True)
            ]
            for id in component_cyls:
                # if cyl_to_drip_node[id]:
                # log.warning(f"Cylinder with id {id} has two or more identified to drip nodes: {cyl_to_drip_node[id]} and {drip_node}")
                cyl_to_drip_node[id].append(drip_node)

            # drip_components.append((drip_node, g_drip.subgraph(component_nodes).copy()))

        for cyl in self.cylinders:
            if cyl_to_drip_node[cyl.cyl_id]:
                cyl.drip_node = cyl_to_drip_node[cyl.cyl_id][-1]

        log.info(
            f"{self.file_name} found to have {len(drip_components)} drip components"
        )

        self.stem_flow_component = stem_flow_component
        # self.drip_flow_components = drip_components
        self.drip_graph = g_drip
        self.divide_nodes = divide_nodes
        self.drip_nodes = drip_nodes
        self.cyl_to_drip_node = cyl_to_drip_node

    # @profile
    def calculate_flows(self, plane: str = "XY"):
        """uses subgraphs from FindFlowComponents to aggregate flow characteristics"""
        cyls = self.cylinders
        flow_chars = []

        log.info(
                "Begining sum flows"
        )
        cyls = self.cylinders
        np_flow_chars = [None]*(len(self.drip_nodes) +1)    

        def numpy_flow_chars(lambda_filter:function, drip_cyl, index:int):
            arr= np.array([ 
                    np.array([
                                1,
                                np.float64(cyl.projected_data[plane]["area"]),
                                cyl.surface_area,cyl.angle,cyl.volume,cyl.sa_to_vol 
                            ])      
                    for cyl in cyls if lambda_filter(cyl) ]
                    )
            flow = np.sum(arr, axis = 0)
            np_flow_chars[index] = Flow(flow[0],
                                    flow[1], 
                                    flow[2], 
                                    flow[3], 
                                    flow[4], 
                                    flow[5], 
                                    drip_cyl.cyl_id, 
                                    (drip_cyl.x[0], drip_cyl.y[0], drip_cyl.z[0]))

        numpy_flow_chars(lambda_filter=lambda x: x.is_stem, drip_cyl=self.cylinders[0], index =0 )
        
        flow_chars.append(np_flow_chars[0])
        # log.info(f"summed stem edges {flow_chars}")
        for idx, drip_node in enumerate(self.drip_nodes):
            # drip_cyl = [cyl for cyl in cyls if cyl.cyl_id == drip_node]
            cyl_before_drip = [cyl for cyl in cyls if cyl.cyl_id == drip_node]
            if len(cyl_before_drip) > 1:
                log.warning(f"Error: More that 1 cyl with id {drip_node} found")
            drip_node_loc = (
                    cyl_before_drip[0].x[0],
                    cyl_before_drip[0].y[0],
                    cyl_before_drip[0].z[0],
                )
            drip_cylinders = [
                cyl
                for cyl in cyls
                if cyl.drip_node == drip_node and cyl.cyl_id != drip_node
            ]
            numpy_flow_chars(lambda_filter=lambda x: x.drip_node == drip_node, drip_cyl=self.cylinders[0], index =0 )
        
            # self.cylinders[cylinders.is_stem]


            # edge_attributes = {
            #     (u, v): {
            #         "dripNode": drip_node,
            #         "flowType": "drip",
            #         "flowID": (idx + 1),
            #     }
            #     for u, v, _ in edges
            # }


            # [bool_[cyl for cyl in drip_cylinders]]
            flow_chars.append(

                Flow(
                    **{
                        "num_cylinders": len(drip_cylinders),
                        "projected_area": sum(
                            [
                                np.float64(cyl.projected_data[plane]["area"])
                                for cyl in drip_cylinders
                            ]
                        ),
                        "surface_area": sum(
                            [cyl.surface_area for cyl in drip_cylinders]
                        ),
                        "angle_sum": sum([cyl.angle for cyl in drip_cylinders]),
                        "volume": sum([cyl.volume for cyl in drip_cylinders]),
                        "sa_to_vol": sum([cyl.sa_to_vol for cyl in drip_cylinders]),
                        "drip_node_id": drip_node,
                        "drip_node_loc": drip_node_loc,
                    }
                )
            )
            # nx.set_edge_attributes(g, edge_attributes, "dripNode")
            # log.info(f"summed drip edges in component {idx}")
        self.flows = flow_chars

    def identify_stem_paths(self, axis: str):
        # Special case tree attributes
        root = 0
        self.stem_paths = nx.all_simple_paths(self.g, 0, root)

    def find_furthest(node, nodes):
        """Finds the furthest node in the given list from the provided n"""
        nodes = np.asarray(nodes)
        dist_2 = np.sum((nodes - node) ** 2, axis=1)
        return np.argmax(dist_2)

    def find_trunk_lean(self):
        """Draws a straight line from base to tip of the trunk
        Finds the angle of that line from the XZ plane"""
        trunk_beg = lam_filter(self.cylinders, a_lambda=lambda: branch_order == 0)
        trunk_points = [(cyl.x[0], cyl.y[0], cyl.z[0]) for cyl in trunk_beg]
        root = trunk_points[0]
        furthest_afeild, _ = furthest_point(root, trunk_points)
        dx, dy, dz = ((a - b) for a, b in zip(furthest_afeild, root))
        run = np.sqrt(dx**2 + dy**2)
        angle = (
            np.arctan(dz / np.sqrt(dx**2 + dy**2))
            if run > 0
            else np.arctan(0)
        )
        self.trunk_lean = angle
        return angle

    def find_overlap_by_percentile(
        self,
        plane: str = "XY",
        percentiles: list[int] = [25, 50, 75],
        metric: str = "z",
    ):
        percentiles.sort()
        if not self.projections[plane]:
            self.project_cylinders(plane)
        # if eval(metric) not in vars(Cylinder):
        #     log.info(f"Provided metric invalid: {eval(metric)} is not a property of Cylinder")
        non_trunk_polys, _ = lam_filter(self.cylinders, lambda: branch_order != 0)
        cyl_metric = [cyl.z[0] for cyl in non_trunk_polys]
        percentiles_by_metric = np.percentile(cyl_metric, percentiles)
        prev_perc = 0
        poly_matrix = []
        for perc in percentiles_by_metric:
            poly_matrix.append(
                [
                    cyl.projected_data[plane]["polygon"]
                    for cyl in non_trunk_polys
                    if cyl.z[0] <= perc and cyl.z[0] > prev_perc
                ]
            )
            prev_perc = perc
        overlaps = get_projected_overlap(poly_matrix, percentiles)

        return overlaps

    def statistics(self, plane: str = "XY", file_ext: str = ""):
        if not self.pSV:
            self.project_cylinders(plane)
        if not self.hull:
            self.watershed_boundary(component=self.graph)
        if not self.stem_hull:
            if not self.stem_flow_component:
                self.find_flow_components()
                self.calculate_flows(plane=plane)
            self.watershed_boundary(
                component=self.stem_flow_component, plane=plane, stem=True
                # ,save = True, draw=True
                ,file_ext = file_ext + '_stem_hull' 
            )
        dbh = self.get_dbh()

        #     endNodePoly = [self._pSV[n-1] for n in g.nodes if g.degree(n)==1 and n!= 0]
        #     centroids = [x.point_on_surface() for x in endNodePoly]
        #     tot_hull, edge_points = concave_hull(centroids,2.2)
        #     totHullGeo =geo.GeoSeries(tot_hull)
        canopy_cover = self.hull.area
        canopy_boundary = self.hull.boundary.length

        canopy_cover_stem = self.stem_hull.area
        canopy_boundary_stem = self.stem_hull.boundary.length

        log.info("Found hull alpha shape stats")

        # calculate projected areas and (therefore) overlaps
        # Unary union gives a single contiguous polygon when fed many overlapping cylinders
        # The area of the union thus differs from the sum of the areas of its components
        # in that the former counts overlaps only once
        # polys = [poly_dict['polygon'] for poly_dict in self.pSV]
        polys=self.pSV
        tot_poly = unary_union(polys)
        projected_union_area = tot_poly.area
        sum_projected_area = np.sum([poly.area for poly in self.pSV])
        stem_polys = [
            cyl.projected_data[plane]["polygon"]
            for cyl in self.cylinders
            if cyl.is_stem
        ]
        union_poly_stem = unary_union(stem_polys)
        projected_union_area_stem = union_poly_stem.area
        sum_projected_area_stem = np.sum([poly.area for poly in stem_polys])

        log.info("found projected areas")

        # this could techically be conmbined with the above by adding a percentile of 100
        #
        overlap_dict = self.find_overlap_by_percentile(percentiles=[25, 50, 75])

        tot_poly = unary_union(self.pSV)
        projected_area_w_o_overlap = tot_poly.area
        projected_area_w_overlap = np.sum([poly.area for poly in self.pSV])

        min_x = self.extent["min"][0]
        max_x = self.extent["max"][0]
        min_y = self.extent["min"][1]
        max_y = self.extent["max"][1]
        min_z = self.extent["min"][2]
        max_z = self.extent["max"][2]

        if self.flows:
            stem_flow = [flow for flow in self.flows if flow.drip_node_id == 0][0]
        else:
            stem_flow = 0

        total_surface_area = np.sum([cyl.surface_area for cyl in self.cylinders])
        total_volume = np.sum([cyl.volume for cyl in self.cylinders])
        max_bo = np.max([cyl.branch_order for cyl in self.cylinders])

        order_zero_cyls, _ = lam_filter(self.cylinders, lambda: branch_order == 0)
        order_one_cyls, _ = lam_filter(self.cylinders, lambda: branch_order == 1)
        order_two_cyls, _ = lam_filter(self.cylinders, lambda: branch_order == 2)
        order_three_cyls, _ = lam_filter(self.cylinders, lambda: branch_order == 3)
        order_four_cyls, _ = lam_filter(self.cylinders, lambda: branch_order == 4)

        statistics = {
            "total_psa": projected_union_area,
            "psa_w_overlap": sum_projected_area,
            "stem_psa": projected_union_area_stem,
            "stem_psa_w_overlap": sum_projected_area_stem,
            "tot_surface_area": total_surface_area,
            "stem_surface_area": stem_flow.surface_area,
            "tot_hull_area": canopy_cover,
            "tot_hull_boundary": canopy_boundary,
            "stem_hull_area": canopy_cover_stem,
            "stem_hull_boundary": canopy_boundary_stem,

            "num_drip_points": len(self.drip_nodes),

            "max_bo": max_bo,
            "topQuarterTotPsa": overlap_dict[75]["sum_area"],
            "topHalfTotPsa": overlap_dict[50]["sum_area"],
            "topThreeQuarterTotPsa": overlap_dict[25]["sum_area"],
            "TotalShade": sum_projected_area - projected_union_area,
            "top_quarter_shade": overlap_dict[75]["overlap_with_previous"],
            "top_half_shade": overlap_dict[50]["overlap_with_previous"],
            "top_three_quarter_shade": overlap_dict[25]["overlap_with_previous"],
            "DBH": dbh,
            "volume": total_volume,
            "X_max": max_x,
            "Y_max": max_y,
            "Z_max": max_z,
            "X_min": min_x,
            "Y_min": min_y,
            "Z_min": min_z,
            "Order_zero_angle_avg": np.average([cyl.angle for cyl in order_zero_cyls]),
            "Order_zero_angle_std": np.std([cyl.angle for cyl in order_zero_cyls]),
            "Order_one_angle_avg": np.average([cyl.angle for cyl in order_one_cyls]),
            "Order_one_angle_std": np.std([cyl.angle for cyl in order_one_cyls]),
            "Order_two_angle_avg": np.average([cyl.angle for cyl in order_two_cyls]),
            "Order_two_angle_std": np.std([cyl.angle for cyl in order_two_cyls]),
            "Order_three_angle_avg": np.average(
                [cyl.angle for cyl in order_three_cyls]
            ),
            "Order_three_angle_std": np.std([cyl.angle for cyl in order_three_cyls]),
            "order_gr_four_angle_avg": np.average(
                [cyl.angle for cyl in order_four_cyls]
            ),
            "order_gr_four_angle_std": np.std([cyl.angle for cyl in order_four_cyls]),
            "file_name": self.file_name + plane,
        }
        stat_file = save_file(
            self.file_name.replace(".csv", f"_{file_ext}"),
            out_file=statistics,
            subdir="statistics",
            method="statistics",
            overwrite=True
        )

        return stat_file
    
    def generate_flow_file(self, file_ext):
        flow_dicts = [flow.__dict__ for flow in self.flows]
        flow_file = save_file(
            self.file_name.replace(".csv", f"_flows_{file_ext}"),
            out_file=flow_dicts,
            subdir="flows",
            method="flows",
            overwrite=True
        )
        return flow_file

    def describe(self, metric: str, a_lambda: function = lambda: True, **args) -> dict:
        """
        Takes in a metric of the tree and a filter and returns summary stats
        for the given input metric
        Potentially outputs a histogram if requested?
        """
        return {True}

    def compare(
        metric: str,
        a_lambda: function = lambda: True,
        b_lambda: function = lambda: True,
    ):
        """
        Compares the two lists of metrics defined by the two input functions
        """

    def get_drip_points(
        self,
        # metric :str = 'projected_area',
        percentile: int = 30,
        **args,
    ):
        """
        Returns the locations of the drip points.
        Drip points are identified as the subset of drip nodes
        with a value for the provided metric in the percentile given
        """
        if not self.flows:
            log.warning("flows are not defined, run identify and calculate flows")
            return
        # excluding trunk node
        drip_nodes = [
            (f.drip_node_id, f.projected_area, f.drip_node_loc)
            for f in self.flows
            if f.drip_node_id != 0
        ]

        distinct_drip_node_ids = {node[0] for node in drip_nodes}

        flow_by_area_dict = {}
        for node in distinct_drip_node_ids:
            area = np.sum([flow[1] for flow in drip_nodes if flow[0] == node])
            flow_by_area_dict[node] = area

        point_cutoff = np.percentile(
            [area for _, area in flow_by_area_dict.items()], percentile
        )
        drip_points = [
            node for node, area in flow_by_area_dict.items() if area > point_cutoff
        ]

        drip_point_locs = [
            [node[2][0], node[2][1], node[2][2]]
            for node in drip_nodes
            if node[0] in drip_points
        ]
        return drip_point_locs

    def drip_map(
        self, a_lambda: function = lambda: True, scale: int = 1,
        interpolate:bool = False,
        file_ext:str = '', **args
    ) -> None:
        """
        Returns a plot showing the locations of the drip points, subject
        to input params

        a_lambda: function to filter drip points displayed (e.g. those with projected area>10m^2 )
        scale: how large of a boundary to draw around drip points
        """
        drip_point_locs = self.get_drip_points()
        drip_point_locs_x = [pt[0] * scale for pt in drip_point_locs]
        drip_point_locs_y = [pt[1] * scale for pt in drip_point_locs]
        # min_xy = np.min(mins)
        # max_xy = np.max(maxs)
        # x_mesh, y_mesh = np.meshgrid(
        #     np.arange(min_xy, max_xy, 0.05), np.arange(min_xy, max_xy, 0.05)
        # )
        # try:
        # fig, ax = plt.subplots()
        # from geopandas import GeoSeries 
        # except Exception as e:
        #     log.warning(f"drip map not drawn: {e}")
        #     return
        
        # if interpolate:
        #     drip_point_locs_xy = [[pt[0] * scale, pt[1] * scale] for pt in drip_point_locs]

        #     math.floor(np.min(drip_point_locs_x))

        #     mins = self.extent["min"]
        #     maxs = self.extent["max"]
        #     extents = [mins[0], maxs[0], mins[1], maxs[1]]
        #     min_xy = np.min(
        #         [
        #             math.floor(np.min(drip_point_locs_x)),
        #             math.floor(np.min(drip_point_locs_y)),
        #         ]
        #     )
        #     max_xy = np.max(
        #         [math.ceil(np.max(drip_point_locs_x)), math.ceil(np.max(drip_point_locs_y))]
        #     )
        #     x_mesh, y_mesh = np.meshgrid(
        #         np.arange(min_xy, max_xy, 0.005), np.arange(min_xy, max_xy, 0.005)
        #     )

        #     def dist_to_drip(a, b):
        #         distances = distance.cdist([[a, b]], drip_point_locs_xy)
        #         min_dist = np.min(distances)
        #         return math.log(1 / min_dist)

        #     distance_matrix = np.zeros((x_mesh.shape[0], x_mesh.shape[0]))

        #     for a in range(x_mesh.shape[0]):
        #         for b in range(x_mesh.shape[0]):
        #             distance_matrix[a][b] = dist_to_drip(x_mesh[b][a], y_mesh[b][a])


        #     ax.contourf(
        #         y_mesh,
        #         x_mesh,
        #         distance_matrix,
        #         levels=15,
        #         max=0.5,
        #         cmap=plt.cm.Blues,
        #         extend="neither",
        #         extent=extents,
        #     )
            
        # ax.scatter(drip_point_locs_x, drip_point_locs_y)

        # filtered_cyls, _ = lam_filter(self.cylinders, a_lambda, return_all=False)
        # polys = [cyl.projected_data["XY"]["polygon"] for cyl in filtered_cyls]
        # # breakpoint()
        # save_dir = "/".join([str(output_dir), 'draw_drip_map', f"{file_ext}"])#.replace("/", "\\")
        # plt.show()
        # plt.savefig(save_dir, dpi = 3000)
        # if len(polys) > 0:
        #     geoPolys = GeoSeries(polys)
        #     geoPolys.plot(ax=ax)
        # else:
        #     log.warning(
        #         "Drip Map: No cylinders returned for lambda function: {a_lambda}"
        #     )
        # plt.show()
