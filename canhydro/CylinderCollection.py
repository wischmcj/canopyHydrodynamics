"""These constitute our way of 'copying' the QSMs"""


from __future__ import annotations

import copy
import sys
from collections import defaultdict

import networkx as nx
import numpy as np
# import settings
from memory_profiler import LogFile, profile

from canhydro.Cylinder import Cylinder
from canhydro.DataClasses import Flow
from canhydro.global_vars import log, qsm_cols
from canhydro.Plotter import draw_cyls
from canhydro.utils import (concave_hull, intermitent_log, lam_filter,
                            unary_union)

sys.stdout = LogFile()

NAME = "CylinderCollection"


# By inheriting the Model class, lambda cyl : cyl.branch_order = br CC gains managed functionality- like lambda searching
class CylinderCollection:
    cylinders = defaultdict(list)

    # initialize our object level variables for cylinder objects
    def __init__(self) -> None:
        self.file = ""

        # self.collection = CylinderList()

        # Aggregate values from file
        self.surface_area = np.nan
        self.filename = ""
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
            "total_psa": -1,
            "tot_hull_area": -1,
            "stem_flow_hull_area": -1,
            "stem_psa": -1,
            "flowStats": -1,
            "dbh": -1,
            "tot_surface_area": -1,
            "stem_surface_area": -1,
        }

        # Projection Attrs
        self.union_poly = None
        self.stem_path_lengths = []
        self.hull = np.nan
        self.stem_hull = np.nan

        # Special case tree attributes
        self.stem_paths = [[]]  # Cyl collection?
        self.trunk = []  # Collection of cylinders? id list?

        # Graph and Attributes
        self.graph = None
        self.object_graph = None
        self.min_graph = None
        self.flows = None
        self.drip_points = {"x": np.nan, "y": np.nan, "z": np.nan, "flow_id": np.nan}
        self.flow_to_drip = {
            0: 1
        }  # A dictionary of flow ids with values equal to their drip node ids
        self.trunk_nodes = []
        self.drip_loc = np.nan
        self.stemFlowComponent = None
        self.dripFlowComponents = None
        # Calculations using graph results
        self.stemTotal = {
            "attributes": {"cyls": 0, "len": 0, "sa": 0, "pa": 0, "as": 0},
            "loc": {"x": np.nan, "y": np.nan, "z": np.nan},
        }
        self.divide_points = []
        self.stemPolys = []
        self.compGraphs = []

    def aggregate_characteristics(self):
        """Calculates the summations, averages etc. of cylinder characterictics
        that might be of interest"""
        return True

    def create_cyl(self, arr: list):
        cols = qsm_cols
        attrs = {k: arr[v] for (k, v) in cols.items()}
        cyl = Cylinder(**attrs)
        cyl.create_from_list(arr, cols)
        return cyl

    def from_csv(self, file, aggregate_cyls=True):
        """Initializes a new Cyl Collection based on the data in a QSM
        with the configured column locations"""
        self.file = file
        self.filename = file.name
        log.info(f"Processing {str(file)}")
        # self.arr = pd.read_csv(file, header=0)
        self.arr = np.genfromtxt(file, delimiter=",", skip_header=True)[0:, :-1]
        cylinders = [self.create_cyl(row) for row in self.arr]
        self.cylinders = cylinders

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

    def project_cylinders(self, plane: str = "XZ"):
        """Projects cylinders onto the specified plane"""
        if plane not in ("XY", "XZ", "YZ"):
            log.info(f"{plane}: invalid value for plane")
        else:
            polys = []
            log.info(f"Projection into {plane} axis begun for file {self.filename}")
            for idx, cyl in enumerate(self.cylinders):
                poly = cyl.get_projection(plane)
                polys.append(poly)
                # print a progress update once every 10 thousand or so cylinders
                intermitent_log(idx, self.no_cylinders, "Cylinder projection: ")
            # Projection Attrs
            self.union_poly = unary_union(polys)
            self.stem_path_lengths = []
            self.pSV = None  # unsure if I want to keep this attr

    def get_collection_data(self):
        cyl_desc = [cyl.__repr__() for cyl in self.cylinders]
        return cyl_desc

    def draw(
        self,
        plane: str = "XZ",
        a_lambda: function = lambda: True,
        highlight: bool = False,
        **args,
    ):
        """Draws cylinders meeting given characteristics onto the specified plane"""
        if plane not in ("XY", "XZ", "YZ"):
            log.info(f"{plane}: invalid value for plane")
        filtered_cyls, matches = lam_filter(
            self.cylinders, a_lambda, return_all=highlight
        )
        to_draw = [cyl.projected_data[plane]["polygon"] for cyl in filtered_cyls]
        log.info(f"{len(to_draw)} cylinders matched criteria")
        self.union_poly = unary_union(to_draw)
        draw_cyls(collection=to_draw, colors=matches, **args)

    def get_dbh(self):
        g = self.graph
        if self.end_nodes:
            return self.end_nodes
        elif len(g.nodes) > 0:
            start_z = self.extent["min"][2]
            higher_than_breast_radius = lam_filter(
                self.cylinders,
                a_lambda=lambda: branch_order == 0 and z[0] <= start_z + 1.3,
            )
            rbh = np.max(higher_than_breast_radius.radius)
            self.treeQualities["dbh"] = 2 * rbh

    def get_end_nodes(self) -> list[int]:
        g = self.graph
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

    def get_trunk_nodes(self) -> list[int]:
        g = self.graph
        if self.trunk_nodes:
            return self.trunk_nodes
        elif len(g.nodes) > 0:
            trunk_cyls, _ = lam_filter(self.cylinders, lambda: branch_order == 0)
            trunk_nodes = [cyl.cyl_id for cyl in trunk_cyls]
            self.trunk_nodes = trunk_nodes
            return trunk_nodes
        else:
            log.warning(
                "Graph not initialized, run <CylinderCollection>.initialize_graph(**args)"
            )
            return list(None)

    def watershed_boundary(self, plane: str = "XZ", stem_only: bool = False):
        # todo -- check data type of returned (do we really need geo pandas?)
        g = self.graph

        # filter for stem only if option selected
        endNodes = [n for n in g.nodes if g.degree(n) == 1 and n != -1]
        endCyls, _ = lam_filter(self.cylinders, lambda: cyl_id in endNodes)
        centroids = [
            cyl.projected_data[plane]["polygon"].point_on_surface() for cyl in endCyls
        ]
        tot_hull, _ = concave_hull(centroids, 2.2)
        # totHullGeo = geo.GeoSeries(tot_hull)
        # canopyCover = totHullGeo.area
        # canopyBoundary = totHullGeo.boundary.length
        self.hull = np.nan
        # self.stem_hull = np.nan

    @profile
    def initialize_graph(self):
        """This function initialized edge attributes as FULL cylinder dicts"""
        gr = nx.Graph()
        for cyl in self.cylinders:
            attr = cyl.__dict__
            child_node = attr["cyl_id"]
            parent_node = attr["parent_id"]
            gr.add_edge(child_node, parent_node, **attr)
        self.graph = gr

    @profile
    def initialize_object_graph(self):
        """This function initialized edge attributes as cylinder objects"""
        gr = nx.Graph()
        for cyl in self.cylinders:
            attr = cyl.__dict__
            child_node = attr["cyl_id"]
            parent_node = attr["parent_id"]
            gr.add_edge(child_node, parent_node, cylinder=cyl)
        self.object_graph = gr

    @profile
    def initialize_minimal_graph(self):
        """This function initialized edge attributes as cylinder objects"""
        gr = nx.Graph()
        for cyl in self.cylinders:
            attr = cyl.weight_dict()
            child_node = attr["cyl_id"]
            parent_node = attr["parent_id"]
            gr.add_edge(child_node, parent_node)
        self.min_graph = gr

    @profile
    def sum_over_graph(self):
        """This function initialized edge attributes as FULL cylinder dicts"""
        gr = self.graph
        cyls = [
            attr["projected_data"]["XZ"]["polygon"].area
            for u, v, attr in gr.edges(data=True)
            if u > 5000
        ]
        projected_area = np.sum(cyls)
        return projected_area

    @profile
    def sum_over_object_graph(self):
        """This function initialized edge attributes as FULL cylinder dicts"""
        gr = self.object_graph
        areas = [
            attr["cylinder"].projected_data["XZ"]["polygon"].area
            for u, v, attr in gr.edges(data=True)
            if u > 5000
        ]
        projected_area = np.sum(areas)
        return projected_area

    @profile
    def sum_over_min_graph(self):
        """This function initialized edge attributes as FULL cylinder dicts"""
        gr = self.min_graph
        cyl_ids = [
            u for u, v, attr in gr.edges(data=True) if u > 5000
        ]  ## a contrived example
        # cyls = lam_filter(self.cylinders, lambda: cyl_id in cyl_ids)
        # areas = [cyl.projected_data['XZ']["polygon"].area for cyl in cyls]
        areas = [
            cyl.projected_data["XZ"]["polygon"].area
            for cyl in self.cylinders
            if cyl.cyl_id in cyl_ids
        ]

        projected_area = np.sum(areas)
        return projected_area

    def find_flow_components(self, inFlowGradeLim=-1 / 6):
        """Finding Stemflow contributing area"""
        g = self.graph
        # identify drip edges in graph
        drip_edges = [
            (u, v)
            for u, v, attr in g.edges(data=True)
            if attr["angle"] < inFlowGradeLim
        ]

        inFlowGraph = copy.deepcopy(g)  # This could probably just be a subgraph...
        inFlowGraph.remove_edges_from(drip_edges)
        log.info(f"{self.filename} found to have {len(drip_edges)} drip edges")

        # separating the stem flow from the drip flows and the drip flows from each other
        G_drip = copy.deepcopy(g)
        root_node = 0
        stem_flow_component = g.subgraph(
            nx.node_connected_component(inFlowGraph, root_node)
        ).copy()  # returns the connected component containing the root
        G_drip.remove_edges_from(stem_flow_component.edges())
        drip_flow_components = nx.connected_components(G_drip)
        component_graphs = [g.subgraph(c).copy() for c in drip_flow_components]
        log.info(
            f"{self.filename} found to have {len(component_graphs)} drip components"
        )
        for cyl in self.cylinders:
            if cyl.cyl_id in stem_flow_component.nodes():
                cyl.is_stem = True
        self.stemFlowComponent = stem_flow_component
        self.dripFlowComponents = component_graphs
        self.get_trunk_nodes()

    def find_flow_components_object(self, inFlowGradeLim=-1 / 6):
        """Finding Stemflow contributing area"""
        g = self.object_graph
        # identify drip edges in graph
        drip_edges = [
            (u, v)
            for u, v, attr in g.edges(data=True)
            if attr["cylinder"].angle < inFlowGradeLim
        ]
        inFlowGraph = copy.deepcopy(g)  # This could probably just be a subgraph...
        inFlowGraph.remove_edges_from(drip_edges)
        log.info(f"{self.filename} found to have {len(drip_edges)} drip edges")

        # separating the stem flow from the drip flows and the drip flows from each other
        G_drip = copy.deepcopy(g)
        root_node = 0
        stem_flow_component = g.subgraph(
            nx.node_connected_component(inFlowGraph, root_node)
        ).copy()  # returns the connected component containing the root
        G_drip.remove_edges_from(stem_flow_component.edges())
        drip_flow_components = nx.connected_components(G_drip)
        component_graphs = [g.subgraph(c).copy() for c in drip_flow_components]
        log.info(
            f"{self.filename} found to have {len(component_graphs)} drip components"
        )

        self.stemFlowComponent = stem_flow_component
        self.dripFlowComponents = component_graphs

    def find_flow_components_minimal(self, inFlowGradeLim=-1 / 6):
        """Finding Stemflow contributing area"""
        g = self._graph
        # identify drip edges in graph
        all_cyls, is_out = lam_filter(
            self.cylinders, lambda: angle < inFlowGradeLim, return_all=True
        )
        out_key = zip(all_cyls, is_out)
        drip_edges = [(cyl.cyl_id, cyl.parent_id) for cyl, is_out in out_key if is_out]

        print(len(g.edges()))
        inFlowGraph = copy.deepcopy(g)  # This could probably just be a subgraph...
        inFlowGraph.remove_edges_from(drip_edges)
        log.info(f"{self.filename} found to have {len(drip_edges)} drip edges")

        # separating the stem flow from the drip flows and the drip flows from each other
        G_drip = copy.deepcopy(g)
        root_node = 0
        stem_flow_component = g.subgraph(
            nx.node_connected_component(inFlowGraph, root_node)
        ).copy()  # returns the connected component containing the root
        G_drip.remove_edges_from(stem_flow_component.edges())
        drip_flow_components = nx.connected_components(G_drip)
        component_graphs = [g.subgraph(c).copy() for c in drip_flow_components]
        print(component_graphs)
        log.info(
            f"{self.filename} found to have {len(component_graphs)} drip components"
        )

        self.stemFlowComponent = stem_flow_component
        self.dripFlowComponents = component_graphs
        # self.flows = [
        #     {
        #         "cyls": [],
        #         "drip_point": id,
        #         "attributes": {"cyls": 0, "len": 0, "sa": 0, "pa": 0, "as": 0},
        #     }
        # ]
        # self.drip_points = {"x": np.nan, "y": np.nan, "z": np.nan, "flow_id": np.nan}
        # self.flow_to_drip = {
        #     0: 1
        # }  # A dictionary of flow ids with values equal to their drip node ids
        # self.trunk_nodes = []
        # self.drip_loc = np.nan

        # # Calculations using graph results
        # self.stemTotal = {
        #     "attributes": {"cyls": 0, "len": 0, "sa": 0, "pa": 0, "as": 0},
        #     "loc": {"x": np.nan, "y": np.nan, "z": np.nan},
        # }
        # self.divide_points = []
        # self.stemPolys = []

    def calculateFlows(self):
        """uses subgraphs from FindFlowComponents to aggregate flow characteristics"""
        stem_flow_component = self.stemFlowComponent
        drip_flow_components = self.dripFlowComponents
        g = self.graph
        edge_attributes = {}
        flow_chars = {}

        for u, v in stem_flow_component.edges():
            edge_attributes[(u, v)] = {"dripNode": 0, "flowType": "stem", "flowID": 0}

        flow_dict = lambda metric, comp: {
            i: {j: {val}} for i, j, val in comp.edges.data(metric)
        }
        stem_edges = len(stem_flow_component.edges())
        flow_chars.append(
            Flow(
                **{
                    "num_cylinders": stem_edges,
                    "projected_area": nx.cost_of_flow(
                        g, flow_dict("projected_area", stem_flow_component)
                    ),
                    "surface_area": nx.cost_of_flow(
                        g, flow_dict("surface_area", stem_flow_component)
                    ),
                    "angle_sum": nx.cost_of_flow(
                        g, flow_dict("angle", stem_flow_component)
                    ),
                    "volume": nx.cost_of_flow(
                        g, flow_dict("volume", stem_flow_component)
                    ),
                    "sa_to_vol": nx.cost_of_flow(
                        g, flow_dict("sa_to_vol", stem_flow_component)
                    ),
                    "drip_node_id": 0,
                }
            )
        )  # stemflow drips to the trunk

        for idx, flow in enumerate(drip_flow_components):
            if self._projection == "XZ":
                heights = self._df.iloc[:, 4].to_numpy()
            elif self._projection == "YZ":
                heights = self._df.iloc[:, 3].to_numpy()
            else:
                self._heights = heights = self._df.iloc[:, 5].to_numpy()
            nodes = flow.nodes()
            num_cyls = len(flow.edges())
            drip_node = nodes[0]

            for node in nodes:
                if heights[node] < heights[drip_node]:
                    drip_node = node

            for u, v in flow.edges():
                edge_attributes[(u, v)] = {
                    "dripNode": drip_node,
                    "flowType": "drip",
                    "flowID": (idx + 1),
                }

            flow_chars.append(
                Flow(
                    **{
                        "num_cylinders": num_cyls,
                        "projected_area": nx.cost_of_flow(
                            g, flow_dict("projected_area", flow)
                        ),
                        "surface_area": nx.cost_of_flow(
                            g, flow_dict("surface_area", flow)
                        ),
                        "angle_sum": nx.cost_of_flow(g, flow_dict("angle", flow)),
                        "volume": nx.cost_of_flow(g, flow_dict("volume", flow)),
                        "sa_to_vol": nx.cost_of_flow(g, flow_dict("sa_to_vol", flow)),
                        "drip_node_id": drip_node,
                    }
                )
            )

        self.flows = flow_chars

        nx.set_edge_attributes(g, edge_attributes, "dripNode")

    def contractedNodes(node_list, next_node, G):
        # finding and removing the listed nodes to find their edges
        incident_edges = [
            (u, v) for u, v in G.edges() if u in node_list or v in node_list
        ]
        non_incident_edges = [e for e in G.edges() if e not in incident_edges]
        contracted = copy.deepcopy(G)
        contracted.remove_nodes_from(self.get_trunk_nodes())

        neighbors = [
            node
            for node in contracted.nodes()
            if (G.degree(node) != contracted.degree(node))
        ]

        contracted.add_node(next_node)
        edges_to_trunk = [(u, next_node) for u in neighbors]
        contracted.add_edges_from(edges_to_trunk)

        return contracted, neighbors

    def find_trunk_distance(self):
        trunk_nodes = self.get_trunk_nodes()

        trunk_contraction, titans = self.contractedNodes(trunk_nodes, 0, self._graph)
        trunk_paths = nx.shortest_path(trunk_contraction, target=0)
        dists = {node: len(path) - 1 for node, path in trunk_paths.items()}
        return dists

    def identify_stem_paths(self, axis: str):
        # Special case tree attributes
        root = 0
        all_simple_paths(self.g, 0, target, cutoff=None)
        self.stem_paths = [[]]  # Cyl collection?
        self.trunk = []  # Collection of cylinders? id list?

    def find_trunk_lean(self):
        # to populate with x,y,z mins and maxs
        self.aggregate_angle = np.nan
        return True

    def statistics():
        statistics = {}
        return statistics

    # def read_csv(self, df=pd.DataFrame(), polys=[], projection="XY"):
    #     # Columns [ID?,ParentID?,x1,y1,z1,x2,y2,z2,radius,?,?,length,? ,? ,? ,? ,? ,? ,? ,BO]
    #     # Colnums [1  ,2        ,3 , 4, 5, 6, 7, 8,9    ,10,11,12   ,13,14,15,16,17,18,19,20]
    #     # x = x2-x1, y =y2-y1, z=z2-z1
    #     # number of cylinders = cnt values for radius
    #     # Theta  = angle made buy the cylinder axis
    #     self.projection = projection
    #     if "partial" not in self.filename:
    #         self.df_full = pd.read_csv(self.filename, header=0)
    #     else:
    #         self.df_full = df
    #     self.maxBO = np.max(self.BO)
    #     # self.df = self.df_full
    #     self.df = self.df_full  # .iloc[:130,:]
    #     # columns 3 and 6 represent our x values
    #     if projection == "XZ":
    #         self.x = np.transpose(self.df.iloc[:, [3, 6]].to_numpy())
    #         self.y = np.transpose(self.df.iloc[:, [5, 8]].to_numpy())
    #         self.z = np.transpose(self.df.iloc[:, [4, 7]].to_numpy())
    #     elif projection == "YZ":
    #         self.x = np.transpose(self.df.iloc[:, [5, 8]].to_numpy())
    #         self.y = np.transpose(self.df.iloc[:, [3, 6]].to_numpy())
    #         self.z = np.transpose(self.df.iloc[:, [4, 7]].to_numpy())
    #     else:  # 'XY'
    #         self.x = np.transpose(self.df.iloc[:, [3, 6]].to_numpy())
    #         self.y = np.transpose(self.df.iloc[:, [4, 7]].to_numpy())
    #         self.z = np.transpose(self.df.iloc[:, [5, 8]].to_numpy())
    #     # for side view
    #     self.cylID = self.df.iloc[:, 1].to_numpy()
    #     self.pID = self.df.iloc[:, 2].to_numpy()
    #     self.radius = self.df.iloc[:, 9].to_numpy()
    #     self.no_cylinders = self.radius.size
    #     self.cLength = self.df.iloc[:, 12].to_numpy()
    #     self.BO = self.df.iloc[:, 20].to_numpy()
    #     self.branchID = self.df.iloc[:, 4].to_numpy()
    #     self.maxBO = np.max(self.BO)
    #     self.bID = self.df.iloc[:, 24].to_numpy()
    #     if projection == "XZ":
    #         self.dx = self.df.iloc[:, 6].to_numpy() - self.df.iloc[:, 3].to_numpy()
    #         self.dy = self.df.iloc[:, 8].to_numpy() - self.df.iloc[:, 5].to_numpy()
    #         self.dz = self.df.iloc[:, 7].to_numpy() - self.df.iloc[:, 4].to_numpy()
    #     elif projection == "YZ":
    #         self.dx = self.df.iloc[:, 8].to_numpy() - self.df.iloc[:, 5].to_numpy()
    #         self.dy = self.df.iloc[:, 6].to_numpy() - self.df.iloc[:, 3].to_numpy()
    #         self.dz = self.df.iloc[:, 7].to_numpy() - self.df.iloc[:, 4].to_numpy()
    #     else:  # 'XY'
    #         self.dx = self.df.iloc[:, 6].to_numpy() - self.df.iloc[:, 3].to_numpy()
    #         self.dy = self.df.iloc[:, 7].to_numpy() - self.df.iloc[:, 4].to_numpy()
    #         self.dz = self.df.iloc[:, 8].to_numpy() - self.df.iloc[:, 5].to_numpy()
    #     if "partial" in self.filename:
    #         self.pSV = polys
    #     self.theta = np.arctan(self.dz / np.sqrt(self.dx**2 + self.dy**2))
    #     self.output_dir = "".join([DIR, "output/"])
    #     log.info(self.filename + " initialized")

    # def find_flows():
    #     #Replace trunk with single node

    #     #Remove edges with out flow
    #     #Find the connected component with root in it, this is stem flow

    #     #Remove all stem flow from orig graph, these are drip paths
    #         #Find node of minimal height in each component, this is the drip point

    # def network_simplex():
    #     #can be used to calculate flows on graphs with demands
    #     #we could set a demand of generates X volume of flow
