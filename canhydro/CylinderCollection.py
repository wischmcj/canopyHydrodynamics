"""These constitute our way of 'copying' the QSMs"""


from __future__ import annotations

import copy
import sys

import networkx as nx
import numpy as np
from shapely.geometry import Polygon

from memory_profiler import LogFile, profile

from canhydro.Cylinder import Cylinder
from canhydro.DataClasses import Flow
from canhydro.global_vars import log, qsm_cols
from canhydro.Plotter import draw_cyls
from canhydro.utils import (concave_hull, intermitent_log, lam_filter,union, furthest_point)

sys.stdout = LogFile()

NAME = "CylinderCollection"


# By inheriting the Model class, lambda cyl : cyl.branch_order = br CC gains managed functionality- like lambda searching
class CylinderCollection:
    cylinders = dict

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
        self.hull = None
        self.stem_hull = None

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
        self.drip_loc = None
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
        self.trunk_lean = None

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
            self.union_poly = union(polys)
            self.stem_path_lengths = []
            self.pSV = polys  # unsure if I want to keep this attr

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
        self.union_poly = union(to_draw)
        draw_cyls(collection=to_draw, colors=matches, **args)

    def get_dbh(self):
        g = self.graph
        if self.end_nodes:
            return self.end_nodes
        elif len(g.nodes) > 0:
            start_z = self.extent["min"][2]
            higher_than_breast_radius = lam_filter(
                self.cylinders,
                a_lambda=lambda: branch_order == 0 and z[0] >= start_z + 1.3,
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

    def watershed_boundary(self, plane: str = "XZ", component = self.graph, curvature_alpha: float() = 2)-> None:
        """Generates tightly fit concave_hull (alpha shape) around the passed component"""
        """Alpha determines the tightness of the fit of the shape. Too low an alpha results in multiple"""
        if cylinders[1].projected_data.get(plane, False)
            # want to get a ring of points around the component
            # however, the root node also has degree one, so we exclude it with n!= -1
            endNodes = [n for n in component.nodes if component.degree(n) == 1 and n != -1]
            endCyls, _ = lam_filter(self.cylinders, lambda: cyl_id in endNodes)
            boundary_centroids = [
                cyl.projected_data[plane]["polygon"].point_on_surface() for cyl in endCyls
            ]
            hull, _ = concave_hull(boundary_centroids, curvature_alpha)

            if stem_only:
                self.stem_hull = hull
            else:    
                self.hull = hull
        else:
            log.info("")
        
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
    
    def contracted_nodes(self, node_list, next_node, G):
        """Returns a graph made from removing nodes from G,
          and connecting the removed nodes' neighbors via an edge """
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
    
    def find_trunk_distance(self):
        trunk_nodes = self.get_trunk_nodes()

        trunk_contraction, titans = self.contracted_nodes(self.trunk_nodes, 0, self.graph)
        trunk_paths = nx.shortest_path(trunk_contraction, target=0)
        dists = {node: len(path) - 1 for node, path in trunk_paths.items()}
        return dists
    
    def find_flow_components(self, inFlowGradeLim:float= -1/ 6,dist:int = 2):
        """Finding Stemflow contributing area"""
        g = self.graph
        # identify drip edges in graph
        trunk_distance= self.find_trunk_distance()
        drip_edges= [(u,v) for u,v,attr in g.edges(data=True) if attr['angle']<inFlowGradeLim and trunk_distance.get(min(u,v),0)>dist]
        # drip_edges = [
        #     (u, v)
        #     for u, v, attr in g.edges(data=True)
        #     if attr["angle"] < inFlowGradeLim
        # ]

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
        g = self.min_graph
        # identify drip edges in graph
        all_cyls, is_out = lam_filter(
            self.cylinders, lambda: angle < -1/6, return_all=True
        )
        out_key = zip(all_cyls, is_out)
        drip_edges = [(cyl.cyl_id, cyl.parent_id) for cyl, is_out in out_key if is_out]

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

    def flowCost(graph, metric):
        #algo does not play well with floats 
        sum(attr[metric]*10000 for u, v, attr in graph.edges(data=True))/10000
    
    def calculate_flows(self):
        """uses subgraphs from FindFlowComponents to aggregate flow characteristics"""
        stem_flow_component = self.stemFlowComponent
        drip_flow_components = self.dripFlowComponents
        g = self.graph
        edge_attributes = {}
        flow_chars = {}
        G_return = copy.deepcopy(self.graph)

        for u, v in stem_flow_component.edges():
            edge_attributes[(u, v)] = {"dripNode": 0, "flowType": "stem", "flowID": 0}

        flow_dict = lambda comp,metric: {
            i: {j: {val}} for i, j, val in comp.edges.data(metric)
        }
        stem_edges = len(stem_flow_component.edges())
        flow_chars.append(
            Flow(
                **{
                    "num_cylinders": stem_edges,
                    "projected_area": self.flowCost(stem_flow_component, "projected_area"),
                    
                    "surface_area": self.flowCost(stem_flow_component, "surface_area"),
                    
                    "angle_sum":  self.flowCost(stem_flow_component, "angle"),
                    "volume":  self.flowCost(stem_flow_component, "volume"),
                    "sa_to_vol": self.flowCost(stem_flow_component, "sa_to_vol"),
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
                    "projected_area": self.flowCost(flow, "projected_area"),
                    
                    "surface_area": self.flowCost(flow, "surface_area"),
                    
                    "angle_sum":  self.flowCost(flow, "angle"),
                    "volume":  self.flowCost(flow, "volume"),
                    "sa_to_vol": self.flowCost(flow, "sa_to_vol"),
                    "drip_node_id": drip_node,
                }
                )
            )

        self.flows = flow_chars

        nx.set_edge_attributes(g, edge_attributes, "dripNode")

    def identify_stem_paths(self, axis: str):
        # Special case tree attributes
        root = 0
        self.stem_paths = nx.all_simple_paths(self.g, 0, root)

    def find_furthest(node, nodes):
        nodes = np.asarray(nodes)
        dist_2 = np.sum((nodes - node)**2, axis=1)
        return np.argmax(dist_2)

    def find_trunk_lean(self):
        start_z = self.extent["min"][2]
        trunk_beg = lam_filter(
            self.cylinders,
            a_lambda=lambda: branch_order == 0 and z[0] <= start_z + 1.3,
        )
        trunk_points = [(cyl.x[0],cyl.y[0],cyl.z[0]) for cyl in trunk_beg]
        root = trunk_points[0]
        furthest_afeild = furthest_point(root,trunk_points)
        dx, dy, dz = ((a-b) for a,b in zip(furthest_afeild, root))
        run = np.sqrt(self.dx**2 + self.dy**2)
        angle = (
            np.arctan(self.dz / np.sqrt(self.dx**2 + self.dy**2))
            if run > 0
            else np.arctan(0)
        )
        self.trunk_lean = angle
        return angle



    def statistics():

        if not self.hull:
            self.watershed_boundary()
        if not self.stem_hull: 
            self.watershed_boundary()
        if not self.pSV:


        endNodePoly = [self._pSV[n-1] for n in g.nodes if g.degree(n)==1 and n!= 0]
        centroids = [x.point_on_surface() for x in endNodePoly]
        tot_hull, edge_points = concave_hull(centroids,2.2)
        totHullGeo =geo.GeoSeries(tot_hull)
        canopyCover = totHullGeo.area
        canopyBoundary = totHullGeo.boundary.length

        print('found hull stats')

        #calculate overlaps
        totPoly = union(geo.GeoSeries(self._pSV))
        projected_union_area = totPoly.area

        area = 0
        for poly in self._pSV:
            area+=poly.area
        print('found total area')

        canopy_heights = pd.DataFrame(self._z[0])[self._BO>0]
        canopy_percentiles = canopy_heights.describe()
        topQuarterPolys     =  geo.GeoSeries((pd.DataFrame(self._pSV)[self._z[0]>=canopy_percentiles.iloc[6][0]])[0])
        topHalfPolys        =  geo.GeoSeries((pd.DataFrame(self._pSV)[self._z[0]>=canopy_percentiles.iloc[5][0]])[0])
        topThreeQuarterPolys= geo.GeoSeries((pd.DataFrame(self._pSV)[self._z[0]>=canopy_percentiles.iloc[4][0]])[0])

        topQuarterArea     =0
        topHalfArea        =0
        topThreeQuarterArea=0

        print(1)
        for area in topQuarterPolys.area:
            topQuarterArea      += area 
        print(2)

        for area in topHalfPolys.area:
            topHalfArea      += area 
        print(3)
        i=0
        for area in topThreeQuarterPolys.area:
            topThreeQuarterArea      += area 
        print(4)

        topQuarterAggArea =    union(geo.GeoSeries(topQuarterPolys)).area    
        topHalfAggArea     =   union(geo.GeoSeries(topHalfPolys)).area   
        topThreeQuarterAggArea=union(geo.GeoSeries(topThreeQuarterPolys)).area

        totPoly = union(geo.GeoSeries(self._pSV))
        projected_union_area = totPoly.area
        area = 0
        for poly in self._pSV:
            area+=poly.area
        print(5)

        if self._projection == 'XZ':
            X_max = np.max(self._y)
            Y_max = np.max(self._x)
            Z_max = np.max(self._z) 
        elif self._projection == 'YZ':
            X_max = np.max(self._z)
            Y_max = np.max(self._x)
            Z_max = np.max(self._y)   
        else: # 'XY'
            X_max = np.max(self._x)
            Y_max = np.max(self._y)
            Z_max = np.max(self._z)

        nodes = self._graph.nodes()
        filt_nodes = list(nodes) 
        StemPolys=[]
        sumStemSurfaceArea = 0
        sumStemProjectedArea = 0
        for idx,n in enumerate(nodes):
            if idx+1 <= len(filt_nodes): parent = self._pID[n-1]+1
            if n in filt_nodes and n>0: 
                if g.edges[parent,n]['flowType'] =='stem': 
                    StemPolys.append(self._pSV[n-1])
                    sumStemProjectedArea +=self._pSV[n-1].area
                    surface_area = 2* np.pi*self._radius[n-1]*(self._radius[n-1] +self._cLength[n-1])-(2*np.pi*self._radius[n-1]*self._radius[n-1])
                    sumStemSurfaceArea+=surface_area
                    
        totPoly = union(geo.GeoSeries(StemPolys))
        stem_projected_union_area = totPoly.area
        
        WholeStemPoly = union(geo.GeoSeries(self._stemPolys))
        projected_stem_area = WholeStemPoly.area
        
        DivideCentroids = [x.point_on_surface() for x in self._dividePolys]
        hull, edge_points = concave_hull(DivideCentroids, 2.2)

        statistics  = {'total_psa':projected_union_area,
                                            'psa_w_overlap':area,
                                            'stem_psa':stem_projected_union_area,
                                            'stem_psa_w_overlap': sumStemProjectedArea,
                                            'tot_surface_area':np.sum(surface_area),
                                            'stem_surface_area':sumStemSurfaceArea,
                                            'tot_hull_area':canopyCover, 
                                            'tot_hull_boundary':canopyBoundary,
                                            'stem_hull_area':hull.area, 
                                            'stem_hull_boundary':totHullGeo.boundary.length,
                                            'max_bo': np.max(self._BO),
                                            'topQuarterTotPsa':topQuarterArea,     
                                            'topHalfTotPsa':topHalfArea,
                                            'topThreeQuarterTotPsa':topThreeQuarterArea,
                                            'TotalShade': area - projected_union_area ,
                                            'topQuarterShade':topQuarterArea-topQuarterAggArea,     
                                            'topHalfShade':topHalfArea-topHalfAggArea,
                                            'topThreeQuarterShade':topThreeQuarterArea-topThreeQuarterAggArea,
                                            'DBH':DBH,
                                            'volume':np.sum(volume),
                                            'X_max':X_max,
                                            'Y_max':Y_max,
                                            'Z_max':Z_max,
                                            'Order_zero_angle_avg':np.average(pd.DataFrame(angles)[self._BO==0][0]),
                                            'Order_zero_angle_std':np.std(pd.DataFrame(angles)[self._BO==0][0]),
                                            'Order_one_angle_avg':np.average(pd.DataFrame(angles)[self._BO==1][0]),
                                            'Order_one_angle_std':np.std(pd.DataFrame(angles)[self._BO==1][0]),
                                            'Order_two_angle_avg':np.average(pd.DataFrame(angles)[self._BO==2][0]),
                                            'Order_two_angle_std':np.std(pd.DataFrame(angles)[self._BO==2][0]),
                                            'Order_three_angle_avg':np.average(pd.DataFrame(angles)[self._BO==3][0]),
                                            'Order_three_angle_std':np.std(pd.DataFrame(angles)[self._BO==3][0]),
                                            'order_gr_four_angle_avg':np.average(pd.DataFrame(angles)[self._BO>=3][0]),
                                            'order_gr_four_angle_std':np.std(pd.DataFrame(angles)[self._BO>=3][0]),
                                            'fileName':self._fileName + self._projection,
                                            }
        

        return statistics

