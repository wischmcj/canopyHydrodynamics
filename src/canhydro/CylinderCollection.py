from __future__ import annotations

import logging
import toml
import copy
from itertools import chain

import networkx as nx
import numpy as np
from shapely.geometry import Point
from shapely.ops import unary_union

from src.canhydro.Cylinder import create_cyl
from src.canhydro.DataClasses import Flow
from src.canhydro.geometry import (concave_hull, draw_cyls, 
                                    furthest_point,
                                    get_projected_overlap)
from src.canhydro.utils import intermitent_log, lam_filter, save_file
from src.canhydro.import_options import _try_import

#Optional imports    
if has_geopandas := _try_import('geopandas'):
    from geopandas import GeoSeries

if has_mem_profiler := _try_import('memory_profiler'):
    from memory_profiler import LogFile
    import sys
    sys.stdout = LogFile()

if has_matplotlib := _try_import('matplotlib'):
    import matplotlib.pyplot as plt

if has_spatial := _try_import('scipy.spatial'):
    from scipy.spatial import distance


log = logging.getLogger("model")


with open("src/canhydro/user_def_config.toml") as f:
    config = toml.load(f)
    in_flow_grade_lim=config['config_vars']["in_flow_grade_lim"]

NAME = "CylinderCollection"


# By inheriting the Model class, lambda cyl : cyl.branch_order = br CC gains managed functionality- like lambda searching
class CylinderCollection:
    cylinders = dict

    # initialize our object level variables for cylinder objects
    def __init__(self) -> None:
        self.file = ""
        # Aggregate values from file
        self.surface_area = np.nan
        self.file_name = ""
        self.volume = np.nan
        self.avg_sa_to_vol = np.nan
        self.max_branch_order = np.nan
        self.max_reverse_branch_order = np.nan
        self.extent = {
            "min": [np.nan, np.nan, np.nan],
            "max": [np.nan, np.nan, np.nan],
        }
        # to populate with x,y,z mins and maxs
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
        self.hull = None
        self.stem_hull = None

        # Special case tree attributes
        self.stem_paths = [[]]  
        self.trunk = [] 

        # Graph and Attributes
        self.graph = None
        self.digraph = None
        self.flows = None
        self.divide_nodes = None
        self.drip_nodes = None
        self.cyl_to_drip = None
        self.drip_points = {"x": np.nan, "y": np.nan, "z": np.nan, "flow_id": np.nan}
        self.flow_to_drip = {
            0: 1 #allows us to bootstrap stemflow 
                    # (cyl 1 drains to 0, so its child cyl 2 drains to 0, etc)
        }  # A dictionary of flow ids with values equal to their drip node ids
        self.trunk_nodes = []
        self.drip_point_loc = None
        self.stem_flow_component = None
        self.drip_graph = None
        self.drip_flow_components = None
        # Calculations using graph results
        self.flow_chars = {}
        self.trunk_lean = None

    def from_csv(self, file, aggregate_cyls=True):
        """Initializes a new Cyl Collection based on the data in a QSM
        with the configured column locations"""
        self.file = file
        self.file_name = file.name
        log.info(f"Processing {str(file)}")
        self.arr = np.genfromtxt(file, delimiter=",", skip_header=True)[0:, :-1]
        cylinders = [create_cyl(row) for row in self.arr]
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
            polys = []
            log.info(f"Projection into {plane} axis begun for file {self.file_name}")
            for idx, cyl in enumerate(self.cylinders):
                poly = cyl.get_projection(plane)
                polys.append(poly)
                # print a progress update once every 10 thousand or so cylinders
                intermitent_log(idx, self.no_cylinders, "Cylinder projection: ")
            # Used by other functions to know what projections have been run
            self.projections[plane] = True
            self.pSV = polys
            log.info(f"Projection into {plane} axis complete for file {self.file_name}")

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
        draw_cyls(collection=to_draw, colors=matches, **args)

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

    def get_trunk_nodes(self) -> list[int]:
        g = self.digraph
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

    def watershed_boundary(
        self,
        component=None,
        plane: str = "XY",
        curvature_alpha: np.float64 = 2,
        stem: bool = False,
        draw: bool = False,
    ) -> None:
        """Generates tightly fit concave_hull (alpha shape) around the passed component"""
        """Alpha determines the tightness of the fit of the shape. Too low an alpha results in multiple"""
        if not self.projections[plane]:
            self.project_cylinders(plane)
        # want to get a ring of points around the component
        # however, the root node also has degree one, so we exclude it with n!= -1
        component = component if component else self.digraph
        end_nodes = [n for n in component.nodes if component.degree(n) == 1 and n != -1]

        endCyls = [
            cyl
            for cyl in self.cylinders
            if cyl.cyl_id in end_nodes or cyl.branch_order == 0
        ]
        boundary_points = [Point(cyl.x[1], cyl.y[1]) for cyl in endCyls]

        hull, _, _ = concave_hull(boundary_points, curvature_alpha)
        if draw:
            draw_cyls([hull])
        if stem:
            self.stem_hull = hull
        else:
            self.hull = hull

    def initialize_digraph_from(
        self, in_flow_grade_lim=config_vars["in_flow_grade_lim"]
    ):
        """This function creates a directed graph and its undirected counterpart.
        Initializes edge attributes as cylinder objects"""
        gr = nx.DiGraph()
        edges = (
            (
                (int(cyl.cyl_id), int(cyl.parent_id), {"cylinder": cyl})
                if cyl.angle >= in_flow_grade_lim
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

    def find_trunk_distances(self):
        """
        Finds the distance in the graph (in number of nodes) between each node and the closest trunk node
        """
        trunk_nodes = self.get_trunk_nodes()

        trunk_contraction, _ = self.contracted_nodes(trunk_nodes, 0, self.graph)
        trunk_paths = nx.shortest_path(trunk_contraction, target=0)
        dists = {node: len(path) - 1 for node, path in trunk_paths.items()}
        return dists

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
        drip_nodes = [
            node
            for node, out_degree in g.out_degree()
            if out_degree == 0 and node != -1
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
            for idx in component_cyls:
                cyl_to_drip_node[idx].append(drip_node)

            drip_components.append((drip_node, g_drip.subgraph(component_nodes).copy()))

        for cyl in self.cylinders:
            if cyl_to_drip_node[cyl.cyl_id]:
                cyl.drip_node = cyl_to_drip_node[cyl.cyl_id][-1]

        log.info(
            f"{self.file_name} found to have {len(drip_components)} drip components"
        )

        self.stem_flow_component = stem_flow_component
        self.drip_flow_components = drip_components
        self.drip_graph = g_drip
        self.divide_nodes = divide_nodes
        self.drip_nodes = drip_nodes
        self.cyl_to_drip = cyl_to_drip_node

    # @profile
    def calculate_flows(self, plane: str = "XY"):
        """uses subgraphs from FindFlowComponents to aggregate flow characteristics"""
        cyls = self.cylinders
        flow_chars = []
        # edge_attributes = {}

        # log.info("attempting to sum stem edges ")

        # this probably doesn't belong here but its efficient to do it now
        # is 'needed' for statistics section
        stem_cyls = [cyl for cyl in cyls if cyl.is_stem]

        num_stem_edges = len(stem_cyls)
        flow_chars.append(
            Flow(
                **{
                    "num_cylinders": num_stem_edges,
                    "projected_area": np.sum(
                        [
                            np.float64(cyl.projected_data[plane]["area"])
                            for cyl in stem_cyls
                        ]
                    ),
                    "surface_area": np.sum([cyl.surface_area for cyl in stem_cyls]),
                    "angle_sum": np.sum([cyl.angle for cyl in stem_cyls]),
                    "volume": np.sum([cyl.volume for cyl in stem_cyls]),
                    "sa_to_vol": np.sum([cyl.sa_to_vol for cyl in stem_cyls]),
                    "drip_node_id": 0,
                    "drip_node_loc": (self.cylinders[0].x[0], self.cylinders[0].y[0]),
                }
            )
        )
        # log.info(f"summed stem edges {flow_chars}")
        for idx, drip_node in enumerate(self.drip_nodes):
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

            flow_chars.append(
                Flow(
                    **{
                        "num_cylinders": len(drip_cylinders),
                        "projected_area": np.sum(
                            [
                                np.float64(cyl.projected_data[plane]["area"])
                                for cyl in drip_cylinders
                            ]
                        ),
                        "surface_area": np.sum(
                            [cyl.surface_area for cyl in drip_cylinders]
                        ),
                        "angle_sum": np.sum([cyl.angle for cyl in drip_cylinders]),
                        "volume": np.sum([cyl.volume for cyl in drip_cylinders]),
                        "sa_to_vol": np.sum([cyl.sa_to_vol for cyl in drip_cylinders]),
                        "drip_node_id": drip_node,
                        "drip_node_loc": drip_node_loc,
                    }
                )
            )
            # nx.set_edge_attributes(g, edge_attributes, "dripNode")
            log.info(f"summed drip edges in component {idx}")
        self.flows = flow_chars

    def identify_stem_paths(
        self,
    ):
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
        run = np.sqrt(self.dx**2 + self.dy**2)
        angle = (
            np.arctan(self.dz / np.sqrt(self.dx**2 + self.dy**2))
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

    def statistics(self, plane: str = "XY"):
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
            )

        canopy_cover = self.hull.area
        canopy_boundary = self.hull.boundary.length

        canopy_cover_stem = self.stem_hull.area
        canopy_boundary_stem = self.stem_hull.boundary.length

        log.info("Found hull alpha shape stats")

        # calculate projected areas and (therefore) overlaps
        # Unary union gives a single contiguous polygon when fed many overlapping cylinders
        # The area of the union thus differs from the sum of the areas of its components
        # in that the former counts overlaps only once
        tot_poly = unary_union(self.pSV)
        projected_union_area = tot_poly.area
        sum_projected_area = np.sum([poly.area for poly in self.pSV])

        stem_polys = [
            cyl.projected_data[plane]["polygon"]
            for cyl in self.cylinders
            if cyl.is_stem
        ]
        union_poly_stem = unary_union(stem_polys[plane])
        projected_union_area_stem = union_poly_stem.area
        sum_projected_area_stem = np.sum([poly.area for poly in stem_polys[plane]])

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

        stem_flow = [flow for flow in self.flows if flow.drip_node_id == 0][0]

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
            "max_bo": max_bo,
            "topQuarterTotPsa": overlap_dict[75]["sum_area"],
            "topHalfTotPsa": overlap_dict[50]["sum_area"],
            "topThreeQuarterTotPsa": overlap_dict[25]["sum_area"],
            "TotalShade": sum_projected_area - projected_union_area,
            "top_quarter_shade": overlap_dict[75]["overlap_with_previous"],
            "top_half_shade": overlap_dict[50]["overlap_with_previous"],
            "top_three_quarter_shade": overlap_dict[25]["overlap_with_previous"],
            "DBH": self.get_dbh(),
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
        save_file(
            self.file_name.replace(".csv", ""),
            out_file=statistics,
            subdir="statistics",
            method="statistics",
        )

        return statistics

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
        self, a_lambda: function = lambda: True, scale: int = 1, **args
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
        drip_point_locs_xy = [[pt[0] * scale, pt[1] * scale] for pt in drip_point_locs]

        math.floor(np.min(drip_point_locs_x))

        mins = self.extent["min"]
        maxs = self.extent["max"]
        extents = [mins[0], maxs[0], mins[1], maxs[1]]
        # min_xy = np.min(mins)
        # max_xy = np.max(maxs)
        # x_mesh, y_mesh = np.meshgrid(
        #     np.arange(min_xy, max_xy, 0.05), np.arange(min_xy, max_xy, 0.05)
        # )

        min_xy = np.min(
            [
                math.floor(np.min(drip_point_locs_x)),
                math.floor(np.min(drip_point_locs_y)),
            ]
        )
        max_xy = np.max(
            [math.ceil(np.max(drip_point_locs_x)), math.ceil(np.max(drip_point_locs_y))]
        )
        x_mesh, y_mesh = np.meshgrid(
            np.arange(min_xy, max_xy, 0.005), np.arange(min_xy, max_xy, 0.005)
        )

        def dist_to_drip(a, b):
            distances = distance.cdist([[a, b]], drip_point_locs_xy)
            min_dist = np.min(distances)
            return math.log(1 / min_dist)

        distance_matrix = np.zeros((x_mesh.shape[0], x_mesh.shape[0]))

        for a in range(x_mesh.shape[0]):
            for b in range(x_mesh.shape[0]):
                distance_matrix[a][b] = dist_to_drip(x_mesh[b][a], y_mesh[b][a])
        
        if has_matplotlib:
            _, ax = plt.subplots()

            ax.contourf(
                y_mesh,
                x_mesh,
                distance_matrix,
                levels=15,
                max=0.5,
                cmap=plt.cm.Blues,
                extend="neither",
                extent=extents,
            )

            ax.scatter(drip_point_locs_x, drip_point_locs_y)

            filtered_cyls, _ = lam_filter(self.cylinders, a_lambda, return_all=False)
            polys = [cyl.projected_data["XY"]["polygon"] for cyl in filtered_cyls]
            
            if len(polys) > 0:
                geoPolys = GeoSeries(polys)
                geoPolys.plot(ax=ax)
            else:
                log.warning(
                    "Drip Map: No cylinders returned for lambda function: {a_lambda}"
                )
        else:
            log.warning(
                "Unable to draw drip map as module matplotlib.pyplot is not available"
            )
        return distance_matrix
