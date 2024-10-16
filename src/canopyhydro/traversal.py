
from dataclasses import dataclass
import rustworkx as rx
import networkx as nx
import numpy as np

import copy
from canopyhydro.CylinderCollection import CylinderCollection
from canopyhydro.configuration import log, in_flow_grade_lim

class rain_drop(rx.visit.DFSVisitor):
    def __init__(self,divide_nodes):
        self.divide_nodes = divide_nodes
        self.edges = []
        self.results = {}
        self.source = -1
    # def __getitem__(self,item):
    # return chain.from_iterable([res[item] for res in self.results])
    def discover_vertex(self, v, _):
        if self.source == -1:
            self.source = v
    def forward_or_cross_edge(self, edge):
        if edge[1] in self.divide_nodes:
            self.edges.append(edge[2])
            # self.cyl_to_drip_node[edge[2].cyl_id] = source
            edge[2].drip_node = self.source
    def tree_edge(self, edge):
        self.edges.append(edge[2])
        edge[2].drip_node = self.source
    def finish_vertex(self, v, t):
        if v == self.source:
            self.results[v] = self.edges
            self.source = -1
            self.edges = []

class RustyCollection(CylinderCollection):

    def initialize_digraph_from_rust(self,
                                        in_flow_grade_lim=in_flow_grade_lim):
        """This function creates a directed graph and its undirected counterpart.
        Initializes edge attributes as cylinder objects"""

        log.info("Begining initializing graphs")
        gr = rx.PyDAG()
        gr.add_nodes_from(self.cyl_map.values()-1)

        is_out_edge = lambda cyl: (
            cyl.angle < in_flow_grade_lim and cyl.branch_order != 0
        )
        non_base_cyls = [(cyl.cyl_id,cyl.parent_id) for cyl in self.cylinders 
                                            if cyl.parent_id != -1]

        cyls_as_undirected_edges = [
            (self.cyl_id_to_node(cyl.cyl_id), self.cyl_id_to_node(cyl.parent_id)) 
                for cyl in non_base_cyls
        ]
        flow_direction = [is_out_edge(cyl) for cyl in non_base_cyls]
        all_edge_data = zip(cyls_as_undirected_edges, flow_direction)
        edges = [
            (v, u, self.cylinders[u]) if out_flow else (u, v, self.cylinders[u])
            for ((u, v), out_flow) in all_edge_data
        ]

        gr.extend_from_weighted_edge_list(edges)
        self.digraph = gr

    # def classify_nodes():
    #     """Classifies a node based on its out degree"""
    #     if degree > 1:
    #         return "divide"
    #     elif degree == 1:
    #         return "drip"
    #     else:
    #         return "stem"

    # def get_divide_nodes(self):
    #     if self.divide_nodes:
    #         return self.divide_nodes
    #     g = self.digraph
    #     divide_nodes = [
    #         node
    #         for node in g.nodes()
    #         if g.out_degree(node) > 1
    #         or (g.out_degree(node) == 1 and g.in_degree(node) == 0)
    #     ]
    #     self.divide_nodes = divide_nodes
    #     return divide_nodes

    def find_flow_components_rust(self):
        g = self.digraph
        cyl_map = self.cyl_map
        if not isinstance(g, rx.PyDAG):
            msg = "Find Flow Digraph invoked for undirected graph"
            log.error(msg)
            raise TypeError(msg)

        root_node = 0

        is_divide_node = lambda x: (
            self.digraph.out_degree(x) > 1
            or (self.digraph.out_degree(x) == 1 and self.digraph.in_degree(x) == 0)
        )

        divide_nodes = g.filter_nodes(is_divide_node)
        # edge_to_cyl_id = {
        #     (u, v): cyl.cyl_id for (_, (u, v, cyl)) in g.edge_index_map().items()
        # }

        def is_drip_node(node):
            return g.out_degree(node) == 0 or node == root_node  # and node != root_node

        drip_nodes = g.filter_nodes(is_drip_node)
        g.reverse()
        # visitors = {}
        visitor = rain_drop(divide_nodes)
        rx.dfs_search(g, drip_nodes, visitor)
        visitor.results[0].append(self.cylinders[0])
        for edge in visitor.results[0]:
            edge.is_stem = True

        flows_raw = [
            ([x.get_flow_arr() for x in res], self.cylinders[dn])
            for dn, res in visitor.results.items()
        ]
        # breakpoint()
        self.flows = [
            Flow(*np.sum(flow, axis=0), cyl.cyl_id, cyl.get_loc())
            for flow, cyl in flows_raw
        ]

        # for res in visitor.results:
        #     source = res.source_node
        #     cyls = [edge[2].cyl_id for edge in res['edges']]
        #     drip_components[source] = cyls
        #     source = res.source_node
        #     for cyl_id in cyls:
        #         cyl_to_drip_node[cyl_id] = source

        # breakpoint()
        # stem_flow_cylinders = drip_components.pop(0)
        # stem_flow_cylinders.append(root_node)
        # #nodes are named after the cylinder directly preceding them
        # # The root node gets -1 as there is no cylinder preceding it
        # stem_flow_nodes = [edge[0] for edge in res['edges'] if edge[2].cyl_id  ]
        # cyl_to_drip_node[root_node].append(0)

        # for node in stem_flow_cylinders:
        #     self.cylinders[node].is_stem = True

        drip_components = {
            source: [edge.cyl_id for edge in edges]
            for source, edges in visitor.results.items()
        }
        drip_components[0].append(0)
        log.info(
            f"{self.file_name} found to have {len(drip_components)} drip components"
        )
        g.reverse()
        print("reached_End of find flows")
        self.stem_flow_nodes = [
            self.cyl_idx(cyl.cyl_id) for cyl in self.cylinders if cyl.is_stem
        ]
        self.drip_flow_components = drip_components
        # self.drip_graph = g_drip
        self.drip_nodes = drip_nodes
        # self.cyl_to_drip = cyl_to_drip_node

    def calculate_flows_matrix(self, plane: str = "XY"):
        """uses subgraphs from FindFlowComponents to aggregate flow characteristics"""
        cyls = self.cylinders
        log.info("attempting to sum stem edges ")

        # np_drip_array = np.array([[1 if cyl_id in drip_cylinder_ids else 0 for cyl_id in self.cyl_map.keys()]
        #                           for (_,drip_cylinder_ids) in self.drip_flow_components.items()])
        all_flow_data = np.array(
            [
                [
                    [
                        1,
                        np.float16(cyl.projected_data[plane]["area"]),
                        cyl.surface_area,
                        cyl.angle,
                        cyl.volume,
                        cyl.sa_to_vol,
                    ]
                    for cyl in self.cylinders
                    if cyl.cyl_id in flow_cyls
                ]
                for flow_cyls in self.drip_flow_components.values()
            ]
        )

        flow_arrs = [np.sum(flow_data, axis=0) for flow_data in all_flow_data]
        drip_cyls = [(node, self.cylinders[node]) for node in self.drip_nodes]
        drip_node_locs = [
            (node, (cyl.x[0], cyl.y[0], cyl.z[0])) for node, cyl in drip_cyls
        ]
        flows = [
            np.array(*data, node, loc)
            for (node, loc), data in zip(
                drip_node_locs,
                flow_arrs,
            )
        ]
        self.flows = flows

    # def calculate_flows_rust(self, plane: str = "XY"):
    #     all_flow_data = np.array([ [[    1,
    #                                     np.float16(cyl.projected_data[plane]["area"]),
    #                                     cyl.surface_area,
    #                                     cyl.angle,
    #                                     cyl.volume,
    #                                 cyl.sa_to_vol ]
    #                                 for cyl in self.cylinders
    # flows = defaultdict(list)
    # for drip_node, cyl_ids in self.drip_flow_components.items():
    #     flows[drip_node] = [cyl for cyl in self.cylinders if cyl.cyl_id in cyl_ids]

    # radius = self.arr[:, qsm_cols["radius"]]
    # length = self.arr[:, qsm_cols["length"]]
    # projected_area = [np.float16(cyl.projected_data[plane]["area"]) for cyl in self.cylinders]
    # angle = [cyl.angle for cyl in self.cylinders]
    # volume = [cyl.volume for cyl in self.cylinders]
    # sa_to_vol = [cyl.sa_to_vol for cyl in self.cylinders]
    # self.metric_dict = {
    #                 'surface_area': 2 * np.pi * radius * (radius + length) - 2 * np.pi * radius * radius
    #                 ,'volume':  arr[:, qsm_cols['volume']]
    #                 ,

    # }
    # self.surface_area = 2*np.pi*
    # self.surface_area = arr[:, qsm_cols["surface_area"]]
    # self.volume = arr[:, qsm_cols["surface_area"]]
    # self.surface_area = arr[:, qsm_cols["surface_area"]]
