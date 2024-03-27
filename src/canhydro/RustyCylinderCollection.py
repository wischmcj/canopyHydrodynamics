

from __future__ import annotations

import multiprocessing as mp
import copy
from itertools import chain

import rustworkx as rx
import numpy as np

from src.canhydro.DataClasses import Flow
from src.canhydro.global_vars import config_vars, log
from src.canhydro.utils import intermitent_log, lam_filter
from src.canhydro.CylinderCollection import CylinderCollection

# sys.stdout = LogFile()

NAME = "CylinderCollection"


# By inheriting the Model class, lambda cyl : cyl.branch_order = br CC gains managed functionality- like lambda searching
class RustyCylinderCollection(CylinderCollection):

    def initialize_digraph_from(
            self, in_flow_grade_lim=config_vars["in_flow_grade_lim"]
        ):
        """This function creates a directed graph and its undirected counterpart.
        Initializes edge attributes as cylinder objects"""
        log.info(
                "Begining initializing graphs"
        )
        gr = rx.PyDiGraph()
       
        edges = [
            (
                (int(cyl.cyl_id +1), int(cyl.parent_id+1), cyl.cyl_id)
                if cyl.angle >= in_flow_grade_lim
                else (int(cyl.parent_id+1), int(cyl.cyl_id+1), cyl.cyl_id)
            )
            for cyl in self.cylinders
        ]
        gr.add_nodes_from([int(cyl.cyl_id + 1) for cyl in self.cylinders])
        gr.add_nodes_from([0])
        gr.extend_from_weighted_edge_list(edges)
        self.digraph = gr

    def get_end_nodes(self) -> list[int]:
        g = self.digraph
        if self.end_nodes:
            return self.end_nodes
        elif len(g.nodes) > 0:
            end_nodes = [n for n in g.nodes if (g.in_degree(n) + g.out_degree(n)) == 1 and n != -1]
            self.end_nodes = end_nodes
            return end_nodes
        else:
            log.warning(
                "Graph not initialized, run <CylinderCollection>.initialize_graph(**args)"
            
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
        


    def find_drip_component(self, idx, pair):
        drip_node, source_divides = pair
        # import pdb; pdb.set_trace()
        paths = [
            rx.digraph_dijkstra_shortest_paths(self.drip_graph, source =div_node, target = drip_node)[drip_node]
            for div_node in source_divides
        ]
        component_nodes = [node for node in set(chain.from_iterable(paths))]
        component_graph = self.drip_graph.subgraph(component_nodes,preserve_attrs=True).copy()
        component_cyls = [
            attr
            for _, _, attr in component_graph.weighted_edge_list()
        ]

        for cyl_id in component_cyls:
           self.cyl_to_drip_node[int(cyl_id)].append(drip_node)
        log.info(
            f"component cyls  {component_cyls}"
        )
        if idx%50 == 0:
            log.info(
                f"completed drip component {idx}"
            )

        return (drip_node,component_cyls)


    def find_flow_components(self, inFlowGradeLim=-1 / 6):
        log.info(
                "Finding flow components"
        )
        self.cyl_to_drip_node = {int(cyl.cyl_id): [] for cyl in self.cylinders}
        g = self.digraph
        root_node = 0


        g_drip = copy.deepcopy(g)
        divide_nodes = [
            node
            for node in g.nodes()
            if node and ( g.out_degree(node) > 1
            or (g.out_degree(node) == 1 and g.in_degree(node) == 0))
        ]
        drip_nodes = [
            node
            for node in g.nodes()
            if node and g.out_degree(node) == 0
        ]
        stem_comp_nodes = list(rx.ancestors(g, root_node))
        stem_component_edges = [(u,v,cyl_id) for (u,v,cyl_id)  
                                    in g.weighted_edge_list() 
                                        if u in stem_comp_nodes and v in stem_comp_nodes]
        stem_flow_component = rx.PyDiGraph()

        stem_flow_component.add_nodes_from(stem_comp_nodes)
        stem_flow_component.add_nodes_from([0])
        stem_flow_component.extend_from_weighted_edge_list(stem_component_edges)

        # [tup for tup in self.digraph.edge_list() if tup not in [tup for tup in stem_flow_component.edge_list()]]
        stem_cylinders = [ cyl_id for _,_,cyl_id in stem_component_edges ]
        log.info(
                "Setting is stem"
        )
        for cyl in self.cylinders:
            if cyl.cyl_id in stem_cylinders:
                cyl.is_stem = True
        log.info(
                        f"creating g drip: divde nodes {len(divide_nodes)} drip nodes {len(drip_nodes)}"
                )
        g_drip = copy.deepcopy(g)

        g_drip.remove_edges_from([(u,v) for (u,v) in stem_flow_component.edge_list()])


        drip_divide_pairings = [
            (
                drip_node,
                [node for node in divide_nodes if rx.has_path(g_drip, node, drip_node)],
            )
            for drip_node in drip_nodes
        ]
        log.info(
                "starting proecessing of drip components "
        )
        drip_components = []

        log.info(
                "Assessing drip points"
        )

        self.drip_graph = g_drip

        num_pairings = len(drip_divide_pairings)
        
        with mp.Pool(5) as p:
                task_pool = [p.apply_async(self.find_drip_component, args=(idx,pair)) 
                                for idx, pair in enumerate(drip_divide_pairings)]
                component_cyl_tuples = [task.get() for task in task_pool]

        self.drip_node_to_cyl = { k:v for k, v in component_cyl_tuples}


        log.info(
            f"{self.file_name} found to have {len(drip_components)} drip components"
        )

        log.info(self.cyl_to_drip_node)

        self.stem_flow_component = stem_flow_component
        # self.drip_flow_components = drip_components
        self.divide_nodes = divide_nodes
        self.drip_nodes = drip_nodes



    # @profile
    def calculate_flows(self, plane: str = "XY"):
        """uses subgraphs from FindFlowComponents to aggregate flow characteristics"""
        log.info(
                "Begining sum flows"
        )
        cyls = self.cylinders
        flow_chars = [None]*(len(self.drip_nodes) +1)    

        def numpy_flow_chars(lambda_filter:function, drip_cyl, index:int):
            arr= np.array([
                    np.array([1,np.float64(cyl.projected_data[plane]["area"]),cyl.surface_area,cyl.angle,cyl.volume,cyl.sa_to_vol ]) 
                        for cyl in cyls if lambda_filter(cyl)
                    ])
            flow = np.sum(arr, axis = 0)
            flow_chars[index] = Flow(flow, 
                                    flow[1], 
                                    flow[2], 
                                    flow[3], 
                                    flow[4], 
                                    flow[5], 
                                    drip_cyl.cyl_id, 
                                    (drip_cyl.x[0], drip_cyl.y[0], drip_cyl.z[0]))
            
        numpy_flow_chars(lambda_filter=lambda x: x.is_stem, drip_cyl=self.cylinders[0], index =0 )
         
        num_drip_nodes = len(self.drip_nodes)
        
        log.info(f"begining drip flow calculation for {num_drip_nodes} drip nodes") 

        for drip_node, upstream_cyls in enumerate(self.drip_node_to_cyl):
            filt_func = lambda cyl: (cyl in upstream_cyls )

            cyl_before_drip = [cyl for cyl in cyls if cyl.cyl_id == drip_node-1]
            if len(cyl_before_drip) > 1:
                log.warning(f"Error: More that 1 cyl with id {drip_node} found")

            filt_func = lambda cyl: (drip_node in self.cyl_to_drip[cyl.cyl_id] )

            intermitent_log(idx, num_drip_nodes, "running numpy_flow_chars for drip nodes: ")
            numpy_flow_chars(lambda_filter=filt_func, drip_cyl=cyl_before_drip[0], index = idx+1)

        self.flows = flow_chars