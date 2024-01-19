

from __future__ import annotations


class AlternativeCylCollection():

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
    def initialize_minimal_graph(self):
        """This function initialized edge attributes as cylinder objects"""
        gr = nx.Graph()
        for cyl in self.cylinders:
            child_node = cyl.cyl_id
            parent_node = cyl.parent_id
            gr.add_edge(child_node, parent_node)
        self.min_graph = gr

    @profile
    def initialize_graph_from(self):
        """This function initialized edge attributes as FULL cylinder dicts"""
        gr = nx.Graph()
        edges = (
            (int(cyl.cyl_id), int(cyl.parent_id), cyl.__dict__)
            for cyl in self.cylinders
        )
        gr.add_edges_from(edges)
        self.graph = gr

    @profile
    def initialize_minimal_graph_from(self):
        """This function initialized edge attributes as cylinder objects"""
        gr = nx.Graph()
        edges = ((int(cyl.cyl_id), int(cyl.parent_id)) for cyl in self.cylinders)
        gr.add_edges_from(edges)
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
    
    
    def find_flow_components(self, inFlowGradeLim: float = -1 / 6, dist: int = 2):
        """Finding Stemflow contributing area"""
        g = self.graph
        # identify drip edges in graph
        trunk_distance = self.find_trunk_distance()
        drip_edges = [
            (u, v)
            for u, v, attr in g.edges(data=True)
            if attr["angle"] < inFlowGradeLim
            and trunk_distance.get(min(u, v), 0) > dist
        ]
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



    def find_flow_components_minimal(self, inFlowGradeLim=-1 / 6):
        """Finding Stemflow contributing area"""
        g = self.min_graph
        # identify drip edges in graph
        all_cyls, is_out = lam_filter(
            self.cylinders, lambda: angle < -1 / 6, return_all=True
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

    def flowCost(self, graph, metric):
        # algo does not play well with floats
        if metric == "projected_area":
            sum(
                attr["projected_data"]["XY"]["area"] * 10000
                for u, v, attr in graph.edges(data=True)
            ) / 10000
        else:
            sum(attr[metric] * 10000 for u, v, attr in graph.edges(data=True)) / 10000


    @profile
    def calculate_flows(self):
        """uses subgraphs from FindFlowComponents to aggregate flow characteristics"""
        stem_flow_component = self.stemFlowComponent
        drip_flow_components = self.dripFlowComponents
        g = self.graph
        edge_attributes = {}
        flow_chars = []
        G_return = copy.deepcopy(self.graph)

        for u, v in stem_flow_component.edges():
            edge_attributes[(u, v)] = {"dripNode": 0, "flowType": "stem", "flowID": 0}

        flow_dict = lambda comp, metric: {
            i: {j: {val}} for i, j, val in comp.edges.data(metric)
        }
        stem_edges = len(stem_flow_component.edges())
        flow_chars.append(
            Flow(
                **{
                    "num_cylinders": stem_edges,
                    "projected_area": self.flowCost(
                        stem_flow_component, "projected_area"
                    ),
                    "surface_area": self.flowCost(stem_flow_component, "surface_area"),
                    "angle_sum": self.flowCost(stem_flow_component, "angle"),
                    "volume": self.flowCost(stem_flow_component, "volume"),
                    "sa_to_vol": self.flowCost(stem_flow_component, "sa_to_vol"),
                    "drip_node_id": 0,
                }
            )
        )  # stemflow drips to the trunk
        for idx, flow in enumerate(drip_flow_components):
            heights = [attr["z"] for _, _, attr in flow.edges(data=True)]
            nodes = [n for n in flow]
            num_cyls = len(flow.edges())
            drip_node = nodes[0]

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
                        "angle_sum": self.flowCost(flow, "angle"),
                        "volume": self.flowCost(flow, "volume"),
                        "sa_to_vol": self.flowCost(flow, "sa_to_vol"),
                        "drip_node_id": 0,
                    }
                )
            )

        self.flows = flow_chars

        nx.set_edge_attributes(g, edge_attributes, "dripNode")
  @profile
    def calculate_flows_min(self):
        """uses subgraphs from FindFlowComponents to aggregate flow characteristics"""
        stem_flow_component = self.stemFlowComponent
        drip_flow_components = self.dripFlowComponents
        g = self.min_graph
        edge_attributes = {}
        flow_chars = []
        G_return = copy.deepcopy(g)
        stem_cyls = [
            cyl
            for cyl in self.cylinders
            if cyl.cyl_id in [n for n in stem_flow_component]
        ]

        for u, v in [(cyl.cyl_id, cyl.parent_id) for cyl in stem_cyls]:
            edge_attributes[(u, v)] = {"dripNode": 0, "flowType": "stem", "flowID": 0}

        stem_edges = len(stem_cyls)
        flow_chars.append(
            Flow(
                **{
                    "num_cylinders": stem_edges,
                    "projected_area": np.sum([cyl.xz_area for cyl in stem_cyls]),
                    "surface_area": np.sum([cyl.surface_area for cyl in stem_cyls]),
                    "angle_sum": np.sum([cyl.angle for cyl in stem_cyls]),
                    "volume": np.sum([cyl.volume for cyl in stem_cyls]),
                    "sa_to_vol": np.sum([cyl.sa_to_vol for cyl in stem_cyls]),
                    "drip_node_id": 0,
                }
            )
        )  # stemflow drips to the trunk
        for idx, flow in enumerate(drip_flow_components):
            nodes = [n for n in flow]
            drip_cyls = [cyl for cyl in self.cylinders if cyl.cyl_id in nodes]
            heights = [(cyl.z[0], cyl.cyl_id) for cyl in drip_cyls]
            num_cyls = len(drip_cyls)
            drip_node = nodes[0]
            # min_height = np.min([h for h,v in heights ])
            # drip_node = [id for (h,id) in heights if h == min_height][0]

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
                        "projected_area": np.sum([cyl.xz_area for cyl in drip_cyls]),
                        "surface_area": np.sum([cyl.surface_area for cyl in drip_cyls]),
                        "angle_sum": np.sum([cyl.angle for cyl in drip_cyls]),
                        "volume": np.sum([cyl.volume for cyl in drip_cyls]),
                        "sa_to_vol": np.sum([cyl.sa_to_vol for cyl in drip_cyls]),
                        "drip_node_id": 0,
                    }
                )
            )

        self.flows = flow_chars

        nx.set_edge_attributes(g, edge_attributes, "dripNode")
class EfficiencyTester():

    @profile
    def min_graph_test():
        forest = Forester()
        forest.get_file_names(dir=test_input_dir)
        forest.qsm_from_file_names(file_name="4_LargeCollection.csv")
        flexible_collection = forest.cylinder_collections[0]
        flexible_collection.project_cylinders("XZ")
        flexible_collection.initialize_minimal_graph_from()
        proj_area = flexible_collection.sum_over_min_graph()
        flexible_collection.find_flow_components_minimal()
        print(proj_area)

    @profile
    def base_graph_test():
        forest = Forester()
        forest.get_file_names(dir=test_input_dir)
        forest.qsm_from_file_names(file_name="4_LargeCollection.csv")
        flexible_collection = forest.cylinder_collections[0]
        flexible_collection.project_cylinders("XZ")
        flexible_collection.initialize_graph_from()
        proj_area = flexible_collection.sum_over_graph()
        flexible_collection.find_flow_components()
        print(proj_area)

    @profile
    def obj_graph_test():
        forest = Forester()
        forest.get_file_names(dir=test_input_dir)
        forest.qsm_from_file_names(file_name="4_LargeCollection.csv")
        flexible_collection = forest.cylinder_collections[0]
        flexible_collection.project_cylinders("XZ")
        flexible_collection.initialize_object_graph_from()
        proj_area = flexible_collection.sum_over_object_graph()
        flexible_collection.find_flow_components_object()
        print(proj_area)

    if __name__ == "__main__":
        current_GMT = time.gmtime()
        time_stamp = str(calendar.timegm(current_GMT))
        log.info(f"base_graph_test started at {time_stamp}")

        base_graph_test()

        current_GMT = time.gmtime()
        time_stamp = str(calendar.timegm(current_GMT))
        log.info(f"min_graph_test started at {time_stamp}")

        min_graph_test()

        current_GMT = time.gmtime()
        time_stamp = str(calendar.timegm(current_GMT))
        log.info(f"obj_graph_test started at {time_stamp}")

        obj_graph_test()

        current_GMT = time.gmtime()
        time_stamp = str(calendar.timegm(current_GMT))
        log.info(f"finished at {time_stamp}")
