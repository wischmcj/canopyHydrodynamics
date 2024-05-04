 

def old_find_flow_components(self, inFlowGradeLim=-1 / 6):
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