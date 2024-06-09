from src.canhydro.CylinderCollection import CylinderCollection
import rustworkx as rx
class RustyCylinderCollection(CylinderCollection): 
    def initialize_digraph_from_rust(
            self, in_flow_grade_lim
            ):
            """This function creates a directed graph and its undirected counterpart.
            Initializes edge attributes as cylinder objects"""
            log.info(
                    "Begining initializing graphs"
            )
            gr = rx.PyDAG()
            trunk_nodes, _ = lam_filter(self.cylinders, lambda: branch_order == 0)
            edges = list((
                (
                    (self.cyl_map[cyl.cyl_id], self.cyl_map[cyl.parent_id], cyl.cyl_id)
                    if (cyl.angle >= in_flow_grade_lim or cyl.cyl_id in trunk_nodes)
                    else ( self.cyl_map[cyl.parent_id],self.cyl_map[cyl.cyl_id], cyl.cyl_id)
                )
                for cyl in self.cylinders
            ))
            # edges = [
            #     (
            #         (int(cyl.cyl_id +1), int(cyl.parent_id+1), cyl.cyl_id)
            #         if cyl.angle >= in_flow_grade_lim
            #         else (int(cyl.parent_id+1), int(cyl.cyl_id+1), cyl.cyl_id)
            #     )
            #     for cyl in self.cylinders
            # ]
            # gr.add_nodes_from([int(cyl.cyl_id + 1) for cyl in self.cylinders])
            # gr.add_nodes_from([0])
            # gr.extend_from_weighted_edge_list(edges)
            gr.add_nodes_from([idc for idc in self.cyl_map.values()])
            try:
                gr.add_edges_from(edges)
            except IndexError as e:
                breakpoint()
            self.digraph = gr

        def find_flow_components_rust(self):
            g = self.digraph
            cyl_map = self.cyl_map
            if not isinstance(g, rx.PyDAG):
                msg = "Find Flow Digraph invoked for undirected graph"
                log.error(msg)
                raise TypeError(msg)
            root_node = 0

            def is_drip_node(node):
                return g.out_degree(node) == 0 and node != root_node
            
            is_divide_node = lambda x : (self.digraph.out_degree(x)>1 or 
                                        (self.digraph.out_degree(x)==1 
                                            and self.digraph.in_degree(x)==0))
            
            divide_nodes = g.filter_nodes(is_divide_node)
            
            #identify drip component
            edge_lookup = g.edge_index_map()

            drip_nodes = g.filter_nodes(is_drip_node)
                
            stem_flow_nodes =rx.ancestors(g, root_node) | {0}
            stem_cylinders = [self.cylinders[node].cyl_id for node in stem_flow_nodes
                                        if node not in divide_nodes]


            # def is_drip_edge(edge):
            #     return edge.cyl_id not in stem_cylinders
            # non_stem_edges = g.filter_edges(is_drip_edge)
            # g_drip =g.edge_subgraph([(edge_lookup[x][0],edge_lookup[x][1]) 
            #                             for x in non_stem_edges])
            g_drip = g.subgraph([n for n in g.nodes if n not in stem_flow_nodes])
            
            drip_components = dict.fromkeys(drip_nodes)
            cyl_to_drip_node = {cyl_id: [] for cyl_id in cyl_map.keys()}
            for drip_node in drip_nodes:
                #get cyl_ids
                edges = g_drip.filter_edges(lambda cyl_id: rx.has_path(g_drip,cyl_map[cyl_id],drip_node))
                cyl_ids = [edge_lookup[x][2].cyl_id for x in edges]
                #store data
                for cyl_id in cyl_ids:
                    cyl_to_drip_node[cyl_id].append(drip_node)
                    self.update_cyl_by_id(cyl_id,'drip_node',drip_node)
                drip_components[drip_node] = cyl_ids


            for cyl in self.cylinders:
                if cyl.cyl_id in stem_cylinders:
                    cyl.is_stem = True

            log.info(
                f"{self.file_name} found to have {len(drip_components)} drip components"
            )
            print('reached_End of find flows')

            self.stem_flow_nodes = stem_flow_nodes
            self.drip_flow_components = drip_components
            self.drip_graph = g_drip
            self.drip_nodes = drip_nodes
            self.cyl_to_drip = cyl_to_drip_node



# Useful snippets 

# drip_components = [rx.dfs_search(g,[dn],visitor) for dn in drip_nodes]
# drip_edges = [(l[0],l[1]) for l in visitor.edges[0]]
# drip_edges.extend([(l[0],l[1]) for l in visitor.edges[1]])
# edge_lookup = g.edge_index_map()
# edges = [(tup[1][0],tup[1][1]) for tup in edge_lookup.items()]

# rev_trans,node_map = rx.transitive_reduction(g)
# drip_components = {}
# for dn in drip_nodes: 
#     trans_nodes = rev_trans.neighbors(dn)
#     orig_nodes = [node_map[n] for n in trans_nodes]
#     drip_components[dn] = orig_nodes




        # drip_components = dict.fromkeys(drip_nodes)
        # cyl_to_drip_node = {cyl_id: [] for cyl_id in cyl_map.keys()}
        # for drip_node in drip_nodes:
        #     #get cyl_ids
        #     edges = g_drip.filter_edges(lambda cyl_id: rx.has_path(g_drip,cyl_map[cyl_id],drip_node))
        #     cyl_ids = [edge_lookup[x][2].cyl_id for x in edges]
        #     #store data
        #     for cyl_id in cyl_ids:
        #         cyl_to_drip_node[cyl_id].append(drip_node)
        #         self.update_cyl_by_id(cyl_id,'drip_node',drip_node)
        #     drip_components[drip_node] = cyl_ids


        # [(idx,(x,y)) for (idx,(x,y,_)) in gr.edge_index_map().items()]