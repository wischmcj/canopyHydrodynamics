# Functionality Deep Dive

# [flow-identification]Flow Identification

CanoPyHydro's hydrological estimates rely on the classification of QSM cylinders as stemflow contributing or throughfall contributing. The precipitation intercepted by each cylinder is added to a theoretical 'flow', and each flow of water is assumed to flow towards the stem. In the model's simplified view, these flows either reach the stem of the tree or drip to the ground after encountering a cylinder that is too steep to traverse-such points are referred to as 'drip-points'. To identify these 'too-steep' portions of the tree, we choose a 'drip cut-off angle' (configurable by the user) and assume water is only able to flow down branches with an angle greater than the cutoff.
The below diagram demonstrates how a graph based model allows us to use these assumptions to identify which cylinders in a QSM are on some drip-path - and are therefore throughfall contributing - and which are stemflow contributing.

A digraph is then constructed with each edge representing one cylinder, with edges oriented in the presumed direction of flow under the above assumptions. This graph is then traversed, and each cylinder is marked as either contributing to stem flow  - if the stem can be reached from its corresponding edge in the graph - or drip flow - if not.
The algorithm above assigns an id to each of the flows found with 'stemflow' always recieving and id of 0. These flow ids are stored by the cylinder collection in the variable 'cyl_to_drip', a dictionary keyed by cylinder ids and can later be used for calculating the 'size' of the flow (see the Metrics section below) and for creating various visualizations of the canopy watershed.
