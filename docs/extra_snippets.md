
s the level of detail represented in computational models continues to increase, localized estimates of these conditions are now of interest in the study community dynamics (Sybil, some year), soil nutrient distribution (Cavelier 1997), even global climatic modeling (Philipp, \*need to find the specific paper).

The re-distribution of precipitation via vegetation interception therefore influences a wide array of terrestrial hydrological processes: contributing to run off (Savenije, 2004), recharging subsurface water pools (Friesen, 2020), and returning moisture to the atmosphere via transpiration (Coenders-Gerrits et al., 2020).

canoPyHydro is inteneded to be a developlent tool for extracting environmental data from Lidar point clouds. As TLS technology is increasingly available, so too is point cloud data. The functions presented here allow researchers to tap into a treasure trove of existing data to empowering new inferencees.

Drip points are locations in the tree canopy that water flowing towards the stem of the tree cannot pass, due to the branch angle becoming too steep.


# Functionality Deep Dive

# [flow-identification]Flow Identification

CanoPyHydro's hydrological estimates rely on the classification of QSM cylinders as stemflow contributing or throughfall contributing. The precipitation intercepted by each cylinder is added to a theoretical 'flow', and each flow of water is assumed to flow towards the stem. In the model's simplified view, these flows either reach the stem of the tree or drip to the ground after encountering a cylinder that is too steep to traverse-such points are referred to as 'drip-points'. To identify these 'too-steep' portions of the tree, we choose a 'drip cut-off angle' (configurable by the user) and assume water is only able to flow down branches with an angle greater than the cutoff.
The below diagram demonstrates how a graph based model allows us to use these assumptions to identify which cylinders in a QSM are on some drip-path - and are therefore throughfall contributing - and which are stemflow contributing.

A digraph is then constructed with each edge representing one cylinder, with edges oriented in the presumed direction of flow under the above assumptions. This graph is then traversed, and each cylinder is marked as either contributing to stem flow  - if the stem can be reached from its corresponding edge in the graph - or drip flow - if not.
The algorithm above assigns an id to each of the flows found with 'stemflow' always recieving and id of 0. These flow ids are stored by the cylinder collection in the variable 'cyl_to_drip', a dictionary keyed by cylinder ids and can later be used for calculating the 'size' of the flow (see the Metrics section below) and for creating various visualizations of the canopy watershed.
