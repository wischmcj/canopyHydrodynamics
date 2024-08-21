
s the level of detail represented in computational models continues to increase, localized estimates of these conditions are now of interest in the study community dynamics (Sybil, some year), soil nutrient distribution (Cavelier 1997), even global climatic modeling (Philipp, \*need to find the specific paper).

The re-distribution of precipitation via vegetation interception therefore influences a wide array of terrestrial hydrological processes: contributing to run off (Savenije, 2004), recharging subsurface water pools (Friesen, 2020), and returning moisture to the atmosphere via transpiration (Coenders-Gerrits et al., 2020).

CanHydro is inteneded to be a developlent tool for extracting environmental data from Lidar point clouds. As TLS technology is increasingly available, so too is point cloud data. The functions presented here allow researchers to tap into a treasure trove of existing data to empowering new inferencees.

Drip points are locations in the tree canopy that water flowing towards the stem of the tree cannot pass, due to the branch angle becoming too steep.



# [flow-identification]Flow Identification
CanoPyHydro's estimations of precipitation partitioning rely on the identification of 'drip-points' in the canopy. Drip points are locations in the tree canopy that water flowing towards the stem of the tree cannot pass, due to the branch angle becoming too steep. Which partition (stemflow,throughfall) each portion of the tree contributes intercepted precipitation to then relies on their being path from that portion of the tree to the stem of the tree, on which the branch angle allows for the flow of water. To identify these 'too-steep' portions of the tree, we choose a 'drip cut-off angle' (configurable by the user) such that angles more negative than this cutoff do not allow for the flow of intercepted percipitation.  Aflow of water is assumed to flow towards the stem unless the angle of an encountered branch is less than the drip cutoff angle.
A digraph is then constructed with each edge representing one cylinder, with edges oriented in the presumed direction of flow  under these assumptions This graph is then traversed, and each cylinder is marked as either contributing to stem flow  - if the stem can be reached from its corresponding edge in the graph - or drip flow - if not.
