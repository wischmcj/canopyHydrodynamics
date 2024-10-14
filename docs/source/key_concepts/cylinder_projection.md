### Aggregating 2D Area

CanoPyHydro can also provide detailed intra-canopy occlusion data for given heights/depths. For a vertical (XZ) projection, this represents shading by higher branches on lower branches, for hoizontal projections (XZ, YZ) this represents wind exposure (or lack thereof). In future versions, arbitrary projection angles may be used to assist in calclating the effect of occlusion on partitioning in various different weather conditions.
Depending on the goals of the user, the projected area of a collection of cylinders can be given as both:

- a simple sum of the projected area of each cylinder

  - Note that this calculation will ignore overlap between cylinder areas
    - Add example

- the total projected area of the collection of cylinders

  - Using this approach, areas in which cylinder projections over lap are only counted once
    - Add Example



# Canopy Coverage
As a part of the canopy coverage analysis, we are interested in the shape of the canopy. Among other reasons, we are interested in this as the area covered by the canopy is the denominator is several important metrics commonly used to compare one tree to another such as woody area index (WAI). Similarly, users leveraging CanoPyHydro's percipitation partitioning estimates will find it beneficial to produce similar metrics using the amount of rain that becomes throughfall/stemflow and the total amount intercepted by the canopy as a numerator.
**Traditional measurements of canopy coverage area will not sufice for the latter use case** as often subsets of a tree's branches are being investigated. In this case, the shape of the canopy diverges from a simple circle

### Cylinder Overlap (Shade)

In order to better understand how the branches of a tree's canopy overlap, more granular overlap information has also been made available via the 'find_overlap_by_percentile' function.
When considered from a birds eye view (projecting onto the XY plane), this concept can be understood as a facsimile for the 'shade' cast by branches at a certain height in the tree canopy.
Consider the below example \* Add example
The percentile list is used to determine the height at which to calculate shade. As such, the function will look at the overlap between cylinders in the 75%ile by height (in red) with the remaining cylinders (in blue). The returned values thus represent the 'shade' case by the red cylinders on the blue cylinders.

    Following this logic if either the 0%ile or 100%ile is requested, then there will be no overlap reported. In the former case, all cylinders are included in the red group and therefore there are no blue cylinders on which to cast shade. In the latter case, all cylinders are in the blue group and so there are no red cylinders to cast shade onto the blue cylinders.

When considered in the XZ or YZ direction, this calculation can be useful in determining the wind exposure at different canopy depths.
