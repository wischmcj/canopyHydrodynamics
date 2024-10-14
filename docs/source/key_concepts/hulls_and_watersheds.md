

# Canopy Coverage
As a part of the canopy coverage analysis, we are interested in the shape of the canopy. Among other reasons, we are interested in this as the area covered by the canopy is the denominator is several important metrics commonly used to compare one tree to another such as woody area index (WAI). Similarly, users leveraging CanoPyHydro's percipitation partitioning estimates will find it beneficial to produce similar metrics using the amount of rain that becomes throughfall/stemflow and the total amount intercepted by the canopy as a numerator.
**Traditional measurements of canopy coverage area will not sufice for the latter use case** as often subsets of a tree's branches are being investigated. In this case, the shape of the canopy diverges from a simple circle

# Delauny Triangulation
Delauny Triangulation is a method of creating a network of triangles from a set of points. This method is used in a variety of applications, including the creation of a convex hull, which is a shape that encloses a set of points. In the context of canopy coverage, Delauny Triangulation can be used to create a network of triangles that represent a set of points in the canopy. This network can then be used to calculate the area of the canopy as a whole, as well as other metrics such as the density of the canopy.

# Alpha Shapes/Concave Hulls
One key application of delauny triangulation is the creation of concave hulls.

The idea of a concave hull is to create a minimal shape that encloses a set of points. More specifically, the shape is convex if for any 2 of these points that fall on the boundary of the shape and are adjacent, one can draw a line between the two points without intersecting the shape itsself. These shapes are also generalized into Alpha Shapes, as the curvature of the shape's edges is defined by a parameter, alpha, which describes its level of concavity. In fact, a convex hull is frequently referred to ass the smallest convex shape that encloses a set of points.


## Usage Examples
To see these objects in action, navigate to the [Examples](../examples/examples_index.md#examples) page and see
