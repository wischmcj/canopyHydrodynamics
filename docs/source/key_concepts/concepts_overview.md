
# Data Storage Objects

## Cylinders
The Cylinder class is used to represent the 3-D cylinders that make up a QSM. The most important function of Cylinder objects is their ability to return data regarding the projections onto planes. Cylinder objects utilize our custom 'geometry' module to calculate their projections onto the XY, XZ and YZ planes.

## Cylinder Collections
Cylinder Collections are just as they sound and, at the most basic level, a Cylinder Collection is defined as a list of 1 or more Cylinder objects. Cylinder Collections almost always represent QSM's (or parts of a QSM), and are meant to help users explore their QSMs. Below, we demonstrate how one might initialize a cylinder collection using cylinder data (e.g. QSM data) stored in a CSV file.

## Foresters
Forester objects allow users to conveniently create and manage Cylinder Collections. In particular, Foresters are useful for reading in and processing QSM files. When a Forester object is created, available file names are read from the default directory, './data/input/' or, if specified, from the directory specified by the user. This list of available files can be accessed through the Forester.file_names attribute.

## Usage Examples
To see these objects in action, navigate to the [Examples](../examples/examples_index.md) page.

# 3D -> 2D Considerations 

### Aggregating 2D Area

Depending on the goals of the user, the projected area of a collection of cylinders can be given as both:

- a simple sum of the projected area of each cylinder
  - Note that this calculation will ignore overlap between cylinder areas

- the total projected area of the collection of cylinders
  - Using this approach, areas in which cylinder projections overlap are only counted once

### Cylinder Overlap (Shade)

In order to better understand how the branches of a tree's canopy overlap, more granular overlap information has also been made available via the 'find_overlap_by_percentile' function. When considered from a birds eye view (projecting onto the XY plane), this concept can be understood as a facsimile for the 'shade' cast by branches at a certain height in the tree canopy.
Given a list of arbitrary percentiles, and an orientation in which the data is desired, canopyHydro can provide detailed information rearding the overlap between cylinders above and below each percentile in regards to canopy height (in a vertical orientation) or depth (in a horizontal orientation). 
The classic use case focuses on the XY direction (birds-eye view) and evaluates the shade at different heights in the canopy. However, when considered in the XZ or YZ direction, this calculation can be useful in determining wind exposure at different canopy depths.


# Canopy Coverage

## Bounding the Canopy
As a part of the canopy coverage analysis, we are interested in the shape of the canopy. Among other reasons, we are interested in this as the area covered by the canopy is the denominator is several important metrics commonly used to compare one tree to another such as woody area index (WAI). Similarly, users leveraging CanoPyHydro's percipitation partitioning estimates will find it beneficial to produce similar metrics using the amount of rain that becomes throughfall/stemflow and the total amount intercepted by the canopy as a numerator.
**Traditional measurements of canopy coverage area will not sufice for the latter use case** as often subsets of a tree's branches are being investigated. In this case, the shape of the canopy diverges from a simple circle

## Delauny Triangulation
Delauny Triangulation is a method of creating a network of triangles from a set of points. This method is used in a variety of applications, including the creation of a convex hull, which is a shape that encloses a set of points. In the context of canopy coverage, Delauny Triangulation can be used to create a network of triangles that represent a set of points in the canopy. This network can then be used to define a polygon encompassing the canopy, allowing for area calculations as well metrics like leafy area index (LAI) that help define the density of the canopy.

## Alpha Shapes/Concave Hulls
One key application of delauny triangulation is the creation of concave hulls.

The idea of a concave hull is to create a minimal shape that encloses a set of points. More specifically, the shape is convex if for any 2 of these points that fall on the boundary of the shape and are adjacent, one can draw a line between the two points without intersecting the shape itsself. These shapes are also generalized into Alpha Shapes, as the curvature of the shape's edges is defined by a parameter, alpha, which describes its level of concavity. In fact, a convex hull is frequently referred to ass the smallest convex shape that encloses a set of points.


## Usage Examples
To see these objects in action, navigate to the [Examples](../examples/examples_index.md) page and see
