
<head>
   <meta charset=utf-8 />
   <title></title>
   <style>
    div.container {
      display:inline-block;
    };
    p {
      text-align:center;
    };
    img {
      display: block;
      margin-left: auto;
      margin-right: auto;
    };
    figcaption {
      font-size: 10px;
      text-align:center;
    }
   </style>
</head>

# Canopy Coverage
As a part of the canopy coverage analysis, we are interested in the shape of the canopy. Among other reasons, we are interested in this as the area covered by the canopy is the denominator is several important metrics commonly used to compare one tree to another such as woody area index (WAI). Similarly, users leveraging CanoPyHydro's percipitation partitioning estimates will find it beneficial to produce similar metrics using the amount of rain that becomes throughfall/stemflow and the total amount intercepted by the canopy as a numerator.
**Traditional measurements of canopy coverage area will not suffuce for the latter use case** as often subsets of a tree's branches are being investigated. In this case, the shape of the canopy diverges from a simple circle

# Delauny Triangulation
Delauny Triangulation is a method of creating a network of triangles from a set of points. This method is used in a variety of applications, including the creation of a convex hull, which is a shape that encloses a set of points. In the context of canopy coverage, Delauny Triangulation can be used to create a network of triangles that represent a set of points in the canopy. This network can then be used to calculate the area of the canopy as a whole, as well as other metrics such as the density of the canopy.

# Alpha Shapes/Concave Hulls
One key application of delauny triangulation is the creation of concave hulls.

The idea of a concave hull is to create a minimal shape that encloses a set of points. More specifically, the shape is convex if for any 2 of these points that fall on the boundary of the shape and are adjacent, one can draw a line between the two points without intersecting the shape itsself. These shapes are also generalized into Alpha Shapes, as the curvature of the shape's edges is defined by a parameter, alpha, which describes its level of concavity. In fact, a convex hull is frequently referred to ass the smallest convex shape that encloses a set of points.

To help understand this, consider the below example. In this example, we will draw an alpha shape around an example tree.

```{python}
# Drawing various branch orders
from matplotlib import pyplot as plt
import os
os.environ["CANOPYHYDRO_CONFIG"] = "./canopyhydro_config.toml"
from src.canopyhydro.CylinderCollection import CylinderCollection
myCollection = CylinderCollection()
myCollection.from_csv("minimal_tree.csv")

myCollection.draw("XY") # birds-eye view
myCollection.draw("XZ") # 'front' view
```

<div class="container">
    <img src="../imgs/minimal_tree_XZ_hull_tutorial.png" height="300" width="300">
    <img src="../imgs/minimal_tree_XY_hull_tutorial.png" height="300" width="300">
</p>

As you can see above, our tree has only a few, short branches. As is common among groups of smaller groups of cylinders (i.e. branch subsets), there are also significant gaps between the tree's branches.

For this example, we will pass the 'include_alpha_shape' parameter to the 'draw' function, which will automatically call the 'watershed_boundary' function to generate an alpha shape around the tree. As the 'watershed_boundary' function defaults to using the tree's branches' endpoints as the input point set, the alpha shape will be drawn with consideration of the below points

```{python}
import matplotlib.pyplot as plt
# First we pull the end nodes of the tree
# which returns the outermost cylinders in the model
end_nodes = myCollection.get_end_nodes()

# Next we pull the x and y coordinates of these nodes
x_y_coords = [(node.x, node.y) for node in end_nodes]

```


 Now lets see what an alpha shape looks like around these points::

```{python}
myCollection.draw("XY", include_alpha_shape=True)
```
<p align="center">
    <img src="./imgs/minimal_tree_alphashape.png" height="390" width="390">
</p>
<p>
As expected, the alpha shape (the outline of which is in grey above) narrowly encloses the tree's branches as compared to a circular region. In the above example, the alpha shape is defined by an alpha value of 1.8. This is the default value, but can be adjusted by the user. The 'draw' function will draw the most recently produced alpha shape, so one can alter this value as shown below:
</p>

```{python}
myCollection.watershed_boundary(curvature_alpha=1)
myCollection.draw("XY", include_alpha_shape=True)
```
<p align="center">
    <img src="./imgs/minimal_tree_low_alpha_XY_hull_tutorial.png" height="390" width="390">
</p>

```{python}
from canopyhydro import geometry

x_y_coords = sorted([(0.1, 0.),(0.0, 0.0),(0.0, 1.0),(.70, 0.0),(1.0, 1.0),(1.0, 0.0),(0.0, 1.0),(1.0, 1.0),(0.25, 0.3),(0.25, 0.8),(0.5, 0.25),(0.5, 0.5),(0.7, 0.35),(0.8, 0.6),(0.5, 0.75),(0.5, 0.5),])


```

<p align="center">
    <img src="./imgs/alpha_shape_scatter.png" height="390" width="390">
</p>



```{python}
```
