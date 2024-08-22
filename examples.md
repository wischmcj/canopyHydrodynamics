---
title: 'CanoPy: Examples'
tags:
  - Python
  - hydrology
  - rainfall
  - tree
  - stemflow
  - throughfall
authors:
  - name: Collin Wischmeyer
    orcid: 0009-0002-8490-0999
---
# Motivation

Vegetation coverage has a marked effect on the spatiotemporal distribution of terrestrial rainfall. This process, refeered to as precipitation partitioning, is well established in the field of hydrology. However, the influence of partitioning on ecological processes is currently under-represented in statistical models due to the inaccessability of the related measurements.
canoPyHydro seeks to empower researchers with partitioning data; allowing them to leverage existing data sets (e.g. increasingly available terrestrial lidar point clouds) to access a wealth of microclimactic data.

# Summary

The 3D models utilized to represent input trees consist of cylinders of various radii, orientations and locations in space (see [QSMs] section below for more information). As such, these cylinders are the discrete units on which we perform our calculations.
The methods provided fall under two major categories of functionality:

1. Methods used to identify cylinders likely to contribute intercepted moisture to various percipitation partitions (i.e. stemflow, throughfall).
2. Methods used to calculate the physical properties (i.e. surface area, 3D orientation) of tree sections represented by groups of cylinders.

Methods in the former category primarily utilize the orientation of the cylinders - in tandem with user defined parameters - to determine the direction in which intercepted water will flow. Once all cylinders have been assigned to a partition, methods in the latter category can then be applied to provide a wide array of statistics describing each partition - most importantly: contributing branch surface area, 2D projected area and canopy coverage area.

// Probably needs overarching desc
// Our focus is on stem flow with a few drip flow features etc

# Cylinder Partitioning

Cylinders in the input QSM are partitioned by traversing a graph representing the flow of intercepted water on the tree being analyzed. To determine the direction in which an intercepted droplet will travel, a 'drip cut-off angle' (see [User Inputs]) is used. That is, water is assumed to flow towards the stem unless the angle of an encountered branch is less than the drip cutoff angle.

```python
import pandas as pd
// add pictures of how a drip cut off angle effects the drip area
// See old power point

```

A digraph is then constructed with each edge representing one cylinder, with edges oriented in the presumed direction of flow under this assumption. This graph is then traversed, and each cylinder is marked as either contributing to stem flow  - if the stem can be reached from its corresponding edge in the graph - or drip flow - if not.

```python
// add prezi figures for demonstration
//
```

# Statistics

  Though a variety of statistics are available, the majority are straight forward; each being a summation of cylinder characteristics. However, custom functions are available for calculating 2D projected area statistics as well as canopy coverage area.

```python
    import pandas as pd
    // add a table of stats for a given tree? Perhaps a list of definitions from email to Nick
```

# Projeted Area

  In the current version of this tool, projection statistics are only available for the coordinate planes: XY, XZ and YZ. For a tree oriented 'right-side-up', these projections represent the tree as seen from above - XY - and from two, perpendicular 'side' views - XZ and YZ - dependent on the orientation of the point cloud data.

```python
    // add 3 pics, tree w/ x/y/z planes plotted as well. Maybe a 3D pic as werll
```

  In a future version of this tool, projection at an arbitrary angle will be available. At that time, it will be possible to incorperate average rain angle as a model parameter.

## Cylinder Projected Area

  Depending on the goals of the user, the projected area of a collection of cylinders can be given as both:

- a simple sum of the projected area of each cylinder

  - Note that this calculation will ignore overlap between cylinder areas
    // Add example
- the total projected area of the collection of cylinders

  - Using this approach, areas in which cylinder projections over lap are only counted once
    // Add Example

## Cylinder Overlap (Shade)

  In order to better understand how the branches of a tree's canopy overlap, more granular overlap information has also been made available via the 'find_overlap_by_percentile' function.0
  When considered from a birds eye view (projecting onto the XY plane), this concept can be understood as a facsimile for the 'shade' cast by branches at a certain height in the tree canopy.
  Consider the below example
    // Add example
    The percentile list is used to determine the height at which to calculate shade. As such, the function will look at the overlap between cylinders in the 75%ile by height (in red) with the remaining cylinders (in blue). The returned values thus represent the 'shade' case by the red cylinders on the blue cylinders.

    Following this logic if either the 0%ile or 100%ile is requested, then there will be no overlap reported. In the former case, all cylinders are included in the red group and therefore there are no blue cylinders on which to cast shade. In the latter case, all cylinders are in the blue group and so there are no red cylinders to cast shade onto the blue cylinders.

  When considered in the XZ or YZ direction, this calculation can be useful in determining the wind exposure at different canopy depths.

## Canopy Coverage Area

  This area might classically be defined by measuring the radius of a trees canopy. As our method focuses on only portions of the tree canopy, it is useful to determine the area spanned by only those portions of the tree canopy. Using this more specific definition of a classic metric, comparisons can be made using related metrics such as woody area index (WAI).
  When considering the coverage area spanned by the stemflow generating portions of the tree, this metric may also be thought of as an analogous concept to a classical 'watershed'.
  This canopy coverage area is calculated by using a concept familiar to the field of point clound analysis - alpha shapes (link).

# Statement of need

* The idea is to describe just the need here, while moving the introduction to the subject to 'Motivation' *
  For any landscape with forest cover, rain must navigate through tree canopies to reach the ground, marking the initial step in terrestrial rainfall pathways. This ‘net rainfall’ influences all subsequent terrestrial hydrological processes, by contributing to runoff (Savenije, 2004), recharging subsurface water pools (Friesen, 2020), or returning to the atmosphere via transpiration (Coenders-Gerrits et al., 2020). Observations reveal substantial spatiotemporal variability in the amounts, timing, and distribution patterns of net rainfall beneath tree canopies, complicating reliable assessments of terrestrial water balances (Van Stan et al., 2020). The complexity of capturing these dynamics leads to statistically-based monitoring designs that are both labor-intensive and costly (e.g., Voss et al., 2016; Zimmermann & Zimmermann, 2014), challenging current practices in modelling and managing terrestrial water interactions (Gutmann, 2020).

The two types of net rainfall fluxes are throughfall and stemflow (Sadeghi et al., 2020). Throughfall describes the portion of rainfall that reaches the ground directly through gaps in the canopy or by dripping off foliage and branches, while stemflow refers to the water that is channeled down tree stems. Despite sophisticated tools available for scanning trees and creating detailed structural models (see references in Wischmeyer et al., 2024), a definitive method to accurately delineate the origins of these fluxes—critical for pinpointing how much rain falls and where—remains elusive.

HydroCanoPy addresses this gap with an innovative graph-theoretic algorithm that refines the process of determining where rainwater travels once it hits a tree canopy. Using terrestrial LiDAR to scan trees and generate Quantitative Structural Models (QSMs) (Hackenberg et al., 2021), HydroCanoPy delineates the drainage areas for stemflow and throughfall drip points. This approach harnesses detailed canopy structure data to map out precise water pathways, transforming how researchers and practitioners predict and analyze rainfall distribution in forested environments.

HydroCanoPy supports a robust suite of analytical tools tailored to varying environmental conditions and research needs. These include adjustable delineation methods based on assumed storm conditions, comprehensive sensitivity analyses to test model robustness, and the ability to isolate specific branch subnetworks for detailed examination. This versatility ensures that HydroCanoPy can adapt to a wide range of ecological scenarios, offering users an unprecedented level of precision in hydrological modeling.

By bridging the gap between advanced canopy scanning technologies and the need for precise hydrological data, HydroCanoPy empowers researchers and environmental managers to enhance their understanding and management of water flows in forested ecosystems, paving the way for more informed conservation and sustainability practices.

# User Inputs

# Input Data

    In the current version of the tool, point cloud data must first be converted to a quantitative structural model - or QSM - using a tool such as SimpleForest (Hackenburg, 2021) prior to being provided as an input.The current version of this tool focuses on canopy structure as a predictor of water availability. To this end, the tool leverages quantitative structural models; a well known format for representing tree branches as cylinders, allowing for point cloud data fron terrestrial lidar scans to be reduced in complexity while preserving the aggregeate characteristics of the tree canopy. Using these "QSM"'s, we are able to provide a variety of statistical an visualization capabilities for our users.

# QSMs

# Upcoming features

  Future iterations will certaily add functionality to integrate additionall real world data (i.e wind speed and direction, rain intensity and angle, etc.).

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }
