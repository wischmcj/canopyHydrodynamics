---
title: 'CanoPy: A Python Package for Delineating Tree Canopy Drainage Areas'
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
    corresponding: true
    affiliation: 1
  - name: Travis E. Swanson
    orcid: 0000-0002-6879-7621
    affiliation: 2
  - name: Kevin E. Mueller
    orcid: 0000-0002-0739-7472
    affiliation: 1
  - name: Nicholas R. Lewis
    affiliation: 1
  - name: Jillian Bastock
    affiliation: 1
  - name: John T. Van Stan II
    orcid: 000-0002-0692-7064
    corresponding: true
    affiliation: 1

affiliations:
  - name: Department of Biological, Geological, and Environmental Sciences, Cleveland State University, Cleveland OH, USA
    index: 1
  - name: The Water Institute of the Gulf, Baton Rouge LA, USA
    index: 2
date: 20 July 2024
bibliography: paper.bib
---
s the level of detail represented in computational models continues to increase, localized estimates of these conditions are now of interest in the study community dynamics (Sybil, some year), soil nutrient distribution (Cavelier 1997), even global climatic modeling (Philipp, *need to find the specific paper).

The re-distribution of precipitation via vegetation interception therefore influences a wide array of terrestrial hydrological processes: contributing to run off (Savenije, 2004), recharging subsurface water pools (Friesen, 2020), and returning moisture to the atmosphere via transpiration (Coenders-Gerrits et al., 2020).

CanHydro is inteneded to be a developlent tool for extracting environmental data from Lidar point clouds. As TLS technology is increasingly available, so too is point cloud data. The functions presented here allow researchers to tap into a treasure trove of existing data to empowering new inferencees.

<h1 align="center">CanHydro</h1>
<p align="center">
    <img src="./canhydro_logo.jpeg" height="390" width="390">
</p>

  CanoPyMap

</h1>
  <p align="center">
    Leveraging remote sensing to map water availability in tree canopies.
    <br />
    </p>
</p>
<p align="center">
  <a href="#background">Motivation</a> •
  <a href="#summary">Summary</a> •
  <a href="#statement-of-need">Statement of Need</a> •
  <a href="#best-practice">Best Practice</a> •
  <a href="#credits">Credits</a> •
  <a href="examples.md">More Examples</a>
</p>

# [background]Background

 <!-- For any landscape with vegetation cover, rain must navigate through tree canopies to reach the ground,  -->
  Vegetation coverage has a marked effect on the spatiotemporal distribution of terrestrial rainfall, marking the initial step in terrestrial rainfall pathways. As we will discuss below, the importance of this 'precipitation partitioning' is well established in the field of hydrology and is of increasing interest in the modeling of ecological processes. Despite this, the complexity of capturing these dynamics necesitates statistically-based monitoring designs that are both labor-intensive and costly (e.g., Voss et al., 2016; Zimmermann & Zimmermann, 2014).
  Canhydro seeks to empower researchers with percipitation partitioning data - leveraging evermore widely available terrestrial lidar scans (TLS) to provide detailed estimates of water distribution throughout (and below) tree canopies


<!-- The below is an alternatie version, I think I like the above better though
Vegetation coverage has a marked effect on the spatiotemporal distribution of terrestrial rainfall. That is, trees have a measurable effect on when and where rain falls. This process, referred to as precipitation partitioning, is well established in the field of hydrology but the influence of partitioning on ecological processes is currently under-represented in statistical models. This is in large part due to the inaccessability of the related measurements.
Canhydro seeks to empower researchers with percipitation partitioning data; allowing them to leverage existing data sets (e.g. increasingly available terrestrial lidar point clouds) to access a wealth of microclimactic data. -->

# [summary]Summary
The main input to CanoPyHydro is a Qantitative structural model (QSM); models that represent trees as a collection of topologically ordered cylinders. At a high level, canoPyHydro can be thought of as having two parts:

  1. Utilities for the conversion and exploration of QSMs. Methods:
    - Convert QSMs to python objects and calculate individual cyliner metrics
      - i.e. surface area, angle(s)
    - Create 2D and 3D visualizations, complete with robust filtering and highlighting functionality.
    - Evaluate more complex spacial characteristics
      - i.e. Alpha shapes/woody area index and intra-canopy occlusion (both disussed below)

  2. Percipitation partitioning utilities
    - Determining where intercepted percipitation is distributed by each part of the tree

  CanoPyHydro takes a novel approach to the latter in particular, approaching the problem by imagining the tree's canopy as a watershed. Just as series of tributaries join to form larger and larger bodies of water in a traditional watershed, water flowing down the branches of a modeled tree constitutes 'flow's that combine in sequence to produce stemflow and dripflow. Through the use of a graph traversal algorithm, all cylinders in a given QSM are categorized as contributing to either stemflow. The results of this algorithm may then be used in tandem with cylinder geometric properties to estimate the size of each percipitation partion, plot the distribution of drip flow and characterize the 'stem flow generating' portion of the canopy.

  CanoPyHydro's spacial utilities may be of particular interest independent of hydrodynamics. The use of [Alpha Shapes](https://en.wikipedia.org/wiki/Alpha_shape) is a useful tool in the estimation of canopy coverage area (from a birds-eye or side view.) Likewise, CanoPyHydro can also provide detailed *intra*-canopy occlusion data for given heights/depths, providing detailed estimates of the shade within the canopy in the vertical direction and protection from wind/rain in the horizontal direction

# [statement-of-need]Statement of Need

This ‘net rainfall’ influences all subsequent terrestrial hydrological processes, by contributing to runoff (Savenije, 2004), recharging subsurface water pools (Friesen, 2020), or returning to the atmosphere via transpiration (Coenders-Gerrits et al., 2020). Observations reveal substantial spatiotemporal variability in the amounts, timing, and distribution patterns of net rainfall beneath tree canopies, complicating reliable assessments of terrestrial water balances (Van Stan et al., 2020). Unfortunately, the costly and laborious techniques required to characterise processes challenge current practices in modelling and managing terrestrial water interactions (Gutmann, 2020).

The two types of net rainfall fluxes are throughfall and stemflow (Sadeghi et al., 2020). Throughfall describes the portion of rainfall that reaches the ground directly through gaps in the canopy or by dripping off foliage and branches (the latter being referred to as dripflow), while stemflow refers to the water that is channeled down tree stems. Attempts to correlate whole canopy characteristics with stemflow measurements have generated inconslusiver results (Referrence, 1900). Despite these attempts and the sophisticated tools available for scanning trees and creating detailed structural models (see references in Wischmeyer et al., 2024), a definitive method to accurately delineate the origins of these fluxes—critical for pinpointing how much rain falls and where—remains elusive.
<!-- The lack of such a method leaves a gap in our understanding of water availability and is thus a detriment to modeling efforts in a variety of environmental and ecological processes (i.e. the impact of blue-gree architectural techniques water availability, community dynamics (Phlipp Reference)) to their detriment (Referrence, 1900) -->
<!-- in the comment above i wanted to list out some of the areas that we know this will be used. We have a broad scope of potential users so it could potentially be its own, short paragraph -->

CanoPyHydro addresses this gap by employing an innovative, bottom-up approach-supplementing QSMs generated with existing tooling (Hackenberg et al., 2021) with complemntary, graph based models. CanoPyHydro's titular algorithm traverses these graph models, using the resulting labels to percisely delineate drainage areas for stemflow and throughfall drip points. This approach harnesses the detailed canopy structure data availible through terrestrial LiDAR scans to map out precise water pathways, transforming how researchers and practitioners predict and analyze rainfall distribution in forested environments. Furthermore, the tool boasts configuration options that allow for the comparison of rainfall distribution under varying environmental conditions.

CanoPyHydro supports the application of model outputs via a robust suite of analytical tools suitible for a variety of use cases. For example, user friendly filtering capbabilities allow for users to isolate branch subnetworks meeting any arbitrary contition(s) (i.e. only branches with a radius > 10cm, branches with a branch order of 0 within 100cm of the ground,...). These filters may be used in tandem with built in visualization functions to remove or simply highlight specified portions of the tree and to generate descriptive statistics.

For example use cases, see the [USE CASES] section below.

By bridging the gap between advanced canopy scanning technologies and the need for precise hydrological data, CanoPyHydro empowers researchers and environmental managers; enhancing their understanding and management of water flows in forested ecosystems and paving the way for more informed conservation and sustainability practices.

# [functionality]Functionality

## [qsms]QSMs
Quantitative Structural Models are 3D models that approximate the structure of trees via cylinders of various radii, orientations and locations in space (See  \autoref{QSM_fig} below). These models are particularly useful in the reduction of point cloud data as they preserve high level structural data, but are much more compact and less computationally intensive. The QSM's used in the creation of CanoPyHydro were generated by passing TLS point cloud data through a program called SimpleForest (Hackenberg et al., 2021) to  generate .csv files.

![An example of a point cloud (left) and the resulting QSM (right) figure.\label{QSM_fig}](PC_next_to_QSM.png)
and referenced from text using \autoref{fig:example}.

Additional information can be found in [our qsm documentation](../docs) and in the documentation for [SimpleForest](https://www.simpleforest.org/pages/tutorials.html) website

## [flow-identification]Flow Identification

CanoPyHydro's estimations of precipitation partitioning rely on the identification of 'drip-points' in the canopy. The intercepted by each cylinder is considerd to contribute to a 'flow', and each flow of water is assumed to flow towards the stem.Drip points are locations in the tree canopy that water flowing towards the stem of the tree cannot pass, due to the branch angle becoming too steep. Which partition (stemflow,throughfall) each portion of the tree contributes intercepted precipitation to then relies on their being path from that portion of the tree to the stem of the tree, on which the branch angle allows for the flow of water.  To identify these 'too-steep' portions of the tree, we choose a 'drip cut-off angle' (configurable by the user) such that angles more negative than this cutoff do not allow for the flow of intercepted percipitation.
A digraph is then constructed with each edge representing one cylinder, with edges oriented in the presumed direction of flow  under these assumptions This graph is then traversed, and each cylinder is marked as either contributing to stem flow  - if the stem can be reached from its corresponding edge in the graph - or drip flow - if not. Projected area surface area and other metrisc are then calculated for each flow as the sum of the the properties of its contained cylinders.

## [metrics] Metrics

  Though a variety of metrics are available the majority are straight forward, each being a summation of cylinder characteristics. However, custom functions are available for calculating a few more complicated metrics, such as 2D projection and canopy coverage details. The below

## [2d-projection]2D Projection

  In the current version of this tool, projection statistics are only available for the coordinate planes: XY, XZ and YZ. For a tree oriented 'right-side-up', these projections represent the tree as seen from above - XY - and from two, perpendicular 'side' views - XZ and YZ - dependent on the orientation of the point cloud data.

  In a future version of this tool, projection at an arbitrary angle will be available. At that time, it will be possible to incorperate average rain angle as a model parameter.

### [occlusion] Occlusion

  The occlusion of portions of the canopy, as well as the ground itsself has a quantifiable impact on light/UV exposue, surface temperature and wind exposure. In turn, these environmental conditions each impact moisutre availablity via processes such as evapotranspiration. As such, robust utilities for calculating this occlusion are provided to assist in data exploration.

  In the calculation of canopy coverage area, we utilize [Alpha Shapes](https://en.wikipedia.org/wiki/Alpha_shape) rather than a circular region. In the vernacular of some popular python packages alpha shapes are referred to as 'hulls', with the tighly fit version used in canoPyHydro considered 'convcave hulls'. This approach to quantifying canopy coverage provides a lower estimate of canopy coverage than would be measured with a smooth circle.with points along the border of the canopy being connected via concave curves. Doing so mitigate the effect of outlier points-those far further from the enter of the encapsulated figure than the average among said boundary points-byincluding less of the unoccluded space between branches. For this reason, this approach is applicable to horizontal 2D projections (XZ, YZ) in addition to vertical/birds eye view (XZ) projections and provides more accurate measures when boundary points are sparse (i.e. for branch subsets0).

###  Cylinder Overlap (Shade)

  In order to better understand how the branches of a tree's canopy overlap, more granular overlap information has also been made available via the 'find_overlap_by_percentile' function.
  When considered from a birds eye view (projecting onto the XY plane), this concept can be understood as a facsimile for the 'shade' cast by branches at a certain height in the tree canopy.
  Consider the below example
    * Add example
    The percentile list is used to determine the height at which to calculate shade. As such, the function will look at the overlap between cylinders in the 75%ile by height (in red) with the remaining cylinders (in blue). The returned values thus represent the 'shade' case by the red cylinders on the blue cylinders.

    Following this logic if either the 0%ile or 100%ile is requested, then there will be no overlap reported. In the former case, all cylinders are included in the red group and therefore there are no blue cylinders on which to cast shade. In the latter case, all cylinders are in the blue group and so there are no red cylinders to cast shade onto the blue cylinders.

  When considered in the XZ or YZ direction, this calculation can be useful in determining the wind exposure at different canopy depths.

### [aggregate-2d-area] Aggregating 2D Area
 CanoPyHydro can also provide detailed intra-canopy occlusion data for given heights/depths. For a vertical (XZ) projection, this represents shading by higher branches on lower branches, for hoizontal projections (XZ, YZ) this represents wind exposure (or lack thereof). In future versions, arbitrary projection angles may be used to assist in calclating the effect of occlusion on partitioning in various different weather conditions.
  Depending on the goals of the user, the projected area of a collection of cylinders can be given as both:

- a simple sum of the projected area of each cylinder

  - Note that this calculation will ignore overlap between cylinder areas
    * Add example
- the total projected area of the collection of cylinders

  - Using this approach, areas in which cylinder projections over lap are only counted once
    * Add Example


## [visualizations]Visualizations

from src.canhydro.CylinderCollection import CylinderCollection
```{python}
  # Initializing a CylinderCollection object
  myCollection = CylinderCollection()

  # Converting a specified file to a CylinderCollection object
  myCollection.from_csv('example_tree.csv')

  # Requesting an plot of the tree projected onto the XZ plane ('front' view)
  myCollection.project_cylinders('XZ')
  myCollection.project_cylinders('XZ')
  myCollection.draw('XZ')

  # Requesting an plot of the tree projected onto the YZ plane ('side' view)
  myCollection.project_cylinders('YZ')
  myCollection.draw('YZ')

  # Requesting an plot of the tree projected onto the XY plane (Birds eye view)
  myCollection.project_cylinders('XY')
  myCollection.draw('XY')
```
### [filtering-and-highlighting]Filtering and Highlighting


## Canopy Coverage Area

  This area might classically be defined by measuring the radius of a trees canopy. As our method focuses on only portions of the tree canopy, it is useful to determine the area spanned by only those portions of the tree canopy. Using this more specific definition of a classic metric, comparisons can be made using related metrics such as woody area index (WAI).
  When considering the coverage area spanned by the stemflow generating portions of the tree, this metric may also be thought of as an analogous concept to a classical 'watershed'.

# Statement of need

* The idea is to describe just the need here, while moving the introduction to the subject to 'Motivation' *

This ‘net rainfall’ influences all subsequent terrestrial hydrological processes, by contributing to runoff (Savenije, 2004), recharging subsurface water pools (Friesen, 2020), or returning to the atmosphere via transpiration (Coenders-Gerrits et al., 2020). Observations reveal substantial spatiotemporal variability in the amounts, timing, and distribution patterns of net rainfall beneath tree canopies, complicating reliable assessments of terrestrial water balances (Van Stan et al., 2020). The complexity of capturing these dynamics leads to statistically-based monitoring designs that are both labor-intensive and costly (e.g., Voss et al., 2016; Zimmermann & Zimmermann, 2014), challenging current practices in modelling and managing terrestrial water interactions (Gutmann, 2020).

The two types of net rainfall fluxes are throughfall and stemflow (Sadeghi et al., 2020). Throughfall describes the portion of rainfall that reaches the ground directly through gaps in the canopy or by dripping off foliage and branches, while stemflow refers to the water that is channeled down tree stems. Despite sophisticated tools available for scanning trees and creating detailed structural models (see references in Wischmeyer et al., 2024), a definitive method to accurately delineate the origins of these fluxes—critical for pinpointing how much rain falls and where—remains elusive.

CanoPyHydro addresses this gap with an innovative graph-theoretic algorithm that refines the process of determining where rainwater travels once it hits a tree canopy. Using terrestrial LiDAR to scan trees and generate Quantitative Structural Models (QSMs) (Hackenberg et al., 2021), CanoPyHydro delineates the drainage areas for stemflow and throughfall drip points. This approach harnesses detailed canopy structure data to map out precise water pathways, transforming how researchers and practitioners predict and analyze rainfall distribution in forested environments.

CanoPyHydro supports a robust suite of analytical tools tailored to varying environmental conditions and research needs. These include adjustable delineation methods based on assumed storm conditions, comprehensive sensitivity analyses to test model robustness, and the ability to isolate specific branch subnetworks for detailed examination. This versatility ensures that CanoPyHydro can adapt to a wide range of ecological scenarios, offering users an unprecedented level of precision in hydrological modeling.

By bridging the gap between advanced canopy scanning technologies and the need for precise hydrological data, CanoPyHydro empowers researchers and environmental managers to enhance their understanding and management of water flows in forested ecosystems, paving the way for more informed conservation and sustainability practices.

# User Inputs

# Input Data

    In the current version of the tool, point cloud data must first be converted to a quantitative structural model - or QSM - using a tool such as SimpleForest (Hackenburg, 2021) prior to being provided as an input.The current version of this tool focuses on canopy structure as a predictor of water availability. To this end, the tool leverages quantitative structural models; a well known format for representing tree branches as cylinders, allowing for point cloud data fron terrestrial lidar scans to be reduced in complexity while preserving the aggregeate characteristics of the tree canopy. Using these "QSM"'s, we are able to provide a variety of statistical an visualization capabilities for our users.

# QSMs

# Future Direction

  Future iterations will certaily add functionality to integrate additionall real world data (i.e wind speed and direction, rain intensity and angle, etc.).

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge the support of US-NSF DEB-2213623.

# References

Coenders-Gerrits, A.M.J., Schilperoort, B., Jiménez-Rodríguez, C., 2020. Evaporative Processes on Vegetation: An Inside Look. Precipitation Partitioning by Vegetation: A Global Synthesis. https://doi.org/10.1007/978-3-030-29702-2_3

Friesen, J., 2020. Flow Pathways of Throughfall and Stemflow Through the Subsurface. Precipitation Partitioning by Vegetation: A Global Synthesis. https://doi.org/10.1007/978-3-030-29702-2_13

Cavelier, J., Jaramillo, M., Solis, D., de León, D., 1997. Water balance and nutrient inputs in bulk precipitation in tropical montane cloud forest in Panama. J Hydrol (Amst) 193, 83–96.

Dunkerley, D., 2020. A review of the effects of throughfall and stemflow on soil properties and soil erosion. Precipitation Partitioning by Vegetation: A Global Synthesis.
https://doi.org/10.1007/978-3-030-29702-2_12

Hackenberg, J., Calders, K., Demol, M., Raumonen, P., Piboule, A., Disney, M., 2021. SimpleForest - a comprehensive tool for 3d reconstruction of trees from forest plot point clouds. bioRxiv. bioRxiv.

Sadeghi, S.M.M., Gordon, A.G., Van Stan, J.T., 2020. A Global Synthesis of Throughfall and Stemflow Hydrometeorology. Precipitation Partitioning by Vegetation: A Global Synthesis. https://doi.org/10.1007/978-3-030-29702-2_4

Savenije, H.H.G., 2004. The importance of interception and why we should delete the term evapotranspiration from our vocabulary. Hydrol Process 18, 1507–1511. https://doi.org/10.1002/hyp.5563

Van Stan, J.T., Hildebrandt, A., Friesen, J., Metzger, J.C., Yankine, S.A., 2020. Spatial variablity and temporal stability of local net precipitation patterns. Precipitation Partitioning by Vegetation: A Global Synthesis. https://doi.org/10.1007/978-3-030-29702-2_6

Voss, S., Zimmermann, B., Zimmermann, A., 2016. Detecting spatial structures in throughfall data: The effect of extent, sample size, sampling design, and variogram estimation method. J Hydrol (Amst) 540, 527–537. https://doi.org/10.1016/j.jhydrol.2016.06.042

Wischmeyer, C., Swanson, T., Mueller, K., Lewis, N., Bastock, J., Van Stan, I.J.T., 2024. A LiDAR-driven pruning algorithm to delineate canopy drainage areas of stemflow and throughfall drip points. Methods Ecol Evol In press. https://doi.org/10.2139/ssrn.4600550

Zimmermann, A., Zimmermann, B., 2014. Requirements for throughfall monitoring: the roles of temporal scale and canopy complexity. Agric For Meteorol 189, 125–139.
