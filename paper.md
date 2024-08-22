---
TEST TITLE
---

<p align="center">
    <img src=".canopyhydro_logo.jpeg" height="390" width="390">
</p>
<h1 align="center">CanoPyHydro</h1>
  <p align="center">
    Leveraging remote sensing to map water availability in tree canopies.
    </p>
</p>
<p align="center">
  <a href="#summary">Summary</a> •
  <a href="#statement-of-need">Statement of Need</a> •
  <a href="#functionality">Functionality</a> •
  <a href="#flow-evaluation">Flow Evaluation</a> •
  <a href="#filtering-and-highlighting">Filtering and Highlighting</a> •
  <a href="#future-direction">Future Direction</a> •
  <a href="#acknowledgements">Acknowledgements</a> •
  <a href="#references">References</a> •
  <a href="examples.md">More Examples</a>
</p>
<p align="center">
  Vegetation coverage has a marked effect on the spatiotemporal distribution of terrestrial rainfall, marking the initial step in terrestrial rainfall pathways. <!-- As we demonstrate below, the importance of this 'precipitation partitioning' is well established in the field of hydrology and is of increasing interest in the modeling of ecological processes.  -->
  Growing interest from hydrologists and ecologists alike has lead to the numerous efforts to characterize these flows. However these typically correlative, regression-based efforts have struggled with the problem's inherit degrees of freedom, drawing consistent scrutiny farom reviewers for their inconsistency and lack of conclusion
  <!-- (can cite levia X2 and  van stan here ) -->
  canoPyHydro has grown from an effort uncover strong mechanistic inferences into the drivers for the variation in these fluxes and drawn interest from fellow researchers as a tool for a variety of related use-cases, largely. No matter the use case, canoPyHydro provides value to researchers by enriching their (evermore widely available) terrestrial lidar scans (TLS) with percipitation partitioning data; giving them an tailored environment to explore canpy water distribution.
</p>

<!-- plant nutrient uptake and leaching (Aubrey, 2020), litter decomposition (Qualls, 2020),  plant microbiome composition and function (Van Stan et al., 2020). -->
<!-- The below is an alternatie version, I think I like the above better, though

Vegetation coverage has a marked effect on the spatiotemporal distribution of terrestrial rainfall. That is, trees have a measurable effect on when and where rain falls. This process, referred to as precipitation partitioning, is well established in the field of hydrology but the influence of partitioning on ecological processes is currently under-represented in statistical models. This is in large part due to the inaccessability of the related measurements.
canoPyHydro seeks to empower researchers with percipitation partitioning data; allowing them to leverage existing data sets (e.g. increasingly available terrestrial lidar point clouds) to access a wealth of microclimactic data. -->

# Summary

The main inputs to canoPyHydro are Qantitative Structural Models (QSMs). These models simplify TLS point clouds to represent trees as collections of topologically ordered cylinders. CanoPyHydro's functionality can broadly categorized into two groups:

1. Utilities for the conversion and exploration of QSMs. Methods in this category:

  - Convert QSMs to python objects and calculate individual cyliner metrics (i.e. surface area, angle(s))
  - Create 2D and 3D visualizations (with robust filtering and highlighting functionality.
  - Surface a variety of spatial metrics (i.e. inter-canopy occlusion)

2. Percipitation partitioning utilities

   - Determining where intercepted percipitation is distributed by each part of the tree
   - Qantify the structure of the canopy watershed

CanoPyHydro takes a novel approach to the latter in particular, reimagining the trees' canopies as watersheds. By identifying the many tributary-like flows flowing down their branches, and areas of run-off where water drips to the forest floor, canoPyHydro unveils previously unexplored structures and characterizes the stemflow and throughfall generating portions of the canopy.

For tree-data explorers, CanoPyHydro's spacial utilities are of particular interest. The use of [Alpha Shapes](https://en.wikipedia.org/wiki/Alpha_shape) empowers users to apply the concept of canopy coverage area in unexpected, new ways. Likewise, CanoPyHydro can provide detailed _intra_-canopy occlusion data for arbitrary cross-sections of the canopy, providing detailed estimates of the shade within the canopy (in the vertical direction) and protection from wind/rain (in the horizontal direction).

# Statement of Need

This ‘net rainfall’ influences all subsequent terrestrial hydrological processes, by contributing to runoff (Savenije, 2004), recharging subsurface water pools (Friesen, 2020), or returning to the atmosphere via transpiration (Coenders-Gerrits et al., 2020). Observations reveal substantial spatiotemporal variability in the amounts, timing, and distribution patterns of net rainfall beneath tree canopies, complicating reliable assessments of terrestrial water balances (Van Stan et al., 2020). Unfortunately, the costly and laborious techniques required to characterise processes challenge current practices in modelling and managing terrestrial water interactions (Gutmann, 2020).
<!--
  Despite this, the complexity of capturing these dynamics necesitates statistically-based monitoring designs that are both labor-intensive and costly (e.g., Voss et al., 2016; Zimmermann & Zimmermann, 2014). -->


The two types of net rainfall fluxes are throughfall and stemflow (Sadeghi et al., 2020). Throughfall describes the portion of rainfall that reaches the ground directly through gaps in the canopy or by dripping off foliage and branches (the latter being referred to as dripflow), while stemflow refers to the water that is channeled down tree stems. Attempts to correlate whole canopy characteristics with stemflow measurements have generated inconslusiver results (Referrence, 1900). Despite these attempts and the sophisticated tools available for scanning trees and creating detailed structural models (see references in Wischmeyer et al., 2024), a definitive method to accurately delineate the origins of these fluxes—critical for pinpointing how much rain falls and where—remains elusive.

CanoPyHydro addresses this gap by employing an innovative, bottom-up approach for estimation precipitation redistribution, supplementing QSMs generated using existing tooling (Hackenberg et al., 2021) with complemntary, graph based models. CanoPyHydro's titular algorithm traverses these graph models, using the resulting labels to percisely delineate drainage areas for stemflow and throughfall drip points. This approach harnesses the detailed canopy structure data availible through terrestrial LiDAR scans to map out precise water pathways, transforming how researchers and practitioners predict and analyze rainfall distribution in forested environments. Furthermore, the tool boasts configuration options that allow for the comparison of rainfall distribution under varying environmental conditions.

CanoPyHydro supports the application of model outputs via a robust suite of analytical tools suitible for a variety of use cases. For example, user friendly filtering capbabilities allow for users to isolate branch subnetworks meeting any arbitrary contition(s) (i.e. only branches with a radius > 10cm, branches with a branch order of 0 within 100cm of the ground,...). These filters may be used in tandem with built in visualization functions to remove or simply highlight specified portions of the tree and to generate descriptive statistics.

By bridging the gap between advanced canopy scanning technologies and the need for precise hydrological data, CanoPyHydro empowers researchers and environmental managers; enhancing their understanding and management of water flows in forested ecosystems and paving the way for more informed conservation and sustainability practices.

# Functionality Call-Outs

<p align="center">
  <a href="#qsms">QSMs</a> •
  <a href="#flow-identification">Flow Identification</a> •
  <a href="#metrics">Metrics</a> •
  <a href="#2d-Projection">2D Projection</a> •
  <a href="#occlusion">Occlusion</a> •
  <a href="#visualization">Visualization</a>
</p>

## QSMs

Quantitative Structural Models are 3D models that approximate the structure of trees via cylinders of various radii, orientations and locations in space.These models are particularly useful in the reduction of point cloud data as they preserve high level structural data, but are much more compact and less computationally intensive. The QSM's used in the creation of CanoPyHydro were generated by passing TLS point cloud data through a program called SimpleForest (Hackenberg et al., 2021) to generate .csv files.

![Point Cloud and QSM](./imgs/PC_QSM_Plot.png)_(Left to right) A point clound redering, a SimpleForest rendering of the related QSM, a canoPyHydro coloring of said QSM_

The below code demonstrates two different ways that canoPyHydro can read in QSMs

```{python}
  # A CylinderCollection object can be initialized directly
  myCollection = CylinderCollection()

  # Using the details of a QSM model stored in example_tree.csv
  # to create a CylinderCollection object
  myCollection.from_csv('example_tree.csv')


  # Alternatively, the 'Forester' class can be used
  myForester = Forester("data/test/")
  print(f"Files available: {list(map(str,myForester.file_names))}")

  ## ... Read in single QSMs
  myForester.qsm_to_collection("example_tree.csv")
  print(len(myForester.cylinder_collections))
```

Additional information can be found in <a href="./docs/examples.md">our qsm documentation</a> and in the documentation for [SimpleForest](https://www.simpleforest.org/pages/tutorials.html).

### 2D Projections

In the current version of this tool, 2D metrics are available for projections onto the coordinate planes: XY, XZ and YZ. For a tree oriented 'right-side-up', these projections represent the tree as seen from above - XY - and from two, perpendicular 'side' views - XZ and YZ - dependent on the orientation of the point cloud data. 2D projections are critical to the functioning of this tool in variety of ways

- Projected 2D cylinder area is of use when calculating the volume of water generated by a given canopy based on rain intensity
- These projections are used to calculate canopy coverage and, by extension, woody area index; Both are key metrics in the study of stemflow
- By comparing 2D projected areas of different branch subsets, canoPyHydro can provide a variety of detailed occlusion data

## Flow Identification

CanoPyHydro's hydrological estimates rely on the classification of QSM cylinders as stemflow contributing or throughfall contributing. The precipitation intercepted by each cylinder is added to a theoretical 'flow', and each flow of water is assumed to flow towards the stem. In the model's simplified view, these flows either reach the stem of the tree or drip to the ground after encountering a cylinder that is too steep to traverse-such points are referred to as 'drip-points'. To identify these 'too-steep' portions of the tree, we choose a 'drip cut-off angle' (configurable by the user) and assume water is only able to flow down branches with an angle greater than the cutoff.
The below diagram demonstrates how a graph based model allows us to use these assumptions to identify which cylinders in a QSM are on some drip-path - and are therefore throughfall contributing - and which are stemflow contributing.

![Flow ID Alogrithm](./imgscanopyhydro_algo_example.png)
_The above diagram shows a minimal example of a QSM to demonstrate the core concepts of canoPyHydro's flow finding algorithm_

The algorithm above assigns an id to each of the flows found with 'stemflow' always recieving and id of 0. These flow ids are stored by the cylinder collection in the variable 'cyl_to_drip', a dictionary keyed by cylinder ids and can later be used for calculating the 'size' of the flow (see the Metrics section below) and for creating various visualizations of the canopy watershed.

The below code demonstrates how the above is done in practice. Details regarding the various objects and functions used can be found in the <a href="./docs/">`docs`</a> section of this repository.

```{python}
  # Initializing a CylinderCollection object
  myCollection = CylinderCollection()
  myCollection.from_csv('example_tree.csv')

  # Setting a cut-off angle (in radians)
  cut_off_angle = -0.166

  # Initializing the graph based model
  myCollection.initialize_digraph_from(in_flow_grade_lim=cut_off_angle)

  # Running the above described algorithm
  myCollection.find_flow_components()

  # Printing the results of the algorithm

  ## Keys are equal to the cylinder ids of the cylinders in our collection
  cyls = myCollection.cylinders
  print(cyl_to_drip)
```

## Flow Quantification

After the flows in a canopy's watershed have been identified, common statistics regarding these flows can be calculated though the use of the 'calculate_flows' function. In this process flows are characterized based off of the aggregate characteristics of the cylinders that contribute intercepted water to them. In this way, flows are discussed as having:

- A number of cylinders
- A projected area, volume and surface area
  - each being the sum of the same for their contained cyliners
- A surface area to volume ratio
- A sum of the angles of their cylinders - This is available to facilitate the calculation of average flow angle for one or many flows
  Most importantly, each non-stem flow also has a unique drip point and drip point location, representing a point in the canopy at which one would expecte water to drip to the ground.
  Utilizing the above metrics, users can glean important information regarding a tree's watershed. For example, the below graphic uses the projected area data for a tree's flows, along with canoPyHydro's visualization capabilities, to mak the location and relative abundance of moisture beneath two tree canopies.

![Tale of Two Trees Drip Map](./imgscanopyhydro_drip_map_tale_of_2_trees.png) Here we see a side by side comparison of two trees identified as having differing canopy hydrodynamics. The circles represent various drip points in the canopy with the shading based on their respective flow's volume

### Visualization

```{python}
    myCollection = CylinderCollection()
    myCollection.from_csv('example_tree.csv')
    myCollection.project_cylinders('XY')
    myCollection.initialize_digraph_from()
    myCollection.find_flow_components()
    myCollection.calculate_flows()
    myCollection.draw('XY', highlight_lambda=lambda:is_stem, save = True, file_name_ext="docs_ex")
    myCollection.draw('XZ', highlight_lambda=lambda:is_stem, save = True, file_name_ext="docs_ex")
```

![Stem Flow Highlight XY](./imgs/example_tree_XY_docs_ex.png) Here we see a side by side comparison of two trees identified as having differing canopy hydrodynamics. The circles repre
![Stem Flow Highlight XZ](./imgs/example_tree_XZ_docs_ex.png) Here we see a side by side comparison of two trees identified as having differing canopy hydrodynamics. The circles represent various drip points in the canopy with the shading based on their respective flow's volumesent various drip points in the canopy with the shading based on their respective flow's volume

## Metrics

Though a variety of metrics are available through this package, the majority are straight forward, summations of cylinder characteristics. Details regarding these metrics and more are available in the <a href="./docs/metrics_definitions.md">metrics definitions</a> in this repository's <a href="./docs/">documentation</a> directory. However, custom functions are available for calculating a few more complicated metrics, which will be highlighted in this section

## Visualizations

from src.canoPyHydro.CylinderCollection import CylinderCollection-

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

### Occlusion

The occlusion of portions of the canopy, as well as the ground itsself has a quantifiable impact on light/UV exposue, surface temperature and wind exposure. In turn, these environmental conditions each impact moisutre availablity via processes such as evapotranspiration. As such, robust utilities for calculating this occlusion are provided to assist in data exploration.

In the calculation of canopy coverage area, we utilize [Alpha Shapes](https://en.wikipedia.org/wiki/Alpha_shape) rather than a circular region. In the vernacular of some popular python packages alpha shapes are referred to as 'hulls', with the tighly fit version used in canoPyHydro considered 'convcave hulls'. This approach to quantifying canopy coverage provides a lower estimate of canopy coverage than would be measured with a smooth circle.with points along the border of the canopy being connected via concave curves. Doing so mitigate the effect of outlier points-those far further from the enter of the encapsulated figure than the average among said boundary points-byincluding less of the unoccluded space between branches. For this reason, this approach is applicable to horizontal 2D projections (XZ, YZ) in addition to vertical/birds eye view (XZ) projections and provides more accurate measures when boundary points are sparse (i.e. for branch subsets0).

### Cylinder Overlap (Shade)

In order to better understand how the branches of a tree's canopy overlap, more granular overlap information has also been made available via the 'find_overlap_by_percentile' function.
When considered from a birds eye view (projecting onto the XY plane), this concept can be understood as a facsimile for the 'shade' cast by branches at a certain height in the tree canopy.
Consider the below example \* Add example
The percentile list is used to determine the height at which to calculate shade. As such, the function will look at the overlap between cylinders in the 75%ile by height (in red) with the remaining cylinders (in blue). The returned values thus represent the 'shade' case by the red cylinders on the blue cylinders.

    Following this logic if either the 0%ile or 100%ile is requested, then there will be no overlap reported. In the former case, all cylinders are included in the red group and therefore there are no blue cylinders on which to cast shade. In the latter case, all cylinders are in the blue group and so there are no red cylinders to cast shade onto the blue cylinders.

When considered in the XZ or YZ direction, this calculation can be useful in determining the wind exposure at different canopy depths.

### Aggregating 2D Area

CanoPyHydro can also provide detailed intra-canopy occlusion data for given heights/depths. For a vertical (XZ) projection, this represents shading by higher branches on lower branches, for hoizontal projections (XZ, YZ) this represents wind exposure (or lack thereof). In future versions, arbitrary projection angles may be used to assist in calclating the effect of occlusion on partitioning in various different weather conditions.
Depending on the goals of the user, the projected area of a collection of cylinders can be given as both:

- a simple sum of the projected area of each cylinder

  - Note that this calculation will ignore overlap between cylinder areas
    - Add example

- the total projected area of the collection of cylinders

  - Using this approach, areas in which cylinder projections over lap are only counted once
    - Add Example

### Filtering and Highlighting

## Canopy Coverage Area

This area might classically be defined by measuring the radius of a trees canopy. As our method focuses on only portions of the tree canopy, it is useful to determine the area spanned by only those portions of the tree canopy. Using this more specific definition of a classic metric, comparisons can be made using related metrics such as woody area index (WAI).
When considering the coverage area spanned by the stemflow generating portions of the tree, this metric may also be thought of as an analogous concept to a classical 'watershed'.

# Future Direction

Future iterations will certaily add functionality to integrate additionall real world data (i.e wind speed and direction, rain intensity and angle, etc.).

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
