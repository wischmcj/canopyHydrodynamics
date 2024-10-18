---
title: "CanoPyHydro: Leveraging LiDAR to Predict Precipitation Partitioning in Tree Canopies"
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
    affiliation: 1
  - name: Travis E. Swanson
    orcid: 0000-0002-6879-7621
    affiliation: 2
  - name: John T. Van Stan II
    orcid: 000-0002-0692-7064
    corresponding: true
    affiliation: 1

affiliations:
  - name: Department of Biological, Geological, and Environmental Sciences, Cleveland State University, Cleveland OH, USA
    index: 1
  - name: The Water Institute of the Gulf, Baton Rogue LA, USA
    index: 2
date: 31 August 2024
bibliography: paper.bib
---


![canoPyHydro logo](paper/source/imgs/canhydro_logo.png){ align=center width=270 height=400 }

[Summary](#summary)
[Functionality Overview](#functionality)
[Publications and Acknowledgements](#publications-and-acknowledgements)
[Future Direction](#future-direction)

Vegetation coverage has a marked effect on the spatiotemporal
distribution of terrestrial rainfall, marking the initial step in
terrestrial rainfall-to-runoff pathways. Growing interest from hydrologists and environmental scientists has led to numerous attempts to characterize these flows. However, these efforts have largely been correlative and regression-based, lacking clear frameworks to guide meaningful inferences [@Van_Stan:2020]. CanoPyHydro empowers researchers to derive mechanistic inferences into the drivers underlying variability in canopy rainfall drainage fluxes and has garnered interest for its versatility across related use-cases. By integrating precipitation partitioning data with increasingly available terrestrial lidar scans (TLS), canoPyHydro offers a tailored environment to explore canopy water distribution, enhancing the precision and depth of hydrological analyses.

# Summary

Data is ingested by canoPyHydro in the form of quantitative structural models (QSMs), which distill TLS point clouds into topologically ordered cylinders representing a tree's canopy structure. Broadly speaking, CanoPyHydro's functionality is encompassed in two major categories:

*QSM Ingestion and Exploration*
 CanoPyHydro offers robust visualization capabilities for subsetting, highlighting and displaying features of interest. On the quantitative side, a variety of industry standard metrics (i.e. woody area index, trunk lean) are available at the subset level, as well as some more novel metrics. The latter including detailed_intra-canopy shading data as well as the use of [Alpha Shapes](https://en.wikipedia.org/wiki/Alpha_shape) to enable a novel, situationally-improved interpretations of canopy coverage area.

*Characterizing Canopy Watersheds*
CanoPyHydro's novel approach to percipitation partitioning treats tree canopies as watersheds and reveals tributary-like flows within branch networks. It precisely distinguishes which flows reach the trunk (thus becoming stemflow) and delineates areas where water drips to the forest floor (as throughfall).


# Statement of Need

Net rainfall (throughfall + stemflow) reaching the surface beneath plant
canopies influences all subsequent terrestrial hydrological processes,
contributing to runoff [@Savenije:2004], recharging subsurface water
[@Friesen:2020], or returning to the atmosphere via evaporation
[@Coenders-Gerrits:2020]. This redistribution of rainfall has marked effect on a host of environmental processes: plant nutrient uptake and leaching [@Aubrey:2020], litter decomposition [@Qualls:2020], surface runoff [@Gotsch:2018]; [@Ji:2022], soil erosion [@Dunkerley:2020]and plant microbiome composition and function [@Van_Stan-Morris:2020].  However, there is substantial
spatiotemporal variability in net rainfall amount, timing, and
distribution, complicating reliable assessments of terrestrial water
balances [@Van_Stan:2020]. The costly, labor-intensive techniques
required to observe throughfall and stemflow [@Voss:2016;
`@Zimmermann:2014] challenge current approaches to
modelling and managing terrestrial water interactions
[@Gutmann:2020].

Attempts to correlate whole canopy characteristics with stemflow
measurements have produced inconclusive results [@Sadeghi:2020],
despite advances in tree scanning and structural modelling tools (see
references in [@Wischmeyer:2024]). A definitive method for
accurately delineating these flux origins, crucial for understanding
rainfall distribution, remains elusive.

CanoPyHydro addresses this gap with an innovative, bottom-up approach to
estimating precipitation redistribution, supplementing QSMs generated
using existing tools [@Hackenburg:2021] with complemntary,
graph-based models. CanoPyHydro's algorithm traverses these graph models
to percisely delineate stemflow and throughfall drippoint drainage
areas. Leveraging detailed canopy structural data from terrestrial LiDAR
scans (TLS) to map out precise water pathways, CanoPyHydro transforms
how rainfall drainage pathways are predicted and analyzed, offering
configuration options to compare rainfall distribution under varying
environmental conditions.

CanoPyHydro supports the application of model outputs with a robust
suite of analytical tools suitible for various use cases. Its
user-friendly filtering capbabilities allow users to isolate branch
subnetworks based on specified criteria (branches with a radius >10cm, branches with a branch order of 0 within 100cm of the ground,
etc.). These filters may be used in tandem with integrated visualization
functions to further enable the exploration and analysis of tree
structures and hydrological processes.

By bridging the gap between advanced canopy scanning technologies and
the need for precise hydrological insights, CanoPyHydro empowers
researchers and environmental managers, enhancing their understanding of
water flows in forested ecosystems and supporting more informed
conservation and sustainability practices.

# Functionality
 [QSMs](#qsms)    -    [2D Projection](#projection)    -    [Flow Mapping](#flow-identification) -    [Shade](#shading-fraction)    -    [Canopy Coverage Area](#canopy-coverage-area)


## QSMs

Quantitative Structural Models are 3D representations of tree canopies consisting of topologically ordered cylinders. Each cylinder has a radius, an orientation, start/end coordinates and a variety of meta data fields corresponding to a given section of the tree. These models effectively function to isolate canopy structural information from the detailed topological data available in point clound data. The QSMs used in the creation of CanoPyHydro were generated by processing TLS point cloud data using the SimpleForest plug-in [@Hackenburg:2021] for the Computree platform[@Computree_core_team:2024]

![A point cloud, a raw QSM and a QSM with hydrologic characteristics highlighted via canoPyHydro](paper/source/imgs/PC_QSM_Plot.png){ align=center width=500 height=200 }

Figure 2 above llustrates the steps in a workflow leveraging caopPyHydro. A TLS scan is used to generate a point cloud rendering of a tree, then the same point cloud is processed using the Computree plugin *SimpleForest* to generate a Quantitative Structural Model (QSM). This QSM is a lower-resolution representation approximating the tree's branch structure using cylinders
canoPyHydro visualization of the same QSM, highlighting the model’s
features with color to distinguish different hydrological contributions
such as stemflow and throughfall areas.

Additional information can be found in the canoPyHydro documentation under [QSMs](https://canopyhydrodynamics.readthedocs.io/main/key_concepts/qsms.html) and in the documentation for [SimpleForest](https://www.simpleforest.org/pages/tutorials.html).

## Projection

Many of the metrics calculated by CanoPyHydro are derived from
projections of QSMs onto the XY, XZ, and YZ coordinate planes. For a
tree oriented upright, these projections represent the tree from
different perspectives: the XY plane provides a top-down view, while the
XZ and YZ planes offer perpendicular side views, depending on the
orientation of the point cloud data. These 2D projections are essential
for several key functions:

   * The projected 2D cylinder area is used to calculate the yield of
    water \[L m-2\] generated by a given canopy drainage area.
   * Projections are also used to calculate canopy coverage and woody
    area index, both key metrics in the study of canopy precipitation
    partitioning.
   * By comparing the 2D projected areas of different branch subsets,
    CanoPyHydro provides detailed data on in-canopy shading and light
    penetration.

![A 3D cylinder and its XY, XZ and YZ projections](paper/source/imgs/Cylinder_projections_3D.png){ align=center width=250 height=200}

Below, you can see how this concept is extended to whole tree canopies to create visualizations and compare canopy structure.

![A comparison of the XY and XZ projections of two urban trees](paper/source/imgs/canhydro_XZ_XY_Tale_of_2_Trees.png){ align=right height=400 width=400}:

## Flow Identification

CanoPyHydro’s hydrological estimates classify QSM cylinders as
contributing either to stemflow or throughfall. Each cylinder is
assigned to a 'flow' object, which represents the precipitation
intercepted by that cylinder. Water is assumed to flow toward the tree’s
stem unless it encounters a cylinder too steep to traverse. These steep
areas, where water drips off rather than continuing to the stem, are
termed 'drip-points.'

To identify such areas, a user-defined 'drip cut-off angle' is applied,
which assumes water can only flow down branches with angles above the
specified threshold. The below diagram illustrates how graph-based
models use these assumptions to idifferentiate between flows that
contain a drip point (throughfall) and those that do not (stemflow).

 ![A minimal example to demonstrate the core concepts of canoPyHydro's flow finding algorithm.](paper/source/imgs/canhydro_prezi_algo_fig.png){align=center height=350 width=450}

The algorithm assigns an ID to each of the identified flows, with
'stemflow' always receiving an ID of 0. These flow IDs are stored in the
cylinder collection within the 'cyl_to_drip' variable, a dictionary
keyed by cylinder IDs. These IDs can be used later to calculate the flow
'size' (see the Metrics section below) and to generate visualizations of
the canopy watershed.

Detai;s regarding the various objects and functions used can be found in the canoPyHydro [documentation](https://canopyhydrodynamics.readthedocs.io/main/), with .ipynb and .py example files in the [docs](https://github.com/wischmcj/canopyHydrodynamics/tree/v0.1.1/docs/source/examples) section of canoPyHydro git repository.


## Flow Quantification

After the flows in a canopy's watershed have been identified, common
statistics regarding these flows can be calculated though the use of the
'calculate_flows' function. In this process, flows are characterized
based off of the aggregate characteristics their constituent cylinders. In this way, flows are discussed as having:

  *  A number of cylinders
  *  A projected area, volume and surface area
    * the sum of the same for its constituent cylinders
  *  An average angle of decline
  * A unique drip point with corresponding (x,y,z)-coordinates

The last of these, the so-called drip points, represents the location in the canopy at which the water intercepted by the flow's cylinders either drips to the ground (becoming throughfall) or reaches the ground by traversing the stem (becoming stemflow). The below code demonstrates how these flows are identified in practice, and how the data generated can be used to enrich QSM visualizations. The resulting are two images showing views of the same tree from two different angles, with the stemflow contributing cylinders highlighted.

![](paper/source/imgs/example_tree_XY_docs.png){align=left height=200 width=50%} ![](paper/source/imgs/example_tree_XZ_docs.png){align=right height=200 width=50%}

The above demonstrate how this approach can be used to evaluate canopy water availability, but utilizes only a single flow's data - the flow that ends at the tree's base. The below graphic demonstrates how the remaining, throughfall, generating flows can be coupled with rainfall data to map and compare the relative abundance of moisture beneath ree canopies.

![Two trees with differing hydrologic characteristics, with drip points indicated and shaded based on their respective flow's volume.](paper/source/imgs/canhydro_drip_map_tale_of_2_trees.png){align=center height=200 width=450}:

## Metrics

Though a variety of metrics are available through this package, the
majority are straight forward, summations of cylinder characteristics.
Details regarding these metrics and more are available in the
[glossary](https://canopyhydrodynamics.readthedocs.io/main/glossary.html) section of the canoPyHydro documentation. However, there are also custom functions are available for the calculation of more novel metrics. These functions will be highlighted in this section.

## Shading Fraction

The internal shading of the canopy, as well as the ground beneath it,
impacts the energy balance, surface temperature, and wind exposure. In
turn, these environmental conditions influence moisture availability via
related processes like evaporation. To support data exploration,
CanoPyHydro includes robust tools for calculating these shade patterns.

## Canopy Coverage Area

In the calculation of canopy coverage area, we compute [alpha shapes](https://en.wikipedia.org/wiki/Alpha_shape) surrounding canopy boundary points [^1], rather than the traditional circular regions. This approach provides a more precise, often lower, estimate of canopy coverage compared to a circular approximation; an important consideration when assessing throughfall distribution in sparse branch networks.

[^1]: A note on terminology: For a given set of points, the alpha shape with the lowest curvature coefficient thus, tightest fit, is referred to as a concave hull.

# Future Direction
  As we continue to develop CanoPyHydro, several key advancements are planned to broaden its functionality and enhance its utility across a range of research applications. These improvements will focus on integrating real-world environmental data, advancing spatial analysis capabilities to include the creation of QSMs, and optimizing computational efficiency. Below are the primary directions for future development:

   * Leverage existing python libraries for spacial analysis (scipy-spacial, open3d) into canoPyHydro to allow for: the projection of objects at an arbitrary viewing angle, the use of 3D meshes to represent tree structures, and the calculation of canopy coverage area using voxel-based models.

   * Integrate environmental data such as wind speed, wind direction, rain intensity and average rain angle into calcluations for canopy saturation. e. These integrations will enhance the precision of canopy water distribution modeling by factoring in dynamic environmental conditions, allowing researchers to simulate and predict hydrological processes under various weather scenarios

   * Optimizing the computational efficiency of this tool will plan a critical role in
    enabling CanoPyHydro to process a large number of trees. In turn, this will enable the application of CanoPyHydro to larger-scale studies, such as watershed management and forest hydrology research.

Through these future developments, CanoPyHydro will continue to evolve
as a powerful and versatile tool for exploring tree hydrology, driving
new insights into how tree canopies interact with environmental
conditions and contribute to water redistribution in forest ecosystems.

## Publications and Acknowledgements:

CanoPyHydro was developed in the process of authoring [A LiDAR-driven pruning algorithm to delineate canopy drainage areas of stemflow and throughfall drip points.](https://doi.org/10.1111/2041-210X.14378),
which has been accepted for publication by the *British Ecological Society's* ['Methods in Ecology and Evolution'](https://www.britishecologicalsociety.org/publications/journals/methods-in-ecology-and-evolution/). Said paper, and the code within this repository, represents a collaboration between non-academic data professional [Collin Wischmeyer](https://www.linkedin.com/in/collin-wischmeyer-b55659a4), ecohydroloy scholar [Professor John Van Stan](https://expertise.csuohio.edu/csufacultyprofile/detail.cfm?FacultyID=j_vanstan) with notable contributions from industry geo-scientist [Travis Swanson](https://thewaterinstitute.org/our-team/travis-swanson). Likewise, this tool could not exist without the data collected and the ideas put forward by several graduate students working in Cleveland State University's ['Wet Plant Lab'](https://www.researchgate.net/lab/Wet-Plant-Lab-John-Toland-Van-Stan).

This paper's titular publication and related efforts were made possible, in part, by the financial support of US-NSF DEB-2213623.


# References
