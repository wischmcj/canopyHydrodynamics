# Getting started

## Installation

`canoPyHydro` can be installed using pip. This installs our latest stable release with fully-supported features. `canoPyHydro` currently supports Python versions 3.9, 3.10 and 3.11.

```bash
$ pip install canoPyHydro
```

You can also install the latest development version by cloning the GitHub repository and using pip to install from the local directory:

```bash
$ pip install git+https://github.com/canoPyHydro/canoPyHydro.git
```

# Configuration

There are many optional configuration options, but there are only a few that are **necessary** to adjust/check
to ensure the code runs as expected.
The default configuration file can be be found at '/CanopyHydrodynamics/src/canopyhydro/user_def_config.toml'. Configuration options can be set by altering the contents of this file in place. At this time functionality changes must be made to this file (e.g. a custom file location cannot be set)

## QSM File structre
The [qsm] section details the column numbers in which each variable is stored in the input file. To read in this file correctly,
the following columns are required, and each row in the QSM must have a corresponding value for each:
        - cyl_id
        - parent_id
        - x (specify two columns as an array, one for x0 another for x1)
        - y (specify two columns as an array, one for y0 another for y1)
        - z (specify two columns as an array, one for z0 another for z1)
        - radius
        - volume
        - length
the below additional columns are required as well but can be left blank in the input file:
        - branch_order
        - reverse_branch_order
        - segment_id
## Model Parameters
More detail on these parameters can be found in the documentation. For now, all you need to ensure is that both
variables have an integer value
- min_len_drip_flow
    - must be > 0
- min_flow_grade_lim
    - must be > -pi/2 and

## Directories
root_dir
input_dir
output_dir
test_input_dir



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