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
The default configuration file can be be found at '/CanopyHydrodynamics/srccanopyhydro/user_def_config.toml'. Configuration options can be set by altering the contents of this file in place. At this time functionality changes must be made to this file (e.g. a custom file location cannot be set)

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
