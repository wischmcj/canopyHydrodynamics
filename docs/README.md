# CanoPyHydro

The goal of this and future versions of CanoPyHydro is to provide a tool set that empowers researchers and practitioners to gain new perspectives on rainfall distribution in forested environments. A list of publications that have utilized this tool-and influenced its development- can be found at the bottom of this page.

CanoPyHydro provides users access to an innovative, bottom-up approach to estimation precipitation redistribution. By enriching QSM data with additional structure via graph based hydrological models, canoPyHydro allows for the percise delineation of:'

  - Stemflow and throughfall generating areas of the canopy
  - The 'drip points' to which throughfall is directed - complete with their relative volumes
  - 'Divides' and 'confluences' within the canopy that dictate the flow of water through the canopy

The current tool set also boasts several different spacial analysis tools, several of which have been utilized in the study of non-hydrological environmental conditions within tree canopies. These include:

  - Functionality for characterizing the level of obsfucation present at given canopy cross sections
  - Tools for identifying, highlighting and isolating branch subnetworks meeting any arbitrary contition(s)
    - i.e. only branches with a radius > 10cm, branches with a branch order of 0 within 100cm of the ground, ...
  - 2D and 3D visualization functionality to interactively to explore the structure of tree canopies


## Contents:

The interactive jupyter notebook under '.\Cylinders\cli.ipynb' displays the code written for the first draft of the above linked paper how it was run and reviewed. The functionality there-in has been formalized and expanded into the

## Getting Started

1. **Create a Virtual Environment**: Below we use the native 'venv' module to create a virtual environment. This is not strictly necessary, but it is a good practice to keep your project dependencies separate from your system dependencies in a virtual environment.
   ```bash
   python -m venv canHydroVenv
   ```
2. **Activate Your Environment**: Activate your virtual environment to install dependencies and run the project. The commands to activate the virtual environment depend on your operating system and shell. Below are the commands for activating the virtual environment in different operating systems and shells.

```bash
  # if using bash (Mac, Unix)
  source canHydroVenv/bin/activate  
  # if using PowerShell (Windows)
  source canHydroVenv\Scripts\activate.ps1
```

3. **Install canoPyHydro**: canoPHydro is publisehd with PyPA (the python packacing authority). You can install the latest stable release of canoPyHydro using pip. This installs our latest stable release as well as several libraries required for the use of the package's features. canoPyHydro currently supports Python versions 3.9 and up.

```bash
   pip install canoPyHydro
```

4. **Set Configuration Options**: The default configuration file can be found at '/CanopyHydrodynamics/srccanopyhydro/user_def_config.toml'. Configuration options can be set by altering the contents of this file in place. Refer to the getting started guide for more information on configuration options. At this time functionality changes must be made to this file (e.g. a custom file location cannot be set)

```bash
   pip install canoPyHydro
   ```

# Contributing

We welcome contributions to this project! Whether it's reporting a bug, proposing a new feature, or contributing code, we appreciate your help. Here's how you can contribute:

1. **Install Aditional Dependencies**: Some features (linting, git actions, etc.) may require additional dependencies. An additional 'requirements-dev.txt' file has been provided to install these dependencies.

   ```bash
   pip install -r requirements-dev.txt
   ```
2. **Install Pre-commit**: This repository utilizes the ruff pre-commit hook to ensure that all code is linted before being committed. To install pre-commit, run the following commands:

   ```bash
   pip3 install pre-commit
   pre-commit install
   ```
3. **Review the contributing Guidelines **: Check out the documentation, where you can find [contributing guidelines](https://canopyhydrodynamics.readthedocs.io/en/latest/contributing.html). Please note that this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

Thank you for your interest in contributing to our project!

## Publications:

This repository houses python utilities described in 'A LiDAR-driven pruning algorithm to delineate canopy drainage areas of stemflow and throughfall drip points.'
A LiDAR-driven pruning algorithm to delineate canopy drainage areas of stemflow and throughfall drip points.
(https://www.researchgate.net/publication/375530854)


## Packaging

This project was initially packaged with Flit using the the instructions found on the offical python website: https://packaging.python.org/en/latest/tutorials/packaging-projects/.

## Development Setup:

- Pre-requisites

  1. Python version 3.9 or higher
  2. A Virtual environment
     - python -m venv
  3. Activate venv with
      - source venv/bin/activate (zsh, terminal)
      - source venv\Scripts\activate.ps1 (PowerShell)
  4. Install requirements
      - pip install -r requirements_dev.txt
  5. Enabling Pre-commit lining
     2. Install pre-commit - pre-commit install

## Packaging

This project was initially packaged with Flit using the the instructions found on the offical python website: https://packaging.python.org/en/latest/tutorials/packaging-projects/.

## Running

  A command line interface (CLI) is planned for this project, to allow researchers to easily point to, read in and process their lidar scans. Today, however, we have prioritized work needed to create a cahnhydro package, that would allow our code to be utilized in any python project by use of the 'pip install canoPyHydro' command.

  For now, you can get a feel for how the program works in two ways

1. By running the test files, which have been conveniently populated with 'breakpoints' that pause execution and allow for the invesigation of variables
   - This is achieved through running 'pytest tests/test_collection_integration.py'
2. By running the interactive jupyter notebook under '.\Cylinders\cli.ipynb'. This is indeed how most of the data for the paper was generated

## Important Commands

- The projection algorithm is an approximation
- e.g. a vector forming a 3,5,6 triangle (vector from (1,1,1)(4,6,7)) has an angle of 45 degrees or 0.785 rad with the XY plane but the algorithm returns .799 rad

## Known Issues

- The projection algorithm is an approximation
  - e.g. a vector forming a 3,5,6 triangle (vector from (1,1,1)(4,6,7)) has an angle of 45 degrees or 0.785 rad with the XY plane but the algorithm returns .799 rad
- Cylinders occasionally have more than one possible drip point. Our algorithm chooses just one of theses and lists that as the drip node
  - Note that these possible drip nodes more often than not fall very close together

## Wishlist
  - Optimizing the alpha value for alphashapes
      - Can be done locally for areas with different point densities
  - Smoothing cylinders to eliminate false drip points
      -polygon.buffer
  - Creating QSMs from point cloud data
    - would almost certainly need to leverage c++
  - Integrate Point cloud processing libraries like Tree tool
    - https://github.com/porteratzo/TreeTool
  - pip install -U pytreedb
  - A more robust meta manager that stores to a cloud based db
  - Local (maybe also remote) caching
  - 3d plotting

## Tutorials
  The below code can be run at the first breakpoint in the test_collection_integration.py file
  ### Displaying, Filtering and Highlighting
    flexible_collection.draw(plane = 'XZ')
    flexible_collection.draw(plane = 'XZ', a_lambda = lambda: cyl_id>100)
    flexible_collection.draw(plane = 'XZ', filter_lambda = lambda: cyl_id>100)
    flexible_collection.draw(plane = 'XZ', filter_lambda = lambda: cyl_id>50)
    flexible_collection.draw(plane = 'XZ', filter_lambda = lambda: cyl_id>75)
    flexible_collection.draw(plane = 'XZ', filter_lambda = lambda: cyl_id>75, highlight_lambda = lambda:branch_order==2)
    flexible_collection.draw(plane = 'XZ', filter_lambda = lambda: cyl_id>100, highlight_lambda = lambda:branch_order==2)
    flexible_collection.draw(plane = 'XZ', filter_lambda = lambda: cyl_id>100, highlight_lambda = lambda:is_stem)

  ### Draw all projections
    import geopandas as geo  # only import what we need
    import matplotlib.pyplot as plt
    happy_path_projection.project_cylinders('XY')
    happy_path_projection.project_cylinders('XZ')
    happy_path_projection.project_cylinders('YZ')
    xz_poly = [cyl.projected_data['XZ']['polygon'] for cyl in happy_path_projection.cylinders[1:20]]
    xy_poly = [cyl.projected_data['XY']['polygon'] for cyl in happy_path_projection.cylinders[1:20]]
    yz_poly = [cyl.projected_data['YZ']['polygon'] for cyl in happy_path_projection.cylinders[1:20]]
    geoPolys_xy = geo.GeoSeries(xy_poly)
    geoPolys_xz = geo.GeoSeries(xz_poly)
    geoPolys_yz = geo.GeoSeries(yz_poly)
    fig, ax = plt.subplots(3)
    geoPolys_xy.plot(ax=ax[0,0])
    geoPolys_xy.plot(ax=ax[0])
    geoPolys_xz.plot(ax=ax[1])
    geoPolys_yz.plot(ax=ax[2])

## Wishlist

- Optimizing the alpha value for alphashapes
  - Can be done locally for areas with different point densities
- Smoothing cylinders to eliminate false drip points
  -polygon.buffer
- Creating QSMs from point cloud data
  - would almost certainly need to leverage c++
- Integrate Point cloud processing libraries like Tree tool
  - https://github.com/porteratzo/TreeTool
- pip install -U pytreedb
- A more robust meta manager that stores to a cloud based db
- Local (maybe also remote) caching
- 3d plotting