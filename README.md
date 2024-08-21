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

- Pre-requisites

  1. Python version 3.9 or higher
  2. A Virtual environment
      - python -m venv
  3. Activate venv with
      - source venv/bin/activate (zsh, terminal)
      - source venv\Scripts\activate.ps1 (PowerShell)
  4. Install requirements
      - pip install -r requirements.txt

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
