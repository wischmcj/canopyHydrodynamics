<<<<<<< HEAD
# canhydro `<sub>`(formerly dripDropFlow)`</sub>`

=======
# canhydro <sub>(formerly dripDropFlow)</sub>
>>>>>>> ab569e7 (Revert "Merge branch 'main' into ingest-validation-data")
This repository houses a bare bones script relating to 'A LiDAR-driven pruning algorithm to delineate canopy drainage areas of stemflow and throughfall drip points.
=======

This repository houses the python code relating to the below paper, and will recieve periodic updates as our research progresses.
>>>>>>> 8aeaf6b (Revert "Revert "Merge branch 'main' into ingest-validation-data"")
A LiDAR-driven pruning algorithm to delineate canopy drainage areas of stemflow and throughfall drip points.
(https://www.researchgate.net/publication/375530854)

## Contents:

The interactive jupyter notebook under '.\Cylinders\cli.ipynb' displays the code written for the first draft of the above linked paper how it was run and reviewed. The remaning (majority) of this repository represents code written in the process of improving and productionalizing that code.

## Improvements since publication of preprint

The code that you see zipped here is representiitive of the ongoing work in our ***productionalizing*** branch, where you will see in progress improvements such as:

+ Functionality refactored into testable methods
+ Flexible, fast visualization functionality added
+ Linter(s) added for formatting and best practices adherence
+ Fully fledged logging functionality added for visibility into results
+ toml configuration added for model parameterization
+ A pytest based testing suite to surface potential disagreements with previous results 

<<<<<<< HEAD
## Setup:
=======
## Development Setup:

>>>>>>> 8aeaf6b (Revert "Revert "Merge branch 'main' into ingest-validation-data"")
- Pre-requisites

<<<<<<< HEAD

### Start-up
  The setup.py that will perform these tasks is not yet complete so setup must be done manually

  1. Activate venv with -  ~\Venvs\CanopyHydroEnv\Scripts\activate.ps1 (PowerShell)
  2. Install requirements - pip install -r requirements.txt
  3. If making edits to the package

     1. Install dev requirements - pip install -r requirements_dev.txt
     2. Install pre-commit - pre-commit install
     
- pre-commit

  1. set-up - https://pre-commit.com/
  2. run on all files - pre-commit run --all-files
=======
  1. Python version 3.9 or higher
  2. A Virtual environment
     - python -m venv ~\Venvs\CanopyHydroEnv
  3. Activate venv with 
      -  ~\Venvs\CanopyHydroEnv\Scripts\activate.ps1 (PowerShell)
  4. Install requirements 
      - pip install -r requirements_dev.txt
  5. Enabling Pre-commit lining
     2. Install pre-commit - pre-commit install
>>>>>>> 8aeaf6b (Revert "Revert "Merge branch 'main' into ingest-validation-data"")

## Known Issues

<<<<<<< HEAD
## Wishlist</h2>
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


# canhydro
This repository houses a bare bones script relating to 'A LiDAR-driven pruning algorithm to delineate canopy drainage areas of stemflow and throughfall drip points.
A LiDAR-driven pruning algorithm to delineate canopy drainage areas of stemflow and throughfall drip points.
(https://www.researchgate.net/publication/375530854)

## Contents:

The interactive jupyter notebook under '.\Cylinders\cli.ipynb' displays the code written for the above linked paper how it was run and reviewed. The remaning (majority) of this repository represents code written in the proess of improving and productionalizing that code (See Upcoming Improvements)

## Upcoming improvements
The code that you see zipped here is representiitive of the ongoing work in our ***productionalizing*** branch, where you will see in progress improvements such as:

  + Functionality refactored into methods
  + Linter(s) added for formatting and best practices adherence 
  + Fully fledged logging functionality
  + toml configuration enabled set up
  + A pytest based testing framework
+ Functionality refactored into methods
+ Linter(s) added for formatting and best practices adherence
+ Fully fledged logging functionality
+ toml configuration enabled set up
+ A pytest based testing framework
=======
- The projection algorithm is an approximation
  - e.g. a vector forming a 3,5,6 triangle (vector from (1,1,1)(4,6,7)) has an angle of 45 degrees or 0.785 rad with the XY plane but the algorithm returns .799 rad
- Cylinders occasionally have more than one possible drip point. Our algorithm chooses just one of theses and lists that as the drip node
  - Note that these possible drip nodes more often than not fall very close together
>>>>>>> 8aeaf6b (Revert "Revert "Merge branch 'main' into ingest-validation-data"")
