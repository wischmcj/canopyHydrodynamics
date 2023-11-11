# canhydro
Houses code relating to 'A LiDAR-driven pruning algorithm to delineate canopy drainage areas of stemflow and throughfall drip points.
A LiDAR-driven pruning algorithm to delineate canopy drainage areas of stemflow and throughfall drip points.

<h1>CanopyHydrodymaics</h1>

<h2>Setup:</h2>
- Pre-requisites
  1.  Python version 3.9 or higher
  2.  A Virtual environment
      - python -m venv ~\Venvs\CanopyHydroEnv

- Start-up
  1. Activate venv with -  ~\Venvs\CanopyHydroEnv\Scripts\activate
  2. Install requirements - pip install requirements.txt
  3. Install pre-commit - pre-commit install

- pre-commit
  1. set-up - https://pre-commit.com/
  2. run on all files - pre-commit run --all-files

- pre-commit
  1. set-up - https://pre-commit.com/
  2. run on all files - pre-commit run --all-files


<h2>Contents:</h2>

- canhydro Directory

- Data

<h2>Commands:</h2>


<h2>Approach</h2>


<h2>Known Issues</h2>
 - The projection algorithm is an approximation
 - e.g. a vector forming a 3,5,6 triangle (vector from (1,1,1)(4,6,7)) has an angle of 45 degrees or 0.785 rad with the XY plane but the algorithm returns .799 rad


 <h2>Useful scripts</h2>
 Draw all projections
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
