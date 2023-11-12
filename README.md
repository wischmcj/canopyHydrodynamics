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

<h2>Design Choices</h2>
- Graph Construction efficency 
  - tests/test_efficency.py compares run times of various different possible graph constructions 
  - As you can see, setting edges litterally equal to cylinder objects was techinically the most efficient
  - However, given the readability concerns with this method, the slightly less efficent dictionary copying method was chosen
443.47s setup    tests/test_efficency.py::test_min_efficency
436.90s setup    tests/test_efficency.py::test_base_efficency
435.01s setup    tests/test_efficency.py::test_object_efficency
121.44s call     tests/test_efficency.py::test_min_efficency
2.14s call     tests/test_efficency.py::test_base_efficency
1.30s call     tests/test_efficency.py::test_object_efficency
0.00s teardown tests/test_efficency.py::test_object_efficency
0.00s teardown tests/test_efficency.py::test_base_efficency
0.00s teardown tests/test_efficency.py::test_min_efficency



INFO:root:Line #    Mem usage    Increment  Occurrences   Line Contents

INFO:root:=============================================================

INFO:root:    40    965.6 MiB    965.6 MiB           1       @profile

INFO:root:    41                                             def min_graph_test(flexible_collection):  

INFO:root:    42    990.8 MiB     25.2 MiB           1           flexible_collection.initialize_minimal_graph()

INFO:root:    43    995.2 MiB      4.4 MiB           1           proj_area = flexible_collection.sum_over_min_graph()

INFO:root:    44    995.2 MiB      0.0 MiB           1           print(proj_area)

INFO:root:

INFO:root:Line #    Mem usage    Increment  Occurrences   Line Contents

INFO:root:=============================================================

INFO:root:    53   1026.8 MiB   1026.8 MiB           1       @profile           

INFO:root:    54                                             def base_graph_test(flexible_collection): 

INFO:root:    55   1094.5 MiB     67.7 MiB           1           flexible_collection.initialize_graph()

INFO:root:    56   1097.3 MiB      2.8 MiB           1           proj_area = flexible_collection.sum_over_graph()

INFO:root:    57   1097.3 MiB      0.0 MiB           1           print(proj_area)

INFO:root:


INFO:root:Line #    Mem usage    Increment  Occurrences   Line Contents

INFO:root:=============================================================

INFO:root:    47    995.2 MiB    995.2 MiB           1       @profile

INFO:root:    48                                             def obj_graph_test(flexible_collection): 

INFO:root:    49   1026.0 MiB     30.9 MiB           1           flexible_collection.initialize_object_graph()

INFO:root:    50   1026.8 MiB      0.8 MiB           1           proj_area = flexible_collection.sum_over_object_graph()

INFO:root:    51   1026.8 MiB      0.0 MiB           1           print(proj_area)  

INFO:root:

**************Test two (run in reverse order to ensure ordering is not germain to memory usage )************


INFO:root:Line #    Mem usage    Increment  Occurrences   Line Contents

INFO:root:=============================================================

INFO:root:    40   1066.4 MiB   1066.4 MiB           1       @profile

INFO:root:    41                                             def min_graph_test(flexible_collection):  

INFO:root:    42   1089.9 MiB     23.5 MiB           1           flexible_collection.initialize_minimal_graph()

INFO:root:    43   1091.1 MiB      1.3 MiB           1           proj_area = flexible_collection.sum_over_min_graph()

INFO:root:    44   1091.1 MiB      0.0 MiB           1           print(proj_area)

INFO:root:


INFO:root:Line #    Mem usage    Increment  Occurrences   Line Contents

INFO:root:=============================================================

INFO:root:    53    959.4 MiB    959.4 MiB           1       @profile           

INFO:root:    54                                             def base_graph_test(flexible_collection): 

INFO:root:    55   1030.2 MiB     70.8 MiB           1           flexible_collection.initialize_graph()

INFO:root:    56   1033.3 MiB      3.1 MiB           1           proj_area = flexible_collection.sum_over_graph()

INFO:root:    57   1033.3 MiB      0.0 MiB           1           print(proj_area)

INFO:root:


INFO:root:Line #    Mem usage    Increment  Occurrences   Line Contents

INFO:root:=============================================================

INFO:root:    47   1033.3 MiB   1033.3 MiB           1       @profile

INFO:root:    48                                             def obj_graph_test(flexible_collection): 

INFO:root:    49   1065.9 MiB     32.6 MiB           1           flexible_collection.initialize_object_graph()

INFO:root:    50   1066.4 MiB      0.5 MiB           1           proj_area = flexible_collection.sum_over_object_graph()

INFO:root:    51   1066.4 MiB      0.0 MiB           1           print(proj_area)  

INFO:root:



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
