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
- Graph Construction efficiency
  - tests/test_efficiency.py compares run times of various different possible graph constructions
  - As you can see, setting edges literally equal to cylinder objects was technically the most efficient
  - However, given the readability concerns with this method, the slightly less efficient dictionary copying method was chosen

--just set up
443.47s setup    tests/test_efficiency.py::test_min_efficiency
436.90s setup    tests/test_efficiency.py::test_base_efficiency
435.01s setup    tests/test_efficiency.py::test_object_efficiency
121.44s call     tests/test_efficiency.py::test_min_efficiency
2.14s call     tests/test_efficiency.py::test_base_efficiency
1.30s call     tests/test_efficiency.py::test_object_efficiency
0.00s teardown tests/test_efficiency.py::test_object_efficiency
0.00s teardown tests/test_efficiency.py::test_base_efficiency
0.00s teardown tests/test_efficiency.py::test_min_efficiency

--setup and finding stem flow component
620.79s setup    tests/test_efficiency.py::test_base_graph_test[4_LargeCollection.csv]
615.27s setup    tests/test_efficiency.py::test_min_graph_test[4_LargeCollection.csv]
613.70s setup    tests/test_efficiency.py::test_obj_graph_test[4_LargeCollection.csv]
168.80s call     tests/test_efficiency.py::test_min_graph_test[4_LargeCollection.csv]
93.96s call     tests/test_efficiency.py::test_base_graph_test[4_LargeCollection.csv]
82.20s call     tests/test_efficiency.py::test_obj_graph_test[4_LargeCollection.csv]
0.00s teardown tests/test_efficiency.py::test_base_graph_test[4_LargeCollection.csv]
0.00s teardown tests/test_efficiency.py::test_obj_graph_test[4_LargeCollection.csv]
0.00s teardown tests/test_efficiency.py::test_min_graph_test[4_LargeCollection.csv]
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


***test three***


INFO:root:Line #    Mem usage    Increment  Occurrences   Line Contents

INFO:root:=============================================================

INFO:root:    33    957.7 MiB    957.7 MiB           1       @profile

INFO:root:    34                                             def min_graph_test(flexible_collection):

INFO:root:    35    982.9 MiB     25.2 MiB           1           flexible_collection.initialize_minimal_graph()

INFO:root:    36    983.9 MiB      1.0 MiB           1           proj_area = flexible_collection.sum_over_min_graph()

INFO:root:    37   1072.9 MiB     89.0 MiB           1           flexible_collection.find_flow_components_minimal()

INFO:root:    38   1072.9 MiB      0.0 MiB           1           print(proj_area)

INFO:root:



INFO:root:Line #    Mem usage    Increment  Occurrences   Line Contents

INFO:root:=============================================================

INFO:root:    40    987.6 MiB    987.6 MiB           1       @profile

INFO:root:    41                                             def base_graph_test(flexible_collection):

INFO:root:    42   1029.2 MiB     41.6 MiB           1           flexible_collection.initialize_graph()

INFO:root:    43   1031.9 MiB      2.8 MiB           1           proj_area = flexible_collection.sum_over_graph()

INFO:root:    44   1287.3 MiB    255.3 MiB           1           flexible_collection.find_flow_components()

INFO:root:    45   1287.3 MiB      0.0 MiB           1           print(proj_area)

INFO:root:


INFO:root:Line #    Mem usage    Increment  Occurrences   Line Contents

INFO:root:=============================================================

INFO:root:    47    956.6 MiB    956.6 MiB           1       @profile

INFO:root:    48                                             def obj_graph_test(flexible_collection):

INFO:root:    49    990.5 MiB     33.9 MiB           1           flexible_collection.initialize_object_graph()

INFO:root:    50    995.1 MiB      4.7 MiB           1           proj_area = flexible_collection.sum_over_object_graph()

INFO:root:    51   1279.0 MiB    283.9 MiB           1           flexible_collection.find_flow_components_object()

INFO:root:    52   1279.0 MiB      0.0 MiB           1           print(proj_area)

INFO:root:



2799.07s call     tests/test_efficiency.py::test_min_graph_test[4_LargeCollection.csv]
721.71s setup    tests/test_efficiency.py::test_min_graph_test[4_LargeCollection.csv]
682.37s setup    tests/test_efficiency.py::test_obj_graph_test[4_LargeCollection.csv]
634.05s setup    tests/test_efficiency.py::test_base_graph_test[4_LargeCollection.csv]
121.80s call     tests/test_efficiency.py::test_base_graph_test[4_LargeCollection.csv]
112.85s call     tests/test_efficiency.py::test_obj_graph_test[4_LargeCollection.csv]



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
