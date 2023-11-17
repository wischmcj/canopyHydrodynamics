# canhydro
Houses code relating to 'A LiDAR-driven pruning algorithm to delineate canopy drainage areas of stemflow and throughfall drip points.
A LiDAR-driven pruning algorithm to delineate canopy drainage areas of stemflow and throughfall drip points.

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

--Graph initialization
443.47s setup    tests/test_efficiency.py::test_min_efficiency
436.90s setup    tests/test_efficiency.py::test_base_efficiency
435.01s setup    tests/test_efficiency.py::test_object_efficiency
121.44s call     tests/test_efficiency.py::test_min_efficiency
2.14s call     tests/test_efficiency.py::test_base_efficiency
1.30s call     tests/test_efficiency.py::test_object_efficiency

--Initialization and finding stem flow component
620.79s setup    tests/test_efficiency.py::test_base_graph_test[4_LargeCollection.csv]
615.27s setup    tests/test_efficiency.py::test_min_graph_test[4_LargeCollection.csv]
613.70s setup    tests/test_efficiency.py::test_obj_graph_test[4_LargeCollection.csv]
93.96s call     tests/test_efficiency.py::test_base_graph_test[4_LargeCollection.csv]
168.80s call     tests/test_efficiency.py::test_min_graph_test[4_LargeCollection.csv]
82.20s call     tests/test_efficiency.py::test_obj_graph_test[4_LargeCollection.csv]

--Initialization through flow calclulation
634.05s setup    tests/test_efficiency.py::test_base_graph_test[4_LargeCollection.csv]
721.71s setup    tests/test_efficiency.py::test_min_graph_test[4_LargeCollection.csv]
682.37s setup    tests/test_efficiency.py::test_obj_graph_test[4_LargeCollection.csv]
121.80s call     tests/test_efficiency.py::test_base_graph_test[4_LargeCollection.csv]
2799.07s call     tests/test_efficiency.py::test_min_graph_test[4_LargeCollection.csv]
112.85s call     tests/test_efficiency.py::test_obj_graph_test[4_LargeCollection.csv]

    
INFO:my-logger:base_graph_test started at 16999 38382

Line #    Mem usage    Increment  Occurrences   Line Contents
=============================================================
    44    127.4 MiB    127.4 MiB           1       @profile
    45                                             def base_graph_test():
    46    127.4 MiB      0.0 MiB           1           forest = Forester()
    47    127.4 MiB      0.0 MiB           1           forest.get_file_names(dir=test_input_dir)
    48    218.8 MiB     91.4 MiB           1           forest.qsm_from_file_names(file_name="4_LargeCollection.csv")
    49    218.8 MiB      0.0 MiB           1           flexible_collection = forest.cylinder_collections[0]
    50    850.4 MiB    631.6 MiB           1           flexible_collection.project_cylinders("XZ")
    51    926.3 MiB     75.9 MiB           1           flexible_collection.initialize_graph_from()
    52   1002.3 MiB     76.0 MiB           1           proj_area = flexible_collection.sum_over_graph()
    53   1320.9 MiB    318.6 MiB           1           flexible_collection.find_flow_components()
    54   1320.9 MiB      0.0 MiB           1           print(proj_area)
Filename: canhydro/local_run.py

INFO:my-logger:min_graph_test started at 16999 39259
0899sec
Line #    Mem usage    Increment  Occurrences   Line Contents
=============================================================
    32    255.3 MiB    255.3 MiB           1       @profile
    33                                             def min_graph_test():
    34    255.3 MiB      0.0 MiB           1           forest = Forester()
    35    255.3 MiB      0.0 MiB           1           forest.get_file_names(dir=test_input_dir)
    36    300.8 MiB     45.5 MiB           1           forest.qsm_from_file_names(file_name="4_LargeCollection.csv")
    37    300.8 MiB      0.0 MiB           1           flexible_collection = forest.cylinder_collections[0]
    38   1010.7 MiB    709.8 MiB           1           flexible_collection.project_cylinders("XZ")
    39   1017.9 MiB      7.2 MiB           1           flexible_collection.initialize_minimal_graph_from()
    40   1020.8 MiB      3.0 MiB           1           proj_area = flexible_collection.sum_over_min_graph()
    41   1078.0 MiB     57.2 MiB           1           flexible_collection.find_flow_components_minimal()
    42   1078.0 MiB      0.0 MiB           1           print(proj_area)

INFO:my-logger:obj_graph_test started at 16999 40158
 807sec
Line #    Mem usage    Increment  Occurrences   Line Contents
=============================================================
    56    182.3 MiB    182.3 MiB           1       @profile
    57                                             def obj_graph_test():
    58    182.3 MiB      0.0 MiB           1           forest = Forester()
    59    182.3 MiB      0.0 MiB           1           forest.get_file_names(dir=test_input_dir)
    60    241.9 MiB     59.6 MiB           1           forest.qsm_from_file_names(file_name="4_LargeCollection.csv")
    61    241.9 MiB      0.0 MiB           1           flexible_collection = forest.cylinder_collections[0]
    62    989.2 MiB    747.3 MiB           1           flexible_collection.project_cylinders("XZ")
    63   1026.4 MiB     37.2 MiB           1           flexible_collection.initialize_object_graph_from()
    64   1030.2 MiB      3.8 MiB           1           proj_area = flexible_collection.sum_over_object_graph()
    65   1252.6 MiB    222.3 MiB           1           flexible_collection.find_flow_components_object()
    66   1252.6 MiB      0.0 MiB           1           print(proj_area)
INFO:my-logger:finished at 16999 40965





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
