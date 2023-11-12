"""These constitute our way of 'copying' the QSMs"""


from __future__ import annotations

import calendar
import copy
import logging
import math
import os
import time
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path
from pickle import dump, load
from random import random
from time import sleep

import networkx as nx
import numpy as np
import openpyxl
import pandas as pd  # only import what we need
# import settings
from descartes import PolygonPatch
from matplotlib.pyplot import cm
from mpl_toolkits import mplot3d
from shapely.geometry import Point, Polygon
from shapely.ops import transform, unary_union

from canhydro.Cylinder import Cylinder
from canhydro.global_vars import log, qsm_cols
from canhydro.Plotter import draw_cyls
from canhydro.utils import concave_hull, intermitent_log, lam_filter

NAME = "CylinderCollection"


# By inheriting the Model class, lambda cyl : cyl.branch_order = br CC gains managed functionality- like lambda searching
class CylinderCollection:
    cylinders = defaultdict(list)

    # initialize our object level variables for cylider objects
    def __init__(self) -> None:
        self.file = ""

        # self.collection = CylinderList()

        # Aggregate values from file
        self.surface_area = np.nan
        self.filename = ""
        self.volume = np.nan
        self.avg_sa_to_vol = np.nan
        self.max_branch_order = np.nan
        self.max_reverse_branch_order = np.nan
        self.canopy_scope = np.nan  # desc of canopy
        self.extent = {
            "min": [np.nan, np.nan, np.nan],
            "max": [np.nan, np.nan, np.nan],
        }
        # to populate with x,y,z mins and maxs
        self.aggregate_angle = np.nan
        self.descriptive_vectors = np.nan  # Average, median, mode vectors
        self.treeQualities = {
            "total_psa": -1,
            "tot_hull_area": -1,
            "stem_flow_hull_area": -1,
            "stem_psa": -1,
            "flowStats": -1,
            "dbh": -1,
            "tot_surface_area": -1,
            "stem_surface_area": -1,
        }

        # Projection Attrs
        self.union_poly = None
        self.stem_path_lengths = []
        self.hull = np.nan
        self.stem_hull = np.nan

        # Special case tree attributes
        self.stem_paths = [[]]  # Cyl collection?
        self.trunk = []  # Collection of cylinders? id list?

        # Graph and Attributes
        self.graph = None
        self.flows = [
            {
                "cyls": [],
                "drip_point": id,
                "attributes": {"cyls": 0, "len": 0, "sa": 0, "pa": 0, "as": 0},
            }
        ]
        self.drip_points = {"x": np.nan, "y": np.nan, "z": np.nan, "flow_id": np.nan}
        self.flow_to_drip = {
            0: 1
        }  # A dictionary of flow ids with values equal to their drip node ids
        self.trunk_nodes = []
        self.drip_loc = np.nan

        # Calculations using graph results
        self.stemTotal = {
            "attributes": {"cyls": 0, "len": 0, "sa": 0, "pa": 0, "as": 0},
            "loc": {"x": np.nan, "y": np.nan, "z": np.nan},
        }
        self.divide_points = []
        self.stemPolys = []
        self.compGraphs = []

    def aggregate_characteristics(self):
        """Calculates the summations, averages etc. of cylinder characterictics
        that might be of interest"""
        return True

    def create_cyl(self, arr: list):
        cols = qsm_cols
        attrs = {k: arr[v] for (k, v) in cols.items()}
        cyl = Cylinder(**attrs)
        cyl.create_from_list(arr, cols)
        return cyl

    def from_csv(self, file, aggregate_cyls=True):
        """Initializes a new Cyl Collection based on the data in a QSM
        with the configured column locations"""
        self.file = file
        self.filename = file.name
        log.info(f"Processing {str(file)}")
        # self.arr = pd.read_csv(file, header=0)
        self.arr = np.genfromtxt(file, delimiter=",", skip_header=True)[0:, :-1]
        cylinders = [self.create_cyl(row) for row in self.arr]
        self.cylinders = cylinders

        if aggregate_cyls:
            min_x = np.min([cyl.x[0] for cyl in cylinders])
            min_y = np.min([cyl.y[0] for cyl in cylinders])
            min_z = np.min([cyl.z[0] for cyl in cylinders])
            max_x = np.max([cyl.x[1] for cyl in cylinders])
            max_y = np.max([cyl.y[1] for cyl in cylinders])
            max_z = np.max([cyl.z[1] for cyl in cylinders])
            # Aggregate values from file
            self.no_cylinders = len(cylinders)
            self.surface_area = np.sum([cyl.surface_area for cyl in cylinders])
            self.volume = np.sum([cyl.volume for cyl in cylinders])
            self.max_branch_order = np.max([cyl.branch_order for cyl in cylinders])
            self.max_reverse_branch_order = np.max(
                [cyl.reverse_branch_order for cyl in cylinders]
            )
            self.avg_sa_to_vol = (
                np.sum([cyl.sa_to_vol for cyl in cylinders]) / self.no_cylinders
            )
            self.extent = {
                "min": [min_x, min_y, min_z],
                "max": [max_x, max_y, max_z],
            }

        self.descriptive_vectors = np.nan  # Average, median, mode vectors

        self.theta = np.nan
        log.info(f"{file.name} initialized with {self.no_cylinders} cylinders")

    def project_cylinders(self, plane: str = "XZ"):
        """Projects cylinders onto the specified plane"""
        if plane not in ("XY", "XZ", "YZ"):
            log.info(f"{plane}: invalid value for plane")
        else:
            polys = []
            log.info(f"Projection into {plane} axis begun for file {self.filename}")
            for idx, cyl in enumerate(self.cylinders):
                poly = cyl.get_projection(plane)
                polys.append(poly)
                # print a progress update once every 10 thousand or so cylinders
                intermitent_log(idx, self.no_cylinders, "Cylinder projection: ")
            # Projection Attrs
            self.union_poly = unary_union(polys)
            self.stem_path_lengths = []
            self.pSV = None  # unsure if I want to keep this attr

    def get_collection_data(self):
        cyl_desc = [cyl.__repr__() for cyl in self.cylinders]
        return cyl_desc

    def draw(
        self,
        plane: str = "XZ",
        a_lambda: function = lambda: True,
        highlight: bool = False,
        **args,
    ):
        """Draws cylinders meeting given characteristics onto the specified plane"""
        if plane not in ("XY", "XZ", "YZ"):
            log.info(f"{plane}: invalid value for plane")
        filtered_cyls, matches = lam_filter(
            self.cylinders, a_lambda, return_all=highlight
        )
        to_draw = [cyl.projected_data[plane]["polygon"] for cyl in filtered_cyls]
        log.info(f"{len(to_draw)} cylinders matched criteria")
        self.union_poly = unary_union(to_draw)
        draw_cyls(collection=to_draw, colors=matches, **args)

    def get_dbh():
        start_z = self.extent["min"][2]
        higher_than_breast_radius = lam_filter(
            self.cylinders, a_lambda=lambda: branch_order == 0 and z[0] <= start_z + 1.3
        )
        rbh = np.max(higher_than_breast_radius.radius)
        self.treeQualities["dbh"] = 2 * rbh

    def watershed_boundary(self, plane: str = "XZ", stem_only: bool = False):
        # todo -- check data type of returned (do we really need geo pandas?)
        g = self.graph

        # filter for stem only if option selected
        endNodes = [n for n in g.nodes if g.degree(n) == 1 and n != -1]
        endCyls, _ = lam_filter(self.cylinders, lambda: cyl_id in endNodes)
        centroids = [
            cyl.projected_data[plane]["polygon"].point_on_surface() for cyl in endCyls
        ]
        tot_hull, _ = concave_hull(centroids, 2.2)
        breakpoint()
        # totHullGeo = geo.GeoSeries(tot_hull)
        # canopyCover = totHullGeo.area
        # canopyBoundary = totHullGeo.boundary.length
        self.hull = np.nan
        # self.stem_hull = np.nan

    def initialize_graph(self):
        # Graph and Attributes
        gr = nx.Graph()
        for cyl in self.cylinders:
            attr = cyl.__dict__
            child_node = attr["cyl_id"]
            parent_node = attr["parent_id"]
            gr.add_edge(child_node, parent_node, **attr)
        self.graph = gr

        # self.flows = [
        #     {
        #         "cyls": [],
        #         "drip_point": id,
        #         "attributes": {"cyls": 0, "len": 0, "sa": 0, "pa": 0, "as": 0},
        #     }
        # ]
        # self.drip_points = {"x": np.nan, "y": np.nan, "z": np.nan, "flow_id": np.nan}
        # self.flow_to_drip = {
        #     0: 1
        # }  # A dictionary of flow ids with values equal to their drip node ids
        # self.trunk_nodes = []
        # self.drip_loc = np.nan

        # # Calculations using graph results
        # self.stemTotal = {
        #     "attributes": {"cyls": 0, "len": 0, "sa": 0, "pa": 0, "as": 0},
        #     "loc": {"x": np.nan, "y": np.nan, "z": np.nan},
        # }
        # self.divide_points = []
        # self.stemPolys = []

    def find_flows(self):
        # Graph and Attributes
        self.graph = nx.Graph()

    def identify_stem_paths(self, axis: str):
        # Special case tree attributes
        self.stem_paths = [[]]  # Cyl collection?
        self.trunk = []  # Collection of cylinders? id list?

    def find_trunk_lean(self):
        # to populate with x,y,z mins and maxs
        self.aggregate_angle = np.nan
        return True

    # def read_csv(self, df=pd.DataFrame(), polys=[], projection="XY"):
    #     # Columns [ID?,ParentID?,x1,y1,z1,x2,y2,z2,radius,?,?,lenght,? ,? ,? ,? ,? ,? ,? ,BO]
    #     # Colnums [1  ,2        ,3 , 4, 5, 6, 7, 8,9    ,10,11,12   ,13,14,15,16,17,18,19,20]
    #     # x = x2-x1, y =y2-y1, z=z2-z1
    #     # number of cyliders = cnt values for radius
    #     # Theta  = angle made buy the cylinder axis
    #     self.projection = projection
    #     if "partial" not in self.filename:
    #         self.df_full = pd.read_csv(self.filename, header=0)
    #     else:
    #         self.df_full = df
    #     self.maxBO = np.max(self.BO)
    #     # self.df = self.df_full
    #     self.df = self.df_full  # .iloc[:130,:]
    #     # columns 3 and 6 represent our x values
    #     if projection == "XZ":
    #         self.x = np.transpose(self.df.iloc[:, [3, 6]].to_numpy())
    #         self.y = np.transpose(self.df.iloc[:, [5, 8]].to_numpy())
    #         self.z = np.transpose(self.df.iloc[:, [4, 7]].to_numpy())
    #     elif projection == "YZ":
    #         self.x = np.transpose(self.df.iloc[:, [5, 8]].to_numpy())
    #         self.y = np.transpose(self.df.iloc[:, [3, 6]].to_numpy())
    #         self.z = np.transpose(self.df.iloc[:, [4, 7]].to_numpy())
    #     else:  # 'XY'
    #         self.x = np.transpose(self.df.iloc[:, [3, 6]].to_numpy())
    #         self.y = np.transpose(self.df.iloc[:, [4, 7]].to_numpy())
    #         self.z = np.transpose(self.df.iloc[:, [5, 8]].to_numpy())
    #     # for side view
    #     self.cylID = self.df.iloc[:, 1].to_numpy()
    #     self.pID = self.df.iloc[:, 2].to_numpy()
    #     self.radius = self.df.iloc[:, 9].to_numpy()
    #     self.no_cylinders = self.radius.size
    #     self.cLength = self.df.iloc[:, 12].to_numpy()
    #     self.BO = self.df.iloc[:, 20].to_numpy()
    #     self.branchID = self.df.iloc[:, 4].to_numpy()
    #     self.maxBO = np.max(self.BO)
    #     self.bID = self.df.iloc[:, 24].to_numpy()
    #     if projection == "XZ":
    #         self.dx = self.df.iloc[:, 6].to_numpy() - self.df.iloc[:, 3].to_numpy()
    #         self.dy = self.df.iloc[:, 8].to_numpy() - self.df.iloc[:, 5].to_numpy()
    #         self.dz = self.df.iloc[:, 7].to_numpy() - self.df.iloc[:, 4].to_numpy()
    #     elif projection == "YZ":
    #         self.dx = self.df.iloc[:, 8].to_numpy() - self.df.iloc[:, 5].to_numpy()
    #         self.dy = self.df.iloc[:, 6].to_numpy() - self.df.iloc[:, 3].to_numpy()
    #         self.dz = self.df.iloc[:, 7].to_numpy() - self.df.iloc[:, 4].to_numpy()
    #     else:  # 'XY'
    #         self.dx = self.df.iloc[:, 6].to_numpy() - self.df.iloc[:, 3].to_numpy()
    #         self.dy = self.df.iloc[:, 7].to_numpy() - self.df.iloc[:, 4].to_numpy()
    #         self.dz = self.df.iloc[:, 8].to_numpy() - self.df.iloc[:, 5].to_numpy()
    #     if "partial" in self.filename:
    #         self.pSV = polys
    #     self.theta = np.arctan(self.dz / np.sqrt(self.dx**2 + self.dy**2))
    #     self.output_dir = "".join([DIR, "output/"])
    #     log.info(self.filename + " initialized")

    # def find_flows():
    #     #Replace trunk with single node

    #     #Remove edges with out flow
    #     #Find the connected component with root in it, this is stem flow

    #     #Remove all stem flow from orig graph, these are drip paths
    #         #Find node of minimal height in each component, this is the drip point

    # def network_simplex():
    #     #can be used to calculate flows on graphs with demands
    #     #we could set a demand of generates X volume of flow
