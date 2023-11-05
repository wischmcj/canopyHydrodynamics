"""Defines the component parts of the ingested QSM"""

from __future__ import annotations

import calendar
import os
from multiprocessing import Pool
from pathlib import Path
from random import random
from time import sleep

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import toml
from matplotlib.pyplot import cm

# from shapely.geometry import Polygon, Point
# from shapely.ops import unary_union, transform
# from descartes import PolygonPatch
# from mpl_toolkits import mplot3d


# import time
# import copy
# import math
# import openpyxl
# import geopandas as geo


import global_vars as vars

DIR = vars.DIR

NAME = "Cylinder"

class Cylinder:
    # A cylinder could arguably be a data class, it also could extend Shapely.polygon
    #
    # Its purpose is to act as a 1-1 copy of the rows in the QSM data 
    # (rows that represent cylinders in 3-D space)
    #
    #

    #initialize our object level variables for cylider objects
    def __init__(self, CylinderCollection) -> None:
        #Base Attributes from file read in
        # self._cylinderCollection = CylinderCollection
        self.x                  = np.nan #len 2 array
        self.y                  = np.nan #len 2 array
        self.z                  = np.nan #len 2 array
        self.radius             = np.nan
        self.length             = np.nan
        self.branch_order       = np.nan
        self.branch_id          = np.nan
        self.volume             = np.nan
        self.parent_id          = np.nan
        self.rev_branch_order   = np.nan
        self.section_id         = np.nan

        #Calculated based off of projection
        self.projected_data= {
                'XY':{'unit_vect':[],
                        'rev_unit_vect':[],
                        'rise_angle_rad':np.nan,
                        'polygon': None},
                'YZ':{'unit_vect':[],
                        'rev_unit_vect':[],
                        'rise_angle_rad':np.nan,
                        'polygon': None}
        }

        #graph attributes, all in 3-D
        self.flow_id                = np.nan
        self.flow_type              = np.nan
        self.begins_at_drip_point   = np.nan #bool
        self.begins_at_divide_point = np.nan #bool

        #others
        self.stem_path_id = np.nan

    def get_flow_data():
        """Returns the flow ID and flow characteristics of the flow the cyl is contained in"""
        print('Get flow data not written')
