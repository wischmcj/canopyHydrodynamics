"""Defines the component parts of the ingested QSM"""

from __future__ import annotations

import calendar
import os
from multiprocessing import Pool
from pathlib import Path
from random import random
from time import sleep
from math import pi, sqrt

import global_vars as vars
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



DIR = vars.DIR

NAME = "Cylinder"


class Cylinder:
    # A cylinder could arguably be a data class, it also could extend Shapely.polygon
    #
    # Its purpose is to act as a 1-1 copy of the rows in the QSM data
    # (rows that represent cylinders in 3-D space)
    #
    #

    # initialize our object level variables for cylider objects
    def __init__(self, CylinderCollection) -> None:
        # Base Attributes from file read in
        # self._cylinderCollection = CylinderCollection
        self.x = np.nan  # len 2 array
        self.y = np.nan  # len 2 array
        self.z = np.nan  # len 2 array
        self.radius = np.nan
        self.length = np.nan
        self.branch_order = np.nan
        self.branch_id = np.nan
        self.volume = np.nan
        self.parent_id = np.nan
        self.rev_branch_order = np.nan
        self.segment_id = np.nan

        self.surface_area = np.nan

        # Calculated based off of projection
        self.projected_data = {
            "XY": {
                "unit_vect": [],
                "rev_unit_vect": [],
                "rise_angle_rad": np.nan,
                "polygon": None,
            },
            "YZ": {
                "unit_vect": [],
                "rev_unit_vect": [],
                "rise_angle_rad": np.nan,
                "polygon": None,
            },
        }

        # graph attributes, all in 3-D
        self.flow_id = np.nan
        self.flow_type = np.nan
        self.begins_at_drip_point = np.nan  # bool
        self.begins_at_divide_point = np.nan  # bool

        # others
        self.stem_path_id = np.nan

            
            run = math.sqrt( a_leg[curr_cyl_id]**2 + b_leg[curr_cyl_id]**2)
            rise = hypo[curr_cyl_id]
            if run==0:
                slope_idx = 1 # straightDown e.g. is in flow 
            else:
                slope_idx = rise/run
            
            sa_idx = 2*np.pi*radius_idx*(radius_idx +len_idx) - 2*np.pi*radius_idx*radius_idx
            pa_idx = poly_idx.area
            vol_idx = 2*np.pi*len_idx*radius_idx*radius_idx
            sa_to_vol_idx = sa_idx/vol_idx
            ang_idx = np.arctan(slope_idx)

    def calc_surface_area(self):
        radius = self.radius
        length = self.length
        sa = 2*np.pi*radius*(radius + length) - 2*np.pi*radius*radius
        return sa

    def create_from_list(self, attrs:list, columns = global_vars.qsm_cols):
        """creates a cylinder corrosponding to that defined by a given row of the qsm (attrs)"""

        extract = lambda attr : attrs[columns[attr]] #pulls a column from the qsm row (attrs) corrosponding to the input attribute

        self.x = [extract('x')[0], extract('x')[1]]
        self.y = [extract('y')[0], extract('y')[1]]
        self.z = [extract('z')[0], extract('z')[1]]
        self.radius = extract('radius')
        self.length =  extract('length')
        self.surface_area = self.calc_surface_area()
        self.branch_order =  extract('branch_order')
        self.branch_id =  extract('branch_id')
        self.volume = extract('volume')
        self.parent_id = extract('parent_id')
        self.rev_branch_order = extract('reverse_branch_order')
        self.segment_id = extract('segment_id')

    def project(self, axis= 'XY')
        # Calculated based off of projection
        self.projected_data = {
            "XY": {
                "unit_vect": [],
                "rev_unit_vect": [],
                "rise_angle_rad": np.nan,
                "polygon": None,
            },
            "YZ": {
                "unit_vect": [],
                "rev_unit_vect": [],
                "rise_angle_rad": np.nan,
                "polygon": None,
            },
        }


    def get_flow_data():
        """Returns the flow ID and flow characteristics of the flow the cyl is contained in"""
        print("Get flow data not written")
